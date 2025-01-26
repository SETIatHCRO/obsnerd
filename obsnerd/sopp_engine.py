
from sopp.sopp import Sopp
from sopp.builder.configuration_builder import ConfigurationBuilder
from sopp.satellites_filter.filterer import Filterer
from sopp.satellites_filter import filters

from sopp.custom_dataclasses.frequency_range.frequency_range import FrequencyRange
import datetime
import matplotlib.pyplot as plt
from tabulate import tabulate
from . import onutil
import pandas as pd
import numpy as np




def main(start, duration, frequency=None, bandwidth=20.0, az_limit=[0, 360],
         el_limit=0.0, ftype='horizon', search_for=False, orbit_type=None, exclude=False, time_resolution=1,
         ra='58h48m54s', dec='23d23m24s', number_of_rows_to_show=10, row_cadence = 60.0,
         tle_file='tle/active.tle', timezone=None, output_file=False, sat2write=None):
    """
    Parameters
    ----------
    start : onutil.make_datetime
        Time to start observation
    duration : float
        Duration in minutes
    frequency : float
        Frequency in MHz
    bandwith : float
        Bandwidth in MHz
    az_limit : list of float
        Az limits in degrees
    search_for : str or bool
        If str, only allow if str in satellite name
    orbit_type : str
        All, GEO, MEO, LEO, other
    exclude : str or False
        If str, only allow if string not in name
    time_resolution : int
        Time resolution in sec
    ra : str or float
        If ftype == 'beam', RA of observation
    dec : str or float
        If ftype == 'beam', declination of observation
    number_of_row_to_show : int
        Number of rows to show in the ephemeris table
    row_cadence : float
        Cadence of table rows in sec
    tle_file : str
        Name of tle file to use
    timezone : float or None or datetime.timezone
        timezone to use (hours to UTC)
    output_file : bool
        If False don't write
    sat2write : str or None
        If str, then write out the file if satellite name contains this string

    """

    # Filters
    filterer = (
        Filterer()
        .add_filter(filters.filter_frequency(FrequencyRange(bandwidth=bandwidth, frequency=frequency)))
        .add_filter(filters.filter_name_regex(search_for))
        .add_filter(filters.filter_name_does_not_contain(exclude))
        .add_filter(filters.filter_orbit_is(orbit_type))
    )
    print("Verify filterer is ok SE63")
    # filterer = Filterer()
    # """
    # if frequency is not None:
    #     filterer.add_filter(filters.filter_frequency(FrequencyRange(bandwidth=bandwidth, frequency=frequency)))
    # """
    # if search_for:
    #     filterer.add_filter(filters.filter_name_contains(search_for))
    # if exclude:
    #     filterer.add_filter(filters.filter_name_does_not_contain(exclude))
    # if orbit_type in ['leo', 'meo', 'geo']:
    #     filterer.add_filter(getattr(filters, f"filter_is_{orbit_type}")())

    # Observation Window
    starttime = onutil.make_datetime(date=start, tz=timezone)
    stoptime = starttime + datetime.timedelta(minutes=duration)
    configuration = (
        ConfigurationBuilder()
        .set_facility(
            latitude=40.8178049,
            longitude=-121.4695413,
            elevation=986,
            name='HCRO',
            beamwidth=3
        )
        .set_time_window(
            begin=starttime,
            end=stoptime
        )
        .set_frequency_range(
            bandwidth=bandwidth,
            frequency=frequency
        )
        .set_observation_target(
            declination=ra,
            right_ascension=dec
        )
        .set_runtime_settings(
            concurrency_level=8,
            time_continuity_resolution=time_resolution,
            min_altitude=el_limit,
        )
        .set_satellites(tle_file=tle_file)
        .set_satellites_filter(filterer)
        .build()
    )

    # Determine Satellite Interference
    sopp = Sopp(configuration=configuration)

    if ftype == 'horizon':
        events = sopp.get_satellites_above_horizon()
    else:
        events = sopp.get_satellites_crossing_main_beam()

    # Display configuration
    print('\nFinding satellite interference events for:\n')
    print(f'Facility: {configuration.reservation.facility.name}')
    print(f'Location: {configuration.reservation.facility.coordinates} at elevation '
          f'{configuration.reservation.facility.elevation}')
    print(f'Reservation start time: {configuration.reservation.time.begin}')
    print(f'Reservation end time: {configuration.reservation.time.end}')
    print(f'Observation frequency: {configuration.reservation.frequency.frequency} MHz')

    if ftype == 'beam':
        print(f'Observing celestial object at: '
              f'Declination: {configuration.observation_target.declination} '
              f'Right Ascension:{configuration.observation_target.right_ascension}')

    ########################################################################

    print('\n==============================================================\n')
    print(f'There are {len(events)} satellite interference events during the reservation\n')
    print('==============================================================\n')

    shownctr = 0
    jcadence = int(row_cadence / time_resolution)

    ### Frequency info incorporated
    freqData = pd.read_csv('SatList.csv')

    for i, window in enumerate(events, start=1):

        # max_alt = max(window.positions, key=lambda pt: pt.position.altitude)
        az, el, tae, dist = [], [], [], []
        table_data = []

        if sat2write is not None and sat2write not in window.satellite.name:
            continue

        for j, pos in enumerate(window.positions):
            if pos.position.azimuth < az_limit[0] or pos.position.azimuth > az_limit[1]:
                continue
            az.append(pos.position.azimuth)
            el.append(pos.position.altitude)
            tae.append(pos.time)
            dist.append(1000.0 * pos.position.distance_km)
            if len(table_data) < number_of_rows_to_show and not (j % jcadence):
                table_row = [pos.time.strftime('%Y-%m-%dT%H:%M:%S.%f'), f"{pos.position.azimuth:0.3f}", f"{pos.position.altitude:0.3f}"]
                table_data.append(table_row)
        if len(table_data):
            # Query for frequency info
            indFreq = [] 
            try:
                indFreq_id = freqData.query("ID=={}".format(str(window.satellite.tle_information.satellite_number)))["Frequency [MHz]"].values
                indFreq_name = freqData.query("Name=='{}'".format(str(window.satellite.name)))["Frequency [MHz]"].values
                indFreq_ur = list(set(list(indFreq_id) + list(indFreq_name)))
                for freq in indFreq_ur:
                    if(type(freq) != float):
                        indFreq += [float(freq)]
            except:
                indFreq = window.satellite.frequency
            if (frequency != 1575.0):
                #print(frequency)
                #print(type(frequency))
                #print(bandwidth)
                if (len(indFreq) == 0):
                    continue
                freqBools = [((x > (frequency - (bandwidth/2.))) and (x < (frequency + (bandwidth/2.)))) for x in indFreq]
                #print(freqBools)
                if (True not in freqBools):
                    continue
		

            if sat2write is None or sat2write in window.satellite.name:
                fnout = f"{window.satellite.name.replace(' ', '')}.txt"
                print(f"Writing {fnout}")
                with open(fnout, 'w') as fpof:
                    for _t, _a, _e, _d in zip(tae, az, el, dist):
                        print(f"{_t.strftime('%Y-%m-%dT%H:%M:%S.%f')},{_a},{_e},{_d}", file=fpof)
            print('Frequency information:  ', window.satellite.frequency, indFreq)

            print('Orbits/day:  ', window.satellite.tle_information.mean_motion.value * 240.0)
            shownctr += 1
            plt.figure('AzEl Trajectory')
            plt.plot(az, el)
            plt.arrow(az[-2], el[-2], az[-1]-az[-2], el[-1]-el[-2], head_width=2)
            plt.figure('Time Trajectory')
            plt.plot(tae, az)
            plt.plot(tae, el, '--')
            print(f'Satellite interference event #{i}:')
            print(f'Satellite: {window.satellite.name}')
            print(tabulate(table_data))
            if sat2write is None and output_file:
                output_file = None  # Just write out the first one if no satellite name
    ps4 = f" and {search_for}" if search_for else ""
    print(f"Showing {shownctr} entries for {orbit_type}{ps4}")
    plt.figure('AzEl Trajectory')
    plt.xlabel('Az [deg]')
    plt.ylabel('El [deg]')
    plt.grid()
    plt.figure('Time Trajectory')
    plt.xlabel('UTC')
    plt.ylabel('Az/El [deg]')
    plt.grid()
    plt.show()
