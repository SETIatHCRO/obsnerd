from sopp.utilities import read_datetime_string_as_utc
from sopp.satellites_loader.satellites_loader_from_files import SatellitesLoaderFromFiles
from sopp.event_finder.event_finder_rhodesmill.event_finder_rhodesmill import EventFinderRhodesmill
from sopp.path_finder.observation_path_finder_rhodesmill import ObservationPathFinderRhodesmill
from sopp.custom_dataclasses.observation_target import ObservationTarget
from sopp.custom_dataclasses.facility import Facility
from sopp.custom_dataclasses.coordinates import Coordinates
from sopp.custom_dataclasses.time_window import TimeWindow
from sopp.custom_dataclasses.reservation import Reservation
from sopp.custom_dataclasses.runtime_settings import RuntimeSettings
from sopp.custom_dataclasses.frequency_range.frequency_range import FrequencyRange
from sopp.frequency_filter.frequency_filter import FrequencyFilter
import datetime
import matplotlib.pyplot as plt
from tabulate import tabulate
from . import onutil


def find_orbit_type(mmdps):
    """
    Parameter
    ---------
    mmdps : float
        Mean motion in degrees/sec
    """
    orbits_per_day = mmdps * 240.0
    if orbits_per_day < 0.85:
        return 'other'
    if orbits_per_day < 1.5:
        return 'geo'
    if orbits_per_day < 5.0:
        return 'meo'
    return 'leo'


def main(start, duration, frequency, bandwidth=20.0, az_limit=[0, 360],
         el_limit=0.0, ftype='horizon', search_for=False, orbit_type='all', time_resolution=1,
         ra='58h48m54s', dec='23d23m24s', number_of_rows_to_show=10, row_cadence = 60.0,
         tle_file='tle/active.tle', timezone=None, output_file=False):
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
    el_limit : float
        El lower limit in degrees
    ftype : str
        Search type 'beam' or 'horizon'
    search_for : str or bool
        If str, only allow if str in satellite name
    orbit_type : str
        All, GEO, MEO, LEO, other
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
    output_file : str or False
        If str, name of ephemerides file
        If False don't write

    """
    # Facility
    facility = Facility(
        Coordinates(
            latitude=40.8178049,
            longitude=-121.4695413,
        ),
        elevation=986,  # meters
        beamwidth=3,    # degrees
        name='HCRO',
    )

    # Observation Window
    starttime = onutil.make_datetime(date=start, tz=timezone)
    stoptime = starttime + datetime.timedelta(minutes=duration)

    time_window = TimeWindow(
        begin=starttime,
        end=stoptime,
    )

    # Frequency Range
    frequency_range = FrequencyRange(bandwidth=bandwidth, frequency=frequency)

    # Reservation
    reservation = Reservation(
        facility=facility,
        time=time_window,
        frequency=frequency_range
    )

    if ftype == 'beam':
        # Specify Observation Target
        observation_target = ObservationTarget(
            right_ascension=ra,
            declination=dec
        )

        # Antenna Direction Path (going to do automatically)
        antenna_direction_path = ObservationPathFinderRhodesmill(
            facility,
            observation_target,
            time_window
        ).calculate_path()
    else:
        antenna_direction_path = None

    # Load Satellites
    all_satellites = SatellitesLoaderFromFiles(
        tle_file=tle_file,
    ).load_satellites()

    # Filter satellites on frequency (optional, going to do automatically)
    filtered_satellites = FrequencyFilter(
        satellites=all_satellites,
        observation_frequency=frequency_range
    ).filter_frequencies()

    # Runtime Settings
    runtime_settings = RuntimeSettings(
        concurrency_level=8,
        time_continuity_resolution=datetime.timedelta(seconds=time_resolution)
    )

    # Display configuration
    print('\nFinding satellite interference events for:\n')
    print(f'Facility: {reservation.facility.name}')
    print(f'Location: {reservation.facility.coordinates} at elevation '
          f'{reservation.facility.elevation}')
    print(f'Reservation start time: {reservation.time.begin}')
    print(f'Reservation end time: {reservation.time.end}')
    print(f'Observation frequency: {reservation.frequency.frequency} MHz')
    if ftype == 'beam':
        print(f'Observing celestial object at: '
              f'Declination: {observation_target.declination} '
              f'Right Ascension:{observation_target.right_ascension}')

    # Determine Satellite Interference
    if ftype == 'horizon':
        interference_events = EventFinderRhodesmill(
            list_of_satellites=filtered_satellites,
            reservation=reservation,
            antenna_direction_path=antenna_direction_path,
            runtime_settings=runtime_settings,
        ).get_satellites_above_horizon()
    else:
        interference_events = EventFinderRhodesmill(
            list_of_satellites=filtered_satellites,
            reservation=reservation,
            antenna_direction_path=antenna_direction_path,
            runtime_settings=runtime_settings,
        ).get_satellites_crossing_main_beam()

    ########################################################################

    print('\n==============================================================\n')
    print(f'There are {len(interference_events)} satellite interference\n'
          f'events during the reservation\n')
    print('==============================================================\n')

    shownctr = 0
    jcadence = int(row_cadence / time_resolution)
    for i, window in enumerate(interference_events, start=1):
        if search_for and search_for.lower() not in window.satellite.name.lower():
            continue
        this_orbit = find_orbit_type(window.satellite.tle_information.mean_motion.value)
        if orbit_type == 'all':
            pass
        elif orbit_type != this_orbit:
            continue

        # max_alt = max(window.positions, key=lambda pt: pt.position.altitude)
        az, el, tae = [], [], []
        table_data = []

        for j, pos in enumerate(window.positions):
            if pos.position.altitude < el_limit:
                continue
            if pos.position.azimuth < az_limit[0] or pos.position.azimuth > az_limit[1]:
                continue
            az.append(pos.position.azimuth)
            el.append(pos.position.altitude)
            tae.append(pos.time)
            if len(table_data) < number_of_rows_to_show and not (j % jcadence):
                table_row = [pos.time.strftime('%Y-%m-%dT%H:%M:%S.%f'), f"{pos.position.azimuth:0.3f}", f"{pos.position.altitude:0.3f}"]
                table_data.append(table_row)
        if len(table_data):
            if output_file:
                fnout = f"{window.satellite.name.replace(' ', '')}.txt"
                print(f"Writing {fnout}")
                with open(fnout, 'w') as fpof:
                    for _t, _a, _e in zip(tae, az, el):
                        print(f"{_t.strftime('%Y-%m-%dT%H:%M:%S.%f')},{_a},{_e},{1.0}", file=fpof)
            print('Frequency information:  ', window.satellite.frequency)
            print('Orbits/day:  ', window.satellite.tle_information.mean_motion.value * 240.0)
            shownctr += 1
            plt.figure('AzEl Trajectory')
            plt.plot(az, el)
            plt.plot(az[0], el[0], 'ko')
            plt.figure('Time Trajectory')
            plt.plot(tae, az)
            plt.plot(tae, el, '--')
            print(f'Satellite interference event #{i}:')
            print(f'Satellite: {window.satellite.name}')
            print(tabulate(table_data))
            if output_file:
                fpof.close()
            output_file = None  # Just write out the first one
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