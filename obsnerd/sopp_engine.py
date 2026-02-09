
from sopp.sopp import Sopp
from sopp.builder.configuration_builder import ConfigurationBuilder
from sopp.satellites_filter.filterer import Filterer
from sopp.satellites_filter import filters
from sopp.custom_dataclasses.frequency_range.frequency_range import FrequencyRange
import datetime
import matplotlib.pyplot as plt
from tabulate import tabulate
from odsutils import ods_timetools as timetools
from obsnerd.on_observation import Observation
from astropy.coordinates import SkyCoord, EarthLocation
import astropy.units as u
from numpy import array


def main(satname, start, duration, frequency=None, bandwidth=20.0, az_limit=[-180, 360],
         el_limit=0.0, DTC_only=True, ftype='horizon', orbit_type=None, exclude=False, time_resolution=10,
         ra='58h48m54s', dec='23d23m24s', number_of_rows_to_show=10, row_cadence = 60.0,
         init_file='parameters.yaml:observations',
         tle_file='tle/active.tle', verbose=True, show_plots=True):
    """
    Parameters
    ----------
    start : ods_timetools.interpret_date
        Time to start observation
    duration : float
        Duration in minutes
    satname : str or None
        If str, then only show this satellite and write it out
    frequency : float
        Frequency in MHz
    bandwith : float
        Bandwidth in MHz
    az_limit : list of float
        Az limits in degrees
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

    """
    if satname == '*':
        satname = ''
    tracks = {}
    # Filters
    filterer = (
        Filterer()
        .add_filter(filters.filter_frequency(FrequencyRange(bandwidth=bandwidth, frequency=frequency)))
        .add_filter(filters.filter_name_regex(satname))
        .add_filter(filters.filter_name_does_not_contain(exclude))
        .add_filter(filters.filter_orbit_is(orbit_type))
    )

    # Observation Window
    starttime = timetools.interpret_date(start, fmt='datetime')
    stoptime = starttime + datetime.timedelta(minutes=duration)
    configuration = (
        ConfigurationBuilder()
        .set_facility(
            latitude=40.8178049,
            longitude=-121.4695413,
            elevation=1019.0,
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
    location = EarthLocation(lat=configuration.reservation.facility.coordinates.latitude*u.deg,
                             lon=configuration.reservation.facility.coordinates.longitude*u.deg,
                             height=configuration.reservation.facility.elevation*u.m)
    # Determine Satellite Interference
    sopp = Sopp(configuration=configuration)
    events = sopp.get_satellites_above_horizon() if ftype == 'horizon' else sopp.get_satellites_crossing_main_beam()

    # Display configuration
    if verbose:
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
    print(f'There are {len(events)} satellite interference events during the reservation')
    jcadence = int(row_cadence / time_resolution)

    ### Frequency info incorporated
    for i, window in enumerate(events, start=1):
        if DTC_only and 'DTC' not in window.satellite.name:
            continue

        # max_alt = max(window.positions, key=lambda pt: pt.position.altitude)
        az, el, time, dist = [], [], [], []
        table_data = []
        for j, pos in enumerate(window.positions):
            if pos.position.azimuth < az_limit[0] or pos.position.azimuth > az_limit[1]:
                continue
            az.append(pos.position.azimuth)
            el.append(pos.position.altitude)
            time.append(pos.time)
            dist.append(1000.0 * pos.position.distance_km)
            if len(table_data) < number_of_rows_to_show and not (j % jcadence):
                table_row = [pos.time.strftime('%Y-%m-%dT%H:%M:%S.%f'), f"{pos.position.azimuth:0.3f}", f"{pos.position.altitude:0.3f}"]
                table_data.append(table_row)
        if len(az) < 3:
            continue
        srcname = window.satellite.name.replace(' ','').replace('[', '').replace(']', '').replace('-', '')
        eventid = f"{srcname}{i}"
        this_obs = Observation(source=eventid, ptinit=init_file)
        time = timetools.Time(time)
        this_obs.update(az=array(az)*u.deg, el=array(el)*u.deg, time=time, dist=array(dist)*u.m)
        sky = SkyCoord(alt=this_obs.el, az=this_obs.az, obstime=this_obs.time, frame='altaz', location=location)
        this_obs.update(ra=sky.gcrs.ra, dec=sky.gcrs.dec)
        this_obs.calc_properties()
        tracks.setdefault(srcname, [])
        tracks[srcname].append(this_obs)
        if verbose:
            print('Orbits/day:  ', window.satellite.tle_information.mean_motion.value * 240.0)
            fnout = f"{this_obs.source.strip()}.txt"
            print(f"Writing {fnout}")
            with open(fnout, 'w') as fpof:
                for _t, _a, _e, _d in zip(time, az, el, dist):
                    print(f"{_t.strftime('%Y-%m-%dT%H:%M:%S.%f')},{_a},{_e},{_d}", file=fpof)
            print(f'Satellite interference event #{i}:')
            print(f'Satellite: {window.satellite.name}')
            print(tabulate(table_data))
        if show_plots:
            plt.figure('AzEl Trajectory')
            plt.plot(az, el)
            plt.arrow(az[-2], el[-2], az[-1]-az[-2], el[-1]-el[-2], head_width=2)
            plt.figure('Time Trajectory')
            plt.plot(time, az)
            plt.plot(time, el, '--')
    if verbose and frequency:
        print('Frequency information:  ', window.satellite.frequency, frequency)
    if show_plots:
        plt.figure('AzEl Trajectory')
        plt.xlabel('Az [deg]')
        plt.ylabel('El [deg]')
        plt.grid()
        plt.figure('Time Trajectory')
        plt.xlabel('UTC')
        plt.ylabel('Az/El [deg]')
        plt.grid()
        plt.show()
    return tracks