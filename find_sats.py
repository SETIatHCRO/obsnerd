#! /usr/bin/env python
from sopp.utilities import read_datetime_string_as_utc
from sopp.satellites_loader.satellites_loader_from_files import \
    SatellitesLoaderFromFiles
from sopp.event_finder.event_finder_rhodesmill.event_finder_rhodesmill import \
    EventFinderRhodesmill
from sopp.path_finder.observation_path_finder_rhodesmill import \
    ObservationPathFinderRhodesmill
from sopp.custom_dataclasses.observation_target import ObservationTarget
from sopp.custom_dataclasses.facility import Facility
from sopp.custom_dataclasses.coordinates import Coordinates
from sopp.custom_dataclasses.time_window import TimeWindow
from sopp.custom_dataclasses.reservation import Reservation
from sopp.custom_dataclasses.runtime_settings import RuntimeSettings
from sopp.custom_dataclasses.frequency_range.frequency_range import \
    FrequencyRange
from sopp.frequency_filter.frequency_filter import FrequencyFilter
from datetime import datetime, timedelta
import matplotlib.pyplot as plt
from tabulate import tabulate



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


def main(starttime, stoptime, offsettime, frequency, bandwidth=20.0, az_limit=[0, 360],
         el_limit=0.0, ftype='horizon', search_for=False, orbit_type='all',
         tle_file='./satellites.tle', timezone=-8.0):
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
    time_window = TimeWindow(
        begin=read_datetime_string_as_utc(starttime.isoformat()),
        end=read_datetime_string_as_utc(stoptime.isoformat()),
    )
    offsettime = read_datetime_string_as_utc(offsettime.isoformat())

    # Frequency Range
    frequency_range = FrequencyRange(bandwidth=bandwidth, frequency=frequency)

    # Reservation
    reservation = Reservation(
        facility=facility,
        time=time_window,
        frequency=frequency_range
    )

    # Specify Observation Target
    observation_target = ObservationTarget(
        declination='7d24m25.426s',
        right_ascension='5h55m10.3s'
    )

    # Antenna Direction Path (going to do automatically)
    antenna_direction_path = ObservationPathFinderRhodesmill(
        facility,
        observation_target,
        time_window
    ).calculate_path()

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
        time_continuity_resolution=timedelta(seconds=1)
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

    fndctr = 0
    for i, window in enumerate(interference_events, start=1):
        if search_for and search_for not in window.satellite.name:
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
            if pos.time > offsettime and len(table_data) < 10 and not (j % 30):
                local_time = pos.time + timedelta(hours=timezone)
                table_row = [local_time.isoformat(), f"{pos.position.azimuth:0.3f}", f"{pos.position.altitude:0.3f}"]
                table_data.append(table_row)
        if len(table_data):
            print('Frequency information:  ', window.satellite.frequency)
            print('Orbits/day:  ', window.satellite.tle_information.mean_motion.value * 240.0)
            fndctr += 1
            plt.figure('Trajectory')
            plt.plot(az, el)
            plt.plot(az[0], el[0], 'ko')
            print(f'Satellite interference event #{i}:')
            print(f'Satellite: {window.satellite.name}')
            print(tabulate(table_data))
    ps4 = f" and {search_for}" if search_for else ""
    print(f"Found {fndctr} entries for {orbit_type}{ps4}")
    plt.figure('Trajectory')
    plt.xlabel('Az [deg]')
    plt.ylabel('El [deg]')
    plt.grid()
    plt.show()

if __name__ == '__main__':
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument('-t', '--start_time', help="Start time as isoformat string, default is now", default=None)
    ap.add_argument('-d', '--duration', help='Duration in minutes [20]', type=float, default=20.0)
    ap.add_argument('-f', '--frequency', help='Frequency to search in MHz', type=float, default=1575.0)
    ap.add_argument('-b', '--bandwidth', help="Bandwidth in MHz [20]", type=float, default=20.0)
    ap.add_argument('-s', '--search', help="String to search for", default=False)
    ap.add_argument('-o', '--orbit', help='Orbit type: (all, geo, meo, leo) [all]', choices=['all', 'geo', 'meo', 'leo'], default='all')
    ap.add_argument('-e', '--el_limit', help="Lower horizon elevation [20.0]", type=float, default=20.0)
    ap.add_argument('-a', '--az_limit', help="Azimuth range [0,360]", default='0,360')
    ap.add_argument('--offset', help="Number of minutes offset to get positions [10.0]", type=float, default=10.0)
    ap.add_argument('--tz', help='Time zone (hours offset from UTC) [-8]', type=float, default=0.0)
    ap.add_argument('--tle_file', help='Name of tle file', default='tle/active.tle')
    ap.add_argument('--ftype', help='search horizon or beam', choices=['horizon', 'beam'], default='horizon')
    args = ap.parse_args()
    if args.start_time is None:
        args.start_time = datetime.now()
    else:
        args.start_time = datetime.strptime(args.start_time, '%Y-%m-%dT%H:%M') + timedelta(microseconds=1)
    args.start_time -= timedelta(hours=args.tz)  # Convert _to_ UTC
    stop_time = args.start_time + timedelta(minutes=args.duration)
    offset_time = args.start_time + timedelta(minutes=args.offset)
    az_limit = [float(x) for x in args.az_limit.split(',')]

    main(starttime=args.start_time,
         stoptime=stop_time,
         offsettime=offset_time,
         frequency=args.frequency,
         bandwidth=args.bandwidth,
         az_limit = az_limit,
         el_limit = args.el_limit,
         ftype=args.ftype,
         search_for=args.search,
         orbit_type=args.orbit,
         tle_file=args.tle_file,
         timezone=args.tz)
