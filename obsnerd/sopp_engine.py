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
         tle_file='./satellites.tle', timezone=-8.0, output_file=None):
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
    if output_file is not None:
        fpof = open(output_file, 'w')
        print("FS141:  get other info like distance and X,Y,Z !!!")
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
            local_time = (pos.time + timedelta(hours=timezone)).strftime('%Y-%m-%dT%H:%M:%S.%f')
            if output_file is not None:
                print(f"{local_time},{pos.position.azimuth},{pos.position.altitude},{1.0}", file=fpof)
            if pos.time > offsettime and len(table_data) < 10 and not (j % 30):
                table_row = [local_time, f"{pos.position.azimuth:0.3f}", f"{pos.position.altitude:0.3f}"]
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