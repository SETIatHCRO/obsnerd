#!/home/obsuser/miniconda3/envs/ATAobs/bin/python
import argparse
import atexit
from ATATools import ata_control, logger_defaults
from SNAPobs import snap_if
import time
import logging

from SNAPobs.snap_hpguppi import record_in as hpguppi_record_in
from SNAPobs.snap_hpguppi import snap_hpguppi_defaults as hpguppi_defaults
from SNAPobs.snap_hpguppi import auxillary as hpguppi_auxillary
from SNAPobs import snap_config

ap = argparse.ArgumentParser()
ap.add_argument('-t', '--target', help='Name of target', default='Jupiter')
ap.add_argument('-c', '--calibrators', help='List of calibrators (1 or 2)', default='3c454.3,3c454.3')
ap.add_argument('-1', '--freq1', help='Freq 1 in MHz', default=1500, type=float)
ap.add_argument('-2', '--freq2', help='Freq 2 in MHz', default=6000, type=float)
ap.add_argument('--int_cal', help="Integration times on calibrators in minutes (match calibrators length)", default='10,10')
ap.add_argument('--int_target', help='Integration time on target in minutes', default=30.0, type=float)
ap.add_argument('-n', '--number_int', help='Number of secondary/target integrations', default=4, type=int)
ap.add_argument('-b', '--known_bad', help='List of known bad antennas', default='2k,2b,5b')
ap.add_argument('--skip_initial_cal', help="Flag to skip initial cal/wait", action='store_true')
args = ap.parse_args()
args.calibrators = args.calibrators.split(',')
args.int_cal = [float(x) for x in args.int_cal.split(',')]
args.known_bad = args.known_bad.split(',')


def main(target, calibrators, freq1, freq2, int_target, int_cal, number_int, known_bad, skip_initial_cal):
    logger = logger_defaults.getProgramLogger("observe", loglevel=logging.INFO)
    obs_start_delay = 10
    obs_fiddle = 5
    lo_list = ["b", "c"]

    ant_list = snap_config.get_rfsoc_active_antlist()
    for badun in known_bad:
        if badun in ant_list:
            print(f"Removing antenna {badun}")
            ant_list.remove(badun)

    antlo_list = [ant+lo.upper() for lo in lo_list for ant in ant_list]

    ata_control.reserve_antennas(ant_list)
    atexit.register(ata_control.release_antennas, ant_list, False)

    # Set integration length
    d = hpguppi_defaults.hashpipe_targets_LoB.copy()
    d.update(hpguppi_defaults.hashpipe_targets_LoC)
    
    keyval_dict = {'XTIMEINT': int_target}
    hpguppi_auxillary.publish_keyval_dict_to_redis(keyval_dict, d, postproc=False)

    # Set LO b and c
    freqs_b = [freq1]*len(ant_list)
    freqs_c = [freq2]*len(ant_list)
    ata_control.set_freq(freqs_b, ant_list, lo='b', nofocus=True)
    ata_control.set_freq(freqs_c, ant_list, lo='c')
    time.sleep(20)


    if not skip_initial_cal:
        # Observe primary calibrator
        source = calibrators[0]
        ata_control.make_and_track_ephems(source, ant_list)

        # RF + IF tune
        ata_control.autotune(ant_list)
        snap_if.tune_if_antslo(antlo_list)

        # Start recording -- record_in does NOT block
        obs_time = 60 * int_cal[0]
        hpguppi_record_in.record_in(obs_start_delay, obs_time, hashpipe_targets = d)
        print(f"Observing {source} for {obs_time:.2f} seconds...", end='')
        time.sleep(obs_time + obs_start_delay + obs_fiddle)
        print("Done")

        # Calibrate + apply delays
        _ = input("Are we done with calibration?")


    for i in range(number_int):
        # Observe secondary calibrator
        source = calibrators[-1]
        ata_control.make_and_track_ephems(source, ant_list)
        # Start recording -- record_in does NOT block
        obs_time = 60 * int_cal[-1]
        hpguppi_record_in.record_in(obs_start_delay, obs_time, hashpipe_targets = d)
        print(f"Observing {source} for {obs_time:.2f} seconds...", end='')
        time.sleep(obs_time + obs_start_delay + obs_fiddle)
        print("Done")


        #Observe target
        source = target
        ata_control.make_and_track_ephems(source, ant_list)
        # Start recording -- record_in does NOT block
        obs_time = 60 * int_target
        hpguppi_record_in.record_in(obs_start_delay, obs_time, hashpipe_targets = d)
        print(f"Observing {source} for {obs_time:.2f} seconds...", end='')
        time.sleep(obs_time + obs_start_delay + obs_fiddle)
        print("Done")


if __name__ == "__main__":
    main(**vars(args))
