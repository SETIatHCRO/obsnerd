#!/usr/bin/env python3
from ATATools import ata_control
import argparse
import atexit
from datetime import datetime


ap = argparse.ArgumentParser()
ap.add_argument('cmd', help="Action [start, move, freq, end, note]", choices=['start', 'freq', 'move', 'end', 'note'])
ap.add_argument('--az', help="Azimuth in deg", type=float, default=121.958)
ap.add_argument('--el', help="Elevation in deg", type=float, default=23.603)
ap.add_argument('-f', '--freq', help="Frequency in MHz", type=float, default=1680.0)
ap.add_argument('-n', '--note', help="Note to add", default=None)
args = ap.parse_args()

ants = ['1a', '1f', '5c']  # List of antennas that we want to use (the USRPs can use 2 antennas at a time only)

# GOES-16
# cfreq = 1680
# az, el = (121.958, 23.603)

fplog = open('onlog.log', 'a')

if args.cmd == 'start':
    ata_control.move_ant_group(ants, 'none', 'atagr')
    print(f"{datetime.now()} -- start {', '.join(ants)}", file=fplog)
elif args.cmd == 'end':
    atexit.register(ata_control.move_ant_group, ants, 'atagr', 'none')
    atexit.register(ata_control.park_antennas, ants)
    print(f"{datetime.now()} -- end {', '.join(ants)}", file=fplog)
elif args.cmd == 'freq':
    ata_control.set_freq(args.freq, [ants[0]], lo='d')
    ata_control.autotune([ants[0]])
    att = 20  # Attenuation in dB
    ata_control.rf_switch_thread(ants)
    ata_control.set_atten_thread([[f'{ant}x', f'{ant}y'] for ant in ants], [[att, att] for ant in ants])
    with open('metadata.yaml', 'w') as fp:
        print(f"fcen: {args.freq}", file=fp)
    print(f"{datetime.now()} -- freq {args.freq}", file=fplog)
elif args.cmd == 'move':
    ata_control.set_az_el([ants[0]], args.az, args.el)
    print(f"{datetime.now()} -- az {args.az}, el {args.el}", file=fplog)
elif args.cmd == 'note':
    print(f"{datetime.now()} -- Note: {args.note}", file=fplog)

fplog.close()