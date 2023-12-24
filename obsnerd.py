#!/usr/bin/python3
from ATATools import ata_control
import argparse
import atexit
import metadata


ap = argparse.ArgumentParser()
ap.add_argument('cmd', help="Action [start, freq,  move, end, note]", choices=['start', 'freq', 'move', 'end', 'note'])
ap.add_argument('payload', help="Argument for command.", nargs='?', default=None)
args = ap.parse_args()

ants=['1a', '1f', '5c']
defaults = argparse.Namespace(az=121.958, el=23.603, freq=1680.0)


if args.cmd == 'start':
    if args.payload is None:
        print("Please include your name or initials.")
    else:
        ata_control.move_ant_group(ants, 'none', 'atagr')
        metadata.onlog(f"session start: {args.payload}")
elif args.cmd == 'end':
    atexit.register(ata_control.move_ant_group, ants, 'atagr', 'none')
    atexit.register(ata_control.park_antennas, ants)
    metadata.onlog(f"end: {', '.join(ants)}")
elif args.cmd == 'freq':
    freq = defaults.freq if args.payload is None else float(args.payload)
    metadata.onlog(f"fcen: {freq}")
    ata_control.set_freq(freq, [ants[0]], lo='d')
    ata_control.autotune([ants[0]])
    att = 20  # Attenuation in dB
    ata_control.rf_switch_thread(ants)
    ata_control.set_atten_thread([[f'{ant}x', f'{ant}y'] for ant in ants], [[att, att] for ant in ants])
elif args.cmd == 'move':
    if args.payload is None:
        az, el = defaults.az, defaults.el
    elif ',' in args.payload:
        az, el = [float(x) for x in args.payload.split(',')]
    else:
        print(f"Invalid move argument - need az,el")
        az = None
    if az is not None:
        ata_control.set_az_el([ants[0]], az, el)
        metadata.onlog(f"az: {az}, el: {el}")
elif args.cmd == 'note':
    if args.payload is None:
        print("Need to include a note.")
    else:
        metadata.onlog(f"note: {args.payload}")
