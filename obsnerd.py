#!/usr/bin/python3
from ATATools import ata_control
import argparse
import atexit
import metadata


ap = argparse.ArgumentParser()
ap.add_argument('cmd', help="Action [start, freq,  move, end, note]", choices=['start', 'freq', 'move', 'end', 'note'])
ap.add_argument('payload', help="Argument for command.", nargs='?', default=None)
ap.add_argument('-c', '--coord-type', dest='coord_type', help="Coordinate type: azel/radec,source [azel]",
                choices=['azel', 'radec', 'source'], default='azel')
args = ap.parse_args()

ants=['1a', '1f', '5c']
defaults = argparse.Namespace(az=121.958, el=23.603, freq=1680.0, source_name='goes16',
                              x=121.958, y=23.603,
                              ra=23.39077787, dec=58.8077786056, note='radec is for casa')


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
        x, y = defaults.x, defaults.y
    elif ',' in args.payload:
        x, y = [float(_v) for _v in args.payload.split(',')]
    else:
        x = args.payload

    if args.coord_type == 'azel':
        ata_control.set_az_el([ants[0]], x, y)
        metadata.onlog(f"az: {x}, el: {y}")
    elif args.coord_type == 'radec':
        source = ata_control.track_source([ants[0]], radec=[x, y])
        metadata.onlog(f"ra: {x}, dec: {y}")
    elif args.coord_type == 'source':
        source = ata_control.track_source([ants[0]], source=x)
        metadata.onlog(f"source {x}")
elif args.cmd == 'note':
    if args.payload is None:
        print("Need to include a note.")
    else:
        metadata.onlog(f"note: {args.payload}")
