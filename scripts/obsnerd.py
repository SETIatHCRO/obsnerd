#!/usr/bin/env python
from obsnerds import obsnerd_engine as oe
import argparse

ap = argparse.ArgumentParser()
ap.add_argument('cmd', help="Action [start, freq,  move, end, note, source, summary]", choices=['start', 'freq', 'move', 'end', 'note', 'source', 'summary'])
ap.add_argument('arg1', help="Argument for command.", nargs='?', default=None)
ap.add_argument('arg2', help="If source, associated datetime (secondary argument).", nargs='?', default=None)
ap.add_argument('-c', '--coord-type', dest='coord_type', help="For <move> coordinate type: azel/radec/source/traj [azel]",
                choices=['azel', 'radec', 'source', 'traj'], default='azel')
args = ap.parse_args()

session = oe.CommandHandler(arg1=args.arg1, arg2=args.arg2, coord_type=args.coord_type)
getattr(session, args.cmd)()