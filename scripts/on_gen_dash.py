#! /usr/bin/env python
import argparse
from obsnerd import obs_look
from obsnerd.on_sys import AXIS_OPTIONS


ap = argparse.ArgumentParser()
ap.add_argument('obsinfo', help="Name of obsinfo file to use")
ap.add_argument('--script', help="Name of script file", default='dash.sh')
ap.add_argument('--lo', help="LO to use", choices=['A', 'B'], default='A')
ap.add_argument('--cnode', help="Correlator nodes to use", default='all')
ap.add_argument('-a', '--ants', help='antennas to generate for', default='2b')
ap.add_argument('-p', '--pols', help="polarization to generate for", default='xx')
ap.add_argument('-t', '--time_axis', help="Type of 'time' axis for dashboard ([d]atetime, [s]econds, [b]oresight)", default='boresight')
ap.add_argument('--show_diff', help="Flag to include show_diff", action='store_true')


args = ap.parse_args()

args.time_axis = AXIS_OPTIONS[args.time_axis[0].lower()]

obs = obs_look.Look(obsinput=None, lo=args.lo, cnode=args.cnode)
obs.dashboard_gen(obsinfo=args.obsinfo, script_fn=args.script, ants=args.ants, pols=args.pols, taxis=args.time_axis, show_diff=args.show_diff)