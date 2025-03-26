#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import argparse
ap = argparse.ArgumentParser()
from obsnerd import onv_look
from obsnerd.on_sys import AXIS_OPTIONS


if __name__ == '__main__':
    ap.add_argument('source', help="Name of data file (.npz or .uvh5) for single band or source/obsid.")
    ap.add_argument('-a', '--ant', help='antenna to use', default='2a')
    ap.add_argument('-p', '--pol', help="polarization to use", default='xx')
    ap.add_argument('-w', '--waterfall', help="Flag to generate all of the waterfalls.", action='store_true')
    ap.add_argument('-s', '--save', help="Save the generated plots", action='store_true')
    ap.add_argument('-t', '--time_axis', help="Type of 'time' axis for dashboard ([d]atetime, [s]econds, [b]oresight)", default='boresight')
    ap.add_argument('--dash', help="Generate the dashboard", action='store_true')
    ap.add_argument('--log', help="Show the SpaceX log info (provide filename)", default=False)
    ap.add_argument('--show_diff', help="Show the difference from the dashboard", action='store_true')
    ap.add_argument('--cnode', help="If 'source' is an obsid, need cnodes to use.", default='352,544,736,928,1120,1312,1504')
    ap.add_argument('--lo', help="If 'source' is an obsid, need LO to use", default='A')
    ap.add_argument('--freq_unit', help='Frequency unit', choices=['MHz'], default='MHz')

    args = ap.parse_args()
    args.time_axis = AXIS_OPTIONS[args.time_axis[0].lower()]

    look = onv_look.Look(oinput=args.source, lo=args.lo, cnode=args.cnode, freq_unit=args.freq_unit)
    if args.dash:
        for ant in args.ant.split(','):
            for pol in args.pol.split(','):
                look.dashboard(ant=ant, pol=pol, save=args.save, time_axis=args.time_axis, show_diff=args.show_diff, log=args.log)
    if not args.save:
        plt.show()
