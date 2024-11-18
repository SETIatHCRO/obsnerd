#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import argparse
ap = argparse.ArgumentParser()
from obsnerds import starlink_look, starlink_wide

if __name__ == '__main__':
    ap.add_argument('source', help="Name of data file (.npz or .uvh5) for single band or <source> for wide (then lo/cnode).")
    ap.add_argument('-a', '--ants', help='csv list of ants', default=None)
    ap.add_argument('-p', '--pols', help="List of polarizations to use for some options...", default='xx,yy,yx,xy')
    ap.add_argument('-w', '--waterfall', help="Flag to generate all of the waterfalls.", action='store_true')
    ap.add_argument('-s', '--save', help="Save the generated plots", action='store_true')
    ap.add_argument('-t', '--time_axis', help="Type of 'time' axis for dashboard", default='boresight')
    ap.add_argument('--feph', help="Name of feph file (default will read from feph_files.json)", default='auto')
    ap.add_argument('--show_feph', help="Show the feph info in dashboard", action='store_true')
    ap.add_argument('--dash', help="Generate the dashboard", action='store_true')
    ap.add_argument('--cnode', help="If 'wide' cnodes to use.", default='352,544,736,928,1120,1312,1504')
    ap.add_argument('--lo', help="If 'wide' LO to use", default='A')

    args = ap.parse_args()
    if args.ants is not None:
        args.ants = args.ants.split(',')
    args.pols = args.pols.split(',')
    args.cnode = [int(x) for x in args.cnode.split(',')]
    source_type = args.source.split('.')[-1]

    if source_type in ['uvh5', 'npz']:  # Single file as source
        look = starlink_look.Look()
        look.read_source(args.source)

        if args.waterfall:
            for pol in args.pols:
                look.all_wf(pol=pol, save=args.save)
            if not args.save:
                plt.show()
        if args.dash:
            look.dashboard(ant=args.ants[0], pol=args.pols[0], save=args.save, time_axis=args.time_axis, feph=args.feph, show_feph=args.show_feph)
            if not args.save:
                plt.show()
    else:
        wide = starlink_wide.WideBand(source=args.source, lo=args.lo, cnode=args.cnode)
        wide.concat(ant=args.ants[0], pol=args.pols[0])
        if args.dash:
            wide.dashboard(save=args.save, time_axis=args.time_axis, show_feph=args.show_feph)
            if not args.save:
                plt.show()