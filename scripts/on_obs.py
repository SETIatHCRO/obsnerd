#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import argparse
ap = argparse.ArgumentParser()
from obsnerds import obs_look

if __name__ == '__main__':
    ap.add_argument('obsid', help="Name of data file (.npz or .uvh5) for single band or <obsid> for obs (then lo/cnode).")
    ap.add_argument('-a', '--ants', help='csv list of ants', default=None)
    ap.add_argument('-p', '--pols', help="List of polarizations to use for some options...", default='xx,yy,yx,xy')
    ap.add_argument('-w', '--waterfall', help="Flag to generate all of the waterfalls.", action='store_true')
    ap.add_argument('-s', '--save', help="Save the generated plots", action='store_true')
    ap.add_argument('-t', '--time_axis', help="Type of 'time' axis for dashboard ('actual_time, dt, boresight)", default='boresight')
    ap.add_argument('--tag', help="File extenstion.", choices=['npz', 'uvh5'], default='npz')
    ap.add_argument('--dir_data', help="directory for the observation data", default='./data')
    ap.add_argument('--show_obsinfo', help="Show the obsinfo in dashboard", action='store_true')
    ap.add_argument('--dash', help="Generate the dashboard", action='store_true')
    ap.add_argument('--cnode', help="If 'source' is an obsid, need cnodes to use.", default='352,544,736,928,1120,1312,1504')
    ap.add_argument('--lo', help="If 'source' is an obsid, need LO to use", default='A')
    ap.add_argument('--freq_unit', help='Frequency unit', choices=['MHz', 'GHz'], default='MHz')

    args = ap.parse_args()
    if args.ants is not None:
        args.ants = args.ants.split(',')
    args.pols = args.pols.split(',')

    look = obs_look.Look(obsid=args.obsid, lo=args.lo, cnode=args.cnode, tag=args.tag, freq_unit=args.freq_unit, dir_data=args.dir_data)
    if args.dash:
        look.dashboard(ant=args.ants[0], pol=args.pols[0], save=args.save, time_axis=args.time_axis, show_obsinfo=args.show_obsinfo)
    if not args.save:
        plt.show()
