#! /usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
from obsnerd import spacex_api
import json

PROC_FILENAME = 'onchose.json'

ap = argparse.ArgumentParser()
ap.add_argument('-s', '--show_sats', help="Flag to just show DTC list of satellites", action='store_true')
ap.add_argument('-t', '--show_tracks', help="Flag to show satellite tracks", action='store_true')
ap.add_argument('-n', '--satname', help="Name (or partial) of satellite to use", default=None)
ap.add_argument('-p', '--summary_par', help="Parameter to plot", choices=['az', 'el'], default='az')
ap.add_argument('-c', '--choose_tracks', help="csv list of tracks/locations to use", default=None)
ap.add_argument('-x', '--proc_tracks', help="Flag to process tracks", action='store_true')
ap.add_argument('-o', '--obslin_min', help="Observation length in minutes", type=int, default=8)
ap.add_argument('-f', '--freqs', help="List of frequencies to use", default='1990.0,5500.0')
ap.add_argument('--freq_unit', help="Frequency unit", default='MHz')
args = ap.parse_args()
args.freqs = [float(f) for f in args.freqs.split(',')]

sx = spacex_api.SpaceX()
if args.show_sats:
    from tabulate import tabulate
    table = []
    tablerow = []
    for i, sat in enumerate(sx.dtc):
        tablerow.append(sat)
        if i % 10 == 9:
            table.append(tablerow)
            tablerow = []
    print(tabulate(table))
    print("\nPipeline is:\n\ton_chose_spacex.py -s\n\ton_chose_spacex.py -t -n <satname>")
    print("\ton_chose_spacex.py -c <tracks>\n\ton_chose_spacex.py -x\n")
    print("Also, summary_par (az/el) obslin_min and freqs are options")
else:
    try:
        with open(PROC_FILENAME, 'r') as fp:
            inp = json.load(fp)
    except FileNotFoundError:
        inp = {}
    if args.satname is None:
        args.satname = inp['satname']
    else:
        with open(PROC_FILENAME, 'w') as fp:
            json.dump({"satname": args.satname}, fp)
    sx.get_satellite(args.satname)
    sx.get_tracks()
    sx.plot_track_summary(args.summary_par)
    if args.choose_tracks is not None:
        sx.choose_tracks(args.choose_tracks)
        with open(PROC_FILENAME, 'w') as fp:
            json.dump({"satname": args.satname, "tracks": args.choose_tracks}, fp)
    elif 'tracks' in inp:
        sx.choose_tracks(inp['tracks'])
    spacex_api.plt.show()
    if args.proc_tracks:
        sx.proc_tracks(args.obslin_min, args.freqs, args.freq_unit)
        with open(PROC_FILENAME, 'w') as fp:
            json.dump({}, fp)