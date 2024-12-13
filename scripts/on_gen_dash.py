#! /usr/bin/env python
import argparse
from obsnerds import obs_look


ap = argparse.ArgumentParser()
ap.add_argument('obsinfo', help="Name of obsinfo file to use")
ap.add_argument('--script', help="Name of script file", default='dash.sh')
ap.add_argument('--lo', help="LO to use", choices=['A', 'B'], default='A')
ap.add_argument('--cnode', help="Correlator nodes to use", default='all')
ap.add_argument('--ants', help="Ants to use - csv list of 'all'", default='2a,2b,2h,4e')
ap.add_argument('--pols', help="Pols to use", default='xx,xy')
ap.add_argument('--axis', help="Type of axis to use", choices=['b', 'a', 'd'], default='b')
args = ap.parse_args()

obs = obs_look.Look(obsid=None, lo=args.lo, cnode=args.cnode)
obs.dashboard_gen(obsinfo=args.obsinfo, ant=args.ants, pol=args.pols, taxis=args.axis)