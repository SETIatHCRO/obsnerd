#! /usr/bin/env python
import argparse
from obsnerds import obs_look


ap = argparse.ArgumentParser()
ap.add_argument('obsinfo', help="Name of obsinfo file to use")
ap.add_argument('--script', help="Name of script file", default='dash.sh')
ap.add_argument('--lo', help="LO to use", choices=['A', 'B'], default='A')
ap.add_argument('--cnode', help="Correlator nodes to use", default='all')
ap.add_argument('--ants', help="Ants to use - csv list", default='2b')
ap.add_argument('--pols', help="Pols to use - csv list", default='xx')
ap.add_argument('--axis', help="Type of axis to use <[b]oresight/[a]ctual/[d]ifference>", choices=['b', 'a', 'd'], default='b')
args = ap.parse_args()

obs = obs_look.Look(obsinput=None, lo=args.lo, cnode=args.cnode)
obs.dashboard_gen(obsinput=args.obsinfo, script_fn=args.script, ants=args.ants, pols=args.pols, taxis=args.axis)