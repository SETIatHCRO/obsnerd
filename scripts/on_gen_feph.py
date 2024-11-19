#! /usr/bin/env python
from obsnerds import starlink_look, starlink_eph
import argparse

ap = argparse.ArgumentParser()
ap.add_argument('sources', help='Filename containing sources to include.')
args = ap.parse_args()

feph = {'Sources': {}}

look  = starlink_look.Look()
sl = starlink_eph.Eph()

with open(args.sources, 'r') as fp:
    for line in fp:
        fn = line.strip()
        look.read_source(fn)
        srcname = fn.split('_')[0]
        feph['Sources'][srcname] = {'tref': look.times[0].datetime.isoformat()}

sl.write_feph(feph_dict=feph)