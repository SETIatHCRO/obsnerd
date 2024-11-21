#! /usr/bin/env python
from obsnerds import starlink_look, starlink_eph
from os import listdir
import argparse

ap = argparse.ArgumentParser()
ap.add_argument('mjd', help='csv list of mjd to include.')
args = ap.parse_args()
args.mjd = args.mjd.split(',')
fephfn = f"feph_{args.mjd[0].split('.')[0]}.json"

feph = {
    "Filters": {
      "A": {
        "r": [
          1690.0,
          1710.0
        ],
        "b": [
          1939.0,
          1945.0
        ],
        "g": [
          1990.0,
          1995.0
        ],
        "c": [
          1960.0,
          1963.0
        ]
      },
      "B": {
        "r": [
          3740.0,
          3760.0
        ],
        "b": [
          3975.0,
          3990.0
        ]
      }
    },
  "Sources":
  {
  }
}

look  = starlink_look.Look()

for xfn in listdir('.'):
    for mjd in args.mjd:
        if xfn.endswith('.npz') and mjd in xfn:
            print(f"Found {xfn}")
            look.read_obsrec(xfn)
            obsid = '_'.join(xfn.split('_')[:2])
            feph['Sources'][obsid] = {'tref': look.times[0].datetime.isoformat()}
            break
for x in feph['Sources']:
    print(x)

look.put_feph(fn=fephfn, fdict=feph)