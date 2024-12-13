#! /usr/bin/env python
from obsnerds import obs_look
from os import listdir
import argparse

ap = argparse.ArgumentParser()
ap.add_argument('mjd', help='csv list of mjd to include.')
args = ap.parse_args()
args.mjd = args.mjd.split(',')
obsinfofn = f"obsinfo_{args.mjd[0].split('.')[0]}.json"

obsinfo = {
    "dir_data": ".",
    "Filters": {
      "A": {
        "r": [1690.0, 1710.0],
        "b": [1939.0, 1945.0],
        "g": [1990.0, 1995.0],
        "c": [1960.0, 1963.0]
      },
      "B": {
        "r": [3740.0, 3760.0],
        "b": [3975.0, 3990.0],
        "g": [5550.0, 5560.0],
        "m": [1975.0, 1985.0]
      }
    },
  "Sources":
  {
  }
}

look = {}

for xfn in listdir('.'):
    for mjd in args.mjd:
        if xfn.endswith('.npz') and mjd in xfn:
            print(f"Found {xfn}")
            tmp = xfn.split('_')
            obsid = '_'.join(xfn.split('_')[:2])
            look[xfn] = obs_look.Look(obsid)
            look[xfn].read_obsrec(xfn)
            obsinfo['Sources'][obsid] = {'tref': look.times[0].datetime.isoformat()}
            break
for x in obsinfo['Sources']:
    print(x)

look.put_obsinfo(fn=obsinfofn, fdict=obsinfo)