#! /usr/bin/env python
from os import listdir, path
import argparse
import json

ap = argparse.ArgumentParser()
ap.add_argument('mjd', help='csv list of mjd to include.')
ap.add_argument('--data_dir', help="Directory where the data live", default='data')
args = ap.parse_args()
args.mjd = args.mjd.split(',')
obsinfofn = f"obsinfo_{args.mjd[0].split('.')[0]}.json"

obsinfo = {
  "dir_data": args.data_dir,
  "Filters": {
    "A": {
      "r": [1690.0, 1710.0],
      "b": [1939.0, 1945.0],
      "g": [1990.0, 1995.0],
      "c": [1960.0, 1963.0],
      "m": [1975.0, 1985.0],
      "orange": [1155.0, 1165.0],
      "maroon": [1055.0, 1070.0]
    },
    "B": {
      "r": [3740.0, 3760.0],
      "b": [3975.0, 3990.0],
      "g": [5550.0, 5560.0]
    }
  },
  "Sources": {
  }
}

# first update obsmeta
with open('obsmeta.json', 'r') as fp:
  obsmeta = json.load(fp)
obsmeta.update({obsinfofn: set()})
found_files = {}
for xfn in listdir(args.data_dir):
    for mjd in args.mjd:
        if xfn.endswith('.npz') and mjd in xfn:
            # print(f"Found {xfn}")
            obsid = '_'.join(xfn.split('_')[:2])
            obsmeta[obsinfofn].add(obsid)
            found_files[path.join(args.data_dir, xfn)] = obsid
obsmeta[obsinfofn] = sorted(obsmeta[obsinfofn])
print("NEW OBSMETA")
for key, val in obsmeta.items():
    print(key)
    if isinstance(val, list):
        print(', '.join(val))
    else:
        print(val)
with open('obsmeta.json', 'w') as fp:
    json.dump(obsmeta, fp, indent=2)

# Now update obsinfo
for obsid in obsmeta[obsinfofn]:
    obsinfo["Sources"][obsid] = {"offset": 0.0}

with open(obsinfofn, 'w') as fp:
    json.dump(obsinfo, fp, indent=2)
