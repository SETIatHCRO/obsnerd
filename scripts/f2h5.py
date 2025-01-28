#! /usr/bin/env python
from obsnerd import convert2hdf5
import argparse

"""
Convert the gnuradio-companion spectrum file to hdf5
"""

RAW_FILENAME = 'nrdz'

ap = argparse.ArgumentParser()
ap.add_argument('tag', help="Name of output file or prefix for date [meta]", nargs='?', default='meta')
ap.add_argument('-d', '--date', help="Optional date to generate the filename with tag as prefix", default=None)
ap.add_argument('-z', '--timezone', help="Timezone to add (-8 is PST) [0].", type=float, default=0.0)
args = ap.parse_args()

convert2hdf5.convert(RAW_FILENAME, convert2hdf5.make_filename(**vars(args)))
