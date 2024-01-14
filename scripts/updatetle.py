#! /usr/bin/env python
from obsnerd import updatetle_engine
import argparse


ap = argparse.ArgumentParser()
ap.add_argument('--base-url', dest='base_url', help="Base url for tles",
                default='http://celestrak.org/NORAD/elements/')
ap.add_argument('--base-path', dest='base_path', help="Base path for tles",
                default='./tle')
args = ap.parse_args()
updatetle_engine.updatetle(base_path=args.base_path, base_url=args.base_url)
updatetle_engine.update_log()