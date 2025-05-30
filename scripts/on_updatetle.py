#! /usr/bin/env python
from obsnerd import updatetle_engine as ue
import argparse

ap = argparse.ArgumentParser()
ap.add_argument('group', help="Group to update TLEs for (default: *)", nargs='?', default='*')
ap.add_argument('--base-url', dest='base_url', help="Base url for tles",
                default='https://celestrak.org/NORAD/elements/')
ap.add_argument('--base-path', dest='base_path', help="Base path for tles",
                default='./tle')
args = ap.parse_args()
ue.updatetle(group=args.group, base_path=args.base_path, base_url=args.base_url)
ue.update_log()