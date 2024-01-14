#! /usr/bin/env python
import argparse
from datetime import datetime, timedelta
from obsnerds import sopp_engine as se


ap = argparse.ArgumentParser()
ap.add_argument('-t', '--start_time', help="Start time as isoformat string, default is now", default=None)
ap.add_argument('-d', '--duration', help='Duration in minutes [20]', type=float, default=20.0)
ap.add_argument('-f', '--frequency', help='Frequency to search in MHz', type=float, default=1575.0)
ap.add_argument('-b', '--bandwidth', help="Bandwidth in MHz [20]", type=float, default=20.0)
ap.add_argument('-s', '--search', help="String to search for", default=False)
ap.add_argument('-o', '--orbit', help='Orbit type: (all, geo, meo, leo) [all]', choices=['all', 'geo', 'meo', 'leo'], default='all')
ap.add_argument('-e', '--el_limit', help="Lower horizon elevation [20.0]", type=float, default=20.0)
ap.add_argument('-a', '--az_limit', help="Azimuth range [0,360]", default='0,360')
ap.add_argument('--offset', help="Number of minutes offset to get positions [10.0]", type=float, default=10.0)
ap.add_argument('--tz', help='Time zone (hours offset from UTC) [0]', type=float, default=0.0)
ap.add_argument('--tle_file', help='Name of tle file', default='tle/active.tle')
ap.add_argument('--ftype', help='search horizon or beam', choices=['horizon', 'beam'], default='horizon')
ap.add_argument('--output_file', help="Name of output file (None if not supplied)", default=None)
args = ap.parse_args()
if args.start_time is None:
    args.start_time = datetime.now()
else:
    args.start_time = datetime.strptime(args.start_time, '%Y-%m-%dT%H:%M') + timedelta(microseconds=1)
args.start_time -= timedelta(hours=args.tz)  # Convert _to_ UTC
stop_time = args.start_time + timedelta(minutes=args.duration)
offset_time = args.start_time + timedelta(minutes=args.offset)
az_limit = [float(x) for x in args.az_limit.split(',')]

se.main(starttime=args.start_time,
        stoptime = stop_time,
        offsettime = offset_time,
        frequency = args.frequency,
        bandwidth = args.bandwidth,
        az_limit = az_limit,
        el_limit = args.el_limit,
        ftype = args.ftype,
        search_for = args.search,
        orbit_type = args.orbit,
        tle_file = args.tle_file,
        timezone = args.tz,
        output_file = args.output_file)
