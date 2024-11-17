#! /usr/bin/env python
import argparse
from obsnerds import sopp_engine as se


ap = argparse.ArgumentParser()
ap.add_argument('-t', '--time', help="Start time as isoformat string or delay in minutes from now.", default=0.0)
ap.add_argument('-d', '--duration', help='Duration in minutes [20]', type=float, default=20.0)
ap.add_argument('-f', '--frequency', help='Frequency to search in MHz', type=float, default=1575.0)
ap.add_argument('-b', '--bandwidth', help="Bandwidth in MHz [20]", type=float, default=20.0)
ap.add_argument('-s', '--search', help="String to search for", default=False)
ap.add_argument('-x', '--exclude', help="If str, exclude if string is present", default=False)
ap.add_argument('-o', '--orbit', help='Orbit type: (all, geo, meo, leo) [all]', choices=['all', 'geo', 'meo', 'leo'], default='all')
ap.add_argument('-e', '--el_limit', help="Lower horizon elevation [20.0]", type=float, default=20.0)
ap.add_argument('-a', '--az_limit', help="Azimuth range [0,360]", default='0,360')
ap.add_argument('--time_resolution', help="Time resolution of track in seconds", type=int, default=1)
ap.add_argument('--tz', help='Time zone (hours offset from UTC) [0]', type=float, default=None)
ap.add_argument('--tle_file', help='Name of tle file', default='tle/active.tle')
ap.add_argument('--ftype', help='search horizon or beam', choices=['horizon', 'beam'], default='horizon')
ap.add_argument('--output_file', help="Flag to write out ephemerides", action='store_true')
ap.add_argument('--sat2write', help="String in satellite name to write out.", default=None)
args = ap.parse_args()

az_limit = [float(x) for x in args.az_limit.split(',')]

if __name__ == '__main__':
    se.main(start=args.time,
            duration = args.duration,
            frequency = args.frequency,
            bandwidth = args.bandwidth,
            az_limit = az_limit,
            el_limit = args.el_limit,
            ftype = args.ftype,
            search_for = args.search,
            orbit_type = args.orbit,
            exclude = args.exclude,
            time_resolution=args.time_resolution,
            tle_file = args.tle_file,
            timezone = args.tz,
            output_file = args.output_file,
            sat2write=args.sat2write
            )
