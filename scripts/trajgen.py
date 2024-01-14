#! /usr/bin/env python
from obsnerd import trajectory_engine as te
import datetime
import argparse

ap = argparse.ArgumentParser()
ap.add_argument('ttype', help="Type of trajectory ([gal]/file)", nargs='?', default='gal')
ap.add_argument('-s', '--start_time', help="Time to start checking ['now+5' min]", default='now+5')
ap.add_argument('-u', '--using', help="Galactic latitude to use if 'gal' or filename to use [0.0 deg]", default=0.0)
ap.add_argument('-e', '--el_horizon', help='Elevation to start at [30.0 deg]', type=float, default=30.0)
ap.add_argument('-t', '--time_to_track', help='Length of track [10.0 min]', type=float, default=10.0)
ap.add_argument('-r', '--rate', help="Galactic longitude step in trajectory [0.1 deg/sec]", type=float, default=0.1)
ap.add_argument('-w', '--write_out', help="Write out the files.", action='store_true')
ap.add_argument('--tstep', help="Time step in trajectory [1 sec]", type=float, default=1.0)
ap.add_argument('--tz', help="Hours from UTC [0.0 hours]", type=float, default=0.0)
args = ap.parse_args()

if args.start_time.startswith('now'):
    offsp = args.start_time.split('+')
    offsm = args.start_time.split('-')
    if len(offsp) == 2:
        args.start_time = datetime.datetime.now() + datetime.timedelta(minutes=float(offsp[1]))
    elif len(offsm) == 2:
        args.start_time = datetime.datetime.now() - datetime.timedelta(minutes=float(offsm[1]))

traj = te.Trajectory(start_time=args.start_time, time_to_track=args.time_to_track, tstep=args.tstep, el_horizon=args.el_horizon, tz=args.tz)

if args.ttype == 'gal':
    traj.galactic(b=float(args.using), rate=args.rate)
else:
    traj.from_file(filename=args.using)
if args.write_out:
    traj.write()
if traj.valid:
    print(traj.track.name)
    traj.plot()