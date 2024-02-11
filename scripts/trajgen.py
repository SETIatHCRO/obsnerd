#! /usr/bin/env python
from obsnerds import trajectory_engine as te
import argparse

ap = argparse.ArgumentParser()
ap.add_argument('type', help="Type of trajectory ([gal]/file)", nargs='?', default='gal')
ap.add_argument('-t', '--time', help="Start time as isoformat string or delay in minutes from now. 'None' to use file start.", default=0.0)
ap.add_argument('-u', '--using', help="Galactic latitude to use if 'gal' or filename to us for ephem [0.0 deg]", default=0.0)
ap.add_argument('-e', '--el_horizon', help='Elevation to start at [30.0 deg]', type=float, default=30.0)
ap.add_argument('-d', '--duration', help='Duration of track [10.0 min]', type=float, default=10.0)
ap.add_argument('-r', '--rate', help="Galactic longitude step in trajectory [0.1 deg/sec]", type=float, default=0.1)
ap.add_argument('-w', '--write_out', help="Write out the files.", action='store_true')
ap.add_argument('--tstep', help="Time step in trajectory [1 sec]", type=float, default=1.0)
ap.add_argument('--tz', help="Hours from UTC [0.0 hours]", type=float, default=0.0)
args = ap.parse_args()

traj = te.Trajectory(start_time=args.time, duration=args.duration, tstep=args.tstep, el_horizon=args.el_horizon, tz=args.tz)

if args.type == 'gal':
    traj.galactic(b=float(args.using), rate=args.rate)
else:
    traj.from_file(filename=args.using)
if traj.valid:
    print(traj.track.name)
    traj.plot()
    if args.write_out:
        traj.write()