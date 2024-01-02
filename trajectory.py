#! /usr/bin/env python

import astropy.units as u
from astropy.coordinates import AltAz, EarthLocation, SkyCoord
from astropy.time import Time
import numpy as np
import matplotlib.pyplot as plt
import datetime, time


hcro = EarthLocation(lat=40.8178049*u.deg, lon=-121.4695413*u.deg, height=986*u.m)


class Track:
    def __init__(self, timestamp=[], obstime=[], az=[], el=[], ra=[], dec=[], l=[], b=[]):
        self.timestamp = timestamp
        self.obstime = obstime
        self.az = az
        self.el = el
        self.ra = ra
        self.dec = dec
        self.l = l
        self.b = b

    def add(self, **kwargs):
        for key, value in kwargs.items():
            getattr(self, key).append(value)



def main(start_time='2023-12-31 23:59:59', b2use=0.0, el_starting=30.0, time_to_track=20.0, lstep=0.1, tstep=1.0, tz=-8.0):
    """
    Parameters
    ----------
    start_time : str etc
        isoformat or anything Time handles in local time used to start checking
    b2use : float
        galactic b value to use [deg]
    el_starting : float
        Elevation to start using track [deg]
    time_to_track : float
        Length of track in minutes
    lstep : float
        galactic l step for track [deg/sec]
    tstep : float
        Trajectory time step in seconds
    tz : float
        Hours from UTC

    """
    filename = 'track.trk'
    tz = tz * u.hour
    start_time = Time(start_time) - tz

    b2use = b2use * u.deg
    gp_l = np.arange(0, 359) * u.deg
    gp_b = np.ones(len(gp_l)) * b2use


    gal = SkyCoord(frame='galactic', l=gp_l, b=gp_b)
    radec = gal.transform_to('icrs')
    azel = radec.transform_to(AltAz(obstime=start_time, location=hcro))
    above_horizon = np.where(azel.alt.value > el_starting)
    start_horizon = above_horizon[0][0]
    starting_l = gal.l.value[start_horizon]

    plt.figure('At T=0')
    plt.subplot(211)
    plt.plot(radec.ra.value[above_horizon], radec.dec.value[above_horizon], '.', label='radec-init')
    plt.plot(azel.az.value[above_horizon], azel.alt.value[above_horizon], '.', label='azel-init')

    print(f"UTC {start_time}")
    print(f"Start at l = {starting_l}, b={gal.b.value[start_horizon]}") 
    print(f"         RA={radec.ra.value[start_horizon]}, Dec={radec.dec.value[start_horizon]}")
    print(f"         Az={azel.az.value[start_horizon]}, El={azel.alt.value[start_horizon]}")

    track_times = start_time + np.arange(0.0, time_to_track * 60.0, tstep) * u.second
    lstep = lstep*u.deg  #deg/sec

    print(f"Writing {filename}")
    ts = time.mktime(start_time.datetime.timetuple())
    this_track = Track()
    dtns = int(tstep * 1E9)
    with open(filename, 'w') as fp:
        for i in range(len(track_times)):
            this_l = starting_l*u.deg  + i*lstep
            if this_l > 360.0*u.deg:
                this_l = this_l - 360.0*u.deg
            gal = SkyCoord(frame='galactic', l=this_l, b=b2use)
            radec = gal.transform_to('icrs')
            azel = radec.transform_to(AltAz(obstime=track_times[i], location=hcro))
            this_ts = int(ts*1E9) + i*dtns
            this_track.add(timestamp=this_ts, obstime=track_times[i].datetime, az=azel.az.value, el=azel.alt.value,
                           ra=radec.ra.value, dec=radec.dec.value, l=this_l.value, b=b2use.value)
            print(f"{this_ts},{azel.az.value},{azel.alt.value}", file=fp)
    plt.plot(this_track.az, this_track.el, '.', lw=4, label='azel-track')
    plt.grid()
    plt.legend()
    plt.xlabel('RA/Az')
    plt.ylabel('Dec/El')

    return this_track


if __name__ == '__main__':
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument('-s', '--start_time', help="Time to start checking ['now+10' min]", default='now+10')
    ap.add_argument('-b', '--b2use', help='Galactic latitude to use [0.0 deg]', type=float, default=0.0)
    ap.add_argument('-e', '--el_starting', help='Elevation to start at [30.0 deg]', type=float, default=30.0)
    ap.add_argument('-t', '--time_to_track', help='Length of track [20.0 min]', type=float, default=20.0)
    ap.add_argument('-l', '--lstep', help="Galactic longitude step in trajectory [0.1 deg/sec]", type=float, default=0.1)
    ap.add_argument('--tstep', help="Time step in trajectory [1 sec]", type=float, default=1.0)
    ap.add_argument('--tz', help="Hours from UTC [-8.0 hours]", type=float, default=-8.0)
    args = ap.parse_args()

    if args.start_time.startswith('now'):
        offsp = args.start_time.split('+')
        offsm = args.start_time.split('-')
        if len(offsp) == 2:
            args.start_time = datetime.datetime.now() + datetime.timedelta(minutes=float(offsp[1]))
        elif len(offsm) == 2:
            args.start_time = datetime.datetime.now() - datetime.timedelta(minutes=float(offsm[1]))

    traj = main(**vars(args))
    plt.subplot(212)
    plt.plot(traj.obstime, traj.az, '.', label='Az')
    plt.plot(traj.obstime, traj.el, '.', label='El')
    plt.plot(traj.obstime, traj.ra, '.', label='RA')
    plt.plot(traj.obstime, traj.dec, '.', label='Dec')
    plt.grid()
    plt.legend()
    plt.show()