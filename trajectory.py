#! /usr/bin/env python

import astropy.units as u
from astropy.coordinates import AltAz, EarthLocation, SkyCoord
from astropy.time import Time
import numpy as np
import matplotlib.pyplot as plt
import datetime, time

EPHEM_FILENAME = 'track.ephem'
TRAJECTORY_FILENAME = 'track'
TRACK_LOG_FILENAME = 'track.log'
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

    def __repr__(self):
        s = ''
        for i, b in zip([0, -1], ['Start at:', 'End at:']):
            s += f"{b} {self.obstime[i]} UTC\n"
            s += f"\tl={self.l[i]}, b={self.b[i]}\n"
            s += f"\tRA={self.ra[i]}, Dec={self.dec[i]}\n"
            s += f"\tAz={self.az[i]}, El={self.el[i]}\n"
        return s

    def add(self, **kwargs):
        for key, value in kwargs.items():
            getattr(self, key).append(value)

    def log_track(self, filename=TRACK_LOG_FILENAME):
        print(self)
        if filename is not None:
            with open(filename, 'w') as fp:
                print(self, file=fp)
            print(f"Writing log: {filename}")

    def ephem_file(self, filename=EPHEM_FILENAME):
        print(f"Writing ephemerides:  {filename}")
        tai = np.array(self.timestamp, dtype=int)
        az = np.array(self.az, dtype=float)
        el = np.array(self.el, dtype=float)
        ir = np.zeros(len(tai)) + 1E-10
        ephem = ((np.array([tai, az, el, ir], dtype=object)))
        self.ephemtxt = np.savetxt(filename, ephem.T, fmt='%i  %.5f  %.5f  %.10E')

    def full_file(self, filename=TRAJECTORY_FILENAME):
        print(F"Writing full trajectory: {filename}")
        np.savez_compressed(filename, obstime=self.obstime, az=self.az, el=self.el, ra=self.ra, dec=self.dec, l=self.l, b=self.b, allow_pickle=True)

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

    track_times = start_time + np.arange(0.0, time_to_track * 60.0, tstep) * u.second
    lstep = lstep*u.deg  #deg/sec

    ts = time.mktime(start_time.datetime.timetuple()) + 37.0
    this_track = Track()
    dtns = int(tstep * 1E9)
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

    this_track.log_track(TRACK_LOG_FILENAME)
    this_track.ephem_file(EPHEM_FILENAME)
    this_track.full_file(TRAJECTORY_FILENAME)

    plt.plot(this_track.az, this_track.el, '.', lw=4, label='azel-track')
    plt.grid()
    plt.legend()
    plt.xlabel('RA/Az')
    plt.ylabel('Dec/El')

    return this_track


if __name__ == '__main__':
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument('-s', '--start_time', help="Time to start checking ['now+5' min]", default='now+5')
    ap.add_argument('-b', '--b2use', help='Galactic latitude to use [0.0 deg]', type=float, default=0.0)
    ap.add_argument('-e', '--el_starting', help='Elevation to start at [30.0 deg]', type=float, default=30.0)
    ap.add_argument('-t', '--time_to_track', help='Length of track [10.0 min]', type=float, default=10.0)
    ap.add_argument('-l', '--lstep', help="Galactic longitude step in trajectory [0.1 deg/sec]", type=float, default=0.1)
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

    traj = main(**vars(args))
    plt.subplot(212)
    plt.plot(traj.obstime, traj.az, '.', label='Az')
    plt.plot(traj.obstime, traj.el, '.', label='El')
    plt.plot(traj.obstime, traj.ra, '.', label='RA')
    plt.plot(traj.obstime, traj.dec, '.', label='Dec')
    plt.grid()
    plt.legend()
    plt.show()