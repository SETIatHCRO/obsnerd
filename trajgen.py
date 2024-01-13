#! /usr/bin/env python

import astropy.units as u
from astropy.coordinates import AltAz, EarthLocation, SkyCoord
from astropy.time import Time
import numpy as np
import matplotlib.pyplot as plt
import datetime

EPHEM_FILENAME = 'track.ephem'
TRAJECTORY_FILENAME = 'track'
TRACK_LOG_FILENAME = 'track.log'
hcro = EarthLocation(lat=40.8178049*u.deg, lon=-121.4695413*u.deg, height=986*u.m)

class Track:
    """
    Parameters that are tracked:
    timestamp : int (ns timestamps)
    obstime : Time
    az/el/ra/dec/l/b : float
    ir : float (inverse radius)
    """
    def __init__(self, name, **kwargs):
        self.parameters = ['timestamp', 'obstime', 'az', 'el', 'ra', 'dec', 'l', 'b', 'ir']
        self.name = name
        for param in self.parameters:
            setattr(self, param, [])
        for key, val in kwargs.items():
            if key not in self.parameters:
                print(f"{key} not allowed parameter")
                continue
            setattr(self, key, val)

    def __repr__(self):
        s = ''
        for i, lbl in zip([0, -1], ['Start at:', 'End at:']):
            s += f"{lbl} {self.obstime[i]} UTC\n"
            s += f"\tl={self.l[i]:.3f}, b={self.b[i]:.3f}\n"
            s += f"\tRA={self.ra[i]:.3f}, Dec={self.dec[i]:.3f}\n"
            s += f"\tAz={self.az[i]:.3f}, El={self.el[i]:.3f}\n"
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
        if not len(self.ir):
            self.ir = np.zeros(len(tai)) + 1E-20
        ephem = ((np.array([tai, az, el, self.ir], dtype=object)))
        self.ephemtxt = np.savetxt(filename, ephem.T, fmt='%i  %.5f  %.5f  %.10E')

    def full_file(self, filename=TRAJECTORY_FILENAME):
        print(F"Writing full trajectory: {filename}")
        np.savez_compressed(filename, obstime=self.obstime.datetime,
                            az=self.az, el=self.el, ra=self.ra, dec=self.dec, 
                            l=self.l, b=self.b, ir=self.ir, allow_pickle=True)

def _get_track_times(start, duration, step):
    """
    Parameters
    ----------
    start : Time
    duration : float
        duration in seconds
    step : float
        step in seconds

    """
    track_Times = start + np.arange(0.0, duration, step) * u.second
    track_tsns = int(start.unix_tai * 1E9) + np.arange(0, int(duration * 1E9), int(step * 1E9))
    return track_Times, track_tsns

def galactic(ttype='gal', start_time='2023-12-31 23:59:59', using=0.0, el_starting=30.0, time_to_track=20.0, anglerate=0.1, tstep=1.0, tz=-8.0):
    """
    Parameters
    ----------
    start_time : str etc
        isoformat or anything Time handles in local time used to start checking
    using : float
        galactic b value to use [deg]
    el_starting : float
        Elevation to start using track [deg]
    time_to_track : float
        Length of track in minutes
    anglerate : float
        galactic l step for track [deg/sec]
    tstep : float
        Trajectory time step in seconds
    tz : float
        Hours from UTC

    """
    if ttype != 'gal':
        print(f"Incorrect ttype: {ttype}")
        return
    tz = tz * u.hour
    start_time = Time(start_time) - tz
    fixed_b = using * u.deg

    # Get starting galactic longitude to nearest 1/2 degree and plot
    gp_l = np.arange(0.0, 359.9, 0.1) * u.deg
    gp_b = np.ones(len(gp_l)) * fixed_b
    gal = SkyCoord(frame='galactic', l=gp_l, b=gp_b)
    radec = gal.transform_to('icrs')
    azel = radec.transform_to(AltAz(obstime=start_time, location=hcro))

    above_horizon = np.where(azel.alt.value > el_starting)
    if not len(above_horizon[0]):
        print("Nothing above the supplied horizon")
        return False
    start_horizon = above_horizon[0][0]
    starting_l = gal.l.value[start_horizon] * u.deg
    plt.figure('At T=0')
    plt.subplot(211)
    plt.plot(radec.ra.value[above_horizon], radec.dec.value[above_horizon], '.', label='radec-init')
    plt.plot(azel.az.value[above_horizon], azel.alt.value[above_horizon], '.', label='azel-init')

    # Get the track
    track_times, track_tsns = _get_track_times(start_time, time_to_track * 60.0, float(tstep))
    this_track = Track(f'Galactic b={using}', timestamp=track_tsns, obstime=track_times)
    tstep = tstep * u.second
    anglerate = anglerate * u.deg / u.second

    for i in range(len(track_times)):
        this_l = starting_l  + i*anglerate*tstep
        if this_l > 360.0*u.deg:
            this_l = this_l - 360.0*u.deg
        gal = SkyCoord(frame='galactic', l=this_l, b=fixed_b)
        radec = gal.transform_to('icrs')
        azel = radec.transform_to(AltAz(obstime=track_times[i], location=hcro))
        this_track.add(az=azel.az.value, el=azel.alt.value, ra=radec.ra.value, dec=radec.dec.value, l=this_l.value, b=fixed_b.value)

    this_track.log_track(TRACK_LOG_FILENAME)
    this_track.ephem_file(EPHEM_FILENAME)
    this_track.full_file(TRAJECTORY_FILENAME)

    plt.plot(this_track.az, this_track.el, '.', lw=4, label='azel-track')
    plt.grid()
    plt.legend()
    plt.xlabel('RA/Az')
    plt.ylabel('Dec/El')

    return this_track


def from_file(ttype='file', start_time='2023-12-31 23:59:59', using='', el_starting=30.0, time_to_track=20.0, anglerate=0.1, tstep=1.0, tz=-8.0):
    """
    Parameters
    ----------
    start_time : str etc
        isoformat or anything Time handles in local time used to start checking
    using : float
        galactic b value to use [deg]
    el_starting : float
        Elevation to start using track [deg]
    time_to_track : float
        Length of track in minutes
    anglerate : float
        galactic l step for track [deg/sec]
    tstep : float
        Trajectory time step in seconds
    tz : float
        Hours from UTC

    """
    from scipy import interpolate
    if ttype != 'file':
        print(f"Incorrect ttype: {ttype}")
        return
    # Set up track
    tz = tz * u.hour
    start_time = Time(start_time) - tz
    track_times, track_tsns = _get_track_times(start_time, time_to_track * 60.0, float(tstep))
    this_track = Track('Out', timestamp=track_tsns, obstime=track_times)

    tstep = tstep * u.second
    anglerate = anglerate * u.deg / u.second

    # Read input file
    overlap = False
    input_track = Track('In')
    with open(using, 'r') as fp:
        for line in fp:
            data = line.split(',')
            obstime = Time(data[0], format='isot')
            if not overlap and obstime >= track_times[0] and obstime <= track_times[-1]:
                overlap = True
            ts = int(obstime.unix_tai * 1E9)
            input_track.add(timestamp=ts, obstime=obstime, az=float(data[1]), el=float(data[2]), ir=1.0 / float(data[3]))
    if not overlap:
        print("This supplied file and desired track don't overlap")
        return False
    input_track.obstime = Time(input_track.obstime)

    # Overlap, so fit az, el, ir
    f = interpolate.interp1d(input_track.timestamp, input_track.az, fill_value=0.0, bounds_error=False)
    this_track.az = f(this_track.timestamp)
    f = interpolate.interp1d(input_track.timestamp, input_track.el, fill_value=0.0, bounds_error=False)
    this_track.el = f(this_track.timestamp)
    f = interpolate.interp1d(input_track.timestamp, input_track.ir, fill_value=0.0, bounds_error=False)
    this_track.ir = f(this_track.timestamp)

    plt.figure('Input')
    plt.subplot(211)
    plt.plot(input_track.obstime.datetime, input_track.az)
    plt.plot(input_track.obstime.datetime, input_track.el)

    return this_track

if __name__ == '__main__':
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument('ttype', help="Type of trajectory ([gal]/file)", nargs='?', default='gal')
    ap.add_argument('-s', '--start_time', help="Time to start checking ['now+5' min]", default='now+5')
    ap.add_argument('-u', '--using', help="Galactic latitude to use if 'gal' or filename to use [0.0 deg]", default=0.0)
    ap.add_argument('-e', '--el_starting', help='Elevation to start at [30.0 deg]', type=float, default=30.0)
    ap.add_argument('-t', '--time_to_track', help='Length of track [10.0 min]', type=float, default=10.0)
    ap.add_argument('-a', '--anglerate', help="Galactic longitude step in trajectory [0.1 deg/sec]", type=float, default=0.1)
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

    if args.ttype == 'gal':
        traj = galactic(**vars(args))
    else:
        traj = from_file(**vars(args))
    if traj:
        print(traj.name)
        plt.subplot(212)
        plt.plot(traj.obstime.datetime, traj.az, '.', label='Az')
        plt.plot(traj.obstime.datetime, traj.el, '.', label='El')
        if len(traj.ra):
            plt.plot(traj.obstime.datetime, traj.ra, '.', label='RA')
            plt.plot(traj.obstime.datetime, traj.dec, '.', label='Dec')
        plt.grid()
        plt.legend()
        plt.show()