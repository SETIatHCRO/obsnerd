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
HCRO = EarthLocation(lat=40.8178049*u.deg, lon=-121.4695413*u.deg, height=986*u.m)
EL_TRACK_LIMIT = 1.5
AZ_TRACK_LIMIT = 2.0

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
        s = f'Label: {self.name}\n'
        for i, lbl in zip([0, -1], ['Start at:', 'End at:']):
            s += f"{lbl} {self.obstime[i]} UTC\n"
            try:
                s += f"\tl={self.l[i]:.3f}, b={self.b[i]:.3f}\n"
            except IndexError:
                pass
            try:
                s += f"\tRA={self.ra[i]:.3f}, Dec={self.dec[i]:.3f}\n"
            except IndexError:
                pass
            try:
                s += f"\tAz={self.az[i]:.3f}, El={self.el[i]:.3f}\n"
            except IndexError:
                pass
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
        if True:
            print("NEED TO GET DISTANCES FROM SOPP")
            self.ir = np.zeros(len(tai)) + 1E-20
        ephem = ((np.array([tai, az, el, self.ir], dtype=object)))
        self.ephemtxt = np.savetxt(filename, ephem.T, fmt='%i  %.5f  %.5f  %.10E')

    def full_file(self, filename=TRAJECTORY_FILENAME):
        print(F"Writing full trajectory: {filename}")
        np.savez_compressed(filename, obstime=self.obstime.datetime,
                            az=self.az, el=self.el, ra=self.ra, dec=self.dec, 
                            l=self.l, b=self.b, ir=self.ir, allow_pickle=True)

    def rates(self):
        dt = np.diff(self.timestamp) / 1E9
        self.dazdt = np.diff(self.az) / dt
        self.deldt = np.diff(self.el) / dt


class Trajectory:
    def __init__(self, start_time='2023-12-31 23:59:59', time_to_track=20.0, tstep=1.0, el_horizon=30.0, tz=-8.0):
        """
        Parameters
        ----------
        start_time : str etc
            isoformat or anything Time handles in local time used to start checking
        time_to_track : float
            Length of track in minutes
        tstep : float
            Trajectory time step in seconds
        el_horizon : float
            Elevation to start using track [deg]
        tz : float
            Hours from UTC

        """
        self.valid = True
        self.tz = tz * u.hour
        self.start_time = Time(start_time) - self.tz
        self.el_horizon = float(el_horizon)
        self.time_to_track = float(time_to_track) * 60.0 * u.second
        self.tstep = float(tstep) * u.second
        self.track_Time = self.start_time + np.arange(0.0, self.time_to_track.value, self.tstep.value) * u.second
        self.track_tsns = int(self.start_time.unix_tai * 1E9) + np.arange(0, int(self.time_to_track.value * 1E9), int(self.tstep.value * 1E9))

    def galactic(self, b=0.0, rate=0.1):
        self.trajectory_type = 'galactic'
        self.b = b * u.deg
        self.rate = rate * u.deg / u.second

        # Get starting galactic longitude to nearest 1/10 degree and plot
        self.T0 = argparse.Namespace()
        gp_l = np.arange(0.0, 359.9, 0.1) * u.deg
        gp_b = np.ones(len(gp_l)) * self.b
        self.T0.gal = SkyCoord(frame='galactic', l=gp_l, b=gp_b)
        self.T0.radec = self.T0.gal.transform_to('icrs')
        self.T0.azel = self.T0.radec.transform_to(AltAz(obstime=self.start_time, location=HCRO))

        self.T0.above_horizon = np.where(self.T0.azel.alt.value > self.el_horizon)
        if not len(self.T0.above_horizon[0]):
            print("Nothing above the supplied horizon")
            self.valid = False
            return
        starting_l = self.T0.gal.l[self.T0.above_horizon[0][0]]

        # Get the track
        self.track = Track(f'Galactic b={self.b}', timestamp=self.track_tsns, obstime=self.track_Time)        
        R = self.rate*self.tstep
        for i in range(len(self.track_Time)):
            this_l = starting_l  + i*R
            if this_l > 360.0*u.deg:
                this_l = this_l - 360.0*u.deg
            gal = SkyCoord(frame='galactic', l=this_l, b=self.b)
            radec = gal.transform_to('icrs')
            azel = radec.transform_to(AltAz(obstime=self.track_Time[i], location=HCRO))
            self.track.add(az=azel.az.value, el=azel.alt.value, ra=radec.ra.value, dec=radec.dec.value, l=this_l.value, b=self.b.value)
        self.track.rates()

    def from_file(self, filename):
        """
        Parameters
        ----------
        filename : str
            File containing the data

        """
        self.trajectory_type = 'from_file'
        from scipy import interpolate
        self.track = Track('Out', timestamp=self.track_tsns, obstime=self.track_Time)
        # Read input file
        overlap = False
        self.input_track = Track('In')
        with open(filename, 'r') as fp:
            for line in fp:
                data = line.split(',')
                obstime = Time(data[0], format='isot')
                if not overlap and obstime >= self.track_Time[0] and obstime <= self.track_Time[-1]:
                    overlap = True
                ts = int(obstime.unix_tai * 1E9)
                self.input_track.add(timestamp=ts, obstime=obstime, az=float(data[1]), el=float(data[2]), ir=1.0 / float(data[3]))
        if not overlap:
            print("This supplied file and desired track don't overlap")
            self.valid = False
        self.input_track.obstime = Time(self.input_track.obstime)

        # Overlap, so fit az, el, ir
        f = interpolate.interp1d(self.input_track.timestamp, self.input_track.az, fill_value=0.0, bounds_error=False)
        self.track.az = f(self.track.timestamp)
        f = interpolate.interp1d(self.input_track.timestamp, self.input_track.el, fill_value=0.0, bounds_error=False)
        self.track.el = f(self.track.timestamp)
        f = interpolate.interp1d(self.input_track.timestamp, self.input_track.ir, fill_value=0.0, bounds_error=False)
        self.track.ir = f(self.track.timestamp)
        self.track.rates()

    def write(self):
        self.track.log_track(TRACK_LOG_FILENAME)
        self.track.ephem_file(EPHEM_FILENAME)
        self.track.full_file(TRAJECTORY_FILENAME)

    def plot(self):
        if self.trajectory_type == 'galactic':
            self._plot_galactic()
        elif self.trajectory_type == 'from_file':
            self._plot_from_file()
        plt.figure('Rates')
        t_extremes = [self.track.obstime.datetime[1], self.track.obstime.datetime[-1]]
        sgn = '--' if np.any(self.track.dazdt < 0.0) else '-'
        plt.plot(self.track.obstime.datetime[1:], abs(self.track.dazdt), 'k'+sgn, label='daz/dt')
        plt.plot(t_extremes, [AZ_TRACK_LIMIT, AZ_TRACK_LIMIT], 'k-.', lw=3)
        sgn = '--' if np.any(self.track.deldt < 0.0) else '-'
        plt.plot(self.track.obstime.datetime[1:], abs(self.track.deldt), 'b'+sgn, label='del/dt')
        plt.plot(t_extremes, [EL_TRACK_LIMIT, EL_TRACK_LIMIT], 'b-.', lw=3)
        plt.legend()
        plt.ylabel('[deg/sec]')
        plt.grid()

    def _plot_galactic(self):
        plt.figure('Galactic')
        plt.subplot(211)
        plt.plot(self.T0.radec.ra.value[self.T0.above_horizon], self.T0.radec.dec.value[self.T0.above_horizon], '.', label='radec-init')
        plt.plot(self.T0.azel.az.value[self.T0.above_horizon], self.T0.azel.alt.value[self.T0.above_horizon], '.', label='azel-init')
        plt.plot(self.track.az, self.track.el, '.', lw=4, label='azel-track')
        plt.grid()
        plt.legend()
        plt.xlabel('RA/Az [deg]')
        plt.ylabel('Dec/El [deg]')
        plt.subplot(212)
        plt.plot(self.track.obstime.datetime, self.track.az, '.', label='Az')
        plt.plot(self.track.obstime.datetime, self.track.el, '.', label='El')
        plt.plot(self.track.obstime.datetime, self.track.ra, '.', label='RA')
        plt.plot(self.track.obstime.datetime, self.track.dec, '.', label='Dec')
        plt.xlabel('Time')
        plt.ylabel('[deg]')
        plt.grid()
        plt.legend()

    def _plot_from_file(self):
        plt.figure('From_file')
        plt.plot(self.input_track.obstime.datetime, self.input_track.az, 'k', label='Az')
        plt.plot(self.input_track.obstime.datetime, self.input_track.el, 'b', label='El')
        plt.plot(self.track.obstime.datetime, self.track.az, 'ko')
        plt.plot(self.track.obstime.datetime, self.track.el, 'bo')
        plt.xlabel('Time')
        plt.ylabel('[deg]')
        plt.legend()
        plt.grid()

if __name__ == '__main__':
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

    traj = Trajectory(start_time=args.start_time, time_to_track=args.time_to_track, tstep=args.tstep, el_horizon=args.el_horizon, tz=args.tz)

    if args.ttype == 'gal':
        traj.galactic(b=float(args.using), rate=args.rate)
    else:
        traj.from_file(filename=args.using)
    if args.write_out:
        traj.write()
    if traj.valid:
        print(traj.track.name)
        traj.plot()
        plt.show()