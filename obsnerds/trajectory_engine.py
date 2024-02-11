import astropy.units as u
from astropy.coordinates import AltAz, EarthLocation, SkyCoord
from astropy.time import Time
import numpy as np
import yaml
import matplotlib.pyplot as plt
from argparse import Namespace
from . import onutil

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
        self.summary_yaml()
        return yaml.safe_dump(self.summary)

    def summary_yaml(self):
        self.summary = {'Name':  self.name, 'Start': {}, 'Stop': {}}
        for i, lbl in zip([0, -1], ['Start', 'Stop']):
            self.summary[lbl]['time'] = self.obstime[i].datetime.isoformat()
            for par in ['l', 'b', 'ra', 'dec', 'az', 'el']:
                try:
                    self.summary[lbl][par] = float(getattr(self, par)[i])
                except IndexError:
                    pass

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
    def __init__(self, start_time, duration=20.0, tstep=1.0, el_horizon=30.0, tz=0.0):
        """
        Parameters
        ----------
        start_time : str etc
            isoformat or anything Time handles in local time used to start checking
        duration : float
            Length of track in minutes
        tstep : float
            Trajectory time step in seconds
        el_horizon : float
            Elevation to start using track [deg]
        tz : float
            Hours from UTC

        """
        self.valid = True
        self.tz = tz
        if isinstance(start_time, str) and start_time.lower() == 'none':
            self.start_time = None
        else:
            self.start_time = Time(onutil.make_datetime(date=start_time, timezone=tz))
        self.el_horizon = float(el_horizon)
        self.duration = float(duration) * 60.0 * u.second
        self.tstep = float(tstep) * u.second

    def _time_arrays(self):
        self.track_Time = self.start_time + np.arange(0.0, self.duration.value, self.tstep.value) * u.second
        self.track_tsns = int(self.start_time.unix_tai * 1E9) + np.arange(0, int(self.duration.value * 1E9), int(self.tstep.value * 1E9))

    def galactic(self, b=0.0, rate=0.1):
        self.trajectory_type = 'galactic'
        self.b = b * u.deg
        self.rate = rate * u.deg / u.second
        self._time_arrays()

        # Get starting galactic longitude to nearest 1/10 degree and plot
        self.T0 = Namespace()
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
        # Read input file
        self.input_track = Track(f"Input: {filename}")
        with open(filename, 'r') as fp:
            for line in fp:
                data = line.split(',')
                obstime = Time(data[0], format='isot')
                ts = int(obstime.unix_tai * 1E9)
                self.input_track.add(timestamp=ts, obstime=obstime, az=float(data[1]), el=float(data[2]), ir=1.0 / float(data[3]))
        self.input_track.obstime = Time(self.input_track.obstime)
        if self.start_time is None:
            self.start_time = self.input_track.obstime[0]
        self._time_arrays()
        self.track = Track(f"Output: {filename}", timestamp=self.track_tsns, obstime=self.track_Time)

        overlap = self.track_Time[-1] > self.input_track.obstime[0] and self.track_Time[0] < self.input_track.obstime[-1]   
        if not overlap:
            print("This supplied file and desired track don't overlap")
            self.valid = False
            return

        # Overlap, so fit az, el, ir
        f = interpolate.interp1d(self.input_track.timestamp, self.input_track.el, fill_value=0.0, bounds_error=False)
        _el = f(self.track.timestamp)
        valid = np.where(_el > 0.0)
        self.track.timestamp = self.track.timestamp[valid]
        self.track.obstime = self.track.obstime[valid]
        self.track.el = _el[valid]
        f = interpolate.interp1d(self.input_track.timestamp, self.input_track.az, fill_value=0.0, bounds_error=False)
        self.track.az = f(self.track.timestamp)
        f = interpolate.interp1d(self.input_track.timestamp, self.input_track.ir, fill_value=0.0, bounds_error=False)
        self.track.ir = f(self.track.timestamp)
        self.track.rates()

    def write(self):
        if not self.valid:
            print("Trajectory not valid.")
            return
        self.track.log_track(TRACK_LOG_FILENAME)
        self.track.ephem_file(EPHEM_FILENAME)
        self.track.full_file(TRAJECTORY_FILENAME)

    def plot(self):
        if not self.valid:
            print("Trajectory not valid.")
            return
        if self.trajectory_type == 'galactic':
            self._plot_galactic()
        elif self.trajectory_type == 'from_file':
            self._plot_from_file()
        plt.figure('Rates')
        t_extremes = [self.track.obstime.datetime[1], self.track.obstime.datetime[-1]]
        plt.plot(self.track.obstime.datetime[1:], self.track.dazdt, 'k', label='daz/dt')
        plt.plot(t_extremes, [AZ_TRACK_LIMIT, AZ_TRACK_LIMIT], 'k-.', lw=3)
        plt.plot(t_extremes, [-AZ_TRACK_LIMIT, -AZ_TRACK_LIMIT], 'k-.', lw=3)
        plt.plot(self.track.obstime.datetime[1:], self.track.deldt, 'b', label='del/dt')
        plt.plot(t_extremes, [EL_TRACK_LIMIT, EL_TRACK_LIMIT], 'b-.', lw=3)
        plt.plot(t_extremes, [-EL_TRACK_LIMIT, -EL_TRACK_LIMIT], 'b-.', lw=3)
        plt.legend()
        plt.ylabel('[deg/sec]')
        plt.grid()
        plt.show()

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
        plt.plot(self.input_track.obstime.datetime, self.input_track.az, 'k.', label='Az')
        plt.plot(self.input_track.obstime.datetime, self.input_track.el, 'b.', label='El')
        #plt.plot(self.track.obstime.datetime, self.track.az, 'ko')
        #plt.plot(self.track.obstime.datetime, self.track.el, 'bo')
        plt.xlabel('Time')
        plt.ylabel('[deg]')
        plt.legend()
        plt.grid()
        plt.figure('Input Track')
        plt.plot(self.input_track.az, self.input_track.el, '.')
        plt.xlabel('Az')
        plt.ylabel('El')
        plt.grid()
