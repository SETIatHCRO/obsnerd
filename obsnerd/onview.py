#! /usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import h5py
from astropy.time import Time
from datetime import timedelta, datetime
from copy import copy
from f2h5 import HDF5HeaderInfo
import onutil


class StateVariable:
    def __init__(self, kwargs={}, defaults={}):
        self.update(kwargs, defaults)

    def __repr__(self):
        from tabulate import tabulate
        table = []
        for a in dir(self):
            if a.startswith('__') or a == 'update':
                continue
            else:
                val = str(getattr(self, a))
                if len(val) > 50:
                    val = val[:50] + '...'
                table.append([a, val])
        return tabulate(table, headers=['Variable', 'Value'])

    def update(self, kwargs={}, defaults={}):
        """
        Process all of the keywords with provided defaults.
        """
        # Copy over defaults if not present
        for key, val in defaults.items():
            if key not in kwargs:
                kwargs[key] = val

        if 'tz' not in kwargs:
            kwargs['tz'] = 0.0
        if 'log' in kwargs and kwargs['log']:
            kwargs['_ylabel'] = 'log'
            kwargs['_mult'] = 1.0
        else:
            kwargs['log'] = False
            kwargs['_ylabel'] = 'linear'
            kwargs['_mult'] = 1.0
        if  'dB' in kwargs and kwargs['dB']:
            kwargs['_mult'] = 10.0
            kwargs['log'] = True
            kwargs['_ylabel'] = 'dB'
        else:
            kwargs['dB'] = False
        if 'freq' in kwargs and kwargs['freq'] is not None:
            if isinstance(kwargs['freq'], str):
                kwargs['freq'] = [float(x) for x in kwargs['freq'].split(',')]
        else:
            kwargs['freq'] = None
        if 'time' in kwargs and kwargs['time'] is not None:
            if isinstance(kwargs['time'], str):
                kwargs['time'] = [onutil.make_datetime(x) for x in kwargs['time'].split(",")]
        else:
            kwargs['time'] = None

        for key, value in kwargs.items():
            setattr(self, key, value)


class Axis:
    """
    Class handling the Axis definitions, conversions, etc.

    """
    def __init__(self, start, stop, length, unit=""):
        """
        Parameters
        ----------
        start : Axis value (datetime or number)
            Initial value on the Axis
        stop : Axis value (datetime or number)
            Final value on the Axis
        length : int
            Number of data points on the Axis
        unit : str
            Unit of Axis

        """
        self.start = start
        self.stop = stop
        self.length = int(length)
        self.unit = unit
        self.type = 'datetime' if isinstance(start, datetime) else 'number'
    
    def __repr__(self):
        return f"Axis:  {self.type:8s} |  {self.start} - {self.stop} {self.unit} (n={self.length})"

    def array(self, shift=0.0, scale=1.0, label=None):
        """
        Return the generated shifted/scaled data array for the Axis.  Time shift will take UTC to local

        Parameters
        ----------
        shift : float
            Value to shift Axis -- if datetime this is hours, otherwise it is self.unit
        scale : float
            Value to scale Axis -- if datetime this is ignored
        label : str (None to ignore)
            Axis label -- use modified self.unit if None

        """
        self.shift = float(shift)
        self.scale = float(scale)
        self.step = (self.stop - self.start) / self.length
        if label is None:
            if abs(self.shift) < 1E-6:
                self.label = self.unit
            else:
                self.label = f"{self.unit}{'+' if self.shift>0 else '-'}{abs(self.shift):.0f}"
            if abs(scale - 1.0) > 1E-6:
                if self.type == 'number':
                    self.label = f"{self.scale} x {self.label}"  # Ideally change unit name, e.g. MHz -> Hz
                elif self.type == 'datetime':
                    self.scale = 1.0
                    print("Can't scale time -- ignoring non-unity value.")
        else:
            self.label = label

        arr = []
        this_val = copy(self.start)
        if self.type == 'datetime':
            this_val -= timedelta(hours=self.shift)
        else:
            this_val -= self.shift
        for i in range(self.length):
            arrval = this_val if self.type == 'datetime' else self.scale * this_val
            arr.append(arrval)
            this_val += self.step
        return arr
            
    def index(self, value):
        if value is None:
            inds = [0, self.length-1]
        elif isinstance(value, list):
            inds = [int( (x - self.start) / self.step) for x in value]
        else:
            inds = int( (value - self.start) / self.step)
        return inds


class Data:
    def __init__(self, fn, timezone=0.0):
        self.filename = fn
        self.timezone = timezone
        self.name, self.filename_datetime = self._parse_fn()  # Parse time out of filename (hopefully)
        self.h5 = HDF5HeaderInfo()
        with h5py.File(fn, 'r') as fp:
            self.data = np.array(fp[self.h5.data])
            for param in self.h5.float64s:
                if param in fp:
                    setattr(self, param, float(np.float64(fp[param])))
                else:
                    setattr(self, param, None)
            for param in self.h5.strings:
                if param in fp:
                    setattr(self, param, fp[param].asstr()[...])
        # Set time axis
        for field in self.h5.from_datetime:
            setattr(self, field, Time(getattr(self, field), format='jd'))
        self.t_info = Axis(self.tstart.datetime, self.tstop.datetime, len(self.data[:, 0]), 'UTC')
        self.t = self.t_info.array(self.timezone)
        print(self.t_info)
        # Set freq axis
        self.fmin = self.fcen - self.bw / 2.0
        self.fmax = self.fcen + self.bw / 2.0
        self.f_info = Axis(self.fmin, self.fmax, len(self.data[0]), 'MHz')
        self.freq = self.f_info.array()
        print(self.f_info)
        # Other info
        missing_data = self.decimation is None or self.bw is None or self.nfft is None
        self.int_time = None if missing_data else (self.decimation / self.bw) * self.nfft
        self.tle = None if self.tle is None else Time(self.tle, format='jd')
        if self.tle is not None:
            print(f"Updated TLEs: {self.tle.datetime}")
        self.sv = StateVariable(defaults={'tz': 0.0})

    def _parse_fn(self):
        X = self.filename.split('_')
        if len(X) == 3:
            _datetime = onutil.make_datetime(date=f"{X[1]}_{X[2]}")
        else:
            _datetime = None
        return X[0], _datetime

    def _fmt_data(self, data):
        if self.sv.log:
            return self.sv._mult * np.log10(data)
        else:
            return data

    def showdata(self, **kwargs):
        """Print information about the data."""

        print(f"Filename: {self.filename}")
        print(f"Data shape: {np.shape(self.data)}")
        for field in self.h5.metadata:
            val = getattr(self, field)
            if field in self.h5.from_datetime:
                print(f"{field}:  {val.datetime}  (jd={val.jd})")
            else:
                print(f"{field}: {val}")


    def wf(self, **kwargs):
        """
        Make a waterfall plot of the data.

        Parameters:
        -----------
        kwargs:
            num_xticks (int), num_yticks (int), colorbar (bool), log (bool), dB

        """
        defaults = {'log': True, 'dB': False, 'colorbar': True, 'xticks': 12, 'yticks': 6}
        self.sv.update(kwargs, defaults)

        plt.imshow(self._fmt_data(self.data))
        if self.sv.colorbar:
            plt.colorbar()

        plt.xticks(np.linspace(0, len(self.data[0]), self.sv.xticks), [f"{x:.2f}" for x in np.linspace(self.fmin, self.fmax, self.sv.xticks)])
        jds = np.linspace(self.tstart.jd, self.tstop.jd, self.sv.yticks)
        apt = Time(jds, format='jd')
        yticks = [(x + timedelta(hours=self.timezone)).strftime("%H:%M:%S") for x in apt.datetime]
        plt.yticks(np.linspace(0, len(self.data), self.sv.yticks), yticks)
        plt.xlabel(self.f_info.label)
        plt.ylabel(self.t_info.label)
        plt.title(f"{self.t_info.unit}:  {self.t_info.start.strftime('%Y-%m-%d')}")
        plt.tight_layout()

    def _get_ft_slices(self):
        inds = self.f_info.index(self.sv.freq)
        self.fslice = slice(inds[0], inds[-1] + 1)
        self.frange = range(self.fslice.start, self.fslice.stop)
        inds = self.t_info.index(self.sv.time)
        self.tslice = slice(inds[0], inds[-1] + 1)
        self.trange = range(self.tslice.start, self.tslice.stop)

    def spectra(self, **kwargs):
        """Make a 2-D plot of the spectra."""
        defaults = {'log': False, 'dB': False, 'freq': None, 'time': None}
        self.sv.update(kwargs, defaults)
        self._get_ft_slices()

        for i in self.trange:
            plt.plot(self.freq[self.fslice], self._fmt_data(self.data[i][self.fslice]))
        plt.grid()
        plt.xlabel(self.f_info.label)
        plt.ylabel(self.sv._ylabel)

    def make_power(self, **kwargs):
        """
        Get integrated power for each timestep

        """
        defaults = {'norm': '/bin', 'freq': None, 'time': None}
        self.sv.update(kwargs, defaults)
        self._get_ft_slices()
        self.power = np.zeros(self.t_info.length, dtype=float)
        for i in self.frange:
            self.power += self.data[:, i]
        if self.sv.norm == '/bin':
            self.power /= len(self.frange)
        elif self.sv.norm == '/full':
            self.power /= (self.freq[self.frange[0]] - self.freq[self.frange[-1]])

    def load_trajectory(self, filename='track.npz'):
        print(f"Loading {filename} to self.traj")
        self.traj = np.load(filename, allow_pickle=True)
        t_initial = self.traj['obstime'][0]
        l_initial = self.traj['l'][0]
        dldt = (self.traj['l'][-1] - l_initial) / (self.traj['obstime'][-1] - t_initial).total_seconds()
        print(f"Starting:  l={l_initial},  t={t_initial}")
        print(f"dl/dt = {dldt}")
        self.l = []
        for this_time in self.t:
            this_time = this_time - timedelta(hours = -8)
            self.l.append(l_initial + dldt * (this_time - t_initial).total_seconds())
        plt.figure(f"Galaxy @ b = {self.traj['b'][0]}")
        plt.plot(self.l, 10.0*np.log10(self.power))
        plt.grid()
        plt.xlabel(' l [deg]')
        plt.ylabel('dB')

    def series(self, **kwargs):
        """Make a 2-D plot of the time series."""
        defaults = {'log': False, 'dB': False, 'total_power': False, 'freq': None, 'time': None}
        self.sv.update(kwargs, defaults)
        self._get_ft_slices()

        for i in self.frange:
            plt.plot(self.t[self.tslice], self._fmt_data(self.data[self.trange, i]))               
        plt.grid()
        plt.xlabel(self.t_info.label)
        plt.ylabel(self.sv._ylabel)

        # Do extra fitting stuff below
        N = 10
        results_prefix = '--Results--\n' if self.sv.beamfit else ''
        if self.sv.beamfit or self.sv.total_power:
            self.make_power(**kwargs)
            plt.plot(self.t[self.tslice], self._fmt_data(self.power[self.tslice]), lw=4, color='k')
        yaxlim = [plt.axis()[2], plt.axis()[3]]
        if self.filename_datetime is not None:
            print(f"{results_prefix}{'Expected:':{N}s}{self.filename_datetime.isoformat()}")
            if self.filename_datetime>=self.t[self.tslice.start] and self.filename_datetime<=self.t[self.tslice.stop-1]:  
                plt.plot([self.filename_datetime, self.filename_datetime], yaxlim, '--', lw=2, color='k')

        if self.sv.beamfit:
            import beamfit
            coeff, data_fit = beamfit.gaussian(self.power[self.tslice], max(self.power[self.tslice]), len(self.trange) / 2.0, len(self.trange)/4.0 )
            fit_time = int(coeff[1]) + self.tslice.start
            fit_range = [int(coeff[1] - coeff[2]), int(coeff[1] + coeff[2])]
            plt.plot(self.t[self.tslice], self._fmt_data(data_fit), '--', lw=2, color='w')
            plt.plot([self.t[fit_time], self.t[fit_time]], yaxlim, '--', lw=2, color='k')
            print(f"{'Found:':{N}s}{self.t[fit_time].isoformat()}")
            if self.expected is not None:
                  offset = self.t[fit_time] - self.expected
                  print(f"{'Offset:':{N}s}{offset.total_seconds():.1f} sec")
            width = self.t[fit_range[1]] - self.t[fit_range[0]]
            print(f"{'Width:':{N}s}{width.total_seconds():.1f} sec")



if __name__ == '__main__':
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument('fn', help="Name of hdf5 datafile")
    ap.add_argument('-o', '--output_type', help='wf, series, spectra [wf]', choices=['wf', 'series', 'spectra', 'showdata'], default='wf')
    ap.add_argument('-x', '--xticks', help="Number of xticks in waterfall [10]", type=int, default=10)
    ap.add_argument('-y', '--yticks', help="Number of yticks to use in waterfall [4]", type=int, default=4)
    ap.add_argument('-c', '--colorbar', help="Flag to hide colorbar", action='store_false')
    ap.add_argument('-b', '--beamfit', help='Flag to fit for the beam', action='store_true')
    ap.add_argument('-f', '--freq', help='Frequency [range] to use.', default=None)
    ap.add_argument('-t', '--time', help='Time [range] to use.', default=None)
    ap.add_argument('-l', '--log', help="Flag to take log10 of data", action='store_true')
    ap.add_argument('-d', '--dB', help="Flag to convert to dB", action='store_true')
    ap.add_argument('-P', '--total_power', help="Show total power (for series)", action='store_true')
    ap.add_argument('-n', '--norm', help="Norm to apply for total power (Total, [/bin], /full)", choices=['Total', '/bin', '/full'], default='/bin')
    ap.add_argument('--tz', help="Timezone offset to UTC in hours [-8.0]", type=float, default=-8.0)
    args = ap.parse_args()
    obs = Data(args.fn, args.tz)
    getattr(obs, args.output_type)(**vars(args))
    plt.show()
