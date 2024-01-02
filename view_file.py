#! /usr/bin/python3

import numpy as np
import matplotlib.pyplot as plt
import h5py
from astropy.time import Time
from datetime import timedelta, datetime
from copy import copy


def proc_plot_kwargs(kwargs, defaults, time_fmt='%Y-%m-%dT%H:%M:%S'):
    """
    Process all of the keywords with provided defaults.

    """
    for key, val in defaults.items():
        if key not in kwargs:
            kwargs[key] = val
    
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
        kwargs['freq'] = [float(x) for x in kwargs['freq'].split(',')]
    else:
        kwargs['freq'] = None
    if 'time' in kwargs and kwargs['time'] is not None:
        kwargs['time'] = [datetime.strptime(x, time_fmt) for x in kwargs['time'].split(",")]
    else:
        kwargs['time'] = None

    return kwargs


def _fmt_data(data, plog, pmult=1.0):
    if plog:
        return pmult * np.log10(data)
    else:
        return data


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
        Return the generated shifted/scaled data array for the Axis.

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
            this_val += timedelta(hours=self.shift)
        else:
            this_val += self.shift
        for i in range(self.length):
            arrval = this_val if self.type == 'datetime' else self.scale * this_val
            arr.append(arrval)
            this_val += self.step
        return arr
            
    def index(self, value):
        if isinstance(value, list):
            inds = [int( (x - self.start) / self.step) for x in value]
        else:
            inds = int( (value - self.start) / self.step)
        return inds


class Data:
    def __init__(self, fn, timezone=-8.0):
        self.filename = fn
        self.timezone = timezone
        self.name, self.filename_datetime = self._parse_fn()  # Parse time out of filename (hopefully)
        with h5py.File(fn, 'r') as fp:
            self.data = np.array(fp['data'])
            self.jdstart = np.float64(fp['tstart'])  # jd
            self.jdstop = np.float64(fp['tstop'])  # jd
            self.fcen = np.float64(fp['fcen'])  # MHz
            self.bw = np.float64(fp['bw'])  # MHz
            try:
                self.decimation = np.float64(fp['decimation'])
                self.nfft = np.float64(fp['nfft'])
                self.int_time = self.decimation /self.bw * self.nfft                
            except KeyError:
                self.decimation = None
                self.nfft = None
                self.int_time = None
            try:
                self.tle_update_jd = np.float64(fp['tle'])  # jd
                self.tle_update = Time(self.tle_update, format='jd')
            except KeyError:
                self.tle_update_jd = None
                self.tle_update = None
        # Set time axis
        self.tstart = Time(self.jdstart, format='jd')
        self.tstop = Time(self.jdstop, format='jd')
        self.t_info = Axis(self.tstart.datetime, self.tstop.datetime, len(self.data[:, 0]), 'UTC')
        self.t = self.t_info.array(self.timezone)
        print(self.t_info)
        # Set freq axis
        self.fmin = self.fcen - self.bw / 2.0
        self.fmax = self.fcen + self.bw / 2.0
        self.f_info = Axis(self.fmin, self.fmax, len(self.data[0]), 'MHz')
        self.freq = self.f_info.array()
        print(self.f_info)

    def _parse_fn(self):
        X = self.filename.split('_')
        if len(X) == 3:
            _date = '20' + '-'.join( [X[1][i:i+2] for i in range(0, len(X[1]), 2)])
            _time = ':'.join([X[2][i:i+2] for i in range(0, len(X[2].split('.')[0]), 2)])
            _datetime = Time(f"{_date}T{_time}").datetime
        else:
            _datetime = None
        return X[0], _datetime

    def header(self):
        """Print information about the data."""

        print(f"Filename: {self.filename}")
        print(f"Data shape: {np.shape(self.data)}")
        print(f"Start: UTC {self.tstart.datetime} (jd={self.jdstart})")
        print(f"\tLocal: {self.tstart.datetime + timedelta(hours=self.timezone)}")
        print(f"Stop: UTC {self.tstop.datetime} (jd={self.jdstop})")
        print(f"\tLocal: {self.tstop.datetime + timedelta(hours=self.timezone)}")
        print(f"Freq:  {self.fmin:.2f} - {self.fmax:.2f} MHz  (cf = {self.fcen}, bw = {self.bw} MHz)")

    def wf(self, **kwargs):
        """
        Make a waterfall plot of the data.

        Parameters:
        -----------
        kwargs:
            num_xticks (int), num_yticks (int), colorbar (bool), log (bool), dB

        """
        defaults = {'log': True, 'dB': False, 'colorbar': True, 'xticks': 12, 'yticks': 6}
        kwargs = proc_plot_kwargs(kwargs, defaults)

        num_xticks = kwargs['xticks']
        num_yticks = kwargs['yticks']

        if kwargs['log']:
            plt.imshow(kwargs['_mult'] * np.log10(self.data))
        else:
            plt.imshow(self.data)
        if kwargs['colorbar']:
            plt.colorbar()

        plt.xticks(np.linspace(0, len(self.data[0]), num_xticks), [f"{x:.2f}" for x in np.linspace(self.fmin, self.fmax, num_xticks)])
        jds = np.linspace(self.jdstart, self.jdstop, num_yticks)
        apt = Time(jds, format='jd')
        yticks = [(x + timedelta(hours=self.timezone)).strftime("%H:%M:%S") for x in apt.datetime]
        plt.yticks(np.linspace(0, len(self.data), num_yticks), yticks)
        plt.xlabel(self.f_info.label)
        plt.ylabel(self.t_info.label)
        plt.title(f"{self.t_info.unit}:  {self.t_info.start.strftime('%Y-%m-%d')}")
        plt.tight_layout()

    def _get_ft_slices(self, kwargs):
        if kwargs['freq'] is not None:
            frange = self.f_info.index(kwargs['freq'])
        else:
            frange = [0, self.f_info.length-1]
        fslice = slice(frange[0], frange[1] + 1)
        if kwargs['time'] is not None:
            trange = self.t_info.index(kwargs['time'])
        else:
            trange = [0, self.t_info.length - 1]
        tslice = slice(trange[0], trange[1] + 1)

        return fslice, tslice

    def spectra(self, **kwargs):
        """Make a 2-D plot of the spectra."""
        defaults = {'log': False, 'dB': False}
        kwargs = proc_plot_kwargs(kwargs, defaults)
        fslice, tslice = self._get_ft_slices(kwargs)
        trange = range(tslice.start, tslice.stop)

        for i in trange:
            data = _fmt_data(self.data[i][fslice], kwargs['log'], kwargs['_mult'])
            plt.plot(self.freq[fslice], data)
        plt.grid()
        plt.xlabel(self.f_info.label)
        plt.ylabel(kwargs['_ylabel'])

    def series(self, **kwargs):
        """Make a 2-D plot of the time series."""
        defaults = {'log': False, 'dB': False}
        kwargs = proc_plot_kwargs(kwargs, defaults)
        fslice, tslice = self._get_ft_slices(kwargs)
        frange = range(fslice.start, fslice.stop)
        trange = range(tslice.start, tslice.stop)

        for i in frange:
            data = _fmt_data(self.data[trange, i], kwargs['log'], kwargs['_mult'])
            plt.plot(self.t[tslice], data)
        plt.grid()
        plt.xlabel(self.t_info.label)
        plt.ylabel(kwargs['_ylabel'])                

        results_prefix = '--Results--\n' if kwargs['beamfit'] else ''
        yaxlim = [plt.axis()[2], plt.axis()[3]]
        N = 10
        if self.filename_datetime is not None:
            print(f"{results_prefix}{'Expected:':{N}s}{self.filename_datetime.isoformat()}")
            if self.filename_datetime>=self.t[tslice.start] and self.filename_datetime<=self.t[tslice.stop-1]:  
                plt.plot([self.filename_datetime, self.filename_datetime], yaxlim, '--', lw=2, color='k')
        if kwargs['beamfit']:
            import beamfit
            power = np.zeros(self.t_info.length, dtype=float)
            for i in frange:
                power += self.data[:, i] / len(frange)
            data = _fmt_data(power, kwargs['log'], kwargs['_mult'])
            plt.plot(self.t[tslice], data[tslice], lw=4, color='k')
            coeff, data_fit = beamfit.gaussian(power[tslice], max(power[tslice]), len(power[tslice]) / 2.0, len(power[tslice])/4.0 )
            fit_time = int(coeff[1]) + tslice.start
            fit_range = [int(coeff[1] - coeff[2]), int(coeff[1] + coeff[2])]
            data = _fmt_data(data_fit, kwargs['log'], kwargs['_mult'])
            plt.plot(self.t[tslice], data, '--', lw=2, color='w')
            plt.plot([self.t[fit_time], self.t[fit_time]], yaxlim, '--', lw=2, color='k')
            print(f"{'Found:':{N}s}{self.t[fit_time].isoformat()}")
            if self.filename_datetime is not None:
                  offset = self.t[fit_time] - self.filename_datetime
                  print(f"{'Offset:':{N}s}{offset.total_seconds():.1f} sec")
            width = self.t[fit_range[1]] - self.t[fit_range[0]]
            print(f"{'Width:':{N}s}{width.total_seconds():.1f} sec")
        plt.grid()
        plt.xlabel(self.t_info.label)
        plt.ylabel(kwargs['_ylabel'])


if __name__ == '__main__':
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument('fn', help="Name of hdf5 datafile")
    ap.add_argument('-p', '--plot_type', help='wf, series, spectra [wf]', choices=['wf', 'series', 'spectra'], default='wf')
    ap.add_argument('-x', '--xticks', help="Number of xticks in waterfall [10]", type=int, default=10)
    ap.add_argument('-y', '--yticks', help="Number of yticks to use in waterfall [4]", type=int, default=4)
    ap.add_argument('-c', '--colorbar', help="Flag to hide colorbar", action='store_false')
    ap.add_argument('-b', '--beamfit', help='Flag to fit for the beam', action='store_true')
    ap.add_argument('-f', '--freq', help='Frequency [range] to use.', default=None)
    ap.add_argument('-t', '--time', help='Time [range] to use.', default=None)
    ap.add_argument('-l', '--log', help="Flag to take log10 of data", action='store_true')
    ap.add_argument('-d', '--dB', help="Flag to convert to dB", action='store_true')
    ap.add_argument('--tz', help="Timezone offset to UTC in hours [-8.0]", type=float, default=-8.0)
    args = ap.parse_args()
    obs = Data(args.fn, args.tz)
    getattr(obs, args.plot_type)(**vars(args))
    plt.show()
