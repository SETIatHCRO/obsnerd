#! /usr/bin/python3

import numpy as np
import matplotlib.pyplot as plt
import h5py
from astropy.time import Time
from datetime import timedelta
from copy import copy


def proc_plot_kwargs(kwargs, defaults):
    for key, val in defaults.items():
        if key not in kwargs:
            kwargs[key] = val
    kwargs['mult'] = 1.0
    kwargs['ylabel'] = 'linear'
    if kwargs['log']:
        kwargs['ylabel'] = 'log'
    if  kwargs['dB']:
        kwargs['mult'] = 10.0
        kwargs['log'] = True
        kwargs['ylabel'] = 'dB'
    return kwargs


class Data:
    def __init__(self, fn, timezone=-8.0):
        self.filename = fn
        self.timezone = timezone
        self.name, self.datetime = self._parse_fn()
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
        self.tstart = Time(self.jdstart, format='jd')
        self.tstop = Time(self.jdstop, format='jd')
        self.fmin = self.fcen - self.bw / 2.0
        self.fmax = self.fcen + self.bw / 2.0
        df = (self.fmax - self.fmin) / len(self.data[0])
        self.freq = []
        this_f = self.fmin + 0.0
        for i in range(len(self.data[0])):
            self.freq.append(this_f)
            this_f += df
        dt = (self.tstop.datetime - self.tstart.datetime) / len(self.data[:,0])
        self.t = []
        this_t = self.tstart.datetime + timedelta(hours=self.timezone)
        for i in range(len(self.data[:, 0])):
            self.t.append(this_t)
            this_t += dt

    def _parse_fn(self):
        pfn = self.filename.split('_')
        if len(pfn) == 3:
            _date = f"20{pfn[1][:2]}-{pfn[1][2:4]}-{pfn[1][4:6]}"
            _x = pfn[2].split('.')[0]
            _time = f"{_x[:2]}:{_x[2:4]}"
            if len(_x) == 6:
                _time += f":{_x[4:6]}"
            _datetime = f"{_date}T{_time}"
        else:
            _datetime = None
        return pfn[0], _datetime

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
        defaults = {'log': True, 'dB': False, 'colorbar': True, 'xticks': 10, 'yticks': 4}
        kwargs = proc_plot_kwargs(kwargs, defaults)

        num_xticks = kwargs['xticks']
        num_yticks = kwargs['yticks']

        if kwargs['log']:
            plt.imshow(kwargs['mult'] * np.log10(self.data))
        else:
            plt.imshow(self.data)
        if kwargs['colorbar']:
            plt.colorbar()

        plt.xticks(np.linspace(0, len(self.data[0]), num_xticks), [f"{x:.2f}" for x in np.linspace(self.fmin, self.fmax, num_xticks)])
        jds = np.linspace(self.jdstart, self.jdstop, num_yticks)
        apt = Time(jds, format='jd')
        yticks = [(x + timedelta(hours=self.timezone)).strftime("%H:%M:%S") for x in apt.datetime]
        plt.yticks(np.linspace(0, len(self.data), num_yticks), yticks)
        plt.xlabel('MHz')
        plt.ylabel(f"UTC{'+' if self.timezone>0 else '-'}{abs(self.timezone):.0f}")
        plt.title(self.tstart.datetime.strftime('%Y-%m-%d'))

    def spectra(self, **kwargs):
        """Make a 2-D plot of the spectra."""
        defaults = {'log': False, 'dB': False}
        kwargs = proc_plot_kwargs(kwargs, defaults)

        for data in self.data:
            if kwargs['log']:
                plt.plot(self.freq, kwargs['mult'] * np.log10(data))
            else:
                plt.plot(self.freq, data)
        plt.grid()
        plt.xlabel('MHz')
        plt.ylabel(kwargs['ylabel'])

    def series(self, **kwargs):
        defaults = {'log': False, 'dB': False}
        kwargs = proc_plot_kwargs(kwargs, defaults)

        self.plot_max = 0.0
        self.plot_min = 1E20

        for i in range(len(self.data[:,0])):
            if kwargs['log']:
                data = np.log10(self.data[:, i])
            else:
                data = self.data[:, i]
            plt.plot(self.t, data)
            this_max = np.max(data)
            this_min = np.min(data)
            if this_max > self.plot_max:
                self.plot_max = copy(this_max)
            if this_min < self.plot_min:
                self.plot_min = copy(this_min)
        if self.datetime is not None:
            self.expected(self.datetime, kwargs['tz'])
        if kwargs['freq'] is not None:
            import beamfit
            ifg = int((float(kwargs['freq']) - self.freq[0]) / (self.freq[1] - self.freq[0]))
            plt.plot(self.t, self.data[:, ifg], lw=4, color='k')
            coeff, data_fit = beamfit.fit_it(self.data[:, ifg], max(self.data[:, ifg]), len(self.data[:, ifg]) / 2.0, len(self.data[:, ifg])/4.0 )
            plt.plot(self.t, data_fit, '--', color='w')
            ctrt = int(coeff[1])
            rngt = [int(coeff[1] - coeff[2]), int(coeff[1] + coeff[2])]
            print(f"Found: {self.t[ctrt]}")
            print(f"Width:  {self.t[rngt[1]] - self.t[rngt[0]]}")
        plt.grid()
        plt.xlabel(f"UTC{'+' if self.timezone>0 else '-'}{abs(self.timezone):.0f}")
        plt.ylabel(kwargs['ylabel'])

    def expected(self, ts, tz_offset=0.0):
        ts = Time(ts).datetime
        plt.plot([ts, ts], [self.plot_min, self.plot_max], color='k', lw=3)

if __name__ == '__main__':
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument('fn', help="Name of hdf5 datafile")
    ap.add_argument('-p', '--plot_type', help='wf, series, spectra [wf]', choices=['wf', 'series', 'spectra'], default='wf')
    ap.add_argument('-x', '--xticks', help="Number of xticks in waterfall [10]", type=int, default=10)
    ap.add_argument('-y', '--yticks', help="Number of yticks to use in waterfall [4]", type=int, default=4)
    ap.add_argument('-c', '--colorbar', help="Flag to hide colorbar", action='store_false')
    ap.add_argument('-f', '--freq', help='Frequency to use for fitting beam crossing', default=None)
    ap.add_argument('-l', '--log', help="Flag to take log10 of data", action='store_true')
    ap.add_argument('-d', '--dB', help="Flag to convert to dB", action='store_true')
    ap.add_argument('--tz', help="Timezone offset to UTC in hours [-8.0]", type=float, default=-8.0)
    args = ap.parse_args()
    obs = Data(args.fn, args.tz)
    getattr(obs, args.plot_type)(**vars(args))
    plt.show()
