#! /usr/bin/python3

import numpy as np
import matplotlib.pyplot as plt
import h5py
from astropy.time import Time
import numpy as np


class Data:
    def __init__(self, fn):
        self.filename = fn
        with h5py.File(fn, 'r') as fp:
            self.data = np.array(fp['data'])
            self.jdstart = np.float64(fp['tstart'])  # jd
            self.jdstop = np.float64(fp['tstop'])  # jd
            self.fcen = np.float64(fp['fcen'])  # MHz
            self.bw = np.float64(fp['bw'])  # MHz
        self.tstart = Time(self.jdstart, format='jd')
        self.tstop = Time(self.jdstop, format='jd')
        self.fmin = self.fcen - self.bw / 2.0
        self.fmax = self.fcen + self.bw / 2.0

    def header(self):
        """Print information about the data."""

        print(f"Filename: {self.filename}")
        print(f"Data shape: {np.shape(self.data)}")
        print(f"Start: {self.tstart.datetime} (jd={self.jdstart})")
        print(f"Stop: {self.tstop.datetime} (jd={self.jdstop})")
        print(f"Freq:  {self.fmin:.2f} - {self.fmax:.2f} MHz  (cf = {self.fcen}, bw = {self.bw} MHz)")

    def wf(self, **kwargs):
        """
        Make a waterfall plot of the data.

        Parameters:
        -----------
        kwargs:
            num_xticks (int), num_yticks (int), colorbar (bool), log (bool)

        """
        if 'log' in kwargs and not kwargs['log']:
            plt.imshow(self.data)
        else:  # Default to using log
            plt.imshow(np.log10(self.data))
        if 'colorbar' in kwargs and not kwargs['colorbar']:
            pass
        else:
            plt.colorbar()  # Default to showing colorbar
        if 'xticks' in kwargs:
            num_xticks = int(kwargs['xticks'])
        else:  # Default to 10
            num_xticks = 10
        if 'yticks' in kwargs:
            num_yticks = int(kwargs['yticks'])
        else:  # Default to 4
            num_yticks = 4

        plt.xticks(np.linspace(0, len(self.data[0]), num_xticks), [f"{x:.2f}" for x in np.linspace(self.fmin, self.fmax, num_xticks)])
        jds = np.linspace(self.jdstart, self.jdstop, num_yticks)
        apt = Time(jds, format='jd')
        yticks = [x.strftime("%H:%M:%S") for x in apt.datetime]
        plt.yticks(np.linspace(0, len(self.data), num_yticks), yticks)
        plt.xlabel('MHz')
        plt.title(self.tstart.datetime.strftime('%Y-%m-%d'))

    def plot(self, **kwargs):
        """Make a 2-D plot of the spectra."""

        if 'log' not in kwargs:
            kwargs['log'] = False

        for data in self.data:
            if kwargs['log']:
                plt.plot(np.log10(data))
            else:
                plt.plot(data)


if __name__ == '__main__':
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument('fn', help="Name of hdf5 datafile")
    ap.add_argument('-x', '--xticks', help="Number of xticks in waterfall [10]", format=int, default=10)
    ap.add_argument('-y', '--yticks', help="Number of yticks to use in waterfall [4]", format=int, default=4)
    ap.add_argument('-c', '--colorbar', help="Flag to use colorbar [True]", format=bool, default=True)
    ap.add_argument('-l', '--log', help="Flag to take log10 of data [True]", format=bool, default=True)
    args = ap.parse_args()
    obs = Data(args.fn)
    obs.wf(**vars(args))
    plt.show()
