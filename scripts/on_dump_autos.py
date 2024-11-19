from pyuvdata import UVData
from astropy.time import Time
import numpy as np
from copy import copy

FREQ_CONVERT = {'MHz': 1E6, 'GHz': 1E9}


class Dump:
    def __init__(self, fn, freq_unit='MHz'):
        self.fnuvh5 = fn
        self.freq_unit = freq_unit
        self.read_uvh5()

    def read_uvh5(self):
        self.uv = UVData()
        self.uv.read(self.fnuvh5)
        self.ant_numbers = self.uv.get_ants()
        self.ant_names = np.array(self.uv.antenna_names)[np.unique(self.uv.ant_1_array) - 1]
        self.ant_map = {}
        for antno, antna in zip(self.ant_numbers, self.ant_names):
            self.ant_map[antna] = antno
        self.freqs = self.uv.freq_array[0] / FREQ_CONVERT[self.freq_unit]

    def dump_autos(self, fnout, ants=None, pols=['xx', 'yy', 'xy', 'yx']):
        self.source = fnout  # For better or worse, these are the same
        if ants is None:
            ants = self.ant_names
            antstr = 'all'
        else:
            antstr = ','.join(ants)
        fn = f"{self.source}.npz"
        outdata = {'ants': ants, 'freqs': self.freqs, 'pols': pols, 'source': self.source, 'uvh5': self.fnuvh5, 'freq_unit': self.freq_unit}
        print(f"Dumping autos in {self.fnuvh5} for {antstr} to {fn} for {pols}")
        for ant in self.ant_names:
            for pol in pols:
                self.get_bl(ant, pol=pol)
                outdata[f"{ant}{pol}"] = copy(self.data)
        outdata['times'] = self.times.jd
        np.savez(fn, **outdata)

    def get_bl(self, a, b=None, pol='xx'):
        self.a = a
        self.b = b
        self.pol = pol
        if self.b is None:
            self.b = self.a
        print(f"Reading ({self.a}, {self.b}){self.pol}", end='')
        self.ano = self.ant_map[self.a]
        self.bno = self.ant_map[self.b]
        self.data = self.uv.get_data(self.ano, self.bno, pol)
        self.times = Time(self.uv.get_times(self.ano, self.bno), format='jd')
        self.datamin = np.min(np.abs(self.data))
        self.datamax = np.max(np.abs(self.data))
        print(f"\tmin={self.datamin}, max={self.datamax}")


if __name__ == '__main__':
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument('filename', help="Name of file to dump.")
    ap.add_argument('-o', '--output', help="Name of npz output file (default is input)", default=None)
    args = ap.parse_args()

    if args.output is None:
        args.output = '.'.join(args.filename.split('.')[:-1]) + '.npz'
    sl = Dump(args.filename)
    sl.dump_autos(args.output)
