from pyuvdata import UVData
from astropy.time import Time
import numpy as np
from copy import copy

FREQ_CONVERT = {'MHz': 1E6, 'GHz': 1E9}


def toMag(x, use_db=True):
    if use_db:
        return 10.0 * np.log10(np.abs(x))
    else:
        return np.abs(x)


class Dump:
    def __init__(self, fn):
        self.fn = fn
        print(f"Dumping autos for {self.fn}")
        self.read_uvh5()
        self.dump_autos(ants=None, pols=['xx', 'yy', 'xy', 'yx'])

    def read_uvh5(self):
        self.source = self.fn.split('_')[4]
        self.uv = UVData()
        self.uv.read(self.fn)
        self.ant_numbers = self.uv.get_ants()
        self.ant_names = np.array(self.uv.antenna_names)[np.unique(self.uv.ant_1_array) - 1]
        self.ant_map = {}
        for antno, antna in zip(self.ant_numbers, self.ant_names):
            self.ant_map[antna] = antno
        self.freqs = self.uv.freq_array[0] / FREQ_CONVERT[self.freq_unit]

    def dump_autos(self, ants=None, pols=['xx', 'yy', 'xy', 'yx']):
        if ants is None:
            ants = self.ant_names
            antstr = 'all'
        else:
            antstr = ','.join(ants)
        fn = f"{self.source}.npz"
        outdata = {'ants': ants, 'freqs': self.freqs, 'pols': pols, 'source': self.source, 'uvh5': self.fn, 'freq_unit': self.freq_unit}
        print(f"Dumping autos in {self.fn} for {antstr} to {fn} for {pols}")
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
        if self.file_type == 'uvh5':
            self.ano = self.ant_map[self.a]
            self.bno = self.ant_map[self.b]
            self.data = self.uv.get_data(self.ano, self.bno, pol)
            self.times = Time(self.uv.get_times(self.ano, self.bno), format='jd')
        elif self.file_type == 'npz':
            self.data = self.npzfile[f"{self.a}{pol}"]
        self.datamin = np.min(np.abs(self.data))
        self.datamax = np.max(np.abs(self.data))
        print(f"\tmin={self.datamin}, max={self.datamax}")


if __name__ == '__main__':
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument('filename', help="Name of file to dump.")

    args = ap.parse_args()
    sl = Dump(args.filename)
