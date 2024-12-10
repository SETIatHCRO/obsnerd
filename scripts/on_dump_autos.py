from pyuvdata import UVData
from astropy.time import Time
import numpy as np
from copy import copy

FREQ_CONVERT = {'MHz': 1E6, 'GHz': 1E9}
UVH5_SRC_IND = 4


class Dump:
    def __init__(self, fn, lo, cnode, freq_unit='MHz'):
        self.fnuvh5 = fn
        self.lo = lo
        self.cnode = cnode
        self.source = fn.split('_')[UVH5_SRC_IND]  # Current standard format
        if self.source[0] == 'S':
            print("AD HOC TO HANDLE THE EXTRA UNDERSCORE")
            self.source += '_' + fn.split('_')[UVH5_SRC_IND+1]
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

    def dump_autos(self, ants=None, pols=['xx', 'yy', 'xy', 'yx']):
        if ants is None:
            ants = self.ant_names
            antstr = 'all'
        else:
            antstr = ','.join(ants)
        outdata = {'ants': ants, 'freqs': self.freqs, 'pols': pols, 'source': self.source, 'uvh5': self.fnuvh5, 'freq_unit': self.freq_unit}
        print(f"Dumping autos in {self.fnuvh5} for {antstr} {pols}", end=' ... ')
        for ant in self.ant_names:
            for pol in pols:
                self.get_bl(ant, pol=pol)
                outdata[f"{ant}{pol}"] = copy(self.data)
        outdata['times'] = self.times.jd  # This assumes that all times in the UVH5 file are the same...
        mjd = outdata['times'][0] - 2400000.5
        obsrec = f"{self.source}_{mjd:.4f}_{self.lo}_{self.cnode}.npz"
        print(f"writing {obsrec}")
        np.savez(obsrec, **outdata)

    def get_bl(self, a, b=None, pol='xx'):
        """
        Pulled from starlink_eph to be stand alone

        """
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
    ap.add_argument('filename', help="Name of UVH5 file to dump.")
    ap.add_argument('--lo', help='LO', choices=['A', 'B'], default='A')
    ap.add_argument('--cnode', help='CNODE', choices=['C0352', 'C0544', 'C0736', 'C0928',  'C1120',  'C1312',  'C1504', 'C0352',
                                                      'C0544', 'C0736', 'C0928', 'C1120', 'C1312', 'C1504'], default='C0928')
    args = ap.parse_args()

    sl = Dump(args.filename, args.lo, args.cnode)
    sl.dump_autos()
