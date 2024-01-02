#! /usr/bin/env python
import numpy as np
import h5py
import yaml
from astropy.time import Time
import metadata


"""
Convert the gnuradio-companion spectrum file to hdf5
"""

RAW_FILENAME = 'nrdz'

class HDF5HeaderInfo:
    def __init__(self):
        self.data = 'data'
        self.float64s = ['tstart', 'tstop', 'fcen', 'bw', 'decimation', 'nfft', 'tle']
        self.from_datetime = ['tstart', 'tstop', 'tle']


def convert(filename, output_file=None, split=4096):
    h5 = HDF5HeaderInfo()
    data = np.fromfile(filename, dtype=float)
    sdata = []
    for i in range(len(data) // split):
        sdata.append(list(data[i*split: (i+1)*split]))
    data = np.array(sdata)

    if output_file is None:
        fsplit = filename.split('.')
        if len(fsplit) == 1:
            output_file = f"{filename}.h5"
        else:
            output_file = f"{'.'.join(fsplit[:-1])}.h5"
    else:
        if not output_file.endswith('h5'):
            output_file = f"{output_file}.h5"
    with open('metadata.yaml', 'r') as fp:
        meta = yaml.safe_load(fp)

    print(f"Writing file {output_file}")
    with h5py.File(output_file, 'w') as fp:
        dset = fp.create_dataset(h5.data, data=data)
        for key, val in meta.items():
            if key in h5.from_datetime:
                val = Time(val, format='datetime').jd
            print(key, val)
            # mset = fp.create_dataset(key, shape=(), dtype=float, data=val)
            mset = fp.create_dataset(key, data=val)
    metadata.onlog(f"Writing {output_file}")

if __name__ == '__main__':
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument('output_file', help="Name of output file")
    args = ap.parse_args()
    convert(RAW_FILENAME, args.output_file)
