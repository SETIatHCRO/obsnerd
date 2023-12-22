#! /usr/bin/env python3
import numpy as np
import h5py
import yaml
from astropy.time import Time


"""
Convert the gnuradio-companion spectrum file to hdf5
"""

def convert(filename, output_file=None, split=4096):
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
    with open('metadata.yaml', 'r') as fp:
        meta = yaml.safe_load(fp)

    with h5py.File(output_file, 'w') as fp:
        dset = fp.create_dataset('data', data=data)
        for key, val in meta.items():
            if key in ['tstart', 'tstop']:
                val = Time(val, format='datetime').jd
            print(key, val)
            # mset = fp.create_dataset(key, shape=(), dtype=float, data=val)
            mset = fp.create_dataset(key, data=val)


if __name__ == '__main__':
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument('fn', help='Name of spectrum file')
    ap.add_argument('-o', '--output_file', help="Name of output file", default=None)
    args = ap.parse_args()
    convert(args.fn, args.output_file)
