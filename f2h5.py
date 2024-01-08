#! /usr/bin/env python
import numpy as np
import h5py
from astropy.time import Time
import metadata
import onutil


"""
Convert the gnuradio-companion spectrum file to hdf5
"""

RAW_FILENAME = 'nrdz'


def make_filename(**kwargs):
    if kwargs['tag'] == 'meta':
        meta = metadata.get_meta()
        kwargs['tag'] = meta['source']
        kwargs['date'] = meta['expected']
    if kwargs['date'] is None:
        return kwargs['tag']
    if 'timezone' not in kwargs:
        kwargs['timezone'] = 0.0
    this_dt = onutil.make_datetime(**kwargs)
    print("F2H5-26:")
    print(this_dt, type(this_dt))

    return f"{kwargs['tag']}_{this_dt.strftime('%y%m%d_%H%M%S')}.h5"


class HDF5HeaderInfo:
    def __init__(self):
        self.data = 'data'
        self.float64s = ['tstart', 'tstop', 'fcen', 'bw', 'decimation', 'nfft', 'tle', 'expected']
        self.strings = ['source', 'move', 'move_data']
        self.from_datetime = ['tstart', 'tstop', 'tle', 'expected']


def convert(input_file, output_file=None, split=4096):
    h5 = HDF5HeaderInfo()
    data = np.fromfile(input_file, dtype=float)
    sdata = []
    for i in range(len(data) // split):
        sdata.append(list(data[i*split: (i+1)*split]))
    data = np.array(sdata)

    if output_file is None:
        fsplit = input_file.split('.')
        if len(fsplit) == 1:
            output_file = f"{input_file}.h5"
        else:
            output_file = f"{'.'.join(fsplit[:-1])}.h5"
    else:
        if not output_file.endswith('h5'):
            output_file = f"{output_file}.h5"
    meta = metadata.get_meta()

    print(f"Writing file {output_file}")
    with h5py.File(output_file, 'w') as fp:
        dset = fp.create_dataset(h5.data, data=data)
        for key, val in meta.items():
            print(f"F2H5-63 -- {key}: {val}   {type(val)}")
            if key in h5.from_datetime:
                val = Time(val, format='datetime').jd
            # if key in h5.from_str:
            fp[key] = val
            print(key, val)
            # mset = fp.create_dataset(key, shape=(), dtype=float, data=val)
            # mset = fp.create_dataset(key, data=val)
    metadata.onlog(f"Writing {output_file}")

if __name__ == '__main__':
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument('tag', help="Name of output file or prefix for date [meta]", nargs='?', default='meta')
    ap.add_argument('-d', '--date', help="Optional date to generate the filename with tag as prefix", default=None)
    ap.add_argument('-z', '--timezone', help="Timezone to add (-8 is PST).", type=float, default=0.0)
    args = ap.parse_args()

    convert(RAW_FILENAME, make_filename(**vars(args)))
