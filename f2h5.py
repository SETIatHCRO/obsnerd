#! /usr/bin/env python
import numpy as np
import h5py
import yaml
from astropy.time import Time
import metadata
from datetime import datetime, timedelta


"""
Convert the gnuradio-companion spectrum file to hdf5
"""

RAW_FILENAME = 'nrdz'

TIME_FORMATS = ['%Y-%m-%dT%H:%M:%S', '%y-%m-%dT%H:%M:%S',
                '%Y-%m-%d %H:%M:%S', '%y-%m-%d %H:%M:%S',
                '%Y/%m/%dT%H:%M:%S', '%y/%m/%dT%H:%M:%S',
                '%Y/%m/%d %H:%M:%S', '%y/%m/%d %H:%M:%S',
                '%d/%m/%YT%H:%M:%S', '%d/%m/%yT%H:%M:%S',
                '%d/%m/%Y %H:%M:%S', '%d/%m/%y %H:%M:%S',
                '%Y%m%dT%H%M%S', '%y%m%dT%H%M%S',
                '%Y%m%d %H%M%S', '%y%m%d %H%M%S',
                '%Y%m%d_%H%M%S', '%y%m%d_%H%M%S',
                '%Y%m%d%H%M%S', '%y%m%d%H%M%S'
                ]

def make_filename(**kwargs):
    if kwargs['date'] is None:
        return kwargs['tag']
    if 'timezone' not in kwargs:
        kwargs['timezone'] = 0.0

    for this_tf in TIME_FORMATS:
        try:
            this_dt = datetime.strptime(kwargs['date'], this_tf) + timedelta(hours=kwargs['timezone'])
            break
        except ValueError:
            continue
    return f"{kwargs['tag']}_{this_dt.strftime('%y%m%d_%H%M%S')}.h5"

class HDF5HeaderInfo:
    def __init__(self):
        self.data = 'data'
        self.float64s = ['tstart', 'tstop', 'fcen', 'bw', 'decimation', 'nfft', 'tle']
        self.from_datetime = ['tstart', 'tstop', 'tle']


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
    ap.add_argument('tag', help="Name of output file or prefix for date.")
    ap.add_argument('-d', '--date', help="Optional date to generate the filename with tag as prefix", default=None)
    ap.add_argument('-z', '--timezone', help="Timezone to add (-8 is PST).", type=float, default=0.0)
    args = ap.parse_args()

    convert(RAW_FILENAME, make_filename(**vars(args)))
