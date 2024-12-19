import numpy as np
import h5py
from astropy.time import Time
from . import metadata, onutil


"""
Convert the gnuradio-companion spectrum file to hdf5
"""


def make_filename(**kwargs):
    if kwargs['tag'] == 'meta':
        meta = metadata.get_meta()
        kwargs['tag'] = meta['source']
        kwargs['date'] = meta['expected']
    if kwargs['date'] is None:
        return kwargs['tag']
    if 'timezone' not in kwargs:
        kwargs['timezone'] = 0.0
    this_dt = onutil.make_datetime(date=kwargs['date'], timezone=kwargs['timezone'])

    return f"{kwargs['tag']}_{this_dt.strftime('%y%m%d_%H%M%S')}.h5"


class HDF5HeaderInfo:
    def __init__(self):
        self.data = 'data'
        self.float64s = ['tstart', 'tstop', 'fcen', 'bw', 'decimation', 'nfft', 'tle', 'expected']
        self.strings = ['source', 'move', 'move_data', 'track']
        self.from_yaml = ['track']
        self.from_datetime = ['tstart', 'tstop', 'tle', 'expected']
        self.metadata = self.float64s + self.strings


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
        _ = fp.create_dataset(h5.data, data=data)
        for key, val in meta.items():
            if key in h5.from_datetime:
                val = Time(val, format='datetime').jd
            # mset = fp.create_dataset(key, shape=(), dtype=float, data=val)
            # mset = fp.create_dataset(key, data=val)
            print(f"\t{key}:  {val}")
            fp[key] = val

    metadata.logger.info(f"Writing {output_file}")

