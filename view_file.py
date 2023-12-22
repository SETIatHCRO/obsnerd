import numpy as np
import matplotlib.pyplot as plt
import h5py
from astropy.time import Time


class Data:
    def __init__(self, fn, datakey='data'):
        with h5py.File(fn, 'r') as fp:
            self.data = np.array(fp[datakey])
            self.jdstart = np.float64(fp['tstart'])  # jd
            self.jdstop = np.float64(fp['tstop'])  # jd
            self.fcen = np.float64(fp['fcen'])  # MHz
            self.bw = np.float64(fp['bw'])  # MHz
        self.tstart = Time(self.jdstart, format='jd')
        self.tsop = Time(self.jdstop, format='jd')
        self.fmin = self.fcen - self.bw / 2.0
        self.fmax = self.fcen + self.bw / 2.0

    def wf(self, num_xticks=10, num_yticks=4, colorbar=False):
        plt.imshow(np.log10(self.data))
        if colorbar:
            plt.colorbar()
        plt.xticks(np.linspace(0, len(self.data[0]), num_xticks), [f"{x:.2f}" for x in np.linspace(self.fmin, self.fmax, num_xticks)])
        jds = np.linspace(self.jdstart, self.jdstop, num_yticks)
        apt = Time(jds, format='jd')
        yticks = [x.strftime("%H:%M:%S") for x in apt.datetime]
        plt.yticks(np.linspace(0, len(self.data), num_yticks), yticks)
        plt.xlabel('MHz')
        plt.title(self.tstart.datetime.strftime('%Y-%m-%d'))

    def plot(self):
        for data in self.data:
            plt.plot(data)


if __name__ == '__main__':
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument('fn', help="Name of hdf5 datafile")
    args = ap.parse_args()
    obs = Data(args.fn)
    obs.wf()
    plt.show()

### IGNORE BELOW
# ATA observing...
def vfDEPR(fn, key='data'):
    if fn.endswith('h5') or fn.endswith('hdf5'):
        ftype = 'real'
    else:
        ftype = 'complex'
    if ftype == 'complex':
        data = np.fromfile(fn, dtype=np.complex64, count=int(1024*1024))
        plt.specgram(data, NFFT=8192, Fs=1)
        plt.show()
    else:
        with h5py.File(fn, 'r') as fp:
            data = np.array(fp[key])
            for i in range(len(data)):
                plt.plot(data[i])
        plt.show()
    return data