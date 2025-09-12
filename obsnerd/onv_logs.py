from datetime import datetime, timedelta
from numpy import ones
from glob import glob
from astropy.time import Time


class TBALog:
    inner = 3.0
    outer = 5.0
    def __init__(self, single_file=None):
        if single_file is not None:
            self.log_files = [single_file]
        else:
            self.log_files = sorted(glob("ATA*.txt") + glob("ATA*.csv"))
        self.times = {'inner': [], 'outer': []}
        self.sats = {'inner': [], 'outer': []}
        for logfn in self.log_files:
            with open(logfn, 'r') as fp:
                for line in fp:
                    data = line.split(',')
                    self.times[data[2]].append(datetime.strptime(data[0], '%Y-%m-%d %H:%M:%S'))
                    self.sats[data[2]].append(data[4].strip())
        for scope in ['inner', 'outer']:
            self.times[scope] = Time(self.times[scope])

    def filter_time(self, central, delta=4):
        """
        Filter the log files to only include those that are within the specified time range.

        Parameters
        ----------
        central : datetime
            Central time to filter around.
        delta : float
            Time range in minutes to filter around the central time.

        """
        dt = timedelta(minutes=delta)
        self.filt_times = {'inner': [], 'outer': []}
        self.filt_sats = {'inner': [], 'outer': []}
        for scope in ['inner', 'outer']:
            for time, sat in zip(self.times[scope], self.sats[scope]):
                if central - dt < time < central + dt:
                    self.filt_times[scope].append(time)
                    self.filt_sats[scope].append(sat)
            if len(self.filt_times[scope]):
                self.filt_times[scope] = Time(self.filt_times[scope])

    def plot(self):
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots()
        ipts = 0.0 * ones(len(self.times['inner']))
        opts = 0.0 * ones(len(self.times['outer']))
        ax.plot(self.times['inner'].datetime, ipts, 'r|', label='Inner')
        ax.plot(self.times['outer'].datetime, opts, 'b|', label='Outer')
        for i in range(len(self.times['inner'])):
            ax.annotate(self.sats['inner'][i], (self.times['inner'].datetime[i], ipts[i]), rotation='vertical', fontsize=8)
        for i in range(len(self.times['outer'])):
            ax.annotate(self.sats['outer'][i], (self.times['outer'].datetime[i], opts[i]), rotation='vertical', fontsize=8)
        plt.show()

    def ods_mon(self, fn='online_ods_mon.csv'):
        from odsutils import ods_engine
        self.ods = ods_engine.ODS(version='latest', output='INFO')
        self.ods.add(fn)