from datetime import datetime
from numpy import ones


class TBALog:
    inner = 3.0
    outer = 5.0
    def __init__(self, logfn):
        self.log_filename = logfn
        self.times = {'inner': [], 'outer': []}
        self.sats = {'inner': [], 'outer': []}
        with open(logfn, 'r') as fp:
            for line in fp:
                data = line.split(',')
                self.times[data[2]].append(datetime.strptime(data[0], '%Y-%m-%d %H:%M:%S'))
                self.sats[data[2]].append(data[4])

    def plot(self):
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots()
        ipts = 0.0 * ones(len(self.times['inner']))
        opts = 0.0 * ones(len(self.times['outer']))
        ax.plot(self.times['inner'], ipts, 'r|', label='Inner')
        ax.plot(self.times['outer'], opts, 'b|', label='Outer')
        for i in range(len(self.times['inner'])):
            ax.annotate(self.sats['inner'][i], (self.times['inner'][i], ipts[i]), rotation='vertical', fontsize=8)
        for i in range(len(self.times['outer'])):
            ax.annotate(self.sats['outer'][i], (self.times['outer'][i], opts[i]), rotation='vertical', fontsize=8)
        plt.show()