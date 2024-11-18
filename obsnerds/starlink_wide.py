from . import starlink_look
import numpy as np

class WideBand:
    def __init__(self, source, lo, cnode=[352, 544, 736, 928, 1120, 1312, 1504]):
        self.filelist = []
        for cn in cnode:
            self.filelist.append(f"{source}_{lo}_C{cn:04d}.npz")
        self.lookall = starlink_look.Look()  # This will use the methods to plot for concatenated
        self.lookall.source = source
        self.lookall.lo = lo
        self.looks = {}
        for fil in self.filelist:
            key = fil.split('.')[0]
            self.looks[key] = starlink_look.Look()
            found_it = self.looks[key].read_npz(fil)
            if not found_it:
                del(self.looks[key])
                continue
            self.looks[key].source = source
            self.looks[key].lo = lo
        ford = {}
        for key, lk in self.looks.items():
            ford[str(lk.freqs[0])] = key
        sford = sorted(ford.keys())
        self.fileorder = []
        for key in sford:
            self.fileorder.append(ford[key])
        self.lookall.times = self.looks[self.fileorder[0]].times
        for key in self.fileorder:
            print(f"{key}:  {self.looks[key].freqs[0]} - {self.looks[key].freqs[-1]}")

    def concat(self, ant, pol='xx'):
        self.ant = ant
        self.pol = pol
        self.lookall.freqs = []
        dataf = []
        for key in self.fileorder:
            self.lookall.freqs += self.looks[key].freqs
            self.looks[key].get_bl(a=ant, pol=pol)
            self.lookall.a = ant
            self.lookall.b = ant
            self.lookall.pol = pol
            dataf.append(self.looks[key].data)
        self.lookall.data = np.concatenate(dataf, axis=1)

    def dashboard(self, use_db=True, save=False, time_axis='diff', feph=False, show_feph='False'):
        """
        Shortcut to the look dashboard with False antenna since self.lookall.data made in concat

        """
        try:
            self.lookall.dashboard(ant=False, pol=self.pol, use_db=use_db, save=save, time_axis=time_axis, feph=feph, show_feph=show_feph)
        except AttributeError:
            print("Run concat first you dummy!")