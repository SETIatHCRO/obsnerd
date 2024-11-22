from . import starlink_look
import numpy as np

class Obs:
    def __init__(self, obsid, lo='A', cnode=[352, 544, 736, 928, 1120, 1312, 1504], tag='npz'):
        self.obsrec_list = []
        self.lookall = starlink_look.Look()  # This will use the methods to plot for concatenated
        self.lookall.obsid = obsid
        self.lookall.lo = lo
        if isinstance(cnode[0], str) and cnode[0][0] == 'C':
            cnodes = cnode
        else:
            cnodes = [f"C{int(x):04d}" for x in cnode]
        self.obsrec_list = [f"{self.lookall.obsid}_{self.lookall.lo}_{x}" for x in cnodes]

        self.looks = {}
        for obsrec, cno in zip(self.obsrec_list, cnodes):
            self.looks[obsrec] = starlink_look.Look()
            found_it = self.looks[obsrec].read_npz(f"{obsrec}.{tag}")
            if not found_it:
                del(self.looks[obsrec])
                continue
            self.looks[obsrec].obsid = obsid
            self.looks[obsrec].obsrec = obsrec
            self.looks[obsrec].lo = lo
            self.looks[obsrec].cnode = cno
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

    def dashboard(self, use_db=True, save=False, time_axis='diff', show_feph=False):
        """
        Shortcut to the look dashboard with False antenna since self.lookall.data made in concat

        """
        try:
            self.lookall.dashboard(ant=False, pol=self.pol, use_db=use_db, save=save, time_axis=time_axis, show_feph=show_feph)
        except AttributeError:
            print("Run concat first you dummy!")