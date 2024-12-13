from pyuvdata import UVData
from astropy.time import Time
import numpy as np
import matplotlib.pyplot as plt
from obsnerds import obs_base
from datetime import datetime
import os.path as path
from copy import copy

FREQ_CONVERT = {'MHz': 1E6, 'GHz': 1E9}


def toMag(x, use_db=True):
    if use_db:
        return 10.0 * np.log10(np.abs(x))
    else:
        return np.abs(x)

def make_cnode(cns):
        """
        Takes a list and makes a list of valid cnodes -> C####

        """
        cns = cns if isinstance(cns, list) else cns.split(',')
        try:
            return [f"C{int(x):04d}" for x in cns]
        except ValueError:
            return cns

def gen_dump_script(date_path, base_path='/mnt/primary/ata/projects/p054/', script_filename='dump_autos.sh', ants='all'):
    from os import walk, listdir, path
    if date_path == '?':
        print(f"Available observation dates in {base_path}:")
        for x in listdir(base_path):
            print(f"\t{x}")
        return
    print(f"Retrieving from {base_path}")
    files = {}
    dbase_path = path.join(base_path, date_path)
    for basedir, _, filelist in walk(dbase_path):
        if base_path in basedir and '/Lo' in basedir:
            lolo, cnode = basedir.split('/')[-1].split('.')
            lo = lolo[2:]
            for fn in filelist:
                if fn.startswith('uvh5_'):
                    fnsplit = fn.split('_')
                    if len(fnsplit) == 6:
                        _, mjd1, mjd2, _, src, _ = fnsplit
                    elif len(fnsplit) == 7:
                        _, mjd1, mjd2, _, src, extra, _ = fnsplit
                        src = src + '_' + extra
                    mjd = float(f"{mjd1}.{mjd2}")
                    obsid = f"{src}_{mjd:.4f}"
                    obsrec = f"{obsid}_{lo}_{cnode}"
                    dfn = path.join(basedir, fn)
                    files[obsrec] = [dfn, lo, cnode]
                
    with open(script_filename, 'w') as fp:
        for obsrec, data in files.items():
            print(f"on_dump_autos.py {data[0]} --lo {data[1]} --cnode {data[2]} --ants {ants}", file=fp)
            print(f"Adding {obsrec}")


ALL_CNODES = ['C0352', 'C0544', 'C0736', 'C0928', 'C1120', 'C1312', 'C1504']


class Look:
    def __init__(self, obsid, lo='A', cnode=ALL_CNODES, tag='npz', freq_unit='MHz', dir_data='.'):
        """
        This initializes, but also reads in all of the data.

        """
        self.obsid = obsid
        self.lo = lo
        self.cnode = make_cnode(cnode)
        self.tag = tag
        self.freq_unit = freq_unit
        self.dir_data = dir_data
        self.obsrec_list = [f"{self.obsid}_{self.lo}_{x}" for x in self.cnode]
        self.get_obsinfo()
        self.npzfile = {}
        self.freqs = []
        for i, obsrec in enumerate(self.obsrec_list):
            if self.tag == 'uvh5':
                if i:
                    print("This will only read in the first obsrec file in the list")
                else:
                    self.read_uvh5(obsrec)
            elif self.tag == 'npz':
                self.read_an_npz(obsrec)
            else:
                print(f"Invalid file {obsrec}.{self.tag}")

    def read_a_uvh5(self, fn):
        print(f"Reading {fn}")
        print("NEED TO MAKE fn AND GET DATA PATH ETC IN HERE -- CURRENTLY NOT LIKELY GOING TO WORK?")
        self.fn = path.join(self.dir_data, fn)
        self.uv = UVData()
        self.uv.read(self.fn)
        self.ant_numbers = self.uv.get_ants()
        self.ant_names = np.array(self.uv.antenna_names)[np.unique(self.uv.ant_1_array) - 1]
        self.ant_map = {}
        for antno, antna in zip(self.ant_numbers, self.ant_names):
            self.ant_map[antna] = antno
        self.freqs = self.uv.freq_array[0] / FREQ_CONVERT[self.freq_unit]

    def read_an_npz(self, obsrec):
        """
        Reads a single obsrec file.
        This will write ant_names, times, into attributes without checking...and appends freqs

        Parameter
        ---------
        obsrec : str
            Obsrec designation for an observation file (without the extension)

        Attributes
        ----------
        npzfile : dictionary
            Contains the npzfile contents, keyed on obsrec
        ant_names : list of str
            List of antenna names - gets overwritten every time.
        freqs : list of float
            List of frequencies - gets appended
        times : list of astropy times
            List of observation time stamps - gets overwritten very time

        """
        fn = path.join(self.obs.obsinfo.dir_data, f"{obsrec}.{self.tag}")
        try:
            self.npzfile[obsrec] = np.load(fn)
        except FileNotFoundError:
            print(f"Couldn't find {fn}")
            return False
        self.ant_names = list(self.npzfile[obsrec]['ants'])
        self.freqs += list(self.npzfile[obsrec]['freqs'])
        self.times = Time(self.npzfile[obsrec]['times'], format='jd')
        return True

    def dump_autos(self, ants=None, pols=['xx', 'yy', 'xy', 'yx']):
        if ants is None:
            ants = self.ant_names
            antstr = 'all'
        else:
            antstr = ','.join(ants)
        outdata = {'ants': ants, 'freqs': self.freqs, 'pols': pols, 'source': self.source, 'uvh5': self.fnuvh5, 'freq_unit': self.freq_unit}
        print(f"Dumping autos in {self.fnuvh5} for {antstr} {pols}", end=' ... ')
        for ant in self.ant_names:
            for pol in pols:
                self.get_bl(ant, pol=pol)
                outdata[f"{ant}{pol}"] = copy(self.data)
        outdata['times'] = self.times.jd  # This assumes that all times in the UVH5 file are the same...
        mjd = outdata['times'][0] - 2400000.5
        obsrec = f"{self.source}_{mjd:.4f}_{self.lo}_{self.cnode}.npz"
        print(f"writing {obsrec}")
        np.savez(obsrec, **outdata)

    def get_bl(self, a, b=None, pol='xx'):
        """
        This reads in a baseline in all of the files read into obsrec_list.

        Parameters
        ----------
        a : str
            Antenna name
        b : str or None
            Another antenna name, same as 'a' if None
        pol : str
            Polarization
        
        Attributes
        ----------
        a : str
            Antenna a
        b : str
            Antenna b
        pol : str
            Polarization
        ano : int
            Antenna number in array (uvh5 only)
        bno : int
            Antenna number in array (uvh5 only)
        data : numpy array
            All of the data
        datamin : float
            Minimum value
        datamax : float
            Maximum value

        """
        self.a = a
        self.b = b
        self.pol = pol
        if self.b is None:
            self.b = self.a
        print(f"Reading ({self.a}, {self.b}){self.pol}", end='')
        dataf = []
        if self.tag == 'uvh5':  # Only one file
            self.ano = self.ant_map[self.a]
            self.bno = self.ant_map[self.b]
            self.data = self.uv.get_data(self.ano, self.bno, pol)
            self.times = Time(self.uv.get_times(self.ano, self.bno), format='jd')
        elif self.tag == 'npz':
            for obsrec in self.obsrec_list:
                dataf.append(self.npzfile[obsrec][f"{self.a}{pol}"])
            self.data = np.concatenate(dataf, axis=1)

        self.datamin = np.min(np.abs(self.data))
        self.datamax = np.max(np.abs(self.data))
        print(f"\tmin={self.datamin}, max={self.datamax}")

    def _axt_xaxis(self):
        if self.time_axis == 'a':  # actual datetime
            return self.times.datetime, 'Time'
        elif self.time_axis == 'd':  # difference
            return (self.times - self.obs.obsinfo.obsid[self.obsid].tref).to_value('sec'), 'sec'
        elif self.time_axis == 'b':  # boresight
            x = (self.times - self.obs.obsinfo.obsid[self.obsid].tref).to_value('sec') - self.obs.obsinfo.obsid[self.obsid].offset
            self.xp = self.obs.obsinfo.obsid[self.obsid].off_times
            self.yp = self.obs.obsinfo.obsid[self.obsid].off_boresight
            # Now extrapolate xp, yp to the limits of x...
            m_lo = (self.yp[1] - self.yp[0]) / (self.xp[1] - self.xp[0])
            b_lo = self.yp[0] - m_lo * self.xp[0]
            new_ylo = m_lo * x[0] + b_lo
            m_hi = (self.yp[-1] - self.yp[-2]) / (self.xp[-1] - self.xp[-2])
            b_hi = self.yp[-1] - m_hi * self.xp[-1]
            new_yhi = m_hi * x[-1] + b_hi
            self.xp = [x[0]] + self.xp + [x[-1]]
            self.yp = [new_ylo] + self.yp + [new_yhi]
            # ...ugly but working
            bbb = np.interp(x, self.xp, self.yp)
            return bbb, 'deg'
        
    def _invert_axis(self, dat, num=8, precision=-1):
        if isinstance(dat[0], datetime):
            dat = (self.times.jd - self.times[0].jd) * 24.0 * 3600.0
        idat = list(np.arange(len(dat)))
        xstart = np.round(np.floor(dat[0]), precision)
        xstep = np.round((dat[-1] - dat[0]) / num, precision)
        x = [int(xx) for xx in np.arange(xstart, dat[-1], xstep)]
        m = np.round(np.interp(x, dat, idat), 0)
        return m, x

    def dashboard_gen(self, obsid, lo='A', pol='xx', ant='2b', taxis='b'):
        """
        Parameters
        ----------
        obsid : str
            Can either be an obsid or a source in the obsinfo, or a json filename
        lo : str
            LO to use [A/B]
        pol : str
            pol to use [xx,yy,xy,yz]
        taxis : str
            t/x axis to use in plot [a/b/d]
    
        """
        self.obsid = obsid
        self.get_obsinfo()
        with open('dash.sh', 'w') as fp:
            for obsid in self.obs.obsinfo.obsid:
                print(f"on_obs.py {obsid} -a {ant} -t {taxis} --lo {lo} -p {pol} --dash -s", file=fp)

    def dashboard(self, ant='2b', pol='xx', use_db=True, save=False, time_axis='diff', show_obsinfo=False):
        """
        Parameters
        ----------
        ant : str
            antenna to use
        pol : str
            pol to use [xx,yy,xy,yz]
        taxis : str
            t/x axis to use in plot [a/b/d]
    
        """
        self.time_axis = time_axis[0].lower()
        x_axis_time, xlabel = self._axt_xaxis()
        x_ticks, x_labels = self._invert_axis(self.freqs)
        y_ticks, y_labels = self._invert_axis(x_axis_time)

        plt.figure('Dashboard', figsize=(16, 9))
        #, gridspec_kw={'width_ratios': [3, 1]}
        suptitle = f"{self.obsid}: {self.obs.obsinfo.obsid[self.obsid].tref.datetime.isoformat()}"
        try:
            suptitle += f"  ({self.obs.obsinfo.obsid[self.obsid].bf_distance} km)"
        except (AttributeError, KeyError):
            pass
        plt.suptitle(suptitle)
        # Water fall plot
        axwf = plt.subplot2grid((2, 2), (0, 0), rowspan=2, colspan=1)
        axt = plt.subplot2grid((2, 2), (0, 1), rowspan=1, colspan=1)
        axf = plt.subplot2grid((2, 2), (1, 1), rowspan=1, colspan=1)
        if ant:
            self.get_bl(ant, pol=pol)
        axwf.imshow(toMag(self.data, use_db))
        axwf.set_aspect('auto')
        axwf.set_xlabel('Freq')
        axwf.set_ylabel(xlabel)
        axwf.set_xticks(x_ticks, x_labels)
        axwf.set_yticks(y_ticks, y_labels)

        # Time plot
        filter = self.obs.obsinfo.filters[self.lo]
        used_filter = {}
        for clr, filt in filter.items():
            if filt[1] < self.freqs[0] or filt[0] > self.freqs[-1]:
                continue
            used_filter[clr] = filt
            tax2 = self.get_sum(over='freq', dmin=filt[0], dmax=filt[1], use_db=use_db)
            axt.plot(x_axis_time, tax2, label=f"{filt[0]}-{filt[1]}", color=clr)
        if show_obsinfo:
            self.plot_obsinfo_times(axt)
        axt.legend()
        axt.set_xlabel(xlabel)
        if use_db:
            axt.set_ylabel('dB')
        # axt.set_xlim(left=-15.0, right=15.0)

        # Frequency plot
        for i in range(len(self.times)):
            axf.plot(self.freqs, toMag(self.data[i], use_db), '0.8')
        axf.set_xlabel(self.freq_unit)
        try:
            dmin = self.obs.obsinfo.obsid[self.obsid].t0.jd
            dmax = self.obs.obsinfo.obsid[self.obsid].t1.jd
            if dmin < self.times.jd[0]:
                dmin = self.times.jd[0]
        except AttributeError:
            dmin = self.times.jd[0]
            dmax = self.times.jd[-1]
        fax1 = self.get_sum(over='time', dmin=dmin, dmax=dmax, norm=True, use_db=use_db)
        axf.plot(self.freqs, fax1, 'k', linewidth=2)
        if use_db:
            axf.set_ylabel('dB')
        axmin, axmax = axf.axis()[2], axf.axis()[3]
        for clr, filt in used_filter.items():
            axf.plot([filt[0], filt[0]], [axmin, axmax], color=clr)
            axf.plot([filt[1], filt[1]], [axmin, axmax], color=clr)
        axf.set_ylim(bottom=axmin, top=axmax)
        if save:
            fn = f"{self.obsid}_{ant}_{pol}.png"
            plt.savefig(fn)

    def get_sum(self, over='freq', dmin=1990.0, dmax=1995.0, norm=False, plotting='amplitude', use_db=True):
        ave = []
        if over[0].lower() == 'f':
            npfr = np.array(self.freqs)
            lp = np.where(npfr >= dmin)
            up = np.where(npfr[lp] <= dmax)
            ind = lp[0][up[0]]
            for i in range(len(self.times)):
                ave.append(np.sum(np.abs(self.data[i][ind])))
            ave = np.array(ave)
            if norm:
                ave /= len(self.times)
        elif over[0].lower() == 't':
            lp = np.where(self.times.jd >= dmin)
            up = np.where(self.times.jd[lp] <= dmax)
            ind = lp[0][up[0]]
            for i in range(len(self.freqs)):
                ave.append(np.sum(np.abs(self.data[:, i][ind])))
            ave = np.array(ave)
            if norm:
                ave /= len(self.freqs)
        ave = toMag(ave, use_db)
        return ave

    def all_wf(self, pol='xx', use_db=True, save=False):
        fig, axs = plt.subplots(nrows=4, ncols=7, figsize=(16, 9), tight_layout=True)
        ctr = 0
        for i in range(4):
            for j in range(7):
                self.get_bl(self.ant_names[ctr], pol=pol)
                axs[i][j].imshow(10.0 * toMag(self.data, use_db))
                axs[i][j].set_aspect('auto')
                axs[i][j].set_title(f"{self.a}")
                axs[i][j].set_xticks([])
                axs[i][j].set_yticks([])
                ctr += 1
        fig.suptitle(f"{self.obsid}:{pol}")
        if save:
            fn = f"{self.obsid}_{pol}.png"
            plt.savefig(fn)

    def plot_wf(self, plotting='amplitude', use_db=True):
        ptitle = (f'WF: {self.obsid} - ({self.a},{self.b}){self.pol}')
        fig, ax = plt.subplots()
        ax.set_title(ptitle)
        if plotting[0].lower() == 'a':
            pdata = 10.0 * toMag(self.data, use_db)
        ax.imshow(pdata)
        ax.set_aspect('auto')

    def plot_freqs(self, ch='all', plotting='amplitude', use_db=True):
        plt.figure(f'Freqs: {self.obsid} - ({self.a}, {self.b})')
        for i in range(len(self.times)):
            if plotting[0].lower() == 'a':
                pdata = toMag(self.data[i], use_db)
            plt.plot(self.freqs, pdata)        
        plt.xlabel(self.freq_unit)
    
    def plot_times(self, trange='all', plotting='amplitude', use_db=True):
        plt.figure(f'Times: {self.obsid} - ({self.a}, {self.b})')
        for i in range(len(self.freqs)):
            if plotting[0].lower() == 'a':
                pdata = toMag(self.data[:, i], use_db)
            plt.plot(self.times.datetime, pdata)

    def plot_obsinfo_times(self, ax):
        if self.obs is None:
            return
        dmin = ax.axis()[2]
        dmax = ax.axis()[3]
        this_src = self.obs.obsinfo.obsid[self.obsid]
        try:
            tref = this_src.tref
        except AttributeError:
            tref = self.times[len(self.times // 2)]
        try:
            t0 = this_src.t0
        except AttributeError:
            t0 = self.times[0]
        try:
            t1 = this_src.t1
        except AttributeError:
            t1 = self.times[-1]
        if self.time_axis == 'a':  # absolute
            sups = [[tref.datetime, 'r--'],
                    [t0.datetime, 'k--'],
                    [t1.datetime, 'k--']]
        elif self.time_axis == 'd':  # diff
            sups = [[0.0, 'r--'],
                    [(t0 - tref).to_value('sec'), 'k--'],
                    [(t1 - tref).to_value('sec'), 'k--']]
        elif self.time_axis == 'b':
            x = [-this_src.offset,
                 (t0 - tref).to_value('sec') - this_src.offset,
                 (t1 - tref).to_value('sec') - this_src.offset
            ]
            b = np.interp(x, self.xp, self.yp)
            sups = [[b[0], 'r--'], [b[1], 'k--'], [b[2], 'k--']]
        for x, lsc in sups:
            ax.plot([x, x], [dmin, dmax], lsc, linewidth=2)
        ax.set_ylim(bottom=dmin, top=dmax)

    def get_obsinfo(self):
        self.obs = obs_base.Base()
        self.obs.read_obsinfo(obs=self.obsid)

    def put_obsinfo(self, fn, fdict):
        self.obs = obs_base.Base()
        self.obs.write_obsinfo(fn, fdict)