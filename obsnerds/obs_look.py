from pyuvdata import UVData
from astropy.time import Time
import numpy as np
import matplotlib.pyplot as plt
from obsnerds import obs_base
from datetime import datetime
import os.path as path
from . import obs_sys as OS


def toMag(x, use_db=True):
    if use_db:
        return 10.0 * np.log10(np.abs(x))
    else:
        return np.abs(x)


class Look:
    def __init__(self, obsinput=None, lo='A', cnode='all', freq_unit='MHz', dir_data='.'):
        """
        This initializes, but also reads in all of the data.

        Parameters
        ----------
        obsinput : None, str
            Either a file (ends with .uvh5 or .npz) or obsid

        """
        self.obsid = None
        self.lo = lo
        self.cnode = OS.make_cnode(cnode)
        self.freq_unit = freq_unit
        self.dir_data = dir_data
        self.npzfile = {}
        self.freqs = []
        if obsinput is not None:
            if obsinput.endswith('.uvh5') or obsinput.endswith('.npz'):
                self.obsrec_files = obsinput.split(',')
            else:
                self.obsid = obsinput
                self.get_obsinfo()
                self.source, self.mjd = OS.split_obsid(self.obsid)
                self.obsrec_files = [f"{self.obsid}_{self.lo}_{x}.npz" for x in self.cnode]
            self.read_in_files()

    def read_in_files(self):
        for obsrec in self.obsrec_files:
            if obsrec.endswith('.uvh5'):
                self.read_a_uvh5(obsrec)
                if len(self.obsrec_files) > 1:
                    print("Only reading in first uvh5 file")
                    break
            elif obsrec.endswith('npz'):
                self.read_an_npz(obsrec)
            else:
                print(f"Invalid file {obsrec}")

    def read_a_uvh5(self, fn):
        print(f"Reading {fn}")
        self.file_type = 'uvh5'
        self.uvh5_pieces = OS.parse_uvh5_filename(fn)
        self.fn = fn
        self.source = self.uvh5_pieces['source']
        self.uv = UVData()
        self.uv.read(self.fn)
        self.ant_numbers = self.uv.get_ants()
        self.ant_names = np.array(self.uv.antenna_names)[np.unique(self.uv.ant_1_array) - 1]
        self.ant_map = {}
        for antno, antna in zip(self.ant_numbers, self.ant_names):
            self.ant_map[antna] = antno
        self.freqs = self.uv.freq_array[0] / OS.FREQ_CONVERT[self.freq_unit]
        return True

    def read_an_npz(self, obsrec_file):
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
        self.file_type = 'npz'
        fn = path.join(self.obs.obsinfo.dir_data, obsrec_file)
        X = OS.split_obsrec(obsrec_file)
        try:
            self.npzfile[obsrec_file] = np.load(fn)
        except FileNotFoundError:
            print(f"Couldn't find {fn}")
            return False
        self.ant_names = list(self.npzfile[obsrec_file]['ants'])
        self.freqs += list(self.npzfile[obsrec_file]['freqs'])
        self.times = Time(self.npzfile[obsrec_file]['times'], format='jd')
        return True

    def get_bl(self, a, b=None, pol='xx'):
        """
        This reads in a baseline in all of the files read into obsrec_files.

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
        if self.file_type == 'uvh5':  # Only one file
            self.ano = self.ant_map[self.a]
            self.bno = self.ant_map[self.b]
            self.data = self.uv.get_data(self.ano, self.bno, pol)
            self.times = Time(self.uv.get_times(self.ano, self.bno), format='jd')
        elif self.file_type == 'npz':
            for obsrec_file in self.obsrec_files:
                dataf.append(self.npzfile[obsrec_file][f"{self.a}{pol}"])
            self.data = np.concatenate(dataf, axis=1)

        self.datamin = np.min(np.abs(self.data))
        self.datamax = np.max(np.abs(self.data))
        print(f"\tmin={self.datamin}, max={self.datamax}")

    def get_time_axes(self):
        """
        Get the values for the various "time" axis options (which is vertical in the waterplot)
        
        Attribute
        ---------
        taxes : dict
            dictionary containing the axis values and axis label

        """

        self.taxes = {
            'datetime': {
                'values': self.times.datetime,
                'label': 'Time'},
            'seconds': {
                'values': (self.times - self.obs.obsinfo.obsid[self.obsid].tref).to_value('sec'),
                'label': 'Seconds',
                'reference': self.obs.obsinfo.obsid[self.obsid].tref}
        }
        try:
            self.xp = self.obs.obsinfo.obsid[self.obsid].off_times
            self.yp = self.obs.obsinfo.obsid[self.obsid].off_boresight
        except AttributeError:
            self.taxes['boresight'] = None
            return
        x = self.taxes['seconds']['values'] - self.obs.obsinfo.obsid[self.obsid].offset
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
        self.taxes['boresight'] = {
            'values': np.interp(x, self.xp, self.yp),
            'label': 'Degrees',
            'offset': self.obs.obsinfo.obsid[self.obsid].offset
        }
        
    def _get_wf_ticks(self, dat, ticks=8, precision=-1):
        if isinstance(dat[0], datetime):
            dat = (self.times.jd - self.times[0].jd) * 24.0 * 3600.0
        idat = list(np.arange(len(dat)))
        if isinstance(ticks, (float, int)):
            xstart = np.round(np.floor(dat[0]), precision)
            xstep = np.round((dat[-1] - dat[0]) / ticks, precision)
            x = [int(xx) for xx in np.arange(xstart, dat[-1], xstep)]
        elif isinstance(ticks, list):
            if len(ticks) == 3:
                x = [int(xx) for xx in np.arange(ticks[0], ticks[1]+1, ticks[2])]
            else:
                x = ticks
        else:
            raise ValueError("Invalid ticks parameters")
        m = np.round(np.interp(x, dat, idat), 0)
        return m, x

    def dashboard_gen(self, obsinfo, script_fn='dash.sh', ants='2b,4e', pols='xx,xy', taxis='b'):
        """
        Parameters
        ----------
        obsinfo : str
            Name of obsinfo file to use
        script_fn : str
            Name of bash script file to write
        ant : str
            ants to use
        pol : str
            pol to use [xx,yy,xy,yz]
        taxis : str
            t/x axis to use in plot [a/b/d]
    
        """
        self.obs = obs_base.Base()
        self.obs.read_obsinfo(obs=obsinfo)
        cnode = ','.join(self.cnode)
        ants = OS.listify(ants)
        pols = OS.listify(pols, {'all': ['xx', 'xy', 'yy', 'yx']})
        with open(script_fn, 'w') as fp:
            for obsid in self.obs.obsinfo.obsid:
                for ant in ants:
                    for pol in pols:
                        print(f"on_obs.py {obsid} -a {ant}  -p {pol} -t {taxis} --lo {self.lo} --cnode {cnode} --dash -s", file=fp)

    def dashboard(self, ant='2b', pol='xx', time_axis='seconds', **kwargs):
        """
        Parameters
        ----------
        ant : str
            antenna to use
        pol : str
            pol to use [xx,yy,xy,yz]
        time_axis : str
            t axis to use in plot [datetime/boresight/seconds]
        kwargs : use_db, save, t_wfticks, f_wfticks, zoom_time, zoom_freq, filter_time
    
        """
        use_db = kwargs['use_db'] if 'use_db' in kwargs else True
        save = kwargs['save'] if 'save' in kwargs else False
        t_wfticks = kwargs['t_wfticks'] if 't_wfticks' in kwargs else [-120, 120, 20]
        f_wfticks = kwargs['f_wfticks'] if 'f_wfticks' in kwargs else 8
        zoom_time = kwargs['zoom_time'] if 'zoom_time' in kwargs else False
        zoom_freq = kwargs['zoom_freq'] if 'zoom_freq' in kwargs else False
        filter_time = kwargs['filt_time'] if 'filt_time' in kwargs else {}
        filter_time.update({'k': [self.times.jd[0], self.times.jd[-1]]})

        self.time_axis = OS.AXIS_OPTIONS[time_axis[0].lower()]
        self.get_time_axes()
        if ant:  # Otherwise use existing
            self.get_bl(ant, pol=pol)
        self.apply_obsinfo_freq_filters(use_db=use_db)

        plt.figure('Dashboard', figsize=(16, 9))
        #, gridspec_kw={'width_ratios': [3, 1]}
        suptitle = f"{self.obsid}: {self.obs.obsinfo.obsid[self.obsid].tref.datetime.isoformat()}"
        try:
            suptitle += f"  ({self.obs.obsinfo.obsid[self.obsid].bf_distance} km)"
        except (AttributeError, KeyError):
            pass
        suptitle += f'  --  ({self.a},{self.b}) {self.pol}'
        plt.suptitle(suptitle)
        axwf = plt.subplot2grid((2, 2), (0, 0), rowspan=2, colspan=1)
        axt = plt.subplot2grid((2, 2), (0, 1), rowspan=1, colspan=1)
        axf = plt.subplot2grid((2, 2), (1, 1), rowspan=1, colspan=1)

        # Water fall plot
        f_wfticks, f_wftick_labels = self._get_wf_ticks(self.freqs, ticks=f_wfticks)
        t_wfticks, t_wftick_labels = self._get_wf_ticks(self.taxes[time_axis]['values'], ticks=t_wfticks)
        axwf.imshow(toMag(self.data, use_db))
        axwf.set_aspect('auto')
        axwf.set_xlabel('Freq')
        axwf.set_ylabel(self.taxes[time_axis]['label'])
        axwf.set_xticks(f_wfticks, f_wftick_labels)
        axwf.set_yticks(t_wfticks, t_wftick_labels)

        # Time plot
        for clr, data in self.filtered_power.items():
            axt.plot(self.taxes[time_axis]['values'], data, label=f"{self.used_filters[clr][0]}-{self.used_filters[clr][1]}", color=clr)
        axt.legend()
        axt.set_xlabel(self.taxes[time_axis]['label'])
        ylabel = 'dB' if use_db else 'linear'
        axt.set_ylabel(ylabel)
        if zoom_time:
            axt.set_xlim(left=zoom_time[0], right=zoom_time[1])

        # Frequency plot
        for i in range(len(self.times)):
            axf.plot(self.freqs, toMag(self.data[i], use_db), '0.8')
        axf.set_xlabel(self.freq_unit)
        for clr, filt in filter_time.items():
            fax1 = self.power_filter(over='time', dmin=filt[0], dmax=filt[1], norm=True, use_db=use_db)
            axf.plot(self.freqs, fax1, clr, linewidth=2 if clr=='k' else 1)
        axf.set_ylabel(ylabel)
        axmin, axmax = axf.axis()[2], axf.axis()[3]
        for clr, filt in self.used_filters.items():
            axf.plot([filt[0], filt[0]], [axmin, axmax], color=clr)
            axf.plot([filt[1], filt[1]], [axmin, axmax], color=clr)
        axf.set_ylim(bottom=axmin, top=axmax)
        if zoom_freq:
            axt.set_xlim(left=zoom_freq[0], right=zoom_freq[1])

        if save:
            fn = f"{self.obsid}_{ant}_{pol}.png"
            plt.savefig(fn)

    def apply_obsinfo_freq_filters(self, use_db=True):
        self.filtered_power = {}
        self.used_filters = {}
        filters = self.obs.obsinfo.filters[self.lo]
        for clr, filt in filters.items():
            if filt[1] < self.freqs[0] or filt[0] > self.freqs[-1]:
                continue
            self.used_filters[clr] = filt
            self.filtered_power[clr] = self.power_filter(over='freq', dmin=filt[0], dmax=filt[1], use_db=use_db)

    def power_filter(self, over='freq', dmin=1990.0, dmax=1995.0, norm=False, use_db=True):
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

    def get_obsinfo(self):
        self.obs = obs_base.Base()
        self.obs.read_obsinfo(obs=self.obsid)
