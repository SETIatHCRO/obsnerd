from pyuvdata import UVData
from astropy.time import Time
import numpy as np
import matplotlib.pyplot as plt
from . import onv_base
from datetime import datetime
import os.path as path
from . import on_sys as OS
from odsutils import ods_tools as tools


def toMag(x, use_dB=True):
    if use_dB:
        return 10.0 * np.log10(np.abs(x))
    else:
        return np.abs(x)


class Filter:
    def __init__(self, ftype=None, unit=None, lo=None, hi=None, norm=False, color='k', shape='rect', invert=False):
        self.use = True
        self.ftype = ftype
        self.axis = OS.FILTER_AXIS[self.ftype]
        self.unit = unit
        self.lo = lo
        self.hi = hi
        self.norm = norm
        self.color = color
        self.shape = 'rect'  # only option for now
        self.invert = invert

    def __str__(self):
        s = '<'
        for p in ['use', 'ftype', 'unit', 'lo', 'hi', 'norm', 'color']:
            s += f"{p}={getattr(self, p)},\t"
        s += '>'
        return s

    def apply(self, x, data, use_dB=None):
        """
        Parameters
        ----------
        x : array/list
            x axis values
        y : array/list
            y axis values
        kwargs : see list under __init__

        """
        if self.hi < x[0] or self.lo > x[-1]:
            self.use = False
            return
        self.use = True
        if not isinstance(x, np.ndarray):
            x = np.array(x)
        lp = np.where(x >= self.lo)
        up = np.where(x[lp] <= self.hi)
        self.inds = lp[0][up[0]]
        filterarr = np.zeros(data.shape)
        self.indicator = [[self.lo, self.hi]]
        if self.invert:
            self._invert(len(x), xmm=[min(x), max(x)])
        if self.axis == 0:
            filterarr[self.inds, :] = 1
        elif self.axis == 1:
            filterarr[:, self.inds] = 1
        power = np.sum(data * filterarr, axis=self.axis)
        self.power = power if not self.norm else power / len(self.inds)
        return self.use

    def _invert(self, N, xmm):
        new_ind = []
        for i in range(N):
            if i not in self.inds:
                new_ind.append(i)
        self.inds = new_ind
        self.indicator = [[xmm[0], self.lo], [self.hi, xmm[1]]]


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
        self.filters = {}
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
            Obsrec designation for an observation file

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

    def dashboard_gen(self, obsinfo, script_fn='dash.sh', ants='2b,4e', pols='xx,xy', taxis='b', show_diff=False):
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
        self.obs = onv_base.Base()
        self.obs.read_obsinfo(obs=obsinfo)
        cnode = ','.join(self.cnode)
        ants = tools.listify(ants)
        pols = tools.listify(pols, {'all': ['xx', 'xy', 'yy', 'yx']})
        show_diff = ' --show_diff ' if show_diff else ''
        with open(script_fn, 'w') as fp:
            for obsid in self.obs.obsinfo.obsid:
                for ant in ants:
                    for pol in pols:
                        print(f"on_look.py {obsid} -a {ant}  -p {pol} -t {taxis} --lo {self.lo} --cnode {cnode} {show_diff} --dash -s", file=fp)

    def dashboard(self, ant='2b', pol='xx', time_axis='seconds', transit_time=4.0, **kwargs):
        """
        Parameters
        ----------
        ant : str
            antenna to use
        pol : str
            pol to use [xx,yy,xy,yz]
        time_axis : str
            t axis to use in plot [datetime/boresight/seconds]
        transit_time : float
            time around zero to average
        kwargs : use_dB, save, t_wfticks, f_wfticks, zoom_time, zoom_freq, filter_time
    
        """
        use_dB = kwargs['use_dB'] if 'use_dB' in kwargs else True
        save = kwargs['save'] if 'save' in kwargs else False
        t_wfticks = kwargs['t_wfticks'] if 't_wfticks' in kwargs else [-120, 120, 20]
        f_wfticks = kwargs['f_wfticks'] if 'f_wfticks' in kwargs else 8
        zoom_time = kwargs['zoom_time'] if 'zoom_time' in kwargs else False
        zoom_freq = kwargs['zoom_freq'] if 'zoom_freq' in kwargs else False
        filter_time = kwargs['filt_time'] if 'filt_time' in kwargs else {}
        show_diff = kwargs['show_diff'] if 'show_diff' in kwargs else False

        self.time_axis = OS.AXIS_OPTIONS[time_axis[0].lower()]
        self.get_time_axes()
        if ant:  # Otherwise use existing
            self.get_bl(ant, pol=pol)

        # Set up filters
        self.filters = {}
        tt = transit_time / 2.0
        self.filters['on:boresight'] = Filter(ftype='time', unit=self.time_axis, lo=-tt, hi=tt, norm=True, color='r')
        self.filters['off:boresight'] = Filter(ftype='time', unit=self.time_axis, lo=-tt, hi=tt, norm=True, color='k', invert=True)
        for clr, filt in filter_time:
            self.filters[f"time:{filt[0]}-{filt[1]}:{clr}"] = Filter(ftype='time', unit=self.time_axis, lo=filt[0], hi=filt[1], norm=True, color=clr)
        for clr, filt in self.obs.obsinfo.filters[self.lo].items():
            self.filters[f"freq:{filt[0]}-{filt[1]}:{clr}"] = Filter(color=clr, ftype='freq', unit='MHz', lo=filt[0], hi=filt[1])

        plt.figure('Dashboard', figsize=(16, 9))
        #, gridspec_kw={'width_ratios': [3, 1]}
        self.suptitle = f"{self.obsid}: {self.obs.obsinfo.obsid[self.obsid].tref.datetime.isoformat()}"
        try:
            self.suptitle += f"  ({self.obs.obsinfo.obsid[self.obsid].bf_distance} km)"
        except (AttributeError, KeyError):
            pass
        self.suptitle += f'  --  ({self.a},{self.b}) {self.pol}'
        plt.suptitle(self.suptitle)
        axwf = plt.subplot2grid((2, 2), (0, 0), rowspan=2, colspan=1)
        axt = plt.subplot2grid((2, 2), (0, 1), rowspan=1, colspan=1)
        axf = plt.subplot2grid((2, 2), (1, 1), rowspan=1, colspan=1)

        # Water fall plot
        f_wfticks, f_wftick_labels = self._get_wf_ticks(self.freqs, ticks=f_wfticks)
        t_wfticks, t_wftick_labels = self._get_wf_ticks(self.taxes[time_axis]['values'], ticks=t_wfticks)
        axwf.imshow(toMag(self.data, use_dB))
        axwf.set_aspect('auto')
        axwf.set_xlabel('Freq')
        axwf.set_ylabel(self.taxes[time_axis]['label'])
        axwf.set_xticks(f_wfticks, f_wftick_labels)
        axwf.set_yticks(t_wfticks, t_wftick_labels)

        # Time plot
        for key in [k for k in self.filters if self.filters[k].ftype=='freq']:
            if self.filters[key].apply(self.freqs, self.data):
                axt.plot(self.taxes[time_axis]['values'], toMag(self.filters[key].power, use_dB), color=self.filters[key].color)
        axt.set_xlabel(self.taxes[time_axis]['label'])
        ylabel = 'dB' if use_dB else 'linear'
        axt.set_ylabel(ylabel)
        axtlim = axt.axis()

        # Frequency plot
        for i in range(len(self.times)):
            axf.plot(self.freqs, toMag(self.data[i], use_dB), '0.8')
        axf.set_xlabel(self.freq_unit)
        for key in [k for k in self.filters if self.filters[k].ftype=='time']:
            if self.filters[key].apply(self.taxes[time_axis]['values'], self.data):
                axf.plot(self.freqs, toMag(self.filters[key].power, use_dB), color=self.filters[key].color)
        axf.set_ylabel(ylabel)
        axflim = axf.axis()

        # Indicate filters and set limits
        # ... time
        axtshade = axtlim[2] + 0.15 * (axtlim[3] - axtlim[2])
        for key in [k for k in self.filters if self.filters[k].ftype=='time' and self.filters[k].use]:
            for xlim in self.filters[key].indicator:
                axt.fill_between(xlim, [axtshade, axtshade], color=self.filters[key].color)
        axt.set_ylim(bottom=axtlim[2], top=axtlim[3])
        if zoom_time:
            axt.set_xlim(left=zoom_time[0], right=zoom_time[1])
        else:
            axt.set_xlim(left=self.taxes[time_axis]['values'][0], right=self.taxes[time_axis]['values'][-1])
        # ... freq
        axfshade = axflim[2] + 0.15 * (axflim[3] - axflim[2])
        for key in [k for k in self.filters if self.filters[k].ftype=='freq' and self.filters[k].use]:
            for xlim in self.filters[key].indicator:
                axf.fill_between(xlim, [axfshade, axfshade], color=self.filters[key].color)
        axf.set_ylim(bottom=axflim[2], top=axflim[3])
        if zoom_freq:
            axf.set_xlim(left=zoom_freq[0], right=zoom_freq[1])
        else:
            axf.set_xlim(left=self.freqs[0], right=self.freqs[-1])

        fn = f"{self.obsid}_{ant}_{pol}.png"
        if save:
            plt.savefig(fn)

        if show_diff:
            self.plot_diff(use_dB=use_dB, save=save, fn=f'Diff_{fn}')

    def plot_diff(self, use_dB=True, save=False, fn=None, num='on:boresight', den='off:boresight'):
        plt.figure("On / Off: "+self.suptitle, figsize=(12,6))
        plt.title("On / Off: "+self.suptitle)
        dp = self.filters[num].power / self.filters[den].power
        plt.plot(self.freqs, toMag(dp, use_dB), color='k')
        plt.fill_between(self.freqs, toMag(dp, use_dB), color='0.8')
        plt.xlabel('Freq [MHz]')
        ylabel = 'dB' if use_dB else ''
        plt.ylabel(ylabel)
        plt.grid()
        plt.axis([1900.0, 2060, -0.35, 3.5])
        if save:
            plt.savefig(fn)

    def all_wf(self, pol='xx', use_dB=True, save=False):
        fig, axs = plt.subplots(nrows=4, ncols=7, figsize=(16, 9), tight_layout=True)
        ctr = 0
        for i in range(4):
            for j in range(7):
                self.get_bl(self.ant_names[ctr], pol=pol)
                axs[i][j].imshow(10.0 * toMag(self.data, use_dB))
                axs[i][j].set_aspect('auto')
                axs[i][j].set_title(f"{self.a}")
                axs[i][j].set_xticks([])
                axs[i][j].set_yticks([])
                ctr += 1
        fig.suptitle(f"{self.obsid}:{pol}")
        if save:
            fn = f"{self.obsid}_{pol}.png"
            plt.savefig(fn)

    def plot_wf(self, plotting='amplitude', use_dB=True):
        ptitle = (f'WF: {self.obsid} - ({self.a},{self.b}){self.pol}')
        fig, ax = plt.subplots()
        ax.set_title(ptitle)
        if plotting[0].lower() == 'a':
            pdata = 10.0 * toMag(self.data, use_dB)
        ax.imshow(pdata)
        ax.set_aspect('auto')

    def plot_freqs(self, use_dB=True):
        plt.figure(f'Freqs: {self.obsid} - ({self.a}, {self.b})')
        for i in range(len(self.times)):
            plt.plot(self.freqs, toMag(self.data[i], use_dB))        
        plt.xlabel(self.freq_unit)
    
    def plot_times(self, use_dB=True):
        plt.figure(f'Times: {self.obsid} - ({self.a}, {self.b})')
        for i in range(len(self.freqs)):
            plt.plot(self.times.datetime, toMag(self.data[:, i], use_dB))

    def get_obsinfo(self):
        self.obs = onv_base.Base()
        self.obs.read_obsinfo(obs=self.obsid)
