from pyuvdata import UVData
from astropy.time import Time
import numpy as np
from copy import copy
import matplotlib.pyplot as plt
from datetime import datetime
import os.path as path
from . import on_sys, on_track, onv_logs
from odsutils import ods_tools as tools
import astropy.units as u


DEFAULT_DATA_DIR = 'data'

def toMag(x, use_dB=True):
    absx = np.abs(x)
    if use_dB:
        bad = np.where(absx < 1E-12)
        absx[bad] = 1E-6
        return 10.0 * np.log10(absx)
    else:
        return absx

class Filter:
    def __init__(self, ftype=None, unit=None, lo=None, hi=None, norm=False, color='k', shape='rect', invert=False):
        self.use = True
        self.ftype = ftype
        self.axis = on_sys.FILTER_AXIS[self.ftype]
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
            self.do_invert(len(x), xmm=[min(x), max(x)])
        if self.axis == 0:
            try:
                filterarr[self.inds, :] = 1
            except IndexError:
                pass
        elif self.axis == 1:
            try:
                filterarr[:, self.inds] = 1
            except IndexError:
                pass
        power = np.sum(data * filterarr, axis=self.axis)
        self.power = power if not self.norm else power / len(self.inds)
        return self.use

    def do_invert(self, N, xmm):
        new_ind = []
        for i in range(N):
            if i not in self.inds:
                new_ind.append(i)
        self.inds = new_ind
        self.indicator = [[xmm[0], self.lo], [self.hi, xmm[1]]]


class Look:
    def __init__(self, oinput=None, lo='A', cnode='all', freq_unit='MHz'):
        """
        This initializes the Look class.  Generally will then use "get" to read in or prep data.

        Parameters
        ----------
        lo : str
            LO designation
        cnode : str
            Cnode designation
        freq_unit : str
            Frequency unit

        """
        self.lo = on_sys.make_lo(lo)
        self.cnode = on_sys.make_cnode(cnode)
        self.freq_unit = freq_unit
        self.npzfile = {}
        self.freqs = []
        self.filters = {}
        self.obs = None
        self.get(oinput)

    def get(self, oinput):
        """
        This reads in the input and sets up the class for further processing.

        Input may be one of the following:
         - uvh5 file (.uvh5)
         - npz file (.npz)
         - obsinfo filename (.json or mjd)
         - obsid (str parseable into source/mjd)
         - source name (str with files in obs.dir_data)

        """
        # Some pre-processing
        if oinput is None:
            return
        try:
            oinput_is_mjd = float(oinput)
        except (ValueError, TypeError):
            oinput_is_mjd = False
        # Now get info
        if oinput.endswith('.uvh5'):
            self.read_a_uvh5(oinput)
        elif oinput.endswith('.npz'):
            self.read_an_npz(oinput)
        elif oinput_is_mjd or oinput.endswith('.json'):  # Read in obsinfo file
            self.obs = on_track.read_obsinfo(oinput)
            print(f"Reading {self.obs.filename} into self.obs")
            self.found_obsids = []
            for source in self.obs.sources:
                self.found_obsids.append(self.obs.sources[source].obsid)
            print("Found obsids:  ", ', '.join(self.found_obsids))
        else:  # This is an obsid or source, and want to end up with obsid, source, mjd
            dir_data = DEFAULT_DATA_DIR if self.obs is None else self.obs.dir_data
            _source, self.mjd = on_sys.split_obsid(oinput)
            self.source = oinput if self.mjd is None else _source
            self.obsid = on_track.get_obsid_from_source(self.source, dir_data) if self.mjd is None else oinput
            if self.obsid is None:
                print(f"Couldn't find obsid from {self.source} in {dir_data}")
                return
            else:
                self.source, self.mjd = on_sys.split_obsid(self.obsid)
            print(f"Setting obsrec_files for {self.obsid}")
            self.obsrec_files = [] if self.obsid is None else [f"{self.obsid}_{self.lo}_{x}.npz" for x in self.cnode]
            if self.obs is None:
                self.get(str(self.mjd))

    def read_in_files(self):
        self.freqs = []
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
        if self.file_type == 'npz':
            self.freqs = [f.to_value(self.freq_unit) for f in self.freqs]

    def read_a_uvh5(self, fn):
        print(f"Reading {fn}")
        self.file_type = 'uvh5'
        self.uvh5_pieces = on_sys.parse_uvh5_filename(fn)
        self.fn = fn
        self.source = self.uvh5_pieces['source']
        self.uv = UVData()
        self.uv.read(self.fn)
        self.ant_numbers = self.uv.get_ants()
        self.ant_names = np.array(self.uv.antenna_names)[np.unique(self.uv.ant_1_array) - 1]
        self.ant_map = {}
        for antno, antna in zip(self.ant_numbers, self.ant_names):
            self.ant_map[antna] = antno
        self.freqs = (self.uv.freq_array[0] * u.Unit('Hz')).to_value(self.freq_unit)
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
        try:
            self.fn = path.join(self.obs.dir_data, obsrec_file)
        except AttributeError:
            self.fn = obsrec_file
        try:
            self.npzfile[obsrec_file] = np.load(self.fn)
        except FileNotFoundError:
            print(f"Couldn't find {self.fn}")
            self.npzfile[obsrec_file] = None
            return False
        print(f"Reading {self.fn}")
        self.ant_names = list(self.npzfile[obsrec_file]['ants'])
        self.freq_unit = str(self.npzfile[obsrec_file]['freq_unit'])
        self.freqs += list(self.npzfile[obsrec_file]['freqs'] * u.Unit(self.freq_unit))
        self.times = Time(self.npzfile[obsrec_file]['times'], format='jd')
        self.pols = list(self.npzfile[obsrec_file]['pols'])
        return True

    def get_bl(self, a, b=None, pol='xx', auto_as_abs=True):
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
        auto_as_abs : bool
            Flag to return autos as absolute value
        
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
        is_auto = (self.a == self.b)
        if self.file_type == 'npz' and not is_auto:
            raise ValueError("For no good reason I've limited npz files to autos only")
        print(f"Reading ({self.a},{self.b}){self.pol}", end='')
        dataf = []
        if self.file_type == 'uvh5':  # Only one file
            self.ano = self.ant_map[self.a]
            self.bno = self.ant_map[self.b]
            self.data = self.uv.get_data(self.ano, self.bno, pol)
            if is_auto and auto_as_abs:
                self.data = np.abs(self.data)
            self.times = Time(self.uv.get_times(self.ano, self.bno), format='jd')
        elif self.file_type == 'npz':
            for obsrec_file in self.obsrec_files:
                if self.npzfile[obsrec_file] is None:
                    continue
                if f"{self.a}{pol}" not in self.npzfile[obsrec_file].keys():
                    continue
                if is_auto and auto_as_abs:
                    dataf.append(np.abs(self.npzfile[obsrec_file][f"{self.a}{pol}"]))
                else:
                    dataf.append(self.npzfile[obsrec_file][f"{self.a}{pol}"])
            if len(dataf):
                self.data = np.concatenate(dataf, axis=1)

        try:
            self.datamin = np.min(np.abs(self.data))
            self.datamax = np.max(np.abs(self.data))
            print(f"\tmin={self.datamin}, max={self.datamax}")
        except AttributeError:
            print("\nNo data found")
            self.data = None

    def get_time_axes(self, log=None):
        """
        Get the values for the various "time" axis options (which is vertical in the waterplot)

        Parameter
        ---------
        logtimes : TBALog
            Log times for SpaceX
        
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
                'values': (self.times - self.this_source.utc).to_value('sec'),
                'label': 'Seconds',
                'reference': self.this_source.utc}
        }
        try:
            self.xp = self.this_source.off_time
            self.yp = self.this_source.off_angle
        except AttributeError:
            self.taxes['boresight'] = None
            return
        x = self.taxes['seconds']['values']
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
        }
        if log:
            for scope in ['inner', 'outer']:
                self.taxes[scope] = {}
                self.taxes[scope]['datetime'] = {
                    'values': log.times[scope],
                    'label': 'Time',
                    'reference': self.this_source.utc
                }
                self.taxes[scope]['seconds'] = {
                    'values': (log.times[scope] - self.this_source.utc).to_value('sec'),
                    'label': 'Seconds',
                    'reference': self.this_source.utc
                }
                self.taxes[scope]['boresight'] = {
                    'values': np.interp(self.taxes[scope]['seconds']['values'], self.xp, self.yp),
                    'label': 'Degrees',
                    'reference': self.this_source.utc
                }
        
    def get_wf_ticks(self, dat, ticks=8, precision=-1, include_0=False):
        if isinstance(dat[0], datetime):
            dat = (self.times.jd - self.times[0].jd) * 24.0 * 3600.0
        idat = list(np.arange(len(dat)))
        if isinstance(ticks, (float, int)):
            xstart = np.round(np.floor(dat[0]), precision)
            xstop = np.round(np.ceil(dat[-1]), precision)
            xstep = np.round((dat[-1] - dat[0]) / ticks, precision)
            x = [int(xx) for xx in np.arange(xstart, xstop, xstep)]
        elif isinstance(ticks, list):
            if len(ticks) == 3:
                x = [int(xx) for xx in np.arange(ticks[0], ticks[1]+1, ticks[2])]
            else:
                x = ticks
        else:
            raise ValueError("Invalid ticks parameters")
        if include_0 and 0 not in x:
            x.append(0)
            x = sorted(x)
        m = np.round(np.interp(x, dat, idat), 0)
        return m, x

    def dashboard_gen(self, script_fn='dash.sh', ants='2b,4e', pols='xx,xy', taxis='b', show_diff=False):
        """
        This generates a bash script to run to generate the dashboard plot for the obsid.
    
        Parameters
        ----------
        oinput : str
            oinput to use to find obsid etc...
        script_fn : str
            Name of bash script file to write
        ant : str
            ants to use
        pol : str
            pol to use [xx,yy,xy,yz]
        taxis : str
            t/x axis to use in plot [a/b/d]
    
        """
        cnode = ','.join(self.cnode)
        ants = tools.listify(ants)
        pols = tools.listify(pols, {'all': ['xx', 'xy', 'yy', 'yx']})
        show_diff = ' --show_diff ' if show_diff else ''
        with open(script_fn, 'w') as fp:
            for src in self.obs.sources:
                for ant in ants:
                    for pol in pols:
                        print(f"on_look.py {src} -a {ant}  -p {pol} -t {taxis} --lo {self.lo} --cnode {cnode} {show_diff} --dash -s", file=fp)

    def dashboard(self, ant='2a', pol='xx', time_axis='seconds', transit_time=4.0, **kwargs):
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
        self.read_in_files()
        if self.obs is None:
            print("No obsinfo found -- limited options.")
        else:
            self.this_source = self.obs.sources[self.source]
        print(f"Dashboard for {self.obsid} - ({ant},{pol})")
        D = {'use_dB': True, 'save': False, 't_wfticks': 8, 'f_wfticks': 8, 'log': False,
             'zoom_time': False, 'zoom_freq': False, 'filter_time': {}, 'show_diff': False}
        D.update(kwargs)

        if D['use_dB'] == 'auto':
            D['use_dB'] = True if pol in ['xx', 'yy'] else False
        if D['log']:
            logs = onv_logs.TBALog(D['log'])
        else:
            logs = None

        self.time_axis = on_sys.AXIS_OPTIONS[time_axis[0].lower()]
        self.get_time_axes(log=logs)
        if ant:  # Otherwise use existing
            self.get_bl(ant, pol=pol)
        if self.data is None:
            return

        # Set up filters
        self.filters = {}
        tt = transit_time / 2.0
        self.filters['on:boresight'] = Filter(ftype='time', unit=self.time_axis, lo=-tt, hi=tt, norm=True, color='r')
        self.filters['off:boresight'] = Filter(ftype='time', unit=self.time_axis, lo=-tt, hi=tt, norm=True, color='k', invert=True)
        for clr, filt in D['filter_time'].items():
            self.filters[f"time:{filt[0]}-{filt[1]}:{clr}"] = Filter(ftype='time', unit=self.time_axis, lo=filt[0], hi=filt[1], norm=True, color=clr)
        for clr, filt in self.obs.filters[self.lo].items():
            self.filters[f"freq:{filt[0]}-{filt[1]}:{clr}"] = Filter(color=clr, ftype='freq', unit='MHz', lo=filt[0], hi=filt[1])

        plt.figure('Dashboard', figsize=(16, 9))
        #, gridspec_kw={'width_ratios': [3, 1]}
        self.suptitle = f"{self.obsid}: {self.this_source.utc.datetime.isoformat(timespec='seconds')}"
        self.suptitle += f'  --  ({self.a},{self.b}) {self.pol}'
        plt.suptitle(self.suptitle)
        axwf = plt.subplot2grid((2, 2), (0, 0), rowspan=2, colspan=1)
        axt = plt.subplot2grid((2, 2), (0, 1), rowspan=1, colspan=1)
        axf = plt.subplot2grid((2, 2), (1, 1), rowspan=1, colspan=1)

        # Water fall plot
        f_wfticks, f_wftick_labels = self.get_wf_ticks(self.freqs, ticks=D['f_wfticks'], include_0=False)
        include_0 = True if self.time_axis in ['seconds', 'boresight'] else False
        t_wfticks, t_wftick_labels = self.get_wf_ticks(self.taxes[self.time_axis]['values'], ticks=D['t_wfticks'], include_0=include_0)
        axwf.imshow(toMag(self.data, D['use_dB']))
        axwf.set_aspect('auto')
        axwf.set_xlabel('Freq')
        axwf.set_ylabel(self.taxes[self.time_axis]['label'])
        axwf.set_xticks(f_wfticks, f_wftick_labels)
        axwf.set_yticks(t_wfticks, t_wftick_labels)

        # Time plot
        for key in [k for k in self.filters if self.filters[k].ftype=='freq']:
            if self.filters[key].apply(self.freqs, self.data):
                axt.plot(self.taxes[self.time_axis]['values'], toMag(self.filters[key].power, D['use_dB']), color=self.filters[key].color)
        axt.set_xlabel(self.taxes[self.time_axis]['label'])
        ylabel = 'dB' if D['use_dB'] else 'linear'
        axt.set_ylabel(ylabel)
        axtlim = copy(axt.axis())
        text_offset = 0.025 * (axtlim[3] - axtlim[2])
        if D['log']:
            lstyle = {'inner': 'r--', 'outer': 'b:'}
            for scope in ['inner', 'outer']:
                for i, t in enumerate(self.taxes[scope][time_axis]['values']):
                    if t > self.taxes[self.time_axis]['values'][0] and t < self.taxes[self.time_axis]['values'][-1]:
                        axt.plot([t, t], [axtlim[2], axtlim[3]], lstyle[scope])
                        axt.text(t, text_offset + axtlim[3], logs.sats[scope][i], rotation='vertical', fontsize=8)
        axt.axes.set_xlim(axtlim[0], axtlim[1])
        axt.axes.set_ylim(axtlim[2], axtlim[3])

        # Frequency plot
        for i in range(len(self.times)):
            axf.plot(self.freqs, toMag(self.data[i], D['use_dB']), '0.8')
        axf.set_xlabel(self.freq_unit)
        for key in [k for k in self.filters if self.filters[k].ftype=='time']:
            if self.filters[key].apply(self.taxes[self.time_axis]['values'], self.data):
                axf.plot(self.freqs, toMag(self.filters[key].power, D['use_dB']), color=self.filters[key].color)
        axf.set_ylabel(ylabel)
        axflim = copy(axf.axis())

        # Indicate filters and set limits
        # ... time
        axtshade = axtlim[2] + 0.15 * (axtlim[3] - axtlim[2])
        for key in [k for k in self.filters if self.filters[k].ftype=='time' and self.filters[k].use]:
            for xlim in self.filters[key].indicator:
                axt.fill_between(xlim, [axtshade, axtshade], color=self.filters[key].color)
        axt.set_ylim(bottom=axtlim[2], top=axtlim[3])
        if D['zoom_time']:
            axt.set_xlim(left=D['zoom_time'][0], right=D['zoom_time'][1])
        else:
            axt.set_xlim(left=self.taxes[self.time_axis]['values'][0], right=self.taxes[self.time_axis]['values'][-1])
        # ... freq
        axfshade = axflim[2] + 0.15 * (axflim[3] - axflim[2])
        for key in [k for k in self.filters if self.filters[k].ftype=='freq' and self.filters[k].use]:
            for xlim in self.filters[key].indicator:
                axf.fill_between(xlim, [axfshade, axfshade], color=self.filters[key].color)
        axf.set_ylim(bottom=axflim[2], top=axflim[3])
        if D['zoom_freq']:
            axf.set_xlim(left=D['zoom_freq'][0], right=D['zoom_freq'][1])
        else:
            axf.set_xlim(left=self.freqs[0], right=self.freqs[-1])

        fn = f"{self.obsid}_{self.lo}_{ant}_{pol}.png"
        if D['save']:
            plt.savefig(fn)

        if D['show_diff']:
            self.plot_diff(use_dB=D['use_dB'], save=D['save'], fn=f'Diff_{fn}')

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