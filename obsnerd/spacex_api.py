from astropy import coordinates as coord
from astropy import units as u
from astropy.time import Time, TimeDelta
from datetime import datetime
import requests
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from argparse import Namespace


BASE_URL = "https://api.starlink.com/public-files/ephemerides/"
MANIFEST = 'MANIFEST.txt'
ATA = coord.EarthLocation(lat=40.817431*u.deg, lon=-121.470736*u.deg, height=1019*u.m)
CELESTRAK = 'http://celestrak.org/NORAD/elements/gp.php?GROUP=active&FORMAT=tle'
FZERO = 1E-9


class Track:
    def __init__(self, satname, elmin):
        self.satname = satname
        self.elmin = elmin

    def update(self, state, ind, azdot, t, az, el):
        if state == 'start':
            self.start = {'ind': ind, 'azdot': azdot, 't': t, 'az': az, 'el': el}
        elif state == 'stop':
            self.stop = {'ind': ind, 'azdot': azdot, 't': t, 'az': az, 'el': el}
        elif state == 'peak':
            self.peak = {'ind': ind, 'azdot': azdot, 't': t, 'az': az, 'el': el}

    def set_duration(self):
        self.duration = self.stop['t'] - self.start['t']
        self.midpoint = self.start['t'] + self.duration / 2.0
        self.maxazdot = self.peak['azdot']


class SpaceX:
    def __init__(self, base_url=BASE_URL, manifest_file=MANIFEST, elmin=15.0, loc=ATA):
        self.base_url = base_url
        self.manifest = manifest_file
        self.loc = loc
        manifest_url = self.base_url + self.manifest
        print(f"Reading {manifest_url}")
        try:
            xxx = requests.get(manifest_url)
        except Exception as e:
            print(f"Error reading {manifest_url}:  {e}")
        self.manifest = xxx.text.splitlines()
        try:
            xxx = requests.get(CELESTRAK)
        except Exception as e:
            print(f"Error reading {CELESTRAK}")
        self.celestrak = xxx.text.splitlines()
        self.elmin = elmin * u.deg
        self.dtc = []
        for line in self.celestrak:
            if 'starlink' in line.lower() and 'dtc' in line.lower():
                self.dtc.append(line.split('[')[0].strip())

    def get_satellite(self, satname):
        self.satname = satname
        for entry in self.manifest:
            if satname in entry:
                ephem_url = self.base_url + entry
                self.satinfo = entry
                break

        print(f"Getting {ephem_url}")
        try:
            xxx = requests.get(ephem_url)
        except Exception as e:
            print(f"Error reading satellite")
        data = xxx.text.splitlines()

        self.meta, self.times, self.r_eci, self.v_eci, self.df = parse_ephemeris(data)

        self.times = Time(self.times)
        self.pst = self.times - TimeDelta(8 * 3600, format='sec')
        satpath = coord.CartesianRepresentation(x=self.df['r_ecix'],
                                                y=self.df['r_eciy'],
                                                z=self.df['r_eciz'])
        self.gcrs = coord.GCRS(satpath, obstime=self.times)
        self.azel = coord.SkyCoord(self.gcrs.ra, self.gcrs.dec).transform_to(coord.AltAz(location=self.loc, obstime=self.times))
        dt = np.diff(self.pst.jd) * 24 * 3600
        self.daz = self.azel.az.diff().to_value('deg')
        wrap = np.where(abs(self.daz) > 180.0)
        dwrap = wrap[0] - 1
        self.daz[wrap] = self.daz[dwrap]
        self.azdot = self.daz / dt
        self.azdot = np.insert(self.azdot, 0, self.azdot[0])
        self.azdot = abs(self.azdot)
        self.eldot = self.azel.alt.diff().to_value('deg') / dt
        self.eldot = np.insert(self.eldot, 0, self.eldot[0])
        #self.eldot = abs(self.eldot)
        notviewable = np.where(self.azel.alt < self.elmin)
        self.azdot[notviewable] = 0.0
        self.eldot[notviewable] = 0.0

    def get_tracks(self):
        self.tracks = []
        in_track = self.azdot[0] > 0.01
        if in_track:
            print("STARTING IN A TRACK -- HANDLE IT!!!")
            print("CAN'T YET SO SKIPPING")
            return
        this_track = Track(satname=self.satname, elmin=self.elmin)
        for i, azdot in enumerate(self.azdot):
            if not in_track and azdot > FZERO:
                this_track.update('start', ind=i, azdot=azdot, t=self.pst[i], az=self.azel.az[i], el=self.azel.alt[i])
                in_track = True
                peak = {'value': 0.0, 'index': i}
            elif in_track and azdot < FZERO:
                this_track.update('stop', ind=i-1, azdot=azdot, t=self.pst[i-1], az=self.azel.az[i-1], el=self.azel.alt[i-1])
                this_track.update('peak', ind=peak['index'], azdot=peak['value'], t=self.pst[peak['index']],
                                  az=self.azel.az[peak['index']], el=self.azel.alt[peak['index']])
                in_track = False
                this_track.set_duration()
                self.tracks.append(this_track)
                this_track = Track(satname=self.satname, elmin=self.elmin)
            if in_track:
                if azdot > peak['value']:
                    peak = {'value': azdot, 'index': i}

    def plot_track_summary(self, par='el'):
        self.par = par
        for i, trk in enumerate(self.tracks):
            # plt.plot([trk.start['t'].datetime, trk.peak['t'].datetime, trk.stop['t'].datetime],
            #          [trk.start[par].value, trk.peak[par].value, trk.stop[par].value], 'r-')
            if par == 'el':
                plt.plot(self.pst[trk.start['ind']:trk.stop['ind']+1].datetime, self.azel.alt[trk.start['ind']:trk.stop['ind']+1], 'k')
            elif par == 'az':
                plt.plot(self.pst[trk.start['ind']:trk.stop['ind']+1].datetime, self.azel.az[trk.start['ind']:trk.stop['ind']+1], 'k')
            plt.plot(trk.peak['t'].datetime, trk.peak[par].value, 'c.')
            plt.text(trk.peak['t'].datetime, trk.peak[par].value, f"{i}")
            # plt.plot(trk.start['t'].datetime, trk.start[par].value, 'g|')
            # plt.plot(trk.stop['t'].datetime, trk.stop[par].value, 'b|')
        axlim = plt.axis()
        now = Time.now() - TimeDelta(8 * 3600, format='sec')
        plt.plot([now.datetime, now.datetime], [axlim[2], axlim[3]], 'g--')
        plt.axis(ymin=axlim[2], ymax=axlim[3])

    def choose_tracks(self, ctrks):
        self.chose = [[int(x[:-1]), x[-1]] for x in ctrks.split(',')]
        for this in self.chose:
            if this[1] == 'l':
                ind = (self.tracks[this[0]].start['ind'] + self.tracks[this[0]].peak['ind']) // 2
            elif this[1] == 'm' or this[1] == 'p':
                ind = self.tracks[this[0]].peak['ind']
            elif this[1] == 'r':
                ind = (self.tracks[this[0]].peak['ind'] + self.tracks[this[0]].stop['ind']) // 2
            this.append(ind)
            par = 'alt' if self.par == 'el' else 'az'
            plt.plot(self.pst[ind].datetime, getattr(self.azel, par)[ind], 'rx')

    def proc_tracks(self, obslen_min=8, freqs=[1990.0, 5500.0], freq_unit='MHz'):
        from odsutils import ods_engine
        self.obslen_min = obslen_min
        self.freqs = [(f * u.Unit(freq_unit)).to_value('Hz') for f in freqs]
        self.new_tracks = []
        for this in self.chose:
            newstart = self.times[this[2]] - TimeDelta(self.obslen_min * 60.0 / 2.0, format='sec')
            newstop = self.times[this[2]] + TimeDelta(self.obslen_min * 60.0 / 2.0, format='sec')
            self.new_tracks.append({'src_id':f"S{self.tracks[this[0]].satname}",
                                    'src_ra_j2000_deg': self.gcrs.ra[this[2]].to_value('deg'),
                                    'src_dec_j2000_deg': self.gcrs.dec[this[2]].to_value('deg'),
                                    'src_start_utc': f"{newstart.datetime.isoformat(timespec='seconds')}",
                                    'src_end_utc': f"{newstop.datetime.isoformat(timespec='seconds')}",
                                    'freq_lower_hz': self.freqs[0]-10.0E6,
                                    'freq_upper_hz': self.freqs[0]+10.0E6})
            self.new_tracks.append({'src_id':f"S{self.tracks[this[0]].satname}",
                                    'src_ra_j2000_deg': self.gcrs.ra[this[2]].to_value('deg'),
                                    'src_dec_j2000_deg': self.gcrs.dec[this[2]].to_value('deg'),
                                    'src_start_utc': f"{newstart.datetime.isoformat(timespec='seconds')}",
                                    'src_end_utc': f"{newstop.datetime.isoformat(timespec='seconds')}",
                                    'freq_lower_hz': self.freqs[1]-10.0E6,
                                    'freq_upper_hz': self.freqs[1]+10.0E6})
        with ods_engine.ODS() as ods:
            ods.pipe(self.new_tracks, intake='ods.json', defaults='defaults.json')

    def horizon(self):
        self.up = np.where(self.azel.alt > self.elmin)
        plt.plot(self.pst[self.up].datetime, self.azel.alt[self.up], 'k.')
        plt.plot(self.pst[self.up].datetime, self.azel.alt[self.up], 'k-')
        plt.axis(ymin=self.elmin.to_value('deg')+5.0)

    def bounds(self, hr_start=11, hr_stop=16, el_min=40.0, el_max=60.0):
        self.view = Namespace(pst=[], az=[], el=[], ra=[], dec=[])
        start_day = self.pst[0].value.day
        dt = np.diff(self.pst.jd) * 24 * 3600
        radotcos = (self.gcrs.ra.diff().to_value('deg') / dt) * np.cos(self.gcrs.dec.to_value('rad')[1:])
        decdot = self.gcrs.dec.diff().to_value('deg') / dt
        omega = np.sqrt(radotcos**2 + decdot**2)
        for i in range(len(self.pst)):
            hr = self.pst[i].value.hour
            el = self.azel.alt[i].value
            if self.pst[i].value.day != start_day:
                break
            if hr_start < hr < hr_stop and el_min < el < el_max:
                self.view.pst.append(self.pst[i].datetime)
                self.view.az.append(self.azel.az[i].value)
                self.view.el.append(el)
                self.view.ra.append(self.gcrs.ra[i].value)
                self.view.dec.append(self.gcrs.dec[i].value)     


def parse_ephemeris_metadata(data_input):
    """Parse an ephemeris file's metadata.

    Parameters
    ----------
    data_input : str
        The name of an ephemeris file.

    Returns
    -------
    meta : dict
        The parsed metadata of the ephemeris file.

    Example
    -------
    This text:
    ::
      created:2019-12-30T12:26:33Z
      ephemeris_start:2019-12-30T12:25:18Z ephemeris_stop:2019-12-31T20:25:18Z step_size:60.000000
      ephemeris_source:jspoc_ephemeris_source_propagation

    Will be parsed to:
    :: 
        {'created': datetime(2019, 12, 30, 12, 26, 33),
         'ephemeris_start': datetime(2019, 12, 30, 12, 25, 18),
         'ephemeris_stop':datetime(2019, 12, 31, 20, 25, 18),
         'step_size': 60,
         'ephemeris_source': 'jspoc_ephemeris_source_propagation'}
    """

    # Loop over each of the header lines and store data.
    meta = {}
    if isinstance(data_input, str):
        f = open(data_input,'r')
        line1 = f.readline()
        line2 = f.readline()
        line3 = f.readline()
    elif isinstance(data_input, list):
        line1 = data_input[0]
        line2 = data_input[1]
        line3 = data_input[2]

    try:
        meta[line1[:line1.find(':')]] = datetime.strptime(line1[line1.find(':')+1:-5], '%Y-%m-%d %H:%M:%S')
    except:
        meta['line1'] = line1
    
    try:
        meta[line2[:line2.find(':')]] = datetime.strptime(line2[line2.find(':')+1:-57], '%Y-%m-%d %H:%M:%S')
        meta[line2[40:][:line2[40:].find(':')]] = datetime.strptime(line2[40:][line2[40:].find(':')+1:-18], '%Y-%m-%d %H:%M:%S')
        meta[line2[-13:-4]] = float(line2[-3:])
    except:
        meta['line2'] = line2
    
    try:
        meta[line3[:line3.find(':')]] = line3[line3.find(':')+1:-1]    
    except:
        meta['line3'] = line3
    return meta

def parse_ephemeris(data_input):
    """Parse an ephemeris file.

    Parameters
    ----------
    data_input : str
        The name of an ephemeris file or data

    Returns
    -------
    meta : dict
        The parsed metadata of the ephemeris file.
    times : (N,) Numpy array
        An array of datetime objects indicating satellite state times.
    r_eci : (N,3) Numpy array [meters]
        Inertial positions from the file converted from km to meters.
    v_eci : (N,3) Numpy array [meters/second]
        Inertial velocities from the file converted from km/second to meters/second.
    p_rtn : (N,6,6) Numpy array [SI units]
        The RTN covariance values from the file converted
        from km to meters (squared, per second, etc.)

    Example
    -------
    This text:
    ::
      created:2019-12-30T12:26:33Z
      ephemeris_start:2019-12-30T12:25:18Z ephemeris_stop:2019-12-31T20:25:18Z step_size:60.000000
      ephemeris_source:jspoc_ephemeris_source_propagation
      UVW
      2019364122518.361 4094.0491381513 2859.6314019628 4763.2190701632 -1.2518571159 6.8576664670 -3.0322647905
      1.8599935273e-06 -6.6214191608e-06 3.5110389419e-05 3.9248988820e-09 -2.0852272834e-08 6.7939209526e-07 6.7989595384e-09
      -3.5029804535e-08 2.1185410095e-11 3.5089242659e-11 -1.4734126601e-09 4.8086876808e-09 -3.8209011936e-12 -4.9580509074e-12
      1.1904290648e-12 9.8948088922e-12 -5.8457640564e-11 4.6761468674e-10 5.2675521395e-14 -6.2503425058e-15 6.1090786999e-12

    Will have the first four lines parsed as metadata. The next line will be parsed to a
    datetime, a position, and a velocity. The next three lines will be converted from lower
    triangle values to a symmetric covariance matrix.
    """

    # Read all of the lines.
    if isinstance(data_input, str):
        with open(data_input, 'r') as f:
            lines = f.readlines()
    elif isinstance(data_input, list):
        lines = data_input

    # Parse metadata.
    data_input = lines[:3]
    meta = parse_ephemeris_metadata(data_input)

    # Initialize state lists.
    times = []
    r_eci = []
    v_eci = []
    p_rtn = []

    # Parse states.
    for ii, line in enumerate(lines[4:]):

        # Split the line into tokens.
        tokens = line.split(' ')

        # The first line in each batch of four lines contains a time, position, and velocity.
        # Convert from kilometers to meters.
        if ii % 4 == 0:
            try:
                times.append(datetime.strptime(tokens[0], '%Y%j%H%M%S.%f'))
                # times.append(datetime.strptime(tokens[0], '%Y%j%H%M%S.%f').replace(tzinfo = timezone.utc))
                r_eci.append(np.array(list(map(float, tokens[1:4]))) * 1e3)
                v_eci.append(np.array(list(map(float, tokens[4:7]))) * 1e3)
                cov = np.zeros((6, 6))
            except:
                z=5

        # The next three lines in each batch contain the lower triangle of the covariance matrix.
        # Convert from kilometers squared to meters squared.
        else:
            values = np.array(list(map(float, tokens))) * 1e6
            linear_inds = (7 * ((ii % 4) - 1)) + np.arange(7)
            row_inds, col_inds = np.array(np.tril_indices(6)).T[linear_inds].T
            cov[row_inds, col_inds] = values
            cov[col_inds, row_inds] = values
            if ii % 4 == 3:
                p_rtn.append(cov)

    # Turn things into arrays.
    times, r_eci, v_eci, p_rtn = map(np.array, (times, r_eci, v_eci, p_rtn))

    df = pd.DataFrame({'t_utc':times})
    df[['r_ecix','r_eciy','r_eciz']] = r_eci
    df[['v_ecix','v_eciy','v_eciz']] = v_eci

    return meta, times, r_eci, v_eci, df
