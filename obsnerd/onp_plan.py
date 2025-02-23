import logging
from aocalendar import aocalendar
from . import LOG_FILENAME, LOG_FORMATS, __version__
from odsutils import logger_setup, locations
from odsutils import ods_timetools as ttools
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.time import TimeDelta
import matplotlib.pyplot as plt
import json
import numpy as np


logger = logging.getLogger(__name__)
logger.setLevel('DEBUG')  # Set to lowest

import colorsys

def generate_colors(n):
    """
    Generate a list of N optimally spaced colors using HSL color space.
    
    Parameters:
    n (int): Number of colors to generate.
    
    Returns:
    list: List of RGB color tuples.
    """    
    colors = []
    for i in range(n):
        hue = i / n  # Evenly spaced hues around the color wheel
        rgb = colorsys.hsv_to_rgb(hue, 1, 1)  # Full saturation and value
        colors.append('#{:02x}{:02x}{:02x}'.format(*(int(c * 255) for c in rgb)))  # Convert to hex
    
    return colors


class Plan:
    def __init__(self, conlog='INFO', filelog=False, path='', loc='ata'):
        self.location = locations.Location(name=loc)
        self.log_settings = logger_setup.Logger(logger, conlog=conlog, filelog=filelog, log_filename=LOG_FILENAME, path=path,
                                                conlog_format=LOG_FORMATS['conlog_format'], filelog_format=LOG_FORMATS['filelog_format'])
        logger.info(f"{__name__} ver. {__version__}")
        self.minimum_duration = TimeDelta(0.0, format='sec')
        self.freqs = None
        self.tracks = {}
        self.bandwidth = None
        self.el_limit = None

    def setupcal(self):
        self.this_cal = aocalendar.Calendar(calfile='now', path=self.log_settings.path, conlog=self.log_settings.conlog,
                                            filelog=self.log_settings.filelog, start_new=True, loc=self.location)
        self.this_cal.get_current_time()
        lst = f"{int(self.this_cal.current_lst.hms.h):02d}:{int(self.this_cal.current_lst.hms.m):02d}:{int(self.this_cal.current_lst.hms.s):02d}"
        logger.info(f"Current LST:  {lst}")

    def test_obs(self, az, el, freq=[1990.0, 5990.0], start_in_sec=10*60, obs_len_sec=8 * 60):
        """
        Schedule a test observation.

        Parameters
        ----------
        az : float
            Azimuth in degrees
        el : float
            Elevation in degrees
        obs_len_min : float, optional
            Observation length in minutes, by default 8 * 60

        Returns
        -------
        bool
            True if observation is possible, False otherwise.
        """
        self.setupcal()
        self.tref = ttools.t_delta('now', start_in_sec + obs_len_sec/2.0, 's')
        self.start = ttools.t_delta('now', start_in_sec, 's')
        self.stop = ttools.t_delta(self.start, obs_len_sec, 's')
        obscoord = SkyCoord(alt = el*u.deg, az = az*u.deg, obstime = self.tref, frame = 'altaz', location = self.this_cal.location.loc)
        radec = obscoord.transform_to('gcrs')
        self.this_cal.ods.get_defaults_dict('defaults.json')
        self.this_cal.ods.add_new_record(src_id='Asource', src_ra_j2000_deg=radec.ra.deg, src_dec_j2000_deg=radec.dec.deg,
                                         src_start_utc=self.start, src_end_utc=self.stop,
                                         freq_lower_hz=freq[0] * 1E6, freq_upper_hz=freq[0] * 1E6)
        self.this_cal.ods.add_new_record(src_id='Asource', src_ra_j2000_deg=radec.ra.deg, src_dec_j2000_deg=radec.dec.deg,
                                         src_start_utc=self.start, src_end_utc=self.stop,
                                         freq_lower_hz=freq[1] * 1E6, freq_upper_hz=freq[1] * 1E6)
        self.this_cal.ods.view_ods()
        self.this_cal.ods.write_ods('test_ods.json')

    def get_tracks(self, satname, start, duration, freqs=[1990.0, 5990.0], bandwidth=100.0, freq_unit='MHz', el_limit=15.0,
                   DTC_only=True, time_resolution=10, source='sopp'):
        """
        Get satellite tracks for satellites regex-found from 'satname'.

        Parameters
        ----------
        satname : str
            Regex string to find satellites.
        start : ods_timetools time input
            Start time of the observation (UTC)
        duration : float
            Duration of the observation (minutes)
        freqs : list, optional
            Frequencies to observe, by default [1990.0, 5990.0] (currently ignored in search)
        bandwidth : float, optional
            Bandwidth of the observation, by default 100.0
        freq_unit : str, optional
            Frequency unit, by default 'MHz'
        el_limit : float, optional
            Minimum elevation limit in degrees, by default 15.0
        DTC_only : bool, optional
            If True, only DTC tracks are returned, by default True
        time_resolution : int, optional
            Time resolution in seconds, by default 10
        source : str, optional
            Source of satellite data, by default 'sopp'

        """
        if self.freqs is None:
            self.freqs = [f * u.Unit(freq_unit) for f in freqs]
        else:
            logger.warning(f"Using freqs: {self.freqs}")
        if self.bandwidth is None:
            self.bandwidth = bandwidth * u.Unit(freq_unit)
        else:
            logger.warning(f"Using bandwidth: {self.bandwidth}")
        if self.el_limit is None:
            self.el_limit = el_limit * u.deg
        else:
            logger.warning(f"Using el_limit: {self.el_limit}")
        if source == 'sopp':
            logger.info(f"SOPP not using frequency.")
            freqs = False
            from . import sopp_engine
            these_tracks = sopp_engine.main(satname=satname, start=start, duration=duration, frequency=freqs, el_limit=el_limit,
                                            DTC_only=DTC_only, time_resolution=time_resolution, verbose=False, show_plots=False)
            logger.info(f"Found the following {len(these_tracks)} satellites above {self.el_limit}: {', '.join(these_tracks.keys())}")
            if not(len(these_tracks)):
                logger.warning(f"No satellites found above {self.el_limit}")
                return
            self.tracks.update(these_tracks)
            self.plot_track_summary(satname='all')
        elif source == 'spacex':
            from . import spacex_api
            print("WORKING ON IT")
        else:
            logger.error(f"Unknown source: {source}")
            return

    def plot_track_summary(self, satname='all'):
        if satname == 'all':
            satname = list(self.tracks.keys())
        elif not isinstance(satname, list):
            satname = [satname]
        plt.figure('Track Summary')
        clr_list = generate_colors(len(satname))
        for j, sat in enumerate(satname):
            this_color = clr_list[j]
            labelled = False
            for i, track in enumerate(self.tracks[sat]):
                if track.duration < self.minimum_duration:
                    continue
                if not labelled:
                    plt.plot(track.utc.datetime, track.az, color=this_color, label=sat)
                    labelled = True
                else:
                    plt.plot(track.utc.datetime, track.az, color=this_color)
                plt.plot(track.utc.datetime, track.el, color=this_color)
                plt.plot(track.utc[track.imax].datetime, track.el[track.imax].value, 'c.')
                plt.text(track.utc[track.imax].datetime, track.el[track.imax].value, f"{i}")
                plt.plot(track.utc[track.imax].datetime, track.az[track.imax].value, 'c.')
                plt.text(track.utc[track.imax].datetime, track.az[track.imax].value, f"{i}")
            axlim = plt.axis()
            now = ttools.interpret_date('now')
        plt.plot([now.datetime, now.datetime], [axlim[2], axlim[3]], 'g--')
        plt.axis(ymin=axlim[2], ymax=axlim[3])
        plt.legend()

    def choose_tracks(self, obslen_min=8, minimum_duration_min=5):
        """
        Interactive chooser of tracks to use -- edits the Track instances.

        Parameters
        ----------
        obslen_min : float, optional
            Observation length in minutes, by default 8
        minimum_duration_min : float, optional
            Minimum duration above el_limit in minutes, by default 5

        """
        self.minimum_duration = TimeDelta(minimum_duration_min * 60.0, format='sec')
        self.obslen = TimeDelta(obslen_min * 60.0, format='sec')
        for sat in self.tracks:
            ctrks = input(f"Choose track/ODS(+/-) for {sat}, or '[s]kip' or '[e]nd' (e.g. 0+,1+,2-,3+): ")
            if ctrks[0].lower() == 'e':
                return
            elif ctrks[0].lower() == 's':
                continue
            this_sat = self.tracks[sat]
            for ch in [[int(x[:-1]), x[-1]] for x in ctrks.split(',')]:
                trk, do_ods = ch
                ooddss = [] if do_ods == '+' else None
                ind = this_sat[trk].imax
                this_sat[trk].set_par(iobs=ind, ods=ooddss)
                plt.plot(this_sat[trk].utc[ind].datetime, this_sat[trk].el[ind].value, 'r.')
                plt.plot(this_sat[trk].utc[ind].datetime, this_sat[trk].az[ind].value, 'r.')

    def proc_tracks(self, filter_file='filters.json'):
        """
        After choosing tracks, write the ODS file and the obsinfo.json file.

        Parameter
        ---------
        filter_file : str, optional
            Filter file to use, by default 'filters.json' (just passed through to write_obsinfo)

        """
        from odsutils import ods_engine
        obslen_TD2 = self.obslen / 2.0
        new_tracks = []
        for satname in self.tracks:
            for track in self.tracks[satname]:
                if track.iobs is None:
                    continue
                tobs =  track.utc[track.iobs]
                tstart = tobs - obslen_TD2
                if tstart < track.utc[0]:
                    istart = 0
                else:
                    istart = np.where(track.utc > tstart)[0][0] - 1
                tstop = tobs + obslen_TD2
                if tstop > track.utc[-1]:
                    istop = len(track.utc) - 1
                else:
                    istop = np.where(track.utc < tstop)[0][-1] + 1
                track.set_par(istart=istart, istop=istop, tobs=tobs, tstart=tstart, tstop=tstop)
                if track.ods is None:
                    continue
                for ff in self.freqs:
                    odict = {'src_id':f"{track.source}",
                            'src_ra_j2000_deg': track.ra[track.iobs].to_value('deg'),
                            'src_dec_j2000_deg': track.dec[track.iobs].to_value('deg'),
                            'src_start_utc': f"{tstart.datetime.isoformat(timespec='seconds')}",
                            'src_end_utc': f"{tstop.datetime.isoformat(timespec='seconds')}",
                            'freq_lower_hz': (ff - self.bandwidth / 2.0).to_value('Hz'),
                            'freq_upper_hz': (ff + self.bandwidth / 2.0).to_value('Hz')}
                    track.ods.append(odict)
                    new_tracks.append(odict)
        with ods_engine.ODS() as ods:
            ods.pipe(new_tracks, intake='ods.json', defaults='defaults.json')
        self.write_obsinfo(filter_file=filter_file)

    def write_obsinfo(self, filter_file='filters.json'):
        from astropy.coordinates import angular_separation
        from copy import copy
        try:
            with open(filter_file, 'r') as fp:
                filter = json.load(fp)['Filters']
        except FileNotFoundError:
            filter = {}
        obsinfo = {'dir_data': 'data', 'Filter': filter, "Sources": []}
        start_mjd = None
        for satname in self.tracks:
            for track in self.tracks[satname]:
                if track.iobs is None:
                    continue
                mjd = track.utc[track.istart].mjd
                if start_mjd is None or mjd < start_mjd:
                    start_mjd = copy(mjd)
                satname = f"{track.source}"
                obsinfo[satname] = {}
                obsinfo[satname]['ra'] = track.ra[track.iobs].to_value('deg')
                obsinfo[satname]['dec'] = track.dec[track.iobs].to_value('deg')
                obsinfo[satname]['utc'] = track.tobs.datetime.isoformat(timespec='seconds')
                obsinfo[satname]['az'] = track.az[track.iobs].to_value('deg')
                obsinfo[satname]['el'] = track.el[track.iobs].to_value('deg')
                obsinfo[satname]['off_time'] = []
                obsinfo[satname]['off_angle'] = []
                for i in range(track.istart, track.istop+1):
                    toff = track.utc[i] - track.tobs
                    aoff = angular_separation(track.ra[i], track.dec[i], track.ra[track.iobs], track.dec[track.iobs]) * np.sign(toff)
                    obsinfo[satname]['off_time'].append(np.round(toff.to_value('sec'), 1))
                    obsinfo[satname]['off_angle'].append(np.round(aoff.to_value('deg'), 2))
        fnout = f"obsinfo_{start_mjd:.0f}.json"
        try:
            with open(fnout, 'r') as fp:
                self.obsinfo = json.load(fp)
            print(f"Updating {fnout}")
        except FileNotFoundError:
            self.obsinfo = {}
            print(f"Writing {fnout}")
        self.obsinfo.update(obsinfo)
        with open(fnout, 'w') as fp:
            json.dump(self.obsinfo, fp, indent=4)