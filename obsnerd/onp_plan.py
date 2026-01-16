import logging
from aocalendar import aocalendar
from . import LOG_FILENAME, LOG_FORMATS, __version__, on_sys, on_config
from odsutils import logger_setup, locations
from odsutils import ods_timetools as ttools
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.time import TimeDelta
import matplotlib.pyplot as plt
import json
from copy import copy
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
    def __init__(self, conlog='INFO', filelog=False, config_file=None, loc='ata'):
        self.location = locations.Location(name=loc)
        self.log_settings = logger_setup.Logger(logger, conlog=conlog, filelog=filelog, log_filename=LOG_FILENAME,
                                                conlog_format=LOG_FORMATS['conlog_format'], filelog_format=LOG_FORMATS['filelog_format'])
        logger.info(f"{__name__} ver. {__version__}")
        self.minimum_duration = TimeDelta(0.0, format='sec')
        self.el_limit = None
        self.tracks = {}
        self.config = on_config.Config(config_file=config_file)

    def setupcal(self):
        self.this_cal = aocalendar.Calendar(calfile='now', path=self.log_settings.path, conlog=self.log_settings.conlog,
                                            filelog=self.log_settings.filelog, start_new=True, loc=self.location)
        self.this_cal.get_current_time()
        lst = f"{int(self.this_cal.current_lst.hms.h):02d}:{int(self.this_cal.current_lst.hms.m):02d}:{int(self.this_cal.current_lst.hms.s):02d}"
        logger.info(f"Current LST:  {lst}")

    def test_obs(self, az, el, freq=[1990.0, 1999.0, 1981.0, 2400.0], start_in_sec=10*60, obs_len_sec=8 * 60):
        """
        Schedule a test observation.  HASN'T BEEN UPDATED YET -- DON't USE.

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
        self.this_cal.ods.get_defaults(self.default_ods_default_file)
        self.this_cal.ods.new_record(src_id='Asource', src_ra_j2000_deg=radec.ra.deg, src_dec_j2000_deg=radec.dec.deg,
                                     src_start_utc=self.start, src_end_utc=self.stop,
                                     freq_lower_hz=freq[0] * 1E6, freq_upper_hz=freq[0] * 1E6)
        self.this_cal.ods.new_record(src_id='Asource', src_ra_j2000_deg=radec.ra.deg, src_dec_j2000_deg=radec.dec.deg,
                                     src_start_utc=self.start, src_end_utc=self.stop,
                                     freq_lower_hz=freq[1] * 1E6, freq_upper_hz=freq[1] * 1E6)
        self.this_cal.ods.view_ods()
        self.this_cal.ods.post_ods('test_ods.json')

    def get_tracks(self, satname, start, duration, el_limit=15.0, DTC_only=True, time_resolution=10, source='sopp'):
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
        print("Now run `plan.choose_tracks(auto=True)` to select tracks for observation.")

    def plot_track_summary(self, satname='all', show_az_transit_track=False):
        line_styles = ['-', '--', '-.', ':']
        if satname == 'all':
            satname = list(self.tracks.keys())
        elif not isinstance(satname, list):
            satname = [satname]
        plt.figure('Track Summary')
        clr_list = generate_colors(len(satname))
        for j, sat in enumerate(satname):
            this_color = clr_list[j]
            this_ls = line_styles[j % len(line_styles)]
            labelled = False
            for i, track in enumerate(self.tracks[sat]):
                if track.duration < self.minimum_duration:
                    continue
                if not labelled:
                    plt.plot(track.utc.datetime, track.el, color=this_color, ls=this_ls, label=sat)
                    labelled = True
                else:
                    plt.plot(track.utc.datetime, track.el, color=this_color, ls=this_ls)
                plt.plot(track.utc[track.imax].datetime, track.el[track.imax].value, 'c.')
                plt.text(track.utc[track.imax].datetime, track.el[track.imax].value, f"{i}")
                if show_az_transit_track:
                    plt.plot(track.utc.datetime, track.az, color=this_color, ls='--')
                    plt.plot(track.utc[track.imax].datetime, track.az[track.imax].value, 'c.')
                    plt.text(track.utc[track.imax].datetime, track.az[track.imax].value, f"{i}")
                else:
                    plt.text(track.utc[track.imax].datetime, self.el_limit.value, f"{track.az[track.imax].value:.0f}")
            axlim = plt.axis()
            now = ttools.interpret_date('now')
        plt.plot([now.datetime, now.datetime], [axlim[2], axlim[3]], 'g--')
        plt.axis(ymin=axlim[2], ymax=axlim[3])
        plt.legend()

    def choose_tracks(self, auto=True, obslen_min=8, obsel_deg=25, keyhole_deg=84.0):
        """
        Interactive chooser of tracks to use -- edits the Track instances.

        Parameters
        ----------
        auto : bool
            If True, automatically choose tracks, by default False
        obslen_min : float
            Observation length in minutes, by default 8
        obsel_deg : float
            Minimum el_limit to show in degrees, by default 25
        keyhole_deg : float
            Keyhole angle in degrees, by default 84.0

        """
        self.auto = auto
        self.obslen = TimeDelta(obslen_min * 60.0, format='sec')
        if auto:
            if isinstance(auto, bool):
                track_separation = 1.75 * self.obslen
            else:
                track_separation = TimeDelta(auto * 60.0, format='sec')
            print(f"Automatically choosing tracks at {track_separation.to_value('sec') / 60.0:.1f} min separation.")
        self.track_list = {}
        now = ttools.interpret_date('now')
        for sat in self.tracks:
            for i, track in enumerate(self.tracks[sat]):
                dt = track.utc[track.imax] - now
                key = int(dt.to_value('sec'))
                self.track_list[key] = {"sat": sat, "track": i, "use": 's'}
                track.ptadd(iobs=None, use='s')
        if not auto:
            print("Choose for the following satellite tracks:")
            print("\ty - use and implement ods record")
            print("\tn - use but don't implement ods record")
            print("\ts - skip track")
            print("\te - end choosing tracks\n")
        last_one = ttools.interpret_date('yesterday')
        for key in sorted(self.track_list.keys()):
            track = self.tracks[self.track_list[key]['sat']][self.track_list[key]['track']]
            if track.el[track.imax].to_value('deg') < obsel_deg or track.el[track.imax].to_value('deg') > keyhole_deg:
                continue
            ind = track.imax
            if auto:
                if track.utc[ind] - last_one > track_separation:
                    self.track_list[key]['use'] = input(f"ODS for {track.source} at {track.utc[ind].datetime.strftime('%m-%d %H:%M')} ({key / 60.0:.0f}m) -- {track.el[ind].to_value('deg'):.0f}\u00b0 (y/n)? ")
            else:
                self.track_list[key]['use'] = input(f"{sat}/{i} -- {track.el[ind].to_value('deg'):.0f}\u00b0 @ {track.utc[ind].datetime.strftime('%m-%d %H:%M')} ({key / 60.0:.0f}m) (y/n/s/e):  ")
            if self.track_list[key]['use'] == 'e':
                print("Ending track selection.")
                self.track_list[key]['use'] = 's'
                break
            elif self.track_list[key]['use'] in ['y', 'n']:
                last_one = copy(track.utc[ind])
                track.ptadd(iobs=ind, use=self.track_list[key]['use'])
                plt.plot(track.utc[ind].datetime, track.el[ind].value, 'ro')
                print(f"\t{track.use_def[track.use]} {track.source} at {track.utc[ind].datetime.strftime('%m-%d %H:%M')} ({key / 60.0:.0f}m) -- {track.el[ind].to_value('deg'):.0f}\u00b0")
        print("Now run plan.proc_tracks() to write the obsinfo.json file.")

    def proc_tracks(self, decimal_places=1):
        """
        After choosing tracks, write the obsinfo.json file.

        Parameter
        ---------
        decimal_places : int, optional
            Number of decimal places to use for filenames, by default 1

        """
        obslen_TD2 = self.obslen / 2.0
        self.start_mjd = None
        for key in self.track_list:
            track = self.tracks[self.track_list[key]['sat']][self.track_list[key]['track']]
            if self.track_list[key]['use'] == 's':
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
            track.ptadd(istart=istart, istop=istop, tobs=tobs, tstart=tstart, tstop=tstop)
            mjd = track.utc[track.istart].mjd
            if self.start_mjd is None or mjd < self.start_mjd:
                self.start_mjd = copy(mjd)
        obsinfofn = on_sys.make_obsinfo_filename(self.start_mjd, decimal_places=decimal_places)
        self.write_obsinfo(obsinfofn)
        print(f"Copy the obsinfo file:  scp {obsinfofn} sonata@obs-node1.hcro.org:rfsoc_obs_scripts/p054/obsinfo_rados.json")

    def write_obsinfo(self, obsinfofn):
        from astropy.coordinates import angular_separation
        obsinfo = {'filename': obsinfofn, 'dir_data': 'data', "Sources": {}}
        for k, v in self.config.items().items():
            obsinfo[k] = v
        for satname in self.tracks:
            for track in self.tracks[satname]:
                if track.iobs is None:
                    continue
                obsinfo['Sources'][track.source] = {'obsid': on_sys.make_obsid(track.source, track.utc[track.istart].mjd)}
                obsinfo['Sources'][track.source]['ra'] = track.ra[track.iobs].to_value('deg')
                obsinfo['Sources'][track.source]['dec'] = track.dec[track.iobs].to_value('deg')
                obsinfo['Sources'][track.source]['utc'] = track.tobs.datetime.isoformat(timespec='seconds')
                obsinfo['Sources'][track.source]['az'] = track.az[track.iobs].to_value('deg')
                obsinfo['Sources'][track.source]['el'] = track.el[track.iobs].to_value('deg')
                obsinfo['Sources'][track.source]['start'] = track.tstart.datetime.isoformat(timespec='seconds')
                obsinfo['Sources'][track.source]['stop'] = track.tstop.datetime.isoformat(timespec='seconds')
                obsinfo['Sources'][track.source]['off_time'] = []
                obsinfo['Sources'][track.source]['off_angle'] = []
                obsinfo['Sources'][track.source]['ods'] = True if track.use == 'y' else False
                for i in range(track.istart, track.istop+1):
                    toff = (track.utc[i] - track.tobs).to_value('sec')
                    aoff = angular_separation(track.ra[i], track.dec[i], track.ra[track.iobs], track.dec[track.iobs]) * np.sign(toff)
                    obsinfo['Sources'][track.source]['off_time'].append(np.round(toff, 1))
                    obsinfo['Sources'][track.source]['off_angle'].append(np.round(aoff.to_value('deg'), 2))
        try:
            with open(obsinfofn, 'r') as fp:
                self.obsinfo = json.load(fp)
            print(f"Updating {obsinfofn}")
        except FileNotFoundError:
            self.obsinfo = {}
            print(f"Writing {obsinfofn}")
        self.obsinfo.update(obsinfo)
        with open(obsinfofn, 'w') as fp:
            json.dump(self.obsinfo, fp, indent=4)