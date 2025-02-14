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


class Track:
    def __init__(self, satname, srcname=None):
        self.satname = satname
        self.srcname = satname if srcname is None else srcname
        self.iobs = None
        self.ods = None

    def set_par(self, **kwargs):
        for par, val in kwargs.items():
            if par not in ['iobs', 'istart', 'istop', 'tobs', 'tstart', 'tstop', 'ods']:
                logger.warning(f"Invalid parameter: {par}")
                continue
            setattr(self, par, val)

    def set_track(self, **kwargs):
        for par, val in kwargs.items():
            if par not in ['utc', 'ra', 'dec', 'az', 'el', 'dist']:
                logger.warning(f"Invalid parameter: {par}")
                continue
            setattr(self, par, val)

    def calc_properties(self):
        self.duration = self.utc[-1] - self.utc[0]
        self.imax = np.argmax(self.el)
        dt = np.diff(self.utc.mjd) * 24 * 3600
        self.daz = self.az.diff().to_value('deg')
        wrap = np.where(abs(self.daz) > 180.0)
        dwrap = wrap[0] - 1
        self.daz[wrap] = self.daz[dwrap]
        self.azdot = self.daz / dt
        self.azdot = np.insert(self.azdot, 0, self.azdot[0])
        self.azdot = abs(self.azdot)
        self.eldot = self.el.diff().to_value('deg') / dt
        self.eldot = np.insert(self.eldot, 0, self.eldot[0])


class Plan:
    def __init__(self, conlog='INFO', filelog=False, path='', loc='ata'):
        self.location = locations.Location(name=loc)
        self.log_settings = logger_setup.Logger(logger, conlog=conlog, filelog=filelog, log_filename=LOG_FILENAME, path=path,
                                                conlog_format=LOG_FORMATS['conlog_format'], filelog_format=LOG_FORMATS['filelog_format'])
        logger.info(f"{__name__} ver. {__version__}")

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

    def get_tracks(self, satname, start, duration, freqs=[1990.0, 5990.0], bandwidth=100.0, freq_unit='MHz', el_limit=15.0, source='sopp'):
        self.freqs = [f * u.Unit(freq_unit) for f in freqs]
        self.bandwidth = bandwidth * u.Unit(freq_unit)
        self.el_limit = el_limit * u.deg
        if source == 'sopp':
            logger.info(f"SOPP not using frequency.")
            freqs = False
            from . import sopp_engine
            self.tracks = sopp_engine.main(satname=satname, start=start, duration=duration, frequency=freqs, el_limit=el_limit,
                                           verbose=False, show_plots=False)
            logger.info(f"Found the following satellites above {el_limit}: {self.tracks.keys()}")
            self.satname = list(self.tracks.keys())[0]
            logger.info(f"Using: {self.satname}")
        elif source == 'spacex':
            from . import spacex_api
            print("WORKING ON IT")
        else:
            logger.error(f"Unknown source: {source}")
            return
        self.plot_track_summary()

    def plot_track_summary(self):
        for i, trk in enumerate(self.tracks[self.satname]):
            plt.plot(trk.utc.datetime, trk.az, 'k')
            plt.plot(trk.utc.datetime, trk.el, 'k')
            plt.plot(trk.utc[trk.imax].datetime, trk.el[trk.imax].value, 'c.')
            plt.text(trk.utc[trk.imax].datetime, trk.el[trk.imax].value, f"{i}")
            plt.plot(trk.utc[trk.imax].datetime, trk.az[trk.imax].value, 'c.')
            plt.text(trk.utc[trk.imax].datetime, trk.az[trk.imax].value, f"{i}")
        axlim = plt.axis()
        now = ttools.interpret_date('now')
        plt.plot([now.datetime, now.datetime], [axlim[2], axlim[3]], 'g--')
        plt.axis(ymin=axlim[2], ymax=axlim[3])

    def choose_tracks(self, ctrks, obslen_min=8):
        self.obslen = TimeDelta(obslen_min * 60.0, format='sec')
        this_sat = self.tracks[self.satname]
        for ch in [[int(x[:-1]), x[-1]] for x in ctrks.split(',')]:
            this_duration = min(self.obslen, this_sat[trk].duration)
            Ni = int((this_duration / 2.0) / (this_sat[0].utc[1] - this_sat[0].utc[0])) + 1
            trk, pos = ch
            if pos == 'l':
                ind = Ni
            if pos == 'm' or pos == 'p':
                ind = this_sat[trk].imax
            elif pos == 'l':
                ind = min(Ni, len(this_sat[trk].utc) - Ni)
            elif pos == 'r':
                ind = max(Ni, len(this_sat[trk].utc) - Ni)  
            this_sat[trk].set_par(iobs=ind, ods=[])
            plt.plot(this_sat[trk].utc[ind].datetime, this_sat[trk].el[ind].value, 'r.')
            plt.plot(this_sat[trk].utc[ind].datetime, this_sat[trk].az[ind].value, 'r.')

    def proc_tracks(self):
        from odsutils import ods_engine
        obslen_TD2 = self.obslen / 2.0
        new_tracks = []
        for track in self.tracks[self]:
            if track.iobs is None:
                continue
            tobs =  track.utc[track.iobs]
            tstart = tobs - obslen_TD2
            tstop = tobs + obslen_TD2
            istart = np.where(track.utc > tstart)[0][0] - 1
            istop = np.where(track.utc < tstop)[0][-1] + 1
            track.set_par(istart=istart, istop=istop, tobs=tobs, tstart=tstart, tstop=tstop)
            for ff in self.freqs:
                odict = {'src_id':f"{track.srcname}",
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

    def write_obsinfo(self, filter_file='filters.json'):
        from astropy.coordinates import angular_separation
        from copy import copy
        this_sat = self.tracks[self.satname]
        with open(filter_file, 'r') as fp:
            filter = json.load(fp)['Filters']
        obsinfo = {}
        start_mjd = None
        for track in this_sat:
            if track.iobs is None:
                continue
            mjd = track.utc[track.istart]
            if start_mjd is None:
                start_mjd = copy(mjd)
            satname = f"{track.srcname}"
            obsinfo[satname] = {}
            obsinfo[satname]['ra'] = track.ra[track.iobs].to_value('deg')
            obsinfo[satname]['dec'] = track.dec[track.iobs].to_value('deg')
            obsinfo[satname]['tref'] = track.tobs.datetime.isoformat(timespec='seconds')
            obsinfo[satname]['az'] = track.az[track.iobs].to_value('deg')
            obsinfo[satname]['el'] = track.el[track.iobs].to_value('deg')
            obsinfo[satname]['off_time'] = []
            obsinfo[satname]['off_angle'] = []
            for i in range(track.istart, track.istop+1):
                dt = track.utc[i] - track.tobs
                da = angular_separation(track.ra[i], track.dec[i], track.ra[track.iobs], track.dec[track.iobs]) * np.sign(dt)
                obsinfo[satname]['off_time'].append(np.round(dt.to_value('sec'), 1))
                obsinfo[satname]['off_angle'].append(np.round(da.to_value('deg'), 2))
        fnout = f"obsinfo_{start_mjd.value:.0f}.json"
        with open(fnout, 'w') as fp:
            json.dump(obsinfo, fp, indent=4)
        self.obsinfo = obsinfo
