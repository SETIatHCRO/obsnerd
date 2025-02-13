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
        self.obs_ind = None
        self.ods = None

    def set_par(self, **kwargs):
        for par, val in kwargs.items():
            if par not in ['obs_ind', 'ods']:
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

    def get_tracks(self, satname, start, duration, freqs=False, freq_unit='MHz', source='sopp'):
        if source == 'sopp':
            from . import sopp_engine
            self.tracks = sopp_engine.main(satname=satname, start=start, duration=duration, frequency=freqs, show_plots=False)
            logger.info(f"Found the following satellites: {self.tracks.keys()}")
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
        now = ttools.interpret_date('now') # - TimeDelta(8 * 3600, format='sec')
        plt.plot([now.datetime, now.datetime], [axlim[2], axlim[3]], 'g--')
        plt.axis(ymin=axlim[2], ymax=axlim[3])

    def choose_tracks(self, ctrks):
        self.chose = []
        this_sat = self.tracks[self.satname]
        for ch in [[int(x[:-1]), x[-1]] for x in ctrks.split(',')]:
            trk, pos = ch
            if pos == 'l':
                ind = this_sat[trk].imax // 2
            elif pos == 'm' or pos == 'p':
                ind = this_sat[trk].imax
            elif pos == 'r':
                ind = this_sat[trk].imax  + this_sat[trk].imax // 2
            this_sat[trk].set_par(obs_ind=ind, ods=[])
            self.chose.append({'track': trk, 'type': pos, 'ind': ind})
            plt.plot(this_sat[trk].utc[ind].datetime, this_sat[trk].el[ind].value, 'r.')
            plt.plot(this_sat[trk].utc[ind].datetime, this_sat[trk].az[ind].value, 'r.')

    def proc_tracks(self, obslen_min=8, freqs=[1990.0, 5500.0], freq_unit='MHz'):
        print("ONP142:  SHOULDN'T HAVE TO RE-SPECIFY FREQS")
        from odsutils import ods_engine
        obslen_min = obslen_min
        obslen_TD2 = TimeDelta(obslen_min * 60.0 / 2.0, format='sec')
        freqs = [(f * u.Unit(freq_unit)).to_value('Hz') for f in freqs]
        new_tracks = []
        this_sat = self.tracks[self.satname]
        for this in self.chose:
            newstart = this_sat[this['track']].utc[this['ind']] - obslen_TD2
            newstop = this_sat[this['track']].utc[this['ind']] - obslen_TD2
            for ff in freqs:
                odict = {'src_id':f"{this_sat[this['track']].srcname}",
                         'src_ra_j2000_deg': this_sat[this['track']].ra[this['ind']].to_value('deg'),
                         'src_dec_j2000_deg': this_sat[this['track']].dec[this['ind']].to_value('deg'),
                         'src_start_utc': f"{newstart.datetime.isoformat(timespec='seconds')}",
                         'src_end_utc': f"{newstop.datetime.isoformat(timespec='seconds')}",
                         'freq_lower_hz': ff-10.0E6,
                         'freq_upper_hz': ff+10.0E6}
                this_sat[this['track']].ods.append(odict)
                new_tracks.append(odict)
        with ods_engine.ODS() as ods:
            ods.pipe(new_tracks, intake='ods.json', defaults='defaults.json')

    def write_obsinfo(self, filter_file='filters.json'):
        from astropy.coordinates import angular_separation
        from copy import copy
        keys = ['ra', 'dec', 'az', 'el', 'tref', 'off_time', 'off_angle', 'ods']
        with open(filter_file, 'r') as fp:
            filter = json.load(fp)['Filters']
        obsinfo = {}
        start_mjd = None
        for this in self.chose:
            mjd = self.times[this['ind']].mjd
            if start_mjd is None:
                start_mjd = copy(mjd)
            satname = f"S{self.tracks[this['track']].satname}_{mjd:.4f}"
            obsinfo[satname] = {}
            obsinfo[satname]['ra'] = self.gcrs.ra[this['ind']].to_value('deg')
            obsinfo[satname]['dec'] = self.gcrs.dec[this['ind']].to_value('deg')
            obsinfo[satname]['tref'] = self.times[this['ind']].datetime.isoformat(timespec='seconds')
            obsinfo[satname]['az'] = self.azel.az[this['ind']].to_value('deg')
            obsinfo[satname]['el'] = self.azel.alt[this['ind']].to_value('deg')
            obsinfo[satname]['off_time'] = []
            obsinfo[satname]['off_angle'] = []
            newstart = self.times[this['ind']] - self.obslen_TD2
            newstop = self.times[this['ind']] + self.obslen_TD2
            istart = np.where(self.times > newstart)[0][0] - 1
            istop = np.where(self.times < newstop)[0][-1] + 1
            for i in range(istart, istop+1):
                dt = self.times[i] - self.times[this['ind']]
                da = angular_separation(self.gcrs.ra[i], self.gcrs.dec[i],
                                        self.gcrs.ra[this['ind']], self.gcrs.dec[this['ind']]) * np.sign(dt)
                obsinfo[satname]['off_time'].append(np.round(dt.to_value('sec'), 1))
                obsinfo[satname]['off_angle'].append(np.round(da.to_value('deg'), 2))
        fnout = f"obsinfo_{start_mjd:.0f}.json"
        with open(fnout, 'w') as fp:
            json.dump(obsinfo, fp, indent=4)
        self.obsinfo = obsinfo
