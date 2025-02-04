import logging
from aocalendar import aocalendar
from . import LOG_FILENAME, LOG_FORMATS, __version__
from odsutils import logger_setup
from odsutils import ods_timetools as ttools
from astropy.coordinates import SkyCoord
import astropy.units as u

logger = logging.getLogger(__name__)
logger.setLevel('DEBUG')  # Set to lowest


class Plan:
    def __init__(self, conlog='INFO', filelog=False, path='', loc='ata'):
        self.log_settings = logger_setup.Logger(logger, conlog=conlog, filelog=filelog, log_filename=LOG_FILENAME, path=path,
                                                conlog_format=LOG_FORMATS['conlog_format'], filelog_format=LOG_FORMATS['filelog_format'])
        logger.info(f"{__name__} ver. {__version__}")
        self.this_cal = aocalendar.Calendar(calfile='now', path=path, conlog=conlog, filelog=filelog, start_new=True, loc=loc)
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
        radec = obscoord.transform_to('icrs')
        self.this_cal.ods.get_defaults_dict('defaults.json')
        self.this_cal.ods.add_new_record(src_id='Asource', src_ra_j2000_deg=radec.ra.deg, src_dec_j2000_deg=radec.dec.deg,
                                         src_start_utc=self.start, src_end_utc=self.stop,
                                         freq_lower_hz=freq[0] * 1E6, freq_upper_hz=freq[0] * 1E6)
        self.this_cal.ods.add_new_record(src_id='Asource', src_ra_j2000_deg=radec.ra.deg, src_dec_j2000_deg=radec.dec.deg,
                                         src_start_utc=self.start, src_end_utc=self.stop,
                                         freq_lower_hz=freq[1] * 1E6, freq_upper_hz=freq[1] * 1E6)
        self.this_cal.ods.view_ods()
        self.this_cal.ods.write_ods('test_ods.json')