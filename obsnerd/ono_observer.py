from . import on_sys as OS
import logging
from copy import copy
from odsutils import ods_timetools as ttools
from odsutils import ods_tools as tools
from odsutils import logger_setup, ods_engine
from obsnerd import ono_engine, ono_record

logger = logging.getLogger(__name__)
logger.setLevel('DEBUG')  # Set to lowest
from . import LOG_FILENAME, __version__

DEFAULTS = {'observer': None, 'project_id': None,
            'conlog': 'WARNING', 'filelog': 'INFO', 'path': '.', 'log_filename': LOG_FILENAME,
            'observer': 'me', 'project_id': 'pid', 'ants': 'rfsoc_active-1k', 'lo': ['A', 'B'],
            'attenuation': '8,8', 'focus': '', 'backend': 'xgpu', 'time_per_int': 0.5}

UNITS = {'Hz': 1.0, 'kHz': 1E3, 'MHz': 1E6, 'GHz': 1E9}
SPACEX_LO = 1990 * UNITS['MHz']
SPACEX_HI = 1995 * UNITS['MHz']


class Observer:
    def __init__(self, **kwargs):
        """
        Parameters
        ----------


        """
        kw = copy(DEFAULTS)
        kw.update(kwargs)
        self.logset = logger_setup.Logger(logger, conlog=kw['conlog'], filelog=kw['filelog'],
                                          log_filename=kw['log_filename'], path=kw['path'])
        logger.debug(f"{__name__} ver. {__version__}")
        self.records = []  # ono_records to be made out of ODS
        rec = ono_record.Record()
        for key, val in kw.items():
            if key in rec.fields:
                setattr(self, key, val)
        self.ods = ods_engine.ODS()
        #'observer', 'project_id', 'ants', 'freq', 'lo', 'attenuation', 'focus', 'backend',
        #'source', 'x', 'y', 'coord',
        #'time_per_int', 'start', 'end'

    def get_ods(self, fn, defaults='defaults.json'):
        if fn.endswith('.json') or fn.startswith('http'):
            self.ods.read_ods(fn)
        else:
            self.ods.get_defaults_dict(defaults)
            self.ods.add_from_file(fn)
        self.ods.ods['primary'].sort()
        self.groups = {}
        for entry in self.ods.ods['primary'].entries:
            key = (entry['src_start_utc'].datetime, entry['src_end_utc'].datetime)
            self.groups.setdefault(key, [])
            self.groups[key].append(entry)

    def get_obs_from_ods(self, lo_offset=10.0, lo_unit='MHz'):
        self.records = []
        self.ods.new_ods_instance('output')
        for entries in self.groups.values():
            rec = ono_record.Record(observer=self.observer, project_id=self.project_id, ants=self.ants,
                                    attenuation=self.attenuation, focus=self.focus, backend=self.backend,
                                    time_per_int=self.time_per_int, coord='radec', lo=self.lo)
            freqs = []
            pars = {'src_id': None, 'src_ra_j2000_deg': None, 'src_dec_j2000_deg': None, 'src_start_utc': None, 'src_end_utc': None}
            for i in range(len(entries)):
                freqs.append((entries[i]['freq_lower_hz'] + entries[i]['freq_upper_hz']) / 2.0 - lo_offset*UNITS[lo_unit])
                if not i:  # Get common parameters in first pass through and make consolidated new ods record
                    for par in list(pars.keys()):
                        pars[par] = entries[i][par]
                    new_entry = copy(entries[i])
                    new_entry.update({'freq_lower_hz': SPACEX_LO, 'freq_upper_hz': SPACEX_HI})
                    self.ods.add_new_record('output', **new_entry)
                else:  # Just check if different
                    for par in list(pars.keys()):
                        pars[par] = entries[i][par]
                        if entries[i][par] != pars[par]: logger.error(f"Field mismatch - {par}")
            rec.update(freq=freqs, source=pars['src_id'], x=pars['src_ra_j2000_deg'], y=pars['src_dec_j2000_deg'],
                       start=pars['src_start_utc'], end=pars['src_end_utc'])
            rec.proc()
            self.records.append(rec)

    def update_ods(self, ods_input, ods_output):
        """
        Update the working ods with the new additions.

        Parameters
        ----------
        ods_input : str
            Location of current working ods
        ods_output : str
            Location of the one to write.  If None use ods_input

        """
        self.ods.pipe('output', intake=ods_input, output=ods_output)

    def observe(self, is_actual=True):
        if not is_actual:
            self.backend = 'test'
        self.obs = ono_engine.CommandHandler(observer=self.observer, project_id=self.project_id)
        self.obs.setants(self.ants)  # Assume that all antennas are the same so setants once
        self.obs.setbackend(self.backend)  # And same backend
        for source in self.records:
            self.obs.setrf(freq=source.freq, lo=source.lo, attenuation=source.attenuation)
            self.obs.move(f"{source.x},{source.y}", source.coord)
            tlength = ttools.wait_until(ttools.t_delta(source.start, -1.0*self.obs.obs_start_delay, 's'))
            if tlength is None:
                continue
            self.obs.observe(source.obs_time, source.time_per_int)