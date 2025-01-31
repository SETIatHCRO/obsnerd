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

    def get_ods(self, fn):
        if fn.endswith('.json'):
            self.ods.read_ods(fn)
        else:
            self.ods.get_defaults_dict('defaults.json')
            self.ods.add_from_file(fn)
        self.ods.ods['primary'].sort()
        self.groups = {}
        for entry in self.ods.ods['primary'].entries:
            key = (entry['src_start_utc'].datetime, entry['src_end_utc'].datetime)
            self.groups.setdefault(key, [])
            self.groups[key].append(entry)
        self.ods.new_ods_instance('output')

    def get_obs_from_ods(self):
        self.records = []
        for entries in self.groups.values():
            rec = ono_record.Record(observer=self.observer, project_id=self.project_id, ants=self.ants,
                                    attenuation=self.attenuation, focus=self.focus, backend=self.backend,
                                    time_per_int=self.time_per_int, coord='radec', lo=self.lo)
            freqs = []
            for i in range(len(entries)):
                freqs.append((entries[i]['freq_lower_hz'] + entries[i]['freq_upper_hz']) / 2.0)
                if not i:
                    source = entries[i]['src_id']
                    x = entries[i]['src_ra_j2000_deg']
                    y = entries[i]['src_dec_j2000_deg']
                    start = entries[i]['src_start_utc']
                    end = entries[i]['src_end_utc']
                    new_entry = copy(entries[i])
                    new_entry.update({'freq_lower_hz': 1990000000, 'freq_upper_hz': 1995000000})
                    self.ods.add_new_record('output', **new_entry)
                else:
                    if entries[i]['src_id'] != source: logger.error(f"Sources don't match.")
                    if entries[i]['src_ra_j2000_deg'] != x: logger.error("RAs don't match")
                    if entries[i]['src_dec_j2000_deg'] != y: logger.error("Decs don't match")
                    if entries[i]['src_start_utc'] != start: logger.error("Start times don't match")
                    if entries[i]['src_end_utc'] != end: logger.error("End times don't match")

            rec.update(freq=freqs, source=source, x=x, y=y, start=start, end=end)
            rec.proc()
            self.records.append(rec)
            self.ods.ods['output'].write('new_ods_output.json')

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