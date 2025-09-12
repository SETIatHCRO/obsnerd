import logging
from copy import copy
from odsutils import ods_timetools as ttools
from odsutils import ods_tools as tools
from odsutils import logger_setup, ods_engine
from . import DATA_PATH, ono_engine, on_track, on_sys
import astropy.units as u
import os.path as op

logger = logging.getLogger(__name__)
logger.setLevel('DEBUG')  # Set to lowest
from . import LOG_FILENAME, LOG_FORMATS, __version__

DEFAULTS = {'conlog': 'INFO', 'filelog': 'INFO', 'path': '.', 'log_filename': LOG_FILENAME,
            'observer': 'me', 'project_name': 'Project', 'project_id': 'pid', 'ants': ':ant_file.dat', 'embargo': [],
            'lo': ['A', 'B'], 'attenuation': '8,8', 'focus': '', 'backend': 'xgpu', 'time_per_int_sec': 0.5}

SPACEX_LO = 1990 * u.MHz
SPACEX_HI = 1995 * u.MHz

ODS_URL = 'https://www.seti.org/sites/default/files/HCRO/ods.json'

class Observer:
    def __init__(self, **kwargs):
        """
        Parameters
        ----------
        **kwargs: conlog, filelog, log_filename, path

        """
        kw = copy(DEFAULTS)
        kw.update(kwargs)
        self.log_settings = logger_setup.Logger(logger, conlog=kw['conlog'], filelog=kw['filelog'],
                                                log_filename=kw['log_filename'], path=kw['path'],
                                                filelog_format=LOG_FORMATS['filelog_format'],
                                                conlog_format=LOG_FORMATS['conlog_format'])
        logger.debug(f"{__name__} ver. {__version__}")
        self.records = []  # on_track to be made out of ODS
        track= on_track.Track()
        self.embargo = tools.listify(kw['embargo'])
        self.default_ods_default_file = op.join(DATA_PATH, 'ods_defaults_B.json')
        # Check for antenna file
        if kw['ants'][0] == ':':
            with open(kw['ants'][1:], 'r') as fp:
                kw['ants'] = fp.read().strip()
        for key, val in kw.items():
            if key in track.fields:
                setattr(self, key, val)

    def get_ods(self, ods_input, defaults=None):
        """
        Read an ODS input and make dictionary based on start/end.

        Parmaeters
        ----------
        ods_input : str
            ODS input, file or URL
        defaults : str or None
            Default values to use

        """
        defaults = defaults if defaults is not None else self.default_ods_default_file
        self.ods = ods_engine.ODS(conlog='WARNING', defaults=defaults)
        self.ods.add(ods_input)
        self.ods.ods['primary'].sort()
        self.groups = {}
        group_keys = {}
        for entry in self.ods.ods['primary'].entries:
            key = (entry['src_start_utc'].datetime, entry['src_end_utc'].datetime)
            self.groups.setdefault(key, [])
            self.groups[key].append(entry)
            group_keys[entry['src_start_utc'].datetime] = key
        self.sorted_group_keys = []  # Sorted by start time -- assumes a fair bit of niceness in the groups list
        for key in sorted(group_keys.keys()):
            self.sorted_group_keys.append(group_keys[key])

    def get_obs_from_ods(self, add_to_calendar=False, lo_offset=10.0, lo_unit='MHz', update_source_database=False):
        """
        Generate the obsnerd records based on an ods (as read in get_ods).

        Parmaeters
        ----------
        lo_offset : float
            LO offset to use for frequency.
        lo_unit : str
            Unit of LO offset

        """
        self.records = []
        self.ods.new_ods_instance('OUTPUT_ALL')
        self.ods.new_ods_instance('OUTPUT_ODS')
        if self.groups is None:
            logger.error("Unable to read ODS file")
            return
        for group_key in self.sorted_group_keys:
            entries = self.groups[group_key]
            rec = on_track.Track(observer=self.observer, project_name=self.project_name, project_id=self.project_id,
                                    ants=self.ants, attenuation=self.attenuation, focus=self.focus, backend=self.backend,
                                    time_per_int_sec=self.time_per_int_sec, coord='name', lo=self.lo)
            freqs = []
            pars = {'src_id': None, 'src_ra_j2000_deg': None, 'src_dec_j2000_deg': None, 'src_start_utc': None, 'src_end_utc': None}
            for i in range(len(entries)):
                freq_ods = ((entries[i]['freq_lower_hz'] + entries[i]['freq_upper_hz']) / 2.0) * u.Hz
                freqs.append(freq_ods - lo_offset*u.Unit(lo_unit))
                if not i:  # Get common parameters in first pass through and make consolidated new ods record
                    for par in list(pars.keys()):
                        pars[par] = entries[i][par]
                    new_entry = copy(entries[i])
                    new_entry.update({'freq_lower_hz': SPACEX_LO.to_value('Hz'), 'freq_upper_hz': SPACEX_HI.to_value('Hz')})
                    notes = on_sys.parse_ods_notes(entries[i])
                    self.ods.add(**new_entry, instance_name="OUTPUT_ALL")
                    if notes['ods'] == 'True':
                        self.ods.add(**new_entry, instance_name="OUTPUT_ODS")
                else:  # Just check if different
                    for par in list(pars.keys()):
                        pars[par] = entries[i][par]
                        if entries[i][par] != pars[par]: logger.error(f"Field mismatch - {par}")
            rec.update(freq=freqs * u.Hz, source=pars['src_id'], x=pars['src_ra_j2000_deg'] * u.deg, y=pars['src_dec_j2000_deg'] * u.deg,
                       start=pars['src_start_utc'], end=pars['src_end_utc'])
            if update_source_database:
                ono_engine.update_source(src_id=rec.source, ra_hr=rec.x.to_value('hourangle'), dec_deg=rec.y.to_value('deg'))
            rec.proc()
            self.records.append(rec)
        self.get_overall()
        if add_to_calendar:
            self.update_calendar()

    def get_overall(self):
        if not len(self.records):
            logger.error("Need to make observer records before you can get the overall.")
            return
        kw = {}
        try:
            t0 = min([rec.start for rec in self.records])
            t1 = max([rec.end for rec in self.records])
        except AttributeError:
            logger.error("Need to make observer records before you can get the overall.")
            return
        self.overall = on_track.Track()
        for fld in self.overall.fields:
            try:  # Tragically assume that the first record has most of the same stuff as the rest...
                kw[fld] = getattr(self.records[0], fld)
            except (AttributeError, IndexError):
                continue
        kw['start'], kw['end'] = t0, t1
        self.overall.update(**kw)

    def update_calendar(self):
        # Get times to 5minutes
        t0 = ttools.interpret_date(ttools.interpret_date(self.overall.start, '%Y-%m-%dT%H:%M'), 'datetime')
        t1 = ttools.interpret_date(ttools.interpret_date(self.overall.end, '%Y-%m-%dT%H:%M'), 'datetime')
        t0 = t0.replace(minute=(t0.minute // 5) * 5)
        t1 = ttools.t_delta(t1.replace(minute=(t1.minute // 5) * 5), 5, 'm')
        cal_day = ttools.interpret_date(self.overall.start, '%Y-%m-%d')
        from aocalendar import google_calendar_sync
        self.google_calendar = google_calendar_sync.SyncCal()
        self.google_calendar.get_aocal(calfile=cal_day, path=self.log_settings.path, conlog=self.log_settings.conlog,
                                       filelog=self.log_settings.filelog, start_new=True)
        self.google_calendar.aocal.add(program=self.overall.project_name, pid=self.overall.project_id, observer=self.overall.observer,
                                       utc_start=t0, utc_stop=t1)
        self.google_calendar.add_event_to_google_calendar(self.google_calendar.aocal.events[cal_day][-1])

    def observe_prep(self, add_to_calendar=False,
                     ods_rados = "/opt/mnt/share/ods_project/ods_rados.json",
                     ods_upload = "/opt/mnt/share/ods_upload/ods.json",
                     update_source_database=True,
                     ods_assembly=True):
        self.get_ods(op.basename(ods_rados))
        self.get_obs_from_ods(add_to_calendar=add_to_calendar, update_source_database=update_source_database)
        self.ods.post_ods(ods_rados, instance_name="OUTPUT_ODS")
        if ods_assembly:
            self.ods.assemble_ods(op.dirname(ods_rados), post_to=ods_upload)

    def observe(self, is_actual=True, ods2use = 'ods_rados.json'):
        if not is_actual:
            self.backend = 'test'
        self.get_ods(ods2use)
        self.get_obs_from_ods(add_to_calendar=False)
        if not len(self.records):
            logger.error("Need to make observer records before you can observe.")
            return
        self.obs = ono_engine.CommandHandler(observer=self.observer, project_id=self.project_id, conlog=self.log_settings.conlog, filelog=self.log_settings.filelog)
        self.obs.setants(self.ants)  # Assume that all antennas are the same so setants once
        self.obs.setbackend(self.backend)  # And same backend
        for i, source in enumerate(self.records):
            if source.coord != 'name':
                logger.error(f"Currently only support 'name' coord type, not '{self.coord}'")
                continue
            if not i: print(source.__repr__(use='header'))
            ts = ttools.interpret_date('now', fmt='%H:%M:%S')
            print(f"{ts} -- {i+1}/{len(self.records)}: {source.__repr__(use='short')}")
            print(f"freqs: {', '.join([str(x) for x in source.freq])}")
            these_freq = [x.to_value('MHz') for x in source.freq]
            self.obs.setrf(freq=these_freq, lo=source.lo, attenuation=source.attenuation)
            self.obs.move(source=source.source, coord_type=source.coord)
            tlength = ttools.wait(ttools.t_delta(source.start, -1.0*self.obs.obs_start_delay, 's'))
            if tlength is None:
                continue
            self.obs.take_data(source.obs_time_sec, source.time_per_int_sec)
        self.obs.release_ants()