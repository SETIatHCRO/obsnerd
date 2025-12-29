import json
import logging
from copy import copy
from odsutils import ods_timetools as ttools
from odsutils import ods_tools as tools
from odsutils import logger_setup, ods_engine
from . import DATA_PATH, ono_engine, on_track
import astropy.units as u
import os.path as op

logger = logging.getLogger(__name__)
logger.setLevel('DEBUG')  # Set to lowest
from . import LOG_FILENAME, LOG_FORMATS, __version__

DEFAULTS = {'conlog': 'INFO', 'filelog': 'INFO', 'path': '.', 'log_filename': LOG_FILENAME,
            'observer': 'me', 'project_name': 'Project', 'project_id': 'p054', 'embargo': [],
            'attenuation': '8,8', 'focus': '', 'backend': 'xgpu', 'time_per_int_sec': 0.5}

SPACEX_LO = 1990 * u.MHz
SPACEX_HI = 1995 * u.MHz

ODS_URL = 'https://ods.hcro.org/ods.json'

class Observer:
    def __init__(self, obsfile='obsfile_rados.json', **kwargs):
        """
        Parameters
        ----------
        obsfile : str
            Obsinfo file to use
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
        self.obsfile = obsfile
        self.obsinfo = json.load(open(self.obsfile, 'r'))
        self.freqs = [x * u.Unit(self.obsinfo['Freq_unit']) for x in self.obsinfo['Tunings'].values()]
        self.lo = list(self.obsinfo['Tunings'].keys())
        self.ants = self.obsinfo['Ants']

        for key, val in kw.items():
            if key in track.fields:
                setattr(self, key, val)

    def get_obs(self, add_to_calendar=False, update_source_database=False):
        """
        Generate the obsnerd records from the obsinfo in obsfile.

        Fields in on_track.Track:
        'observer', 'project_name', 'project_id', 'ants', 'freq', 'lo', 'attenuation', 'focus', 'backend',
        'source', 'x', 'y', 'coord',
        'start', 'end', 'obs_time_sec', 'time_per_int_sec'

        Parmaeters
        ----------
        add_to_calendar : bool
            Whether to add the observation to the calendar
        lo_offset : float
            LO offset to use for frequency.
        lo_unit : str
            Unit of LO offset
        update_source_database : bool
            Whether to update the source database with new sources.
        'observer', 'project_name', 'project_id', 'ants', 'freq', 'lo', 'attenuation', 'focus', 'backend',
        'source', 'x', 'y', 'coord',
        'start', 'end', 'obs_time_sec', 'time_per_int_sec'
            
        """
        self.records = []
        ctr = 0
        for source_name, info in self.obsinfo['Sources'].items():
            rec = on_track.Track(observer=self.observer, project_name=self.project_name, project_id=self.project_id,
                                 ants=self.ants, freq=self.freqs, lo=self.lo,
                                 start=ttools.interpret_date(info['start'], 'Time'), end=ttools.interpret_date(info['stop'], 'Time'),
                                 source=source_name, x=info['ra']*u.deg, y=info['dec']*u.deg,
                                 attenuation=self.attenuation, focus=self.focus, backend=self.backend,
                                 time_per_int_sec=self.time_per_int_sec, coord='name')
            if update_source_database:
                ono_engine.update_source(src_id=rec.source, ra_hr=rec.x.to_value('hourangle'), dec_deg=rec.y.to_value('deg'))
            rec.proc()
            self.records.append(rec)
            ctr += 1
            print(f"Generated observer record: {ctr}")
            print("-------------------------------")
            print(rec)
            print("===============================================================================")
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
        """
        Prepare for an observation by making and posting the ODS online and adding the source to the source database.

        """ 
        from odsutils import ods_engine
        ods = ods_engine.ODS(engine_url=ODS_URL, default_ods_file=self.default_ods_default_file,
                             log_settings=self.log_settings)
        self.get_obs(add_to_calendar=add_to_calendar, update_source_database=update_source_database)
        ods.post_ods(ods_rados, instance_name="OUTPUT_ODS")
        if ods_assembly:
            ods.assemble_ods(ods_rados, post_to=ods_upload)

    def observe(self, is_actual=True, obsfile2use = 'obsinfo_rados.json'):
        if not is_actual:
            self.backend = 'test'
        self.get_obs(add_to_calendar=False)
        if not len(self.records):
            logger.error("Need to make observer records before you can observe.")
            return
        self.obs = ono_engine.CommandHandler(observer=self.observer, project_id=self.project_id, conlog=self.log_settings.conlog, filelog=self.log_settings.filelog)
        ant_list = self.obs.setants(self.ants)  # Assume that all antennas are the same so setants once...
        self.obs.setbackend(self.backend)  # ...and same backend...
        these_freq = [x.to_value('MHz') for x in self.records[0].freq]
        self.obs.setrf(freq=these_freq, lo=self.records[0].lo, attenuation=self.records[0].attenuation)  # ...and same rf setup
        obsrec = []
        for i, source in enumerate(self.records):
            if source.coord != 'name':
                logger.error(f"Currently only support 'name' coord type, not '{self.coord}'")
                continue
            obsrec.append(source.to_dict())
            if not i: print(source.__repr__(fprnt='header'))
            ts = ttools.interpret_date('now', fmt='%H:%M:%S')
            print(f"{ts} -- {i+1}/{len(self.records)}: {source.__repr__(fprnt='short')}")
            self.obs.move(source=source.source, coord_type=source.coord)
            tlength = ttools.wait(ttools.t_delta(source.start, -1.0*self.obs.obs_start_delay, 's'))
            if tlength is None:
                continue
            self.obs.take_data(source.obs_time_sec, source.time_per_int_sec)
        self.obs.release_ants()
        print("Observation complete -- summary to obsout.yaml")
        import yaml
        yaml.dump(obsrec, open('obsout.yaml', 'w'))