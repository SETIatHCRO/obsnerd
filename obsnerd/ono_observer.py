import json
import logging
from param_track import param_track_timetools as ttools
from param_track import Parameters
from odsutils import logger_setup
from . import DATA_PATH, ono_engine
from .on_observation import Observation
import astropy.units as u
import os.path as op

logger = logging.getLogger(__name__)
logger.setLevel('DEBUG')  # Set to lowest
from . import LOG_FILENAME, LOG_FORMATS, __version__

DEFAULTS = {'conlog': 'INFO', 'filelog': 'INFO', 'path': '.', 'log_filename': LOG_FILENAME, 'embargo': []}

SPACEX_LO = 1990 * u.MHz
SPACEX_HI = 1995 * u.MHz

ODS_URL = 'https://ods.hcro.org/ods.json'

class Observer(Parameters):
    def __init__(self, obsfile='obsinfo_rados.npz:data', **kwargs):
        """
        Parameters
        ----------
        obsfile : str
            Obsinfo file to use
        **kwargs: conlog, filelog, log_filename, path

        """
        super().__init__(ptnote='Observer parameters.', pttype=False, ptverbose=False)
        self.ptadd(**DEFAULTS)
        self.ptset(**kwargs)
        self.log_settings = logger_setup.Logger(logger, conlog=self.conlog, filelog=self.filelog,
                                                log_filename=self.log_filename, path=self.path,
                                                filelog_format=LOG_FORMATS['filelog_format'],
                                                conlog_format=LOG_FORMATS['conlog_format'])
        logger.debug(f"{__name__} ver. {__version__}")
        self.default_ods_default_file = op.join(DATA_PATH, 'ods_defaults_B.json')
        self.obsfile = obsfile
        

    def get_obs(self, add_to_calendar=False, update_source_database=False):
        """
        Generate the obsnerd records from the obsinfo in obsfile.

        Parameters
        ----------
        add_to_calendar : bool
            Whether to add the observation to the calendar
        update_source_database : bool
            Whether to update the source database with new sources.
            
        """
        self.obsinfo = Parameters(ptnote='Obsinfo file', ptinit=self.obsfile, pttype=False, ptverbose=False)
        if not len(self.obsinfo.observations):
            logger.error("Need to make observer records before you can get the overall.")
            return
        ctr = 0
        for source_name, info in self.obsinfo.observations.items():
            if update_source_database:
                ono_engine.update_source(src_id=info.source, ra_hr=info.ra.to_value('hourangle'), dec_deg=info.dec.to_value('deg'))
            ctr += 1
            print(f"Generated observer record {ctr}:  {info.source} -- {info.start} to {info.stop}")
            if source_name != info.source:
                print(f"{source_name} and {info.source} are not the same.")
        self.session = Parameters(ptnote="Overall session", ptstrict=False, ptverbose=False)
        self.session.ptset(sources=list(self.obsinfo.observations.keys()),
                           start=min([rec.start for rec in self.obsinfo.observations.values()]),
                           stop=max([rec.stop for rec in self.obsinfo.observations.values()]))
        self.session.ptset(observer=self.obsinfo.observer, project_name=self.obsinfo.project_name, project_id=self.obsinfo.project_id)
        if add_to_calendar:
            self.update_calendar()

    def update_calendar(self):
        # Get times to 5minutes
        t0 = ttools.interpret_date(ttools.interpret_date(self.session_start, '%Y-%m-%dT%H:%M'), 'datetime')
        t1 = ttools.interpret_date(ttools.interpret_date(self.session_stop, '%Y-%m-%dT%H:%M'), 'datetime')
        t0 = t0.replace(minute=(t0.minute // 5) * 5)
        t1 = ttools.t_delta(t1.replace(minute=(t1.minute // 5) * 5), 5, 'm')
        cal_day = ttools.interpret_date(self.session.start, '%Y-%m-%d')
        from aocalendar import google_calendar_sync
        self.google_calendar = google_calendar_sync.SyncCal()
        self.google_calendar.get_aocal(calfile=cal_day, path=self.log_settings.path, conlog=self.log_settings.conlog,
                                       filelog=self.log_settings.filelog, start_new=True)
        self.google_calendar.aocal.add(program=self.session.project_name, pid=self.session.project_id, observer=self.session.observer,
                                       utc_start=t0, utc_stop=t1)
        self.google_calendar.add_event_to_google_calendar(self.google_calendar.aocal.events[cal_day][-1])

    def rec2ods(self, record, ignore_freq=True):
        """
        Convert an observation to an ODS dictionary.

        Parameters
        ----------
        record : on_track.Track
            Record to convert

        Returns
        -------
        dict
            ODS dictionary

        """
        ods_rec = {'freq_lower_hz': SPACEX_LO.to_value('Hz'), 'freq_upper_hz': SPACEX_HI.to_value('Hz')}
        for a, b in record.ods_mapping.items():
            val = getattr(record, a, None)
            if a.startswith('freq'):
                if ignore_freq:
                    continue
                else:
                    print("FREQUENCY NOTE YET IMPLEMENTED IN rec2ods")
                    continue
            if val is None:
                continue
            if a in ['x', 'y']:
                ods_rec[b] = val.to_value('deg')
            else:
                ods_rec[b] = val
        return ods_rec

    def observe_prep(self, add_to_calendar=False, defaults=None,
                     ods_rados = "/opt/mnt/share/ods_project/ods_rados.json",
                     ods_upload = "/opt/mnt/share/ods_upload/ods.json",
                     update_source_database=True,
                     ods_assembly=True):
        """
        Prepare for an observation by making and posting the ODS online and adding the source to the source database.

        """ 
        from odsutils import ods_engine
        self.get_obs(add_to_calendar=add_to_calendar, update_source_database=update_source_database)
        defaults = defaults if defaults is not None else self.default_ods_default_file
        ods = ods_engine.ODS(defaults = defaults)
        ods_update = []
        for rec in self.records:
            ods_rec = self.rec2ods(rec)
            ods_update.append(ods_rec)
        ods.add(ods_update)
        ods.post_ods(ods_rados)
        if ods_assembly:
            ods.assemble_ods(ods_rados, post_to=ods_upload)

    def observe(self, is_actual=True):
        if not is_actual:
            self.backend = 'test'
        self.get_obs(add_to_calendar=False)
        if not len(self.records):
            logger.error("Need to make observer records before you can observe.")
            return
        self.obs = ono_engine.CommandHandler(observer=self.observer, project_id=self.project_id, conlog=self.log_settings.conlog, filelog=self.log_settings.filelog)
        self.obsinfo['Ants'] = self.obs.setants(self.ants)  # Assume that all antennas are the same so setants once...
        self.obs.setbackend(self.backend)  # ...and same backend...
        these_freq = [x.to_value('MHz') for x in self.records[0].freq]
        self.obs.setrf(freq=these_freq, lo=self.records[0].lo, attenuation=self.records[0].attenuation)  # ...and same rf setup
        obsrec = []
        for i, source in enumerate(self.records):
            if source.coord != 'name':
                logger.error(f"Currently only support 'name' coord type, not '{source.coord}'")
                continue
            obsrec.append(source.pt_to_dict())
            if not i: print(source.ptshow(vals_only=True))
            ts = ttools.interpret_date('now', fmt='%H:%M:%S')
            print(f"{ts} -- {i+1}/{len(self.records)}: {source.ptshow(vals_only=True)}")
            print(source.source, source.coord, source.x, source.y)
            self.obs.move(source=source.source, coord_type=source.coord)
            tlength = ttools.wait(ttools.t_delta(source.start, -1.0*self.obs.obs_start_delay, 's'))
            if tlength is None:
                continue
            self.obs.take_data(source.start, source.obs_time_sec, source.time_per_int_sec)
        self.obs.release_ants()
        print("Observation complete -- exit calendar.")
        json.dump(self.obsinfo, open('obsout.json', 'w'), indent=2)