import json
import logging
from param_track import param_track_timetools as ttools
from param_track import Parameters
from odsutils import logger_setup
from . import ono_engine
import astropy.units as u

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
        super().__init__(ptnote='Observer parameters.', ptstrict=False, pttype=False, ptverbose=False)
        self.ptadd(**DEFAULTS)
        self.ptset(**kwargs)
        self.log_settings = logger_setup.Logger(logger, conlog=self.conlog, filelog=self.filelog,
                                                log_filename=self.log_filename, path=self.path,
                                                filelog_format=LOG_FORMATS['filelog_format'],
                                                conlog_format=LOG_FORMATS['conlog_format'])
        logger.debug(f"{__name__} ver. {__version__}")
        self.obsinfo = Parameters(ptnote='Obsinfo', pttype=False, ptverbose=False, filename=obsfile)

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
        self.obsinfo.ptinit(self.obsinfo.filename)  # actually get the records.
        if not len(self.obsinfo.observations):
            logger.error("Need to make observer records before you can get the overall.  'See onp_plan.py'")
            return
        ctr = 0
        for track in self.obsinfo.observations:
            if update_source_database:
                ono_engine.update_source(src_id=track.source, ra_hr=track.ra.to_value('hourangle'), dec_deg=track.dec.to_value('deg'))
            ctr += 1
            print(f"Generated observer record {ctr}:  {track.source} -- {track.start} to {track.stop}")
        self.session = Parameters(ptnote="Overall session", ptstrict=False, ptverbose=False)
        self.session.ptset(sources=[track.source for track in self.obsinfo.observations],
                           start=min([rec.start for rec in self.obsinfo.observations]),
                           stop=max([rec.stop for rec in self.obsinfo.observations]))
        self.session.ptset(observer=self.obsinfo.observer, project_name=self.obsinfo.project_name, project_id=self.obsinfo.project_id)
        if add_to_calendar:
            self.update_calendar()

    def update_calendar(self):
        # Get times to 5minutes
        t0 = ttools.interpret_date(ttools.interpret_date(self.session.start, '%Y-%m-%dT%H:%M'), 'datetime')
        t1 = ttools.interpret_date(ttools.interpret_date(self.session.stop, '%Y-%m-%dT%H:%M'), 'datetime')
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
        ods_rec = {}
        for a, b in self.obsinfo.ods_mapping.items():
            val = getattr(record, a)
            if b == 'src_id':
                ods_rec[b] = val
            elif b.startswith('freq'):
                ods_rec[b] = val.to_value('Hz').item()
            elif b == 'src_ra_j2000_deg' or b == 'src_dec_j2000_deg':
                ods_rec[b] = val.to_value('deg').item()
            elif b == 'src_start_utc' or b == 'src_end_utc':
                ods_rec[b] = val.isot
            elif b == 'corr_integ_time_sec':
                ods_rec[b] = val.to_value('sec')
        return ods_rec

    def observe_prep(self, add_to_calendar=False, defaults=None,
                     ods_rados = "/opt/mnt/share/ods_project/ods_rados.json",
                     ods_upload = "/opt/mnt/share/ods_upload/ods.json",
                     update_source_database=True,
                     ods_assembly=True):
        """
        Prepare for an observation by making and posting the ODS online and adding the source to the source database.

        """ 
        from odsutils import ods_engine, DATA_PATH
        import os.path as op
        self.get_obs(add_to_calendar=add_to_calendar, update_source_database=update_source_database)
        defaults = defaults if defaults is not None else op.join(DATA_PATH, 'ods_defaults_ata_B.json')
        ods = ods_engine.ODS(defaults = defaults)
        ods_update = []
        for rec in self.obsinfo.observations:
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
        if not len(self.obsinfo.observations):
            logger.error("Need to make observer records before you can observe.")
            return
        self.obs = ono_engine.CommandHandler(observer=self.obsinfo.observer, project_id=self.obsinfo.project_id,
                                             conlog=self.log_settings.conlog, filelog=self.log_settings.filelog)
        self.obsinfo.ants = self.obs.setants(self.obsinfo.ants)  # Assume that all antennas are the same so setants once...
        self.obs.setbackend(self.obsinfo.backend)  # ...and same backend...etc...
        these_lo = self.obsinfo.observations[0].lo
        these_freq = [getattr(self.obsinfo.observations[0], f'LO{x.upper()}').to_value('MHz') for x in these_lo]
        these_attenuation = self.obsinfo.observations[0].attenutation
        self.obs.setrf(freq=these_freq, lo=these_lo, attenuation=these_attenuation)  # ...and same rf setup
        num_obs = len(self.obsinfo.observations)
        obsrec = []
        for i, source in enumerate(self.obsinfo.observations):
            if source.coord != 'name':
                logger.error(f"Currently only support 'name' coord type, not '{source.coord}'")
                continue
            obsrec.append(source.pt_to_dict())
            if not i: print(source.ptshow(vals_only=True))
            ts = ttools.interpret_date('now', fmt='%H:%M:%S')
            print(f"{ts} -- {i+1}/{num_obs}: {source.ptshow(vals_only=True)}")
            print(source.source, source.coord, source.ra, source.dec)
            self.obs.move(source=source.source, coord_type=source.coord)
            tlength = ttools.wait(ttools.t_delta(source.start, -1.0*self.obs.obs_start_delay, 's'))
            if tlength is None:
                continue
            self.obs.take_data(source.start, source.obs_time_sec, source.time_per_int_sec)
        self.obs.release_ants()
        print("Observation complete -- exit calendar.")
        json.dump(self.obsinfo, open('obsout.json', 'w'), indent=2)