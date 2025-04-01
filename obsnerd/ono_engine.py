try:
    from ATATools import ata_control  # type: ignore
    from ATATools import ata_sources  # type: ignore
    from SNAPobs.snap_hpguppi import record_in as hpguppi_record_in  # type: ignore
    from SNAPobs.snap_hpguppi import snap_hpguppi_defaults as hpguppi_defaults  # type: ignore
    from SNAPobs.snap_hpguppi import auxillary as hpguppi_auxillary  # type: ignore
    from SNAPobs import snap_config, snap_if  # type: ignore
except ImportError:
    from .ono_debug import Empty
    ata_control = Empty('ata_control')
    ata_sources = Empty('ata_sources')
    hpguppi_record_in = Empty('hpguppi_record_in')
    hpguppi_defaults = Empty('hpguppi_defaults')
    hpguppi_auxillary = Empty('hpguppi_auxillary')
    snap_config = Empty('snap_config')
    snap_if = Empty('snap_if')

import atexit, os
import astropy.units as u
from . import __version__, LOG_FILENAME, ono_record
from . import on_sys
from odsutils import logger_setup
from odsutils import ods_tools as tools
from time import sleep
from copy import copy
import logging
logger = logging.getLogger(__name__)
logger.setLevel('DEBUG')  # Set to lowest


OBS_START_DELAY = 10  # time to prep data until collecting
OBS_DAWDLE = 5  # extra time to "sleep" to make sure things are done
LO_LIST = ['A', 'B']
DEFAULTS = {'observer': None, 'project_id': None, 'project_name': None,
            'conlog': 'WARNING', 'filelog': 'INFO', 'path': '.', 'log_filename': LOG_FILENAME,
            'obs_start_delay': OBS_START_DELAY, 'obs_dawdle': OBS_DAWDLE}


def get_LO_hpguppi(LOs=['A', 'B']):
    # d = hpguppi_defaults.hashpipe_targets_LoA.copy()
    # d.update(hpguppi_defaults.hashpipe_targets_LoB)
    return {'seti-node%i'%i: [0,1] for i in range(1,8)}


def update_source(src_id, ra_hr, dec_deg, owner='ddeboer', category='starlink'):
    """
    Update a source entry -- if exists will delete it and re-add.

    Parameters
    ----------
    src_id : str
        Source ID to update
    ra_hr : float
        Right ascension in hours
    dec_deg : float
        Declination in degrees
    owner : str
        Owner of the source
    category : str
        Category of the source

    """
    sources = [x['Source'] for x in  ata_sources.list_catalog(owner=owner, category=category)]
    if src_id in sources:
        ata_sources.delete_catalog_entry(owner=owner, category=category, source=src_id)
    ata_sources.add_catalog_entry(owner=owner, category=category, source=src_id, ra=ra_hr, dec=dec_deg)


class CommandHandler:
    def __init__(self, **kwargs):
        kw = copy(DEFAULTS)
        kw.update(kwargs)
        self.logset = logger_setup.Logger(logger, conlog=kw['conlog'], filelog=kw['filelog'],
                                          log_filename=kw['log_filename'], path=kw['path'])
        logger.debug(f"{__name__} ver. {__version__}")

        for this_attr in ['observer', 'project_id', 'project_name', 'obs_start_delay', 'obs_dawdle']:
            setattr(self, this_attr, kw[this_attr])
        logger.info(f"session start: {self.observer} -- {self.project_id}")
        self.rec = ono_record.Record(observer=self.observer, project_id=self.project_id)

    def setants(self, ant_list='rfsoc_active', remove_ants=[], park_when_done=True):
        if ant_list.startswith('rfsoc_active'):
            self.ant_list = snap_config.get_rfsoc_active_antlist()
            xxx = ant_list.split('-')
            if len(xxx) == 2:
                remove_ants += tools.listify(xxx[1])
        elif ant_list in on_sys.ANT_LISTS:
            self.ant_list = copy(on_sys.ANT_LIST[ant_list])
        elif isinstance(ant_list, str):
            self.ant_list = tools.listify(ant_list)
        if not len(self.ant_list):
            raise ValueError("No antennas specified.")
        for badun in remove_ants:
            if badun in self.ant_list:
                logger.info(f"Removing antenna {badun}")
                self.ant_list.remove(badun)
        ata_control.move_ant_group(self.ant_list, 'none', 'bfa')
        logger.info(f"antennas:  {(', ').join(self.ant_list)}")
        atexit.register(ata_control.release_antennas, self.ant_list, park_when_done)
        self.rec.update(ants=self.ant_list)

    def release_ants(self, park_when_done=True):
        ata_control.release_antennas(self.ant_list, park_when_done)

    def setrf(self, freq, lo=['A', 'B'], attenuation=[8, 8], focus=False, freq_unit='MHz'):
        """
        Parameters (kwargs)
        -------------------
        freq : *
            Frequency of observation [MHz], must match length of lo
        lo : list
            LO to use A/B/C/D
        attenuation : list of int
            Attenuation setting
        focus : bool
            If True, focus the on max frequ
        freq_unit : astropy.units
            Unit of frequency, default MHz

        """
        if isinstance(freq[0], u.Quantity):
            freq_unit = freq[0].unit
            logger.info(f"freq_unit: {freq_unit}")
            if not isinstance(freq, u.Quantity):
                freq = freq * u.Unit(freq_unit)
        else:
            freq = tools.listify(freq) * u.Unit(freq_unit)
        self.freq = freq
        freq_MHz = [x.to_value('MHz') for x in self.freq][0]
        self.lo = tools.listify(lo)
        self.attenuation = tools.listify(attenuation)
        # antlo_list = [ant+lo.upper() for lo in lo for ant in self.ant_list]
        # snap_if.tune_if_antslo(antlo_list)
        ffoc = 1E9 if not focus else max(freq_MHz)
        logger.info(f"freq: {', '.join([str(x) for x in freq_MHz])} MHz")
        logger.info(f"lo: {', '.join(self.lo)}")
        logger.info(f"attenuation: {', '.join([str(x) for x in self.attenuation])}")
        for fMHz, llo in zip(freq_MHz, self.lo):
            this_freq = [fMHz] * len(self.ant_list)
            ata_control.set_freq(this_freq, self.ant_list, lo=llo.lower(), nofocus=fMHz<ffoc)
        if focus:
            import time
            time.sleep(20)
        ata_control.autotune(self.ant_list)
        logger.warning("REMOVING ATTENTUATION")
        # ata_control.set_atten_thread([[f'{ant}x', f'{ant}y'] for ant in self.ant_list],
        #                              [[self.attenuation[0], self.attenuation[1]] for ant in self.ant_list])
        self.rec.update(freq=self.freq, lo=lo, attenuation=attenuation, focus=ffoc)

    def setbackend(self, backend='xpgu'):
        """
        Parameters (kwargs)
        -------------------
        backend : str
            Name of backend to use

        """
        self.backend = backend
        if self.project_id is None:
            self.project_id = input("Need a project_id:  ")
        logger.info(f"Backend {self.backend} for project {self.project_id}")
        if self.backend == 'xgpu':
            try:
                os.system("ansible-playbook /home/sonata/src/ansible_playbooks/hashpipe/xgpu_record.yml")
                # subprocess.run("ansible-playbook /home/sonata/src/ansible_playbooks/hashpipe/xgpu_record.yml")
            except FileNotFoundError:
                logger.error("Ansible playbook not found, cannot start backend")
        else:
            logger.error(f"Invalid backend: {self.backend} -- no action")
            return
        try:
            # subprocess.run(f"/home/sonata/src/observing_campaign/backend_setup_scripts/set_keys_uvh5_mv_{self.project_id}.py")
            os.system(f"/home/sonata/src/observing_campaign/backend_setup_scripts/set_keys_uvh5_mv_{self.project_id}.py")
        except FileNotFoundError:
            logger.error("Script not found, cannot start backend")
            return
        self.rec.update(backend=backend)

    def move(self, source, coord_type='name', use_ants=None, x_unit='deg', y_unit='deg'):
        """
        Parameters (kwargs)
        -------------------
        source : str
            target source (depends on coord_type) - source_name; x,y; filename
        coord_type : str
            coordination type (azel, radec, name, traj)
        use_ants : list
            List of antennas to move, default to self.ant_list
        x_unit : astropy.units
            Unit of x
        y_unit : astropy.units 
            Unit of y

        """
        if not isinstance(source, str) or len(source) == 0:
            logger.error("Need a valid source")
            return
        if coord_type not in ['azel', 'radec', 'name', 'traj']:
            logger.error(f"Invalid coord_type: {coord_type}")
            return
        use_ants = self.ant_list if use_ants is None else tools.listify(use_ants)
        if not len(use_ants):
            logger.error("Need antennas to move")
            return

        if ',' in source:
            x = self.source.split(',')[0] * u.Unit(x_unit)
            y = self.source.split(',')[1] * u.Unit(y_unit)
            source = coord_type
            logger.info(f'move to: {coord_type}:  {x}, {y}')
        else:
            x, y = None, None
            logger.info(f'move to: {source}')
        self.rec.update(source_name=source, x=x, y=y, coord=coord_type)
        if coord_type == 'azel':
            logger.info(f"azel: {x},{y}")
            sr = ata_control.set_az_el(use_ants, x.to_value('deg'), y.to_value('deg'))
        elif coord_type == 'radec':
            logger.info(f"radec: {x},{y}")
            sr = ata_control.track_source(use_ants, radec=[x.to_value('hourangle'), y.to_value('deg')])
        elif coord_type == 'name':
            # source = ata_control.track_source(use_ants, source=self.source)
            logger.info(f"source: {source}")
            sr = ata_control.make_and_track_ephems(source, use_ants)
        elif coord_type == 'traj':
            logger.info(f"traj: {self.source}")
            from obsnerd.trajectory_engine import TRACK_YAML_FILENAME
            ephem = ata_control.upload_ephemeris(self.source)
            sr = ata_control.track_ephemeris(ephem, use_ants, wait=True)
            try:
                with open(TRACK_YAML_FILENAME, 'r') as fp:
                    for line in fp:
                        if len(line) > 2:
                            logger.info(f"track: {line.strip()}")
            except FileNotFoundError:
                pass

    def take_data(self, obs_time_sec, time_per_int_sec):
        self.rec.update(obs_time_sec=obs_time_sec, time_per_int_sec=time_per_int_sec)
        d = get_LO_hpguppi(self.lo)
        keyval_dict = {'XTIMEINT': time_per_int_sec}
        hpguppi_auxillary.publish_keyval_dict_to_redis(keyval_dict, d, postproc=False)
        hpguppi_record_in.record_in(self.obs_start_delay, obs_time_sec, hashpipe_targets = d)
        sleepy = float(self.obs_start_delay + obs_time_sec + self.obs_dawdle)
        if sleepy > 0.0:
            sleep(sleepy)

    def note(self, note):
        """
        Parameter (kwargs)
        ------------------
        note : str
            Note to log

        """
        logger.info(f"note: {note}")

