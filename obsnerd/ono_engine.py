try:
    from ATATools import ata_control  # type: ignore
    from SNAPobs.snap_hpguppi import record_in as hpguppi_record_in  # type: ignore
    from SNAPobs.snap_hpguppi import snap_hpguppi_defaults as hpguppi_defaults  # type: ignore
    from SNAPobs.snap_hpguppi import auxillary as hpguppi_auxillary  # type: ignore
    from SNAPobs import snap_config, snap_if  # type: ignore
except ImportError:
    from .ono_debug import Empty
    ata_control = Empty('ata_control')
    hpguppi_record_in = Empty('hpguppi_record_in')
    hpguppi_defaults = Empty('hpguppi_defaults')
    hpguppi_auxillary = Empty('hpguppi_auxillary')
    snap_config = Empty('snap_config')
    snap_if = Empty('snap_if')

import atexit
from . import __version__, LOG_FILENAME, ono_record
from . import on_sys as OS
from odsutils import logger_setup
from odsutils import ods_tools as tools
import subprocess
from time import sleep
from copy import copy
import logging
logger = logging.getLogger(__name__)
logger.setLevel('DEBUG')  # Set to lowest


OBS_START_DELAY = 10  # time to prep data until collecting
OBS_DAWDLE = 5  # extra time to "sleep" to make sure things are done
LO_LIST = ['A', 'B']
DEFAULTS = {'observer': None, 'project_id': None,
            'conlog': 'WARNING', 'filelog': 'INFO', 'path': '.', 'log_filename': LOG_FILENAME,
            'obs_start_delay': OBS_START_DELAY, 'obs_dawdle': OBS_DAWDLE}

class CommandHandler:
    def __init__(self, **kwargs):
        kw = copy(DEFAULTS)
        kw.update(kwargs)
        self.logset = logger_setup.Logger(logger, conlog=kw['conlog'], filelog=kw['filelog'],
                                          log_filename=kw['log_filename'], path=kw['path'])
        logger.debug(f"{__name__} ver. {__version__}")

        for this_attr in ['observer', 'project_id', 'obs_start_delay', 'obs_dawdle']:
            setattr(self, this_attr, kw[this_attr])
        logger.info(f"session start: {self.observer} -- {self.project_id}")
        self.rec = ono_record.Record(observer=self.observer, project_id=self.project_id)

    def setants(self, ant_list='rfsoc_active', remove_ants=[], park_when_done=True):
        if ant_list.startswith('rfsoc_active'):
            self.ant_list = snap_config.get_rfsoc_active_antlist()
            xxx = ant_list.split('-')
            if len(xxx) == 2:
                remove_ants += tools.listify(xxx[1])
        elif ant_list in OS.ANT_LISTS:
            self.ant_list = copy(OS.ANT_LIST[ant_list])
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

    def setrf(self, freq, lo=['A', 'B'], attenuation=[8, 8], focus_on=None):
        """
        Parameters (kwargs)
        -------------------
        frequency : list of floats
            Frequency of observation [GHz], must match length of lo
        lo : list
            LO to use A/B/C/D
        attenuation : list of int
            Attenuation setting

        """
        self.freq = tools.listify(freq)
        self.lo = tools.listify(lo)
        self.attenuation = tools.listify(attenuation)
        antlo_list = [ant+lo.upper() for lo in lo for ant in self.ant_list]
        flim = 10000.0 if focus_on is None else max(self.freq)
        logger.info(f"freq: {', '.join([str(x) for x in self.freq])}")
        logger.info(f"lo: {', '.join(self.lo)}")
        logger.info(f"attenuation: {', '.join([str(x) for x in self.attenuation])}")
        need_to_focus = False
        focus_freq = None
        for freq, lo in zip(self.freq, self.lo):
            this_freq = [freq] * len(self.ant_list)
            need_to_focus = freq > (0.99 * flim ) if not need_to_focus  else False
            if need_to_focus: focus_freq = freq
            ata_control.set_freq(this_freq, self.ant_list, lo=lo.lower(), nofocus=freq<flim)
        if need_to_focus:
            import time
            time.sleep(20)
        ata_control.autotune(self.ant_list)
        ata_control.set_atten_thread([[f'{ant}x', f'{ant}y'] for ant in self.ant_list],
                                     [[self.attenuation[0], self.attenuation[1]] for ant in self.ant_list])
        snap_if.tune_if_antslo(antlo_list)
        self.rec.update(freq=freq, lo=lo, attenuation=attenuation, focus=focus_freq)

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
            subprocess.run("ansible-playbook /home/sonata/src/ansible_playbooks/hashpipe/xgpu_record.yml")
        else:
            logger.error(f"Invalid backend: {self.backend} -- no action")
            return
        subprocess.run(f"/home/sonata/src/observing_campaign/backend_setup_scripts/set_keys_uvh5_{self.project_id}.py")
        self.rec.update(backend=backend)

    def move(self, location=None, coord_type='azel', use_ants=None):
        """
        Parameters (kwargs)
        -------------------
        location : str
            target location (depends on coord_type)
        coord_type : str
            coordination type (azel, radec, source, traj)
        use_ants : list
            List of antennas to move, default to group_ants

        """
        self.location = location
        self.coord_type = coord_type
        use_ants = self.ant_list if use_ants is None else tools.listify(use_ants)

        if ',' in self.location:
            x, y = [float(_v) for _v in self.location.split(',')]
        logger.info(f'move to: {self.location}  {self.coord_type}')
        self.rec.update(x=x, y=y, coord=coord_type)

        if self.coord_type == 'azel':
            ata_control.set_az_el(use_ants, x, y)
            logger.info(f"azel: {x},{y}")
        elif self.coord_type == 'radec':
            source = ata_control.track_source(use_ants, radec=[x, y])
            logger.info(f"radec: {x},{y}")
        elif self.coord_type == 'source':
            source = ata_control.track_source(use_ants, source=self.location)
            logger.info(f"source: {self.location}")
        elif self.coord_type == 'traj':
            from obsnerd.trajectory_engine import TRACK_YAML_FILENAME
            ephem = ata_control.upload_ephemeris(self.location)
            ata_control.track_ephemeris(ephem, use_ants, wait=True)
            logger.info(f"traj: {self.location}")
            try:
                with open(TRACK_YAML_FILENAME, 'r') as fp:
                    for line in fp:
                        if len(line) > 2:
                            logger.info(f"track: {line.strip()}")
            except FileNotFoundError:
                pass

    def observe(self, obs_time, time_per_int):
        d = hpguppi_defaults.hashpipe_targets_LoA.copy()
        d.update(hpguppi_defaults.hashpipe_targets_LoB)
        #d = {'seti-node3': [0], 'seti-node6': [1]}
        keyval_dict = {'XTIMEINT': time_per_int}
        hpguppi_auxillary.publish_keyval_dict_to_redis(keyval_dict, d, postproc=True)
        hpguppi_record_in.record_in(self.obs_start_delay, obs_time, hashpipe_targets = d)
        sleepy = float(self.obs_start_delay + obs_time + self.obs_dawdle)
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

