try:
    from ATATools import ata_control  # type: ignore
except ImportError:
    print("ATATools not found.")
    ata_control = None
try:
    from SNAPobs import snap_config  # type: ignore
except ImportError:
    print("SNAPobs not found.")
    snap_config = None
import atexit
from . import metadata, onutil
import subprocess
from copy import copy


ANT_LISTS = {'old_feeds': ['']}


class CommandHandler:
    def __init__(self, **kwargs):
        for key, val in kwargs.items():
            setattr(self, key, val)
        self._setantlists()

    def _setvar(self, initargs, kwargs):
        for key, val in initargs.items():
            if key in kwargs:  # Definitely change if kwarg
                setattr(self, key, kwargs[key])
            elif not hasattr(self, key):  # Only set to initialize value if not already set
                setattr(self, key, val)
        self._setantlists()

    def _setantlists(self, known_lists=['rfsoc_active']):
        known_lists += list(ANT_LISTS.keys())
        for ant_list in ['group_ants', 'use_ants']:
            if hasattr(self, ant_list):
                this_list = copy(getattr(self, ant_list))
                if isinstance(this_list, list):
                    pass
                elif this_list in known_lists:
                    if this_list in ANT_LISTS.keys():
                        setattr(self, ant_list, ANT_LISTS[this_list])
                    elif this_list == 'rfsoc_active':
                        if snap_config is None:
                            print(f"snap_config not loaded to get {this_list}")
                            setattr(self, ant_list, None)
                        else:
                            setattr(self, ant_list, snap_config.get_rfsoc_active_antlist())
                elif isinstance(this_list, str):
                    setattr(self, ant_list, getattr(self, ant_list).split(','))
                else:
                    setattr(self, ant_list, None)
        if hasattr(self, 'use_ants') and self.use_ants is None:
            self.use_ants = self.group_ants

    def start(self, **kwargs):
        """
        Parameters (kwargs)
        -------------------
        initials : str
            Initials of observer
        project_id : str
            Project identifier
        group_ants : list
            list of ants to move into the reserved group

        """
        self._setvar({'initials': None, 'project_id': None, 'group_ants': ['1a']}, kwargs)

        if self.initials is None:
            print("Please include your name or initials.")
        elif ata_control is None:
            self.test(f'test start:  {self.initials} for project {self.project_id}')
        else:
            ata_control.move_ant_group(self.group_ants, 'none', 'atagr')
        metadata.logger.info(f"session start: {self.initials} -- reserving {', '.join(self.group_ants)}")
        metadata.logger.info(f"project_id: {self.project_id}")

    def end(self, **kwargs):
        """
        Parameters (kwargs)
        -------------------
        group_ants : list
            Antenna list in the group to unreserve
        use_ants : list
            Antenna list to park, default use group_ants

        """
        self._setvar({'group_ants': ['1a'], 'use_ants': None}, kwargs)

        if ata_control is None:
            self.test("test end")
        else:
            atexit.register(ata_control.park_antennas, self.use_ants)
            atexit.register(ata_control.move_ant_group, self.group_ants, 'atagr', 'none')
        metadata.logger.info(f"parking: {', '.join(self.use_ants)}")
        metadata.logger.info(f"end: {', '.join(self.group_ants)}")

    def freq(self, **kwargs):
        """
        Parameters (kwargs)
        -------------------
        frequency : float
            Frequency of observation [GHz]
        lo : str
            LO to use A/B/C/D
        use_ants : list or None
            List of antennas to use, if None, assume group_ants
        attenuation : int
            Attenuation setting

        """
        self._setvar({'frequency': None, 'lo': 'A', 'use_ants': None, 'attenuation': 20, 'group_ants': ['1a']}, kwargs)
        if self.frequency is None:
            print("Need a frequency in GHz -- no action")
            return

        metadata.logger.info(f"fcen: {self.frequency}")
        metadata.logger.info(f"lo: {self.lo}")
        metadata.logger.info(f"antennas:  {(', ').join(self.use_ants)}")
        if ata_control is None:
            self.test('Test freq')
        else:
            ata_control.set_freq(self.frequency, self.use_ants, lo=self.lo.lower())
            ata_control.autotune(self.use_ants)
            ata_control.rf_switch_thread(self.use_ants)
            ata_control.set_atten_thread([[f'{ant}x', f'{ant}y'] for ant in self.use_ants], [[self.attenuation, self.attenuation] for ant in self.use_ants])
    
    def backend(self, **kwargs):
        """
        Parameters (kwargs)
        -------------------
        backend : str
            Name of backend to use
        project_id : str
            ID of project

        """
        self._setvar({'backend': 'xgpu', 'project_id': None}, kwargs)
        if self.project_id is None:
            print("Need a project_id")
            return
        metadata.logger.info(f"Backend {self.backend} for project {self.project_id}")
        if ata_control is None:
            print("Test backend")
        else:
            if self.backend == 'xgpu':
                subprocess.run("ansible-playbook /home/sonata/src/ansible_playbooks/hashpipe/xgpu_record.yml")
            else:
                print(f"Invalid backend: {self.backend} -- no action")
                return
            subprocess.run(f"/home/sonata/src/observing_campaign/backend_setup_scripts/set_keys_uvh5_{self.project_id}.py")

    def observe(self, **kwargs):
        print("GET THE GUPPI?ETC STUFF HERE")

    def move(self, **kwargs):
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
        self._setvar({'location': None, 'coord_type': 'azel', 'use_ants': None, 'group_ants': ['1a']}, kwargs)
        if self.location is None:
            print("Need target")
            return

        if ',' in self.location:
            x, y = [float(_v) for _v in self.location.split(',')]
        metadata.logger.info(f'move to: {self.location}  {self.coord_type}')
        if ata_control is None:
            self.test('test move')
            return

        if self.coord_type == 'azel':
            ata_control.set_az_el(self.use_ants, x, y)
            metadata.logger.info(f"azel: {x},{y}")
        elif self.coord_type == 'radec':
            source = ata_control.track_source(self.use_ants, radec=[x, y])
            metadata.logger.info(f"radec: {x},{y}")
        elif self.coord_type == 'source':
            source = ata_control.track_source(self.use_ants, source=self.location)
            metadata.logger.info(f"source: {self.location}")
        elif self.coord_type == 'traj':
            from obsnerds.trajectory_engine import TRACK_YAML_FILENAME
            ephem = ata_control.upload_ephemeris(self.location)
            ata_control.track_ephemeris(ephem, self.use_ants, wait=True)
            metadata.logger.info(f"traj: {self.location}")
            try:
                with open(TRACK_YAML_FILENAME, 'r') as fp:
                    for line in fp:
                        if len(line) > 2:
                            metadata.logger.info(f"track: {line.strip()}")
            except FileNotFoundError:
                pass
    
    def note(self, **kwargs):
        """
        Parameter (kwargs)
        ------------------
        notation : str
            Note to log

        """
        self._setvar({'notation': None}, kwargs)
        if self.notation is None:
            print("Need to include a note.")
        else:
            metadata.logger.info(f"note: {self.notation}")

    def source(self, **kwargs):
        """
        Parameters (kwargs)
        -------------------
        source : str
            Source name
        datestamp : str
            Datestamp

        """
        self._setvar({'name': None, 'datestamp': None}, kwargs)
        self.datestamp = onutil.make_datetime(date=self.datestamp)
        metadata.logger.info([f'source: {self.name}', f'expected: {self.datestamp.isoformat()}'])

    def summary(self, **kwargs):
        summary = metadata.get_summary()

    def test(self, msg='test'):
        print(msg)
        metadata.logger.info(msg)