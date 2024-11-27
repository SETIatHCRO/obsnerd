try:
    from ATATools import ata_control
except ImportError:
    ata_control = None
import atexit
from . import metadata, onutil
import subprocess


class CommandHandler:
    def __init__(self, arg1=None, arg2=None, coord_type=None, group_ants=['1a'], use_ants=['1a']):
        """
        arg1 : str or None
            "arg1" for the empty command call
        arg2 : str or None
            "arg2" for the empty command call
        coord_type : str or None
            supplied coordinate type (azel, radec, source, traj)
        group_ants : list
            list of ants to move into the reserved group
        use_ants : list
            list of ants to point with command

        """
        self.arg1 = arg1
        self.arg2 = arg2
        self.coord_type = coord_type
        self.group_ants = group_ants
        self.use_ants = use_ants
    
    def start(self, initials=None):
        self.initials = self.arg1 if initials is None else initials
        if self.initials is None:
            print("Please include your name or initials.")
        elif ata_control is None:
            self.test(f'test start:  {self.initials}')
        else:
            ata_control.move_ant_group(self.group_ants, 'none', 'atagr')
            metadata.onlog(f"session start: {self.initials} -- {', '.join(self.use_ants)} / {', '.join(self.group_ants)}")

    def end(self):
        if ata_control is None:
            self.test("test end")
        else:
            atexit.register(ata_control.move_ant_group, self.group_ants, 'atagr', 'none')
            atexit.register(ata_control.park_antennas, self.use_ants)
            metadata.onlog(f"end: {', '.join(self.use_ants)} / {', '.join(self.group_ants)}")

    def freq(self, freq=None, att=20):
        self.freq = float(self.arg1) if freq is None else freq
        self.att = att
        metadata.onlog(f"fcen: {self.freq}")
        if ata_control is not None:
            ata_control.set_freq(self.freq, self.use_ants, lo='d')
            ata_control.autotune(self.use_ants)
            ata_control.rf_switch_thread(self.use_ants)
            ata_control.set_atten_thread([[f'{ant}x', f'{ant}y'] for ant in self.use_ants], [[self.att, self.att] for ant in self.use_ants])
    
    def backend(self, betype='xgpu', project_id='p054'):
        metadata.onlog(f"Backend {betype} for project {project_id}")
        if betype == 'xgpu':
            subprocess.run("ansible-playbook /home/sonata/src/ansible_playbooks/hashpipe/xgpu_record.yml")
        else:
            print(f"Invalid backend: {betype} -- no action")
            return
        subprocess.run(f"/home/sonata/src/observing_campaign/backend_setup_scripts/set_keys_uvh5_{project_id}.py")

    def move(self, location=None, coord_type=None):
        self.location = self.arg1 if location is None else location
        self.coord_type = self.coord_type if coord_type is None else coord_type

        if ',' in self.location:
            x, y = [float(_v) for _v in self.location.split(',')]
        if ata_control is None:
            self.test('test move')
            return

        metadata.onlog(f'move to: {self.coord_type}')
        if self.coord_type == 'azel':
            ata_control.set_az_el(self.use_ants, x, y)
            metadata.onlog(f"azel: {x},{y}")
        elif self.coord_type == 'radec':
            source = ata_control.track_source(self.use_ants, radec=[x, y])
            metadata.onlog(f"radec: {x},{y}")
        elif self.coord_type == 'source':
            source = ata_control.track_source(self.use_ants, source=self.location)
            metadata.onlog(f"source: {self.location}")
        elif self.coord_type == 'traj':
            from obsnerds.trajectory_engine import TRACK_YAML_FILENAME
            ephem = ata_control.upload_ephemeris(self.location)
            ata_control.track_ephemeris(ephem, self.use_ants, wait=True)
            metadata.onlog(f"traj: {self.location}")
            try:
                with open(TRACK_YAML_FILENAME, 'r') as fp:
                    for line in fp:
                        if len(line) > 2:
                            metadata.onlog(f"track: {line.strip()}")
            except FileNotFoundError:
                pass
    
    def note(self, notation=None):
        self.notation = self.arg1 if notation is None else notation
        if self.notation is None:
            print("Need to include a note.")
        else:
            metadata.onlog(f"note: {self.notation}")

    def source(self, name=None, datestamp=None):
        self.name = self.arg1 if name is None else name
        self.datestamp = self.arg2 if datestamp is None else datestamp
        self.datestamp = onutil.make_datetime(date=self.datestamp)
        metadata.onlog([f'source: {self.name}', f'expected: {self.datestamp.isoformat()}'])

    def summary(self, **kwargs):
        summary = metadata.get_summary()

    def test(self, msg='test'):
        print(msg)
        metadata.onlog(msg)