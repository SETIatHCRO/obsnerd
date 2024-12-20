import atexit
from . import obs_sys as OS
import time
import logging
from datetime import datetime
from sys import stdout
from astropy.time import Time


try:
    from ATATools import ata_control, logger_defaults
    from SNAPobs import snap_if
    from SNAPobs.snap_hpguppi import record_in as hpguppi_record_in
    from SNAPobs.snap_hpguppi import snap_hpguppi_defaults as hpguppi_defaults
    from SNAPobs.snap_hpguppi import auxillary as hpguppi_auxillary
    from SNAPobs import snap_config
except ImportError:
    from .obs_debug import Empty
    ata_control = Empty()
    logger_defaults = Empty()
    snap_if = Empty()
    hpguppi_record_in = Empty()
    hpguppi_defaults = Empty()
    hpguppi_auxillary = Empty()
    snap_config = Empty()


def augment(arr, N, msg="Lengths must be equal"):
    if not isinstance(arr, list):
        return [arr] * N
    if len(arr) == 1:
        return arr * N
    if len(arr) != N:
        raise ValueError(msg)
    return arr


class Observer:
    def __init__(self, sources, integrations, start_times=None, freqs={'a': 1410, 'b': 5500}, ant_list='rfsoc_active', known_bad='', focus_on='max', obs_fiddle=5, data_record='hpguppi'):
        """
        Parameters
        ----------
        sources : list of str, str
            Names of sources
        integrations : list of float/int
        start_time : list of str of times in isoformat or delays in seconds
        freqs: 
        focus_on : str
            lo ('a', 'b') or 'max'

        """
        self.sources = OS.listify(sources)
        self.integrations = augment([float(x) for x in OS.listify(integrations)], len(self.sources), msg="Number of integrations and sources should match.")
        self.start_times = self.sched(augment(OS.listify(start_times), len(self.sources), msg="Number of start_times and sources should match."))
        self.focus_on = focus_on
        self.logger = logger_defaults.getProgramLogger("observe", loglevel=logging.INFO)
        self.obs_delay = OS.OBS_START_DELAY + obs_fiddle
        self.data_record = data_record
        
        self.ant_list = OS.listify(ant_list, {'rfsoc_active': snap_config.get_rfsoc_active_antlist()})
        if not len(self.ant_list):
            raise ValueError("No antennas specified.")
        for badun in OS.listify(known_bad):
            if badun in self.ant_list:
                print(f"Removing antenna {badun}")
                ant_list.remove(badun)
        self.freqs = OS.dictionify(freqs)
        # Check/set freqs same for all antennas
        max_freq = {'val': 0.0, 'lo': ''}
        for lo in self.freqs:
            self.freqs[lo] = OS.listify(self.freqs[lo])
            if len(self.freqs[lo]) == 1:
                self.freqs[lo] = self.freqs[lo] * len(self.ant_list)
            if self.freqs[lo][0] > max_freq['val']:
                max_freq['val'] = self.freqs[lo][0]
                max_freq['lo'] = lo
        if self.focus_on == 'max':
            self.focus_on = max_freq['lo']

        self.antlo_list = [ant+lo.upper() for lo in self.freqs for ant in self.ant_list]

    def sched(self, tlist):
        """
        Step through self.start_times to make delays and datetimes

        """
        starts = []
        for t in tlist:
            try:
                d = float(t)
            except ValueError:
                d = Time(t).datetime
            starts.append(d)
        return starts

    def setup_session(self):
        yn = input("Have you set up the backend? ")
        if yn[0].lower() == 'n':
            raise RuntimeError("You need to set up the backend.")
        ata_control.reserve_antennas(self.ant_list)
        atexit.register(ata_control.release_antennas, self.ant_list, False)

        self.dnodes = {}
        for lo in self.freqs:
            self.dnodes.update(getattr(hpguppi_defaults, f"hashpipe_targets_Lo{lo.upper()}").copy())

        for lo in self.freqs:
            nofocus = True if lo == self.focus_on else False
            ata_control.set_freq(self.freqs[lo], self.ant_list, lo=lo, nofocus=nofocus)
        print("Focusing feed.")
        time.sleep(OS.FEED_FOCUS_SLEEP)

        ata_control.autotune(self.ant_list)
        snap_if.tune_if_antslo(self.antlo_list)

    def step_obs(self, skip_cal_check=False):
        cal_done = skip_cal_check
        for source, integration, start in zip(self.sources, self.integrations, self.start_times):
            if isinstance(start, datetime):
                start = (start - datetime.now()).total_seconds()
            if start > 1.0:
                print(f"Waiting {start} seconds")
                time.sleep(start)
            self.take_obs(source, integration)
            if not cal_done:
                cal_done = input("Have you calibrated? ")

    def take_obs(self, source, integration):
        keyval_dict = {'XTIMEINT': integration}
        hpguppi_auxillary.publish_keyval_dict_to_redis(keyval_dict, self.dnodes, postproc=False)
        ata_control.make_and_track_ephems(source, self.ant_list)

        hpguppi_record_in.record_in(OS.OBS_START_DELAY, integration, hashpipe_targets = self.dnodes)
        print(f"Observing {source} for {integration:.2f} seconds...", end='')
        stdout.flush()
        time.sleep(integration + self.obs_delay)
        if self.data_record == 'hpguppi':
            hpguppi_record_in.record_in(OS.OBS_START_DELAY, integration, hashpipe_targets = self.dnodes)
        print("Done")
