from datetime import datetime
from time import sleep


class Empty:
    def __init__(self, clsname, *args, **kwargs):
        self.clsname = clsname
        self.hashpipe_targets_LoA = {}
        self.hashpipe_targets_LoB = {}
    def show(self, method):
        ts = datetime.now().strftime('%H:%M:%S')
        print(f"\t{ts} -- {self.clsname}.{method} not loaded.")
    def reserve_antennas(self, a):
        self.show('reserve_antennas')
    def release_antennas(self, a, b):
        self.show('release_antennas')
    def set_freq(self, *args, **kwargs):
        self.show('set_freq')
    def autotune(self, a):
        self.show('autotune')
    def make_and_track_ephems(self, a, b):
        self.show('make_and_track_ephems')
    def record_in(self, *args, **kwargs):
        self.show('record_in')
        print("\t\tPretending to start up ...")
        sleep(args[0])
        ts = datetime.now().strftime('%H:%M:%S')
        print(f"\t\t{ts} Pretending to record")
    def publish_keyval_dict_to_redis(self, *args, **kwargs):
        self.show('publish_keyval_dict_to_redis')
    def move_ant_group(self, a, b, c):
        self.show('move_ant_group')
    def get_rfsoc_active_antlist(self):
        self.show('get_rfsoc_active_antlist')
        return ['1k', '2a', '2b']
    def copy(self):
        self.show('copy')
    def update(self, x):
        self.show('update')
    def tune_if_antslo(self, antlo_list):
        self.show('tune_if_antslo')
    def set_atten_thread(self, a, b):
        self.show('set_atten_thread')
    def track_source(self, a, **kwargs):
        self.show('track_source')
    def delete_catalog_entry(self, **kwargs):
        self.show('delete_catalog')
    def add_catalog_entry(self, **kwargs):
        self.show('add_catalog')