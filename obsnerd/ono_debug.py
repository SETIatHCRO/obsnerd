class Empty:
    def __init__(self, clsname, *args, **kwargs):
        self.clsname = clsname
        self.hashpipe_targets_LoA = {}
        self.hashpipe_targets_LoB = {}
    def reserve_antennas(self, a):
        print(f"{self.clsname} not loaded.")
    def release_antennas(self, a, b):
        print(f"{self.clsname} not loaded.")
    def set_freq(self, *args, **kwargs):
        print(f"{self.clsname} not loaded.")
    def autotune(self, a):
        print(f"{self.clsname} not loaded.")
    def make_and_track_ephems(self, a, b):
        print(f"{self.clsname} not loaded.")
    def record_in(self, *args, **kwargs):
        print(f"{self.clsname} not loaded.")
    def publish_keyval_dict_to_redis(self, *args, **kwargs):
        print(f"{self.clsname} not loaded.")
    def get_rfsoc_active_antlist(self):
        print(f"{self.clsname} not loaded.")
        return ['2b']
    def copy(self):
        print(f"{self.clsname} not loaded.")
    def update(self, x):
        print(f"{self.clsname} not loaded.")
    def tune_if_antslo(self, antlo_list):
        print(f"{self.clsname} not loaded.")