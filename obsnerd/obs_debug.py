class Empty:
    def __init__(self, *args, **kwargs):
        self.hashpipe_targets_LoA = {}
        self.hashpipe_targets_LoB = {}
    def reserve_antennas(self, a):
        pass
    def release_antennas(self, a, b):
        pass
    def set_freq(self, *args, **kwargs):
        pass
    def autotune(self, a):
        pass
    def make_and_track_ephems(self, a, b):
        pass
    def getProgramLogger(self, *args, **kwargs):
        pass
    def record_in(self, *args, **kwargs):
        pass
    def publish_keyval_dict_to_redis(self, *args, **kwargs):
        pass
    def get_rfsoc_active_antlist(self):
        return ['2b']
    def tune_if_antslo(self, a):
        pass
