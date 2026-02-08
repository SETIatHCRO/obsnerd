import numpy as np
from param_track import Parameters
from param_track.param_track_support import listify
from param_track import param_track_timetools as ttools


class Track(Parameters):
    header = ['observer', 'project_name', 'project_id', 'ants', 'focus', 'time_per_int_sec', 'backend', 'focus', 'attenuation', 'coord']
    short = ['freq', 'source', 'x', 'y', 'start', 'stop', 'obs_time']
    use_def = {'y': "use with ods",
               'n': "use without ods",
               's': "skip (don't use)"}

    def __init__(self, **kwargs):
        super().__init__(ptnote='Track parameters', ptinit='track_parameters.yaml', pttype=False, ptverbose=False)
        self.update(**kwargs)

    def view(self, fields_to_show=None):
        self.ptshow(vals_only=True, include_par=fields_to_show)

    def update(self, **kwargs):
        """
        Update parameters, applying listify to some known dtypes.
        If changing time, must include 2 of 3 (and only 2) of [start, stop, obs_time] to keep them consistent.

        """
        timekeys = set(kwargs.keys()).intersection({'start', 'stop', 'obs_time'})
        checked_timekeys = timekeys.copy()
        for key in timekeys:
            if kwargs[key] is None or not kwargs[key]:
                checked_timekeys.remove(key)
        if checked_timekeys == {'start', 'stop'}:
            kwargs['obs_time'] = kwargs['stop'] - kwargs['start']
        elif checked_timekeys == {'start', 'obs_time'}:
            kwargs['stop'] = ttools.t_delta(kwargs['start'], kwargs['obs_time'], 's')
        elif checked_timekeys == {'stop', 'obs_time'}:
            kwargs['start'] = ttools.t_delta(kwargs['stop'], -1.0 * kwargs['obs_time'], 's')
        elif len(checked_timekeys) > 0:
            raise ValueError("When updating time parameters, must include 2 of 3 (and only 2) of 'start', 'stop', 'obs_time' to keep them consistent.")
        
        self._pt_set(**kwargs)

    def calc_properties(self):
        self.duration = self.time[-1] - self.time[0]
        self.imax = np.argmax(self.el)
        dt = np.diff(self.time.mjd) * 24 * 3600
        self.daz = self.az.diff().to_value('deg')
        wrap = np.where(abs(self.daz) > 180.0)
        dwrap = wrap[0] - 1
        self.daz[wrap] = self.daz[dwrap]
        self.azdot = self.daz / dt
        self.azdot = np.insert(self.azdot, 0, self.azdot[0])
        self.azdot = abs(self.azdot)
        self.eldot = self.el.diff().to_value('deg') / dt
        self.eldot = np.insert(self.eldot, 0, self.eldot[0])