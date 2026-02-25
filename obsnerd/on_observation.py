import numpy as np
from param_track import Parameters
from param_track import param_track_timetools as ttools
from .on_sys import make_obsid


header = ['observer', 'project_name', 'project_id', 'ants', 'focus', 'time_per_int_sec', 'backend', 'focus', 'attenuation', 'coord']
short = ['freq', 'source', 'x', 'y', 'start', 'stop', 'obs_time']


class Observation(Parameters):
    def __init__(self, **kwargs):
        ptinit = kwargs.pop('ptinit', None)
        ptnote = kwargs.pop('ptnote', 'Observational parameters')
        super().__init__(ptnote=ptnote, ptinit=ptinit, pttype=False, ptverbose=False)
        # Do it again, to check for the 2 outta 3 thing
        self.update(**kwargs)

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
        self.ptadd(**kwargs)

    def set_obsid(self, **kwargs):
        time = kwargs.get('time', None)
        if time is None:
            try:
                time = self.time[self.istart].mjd
            except (AttributeError, IndexError, KeyError):
                time = None
        self._pt_set(obsid=make_obsid(self.source, time))

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