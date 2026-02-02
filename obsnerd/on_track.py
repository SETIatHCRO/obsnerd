from odsutils import ods_timetools as ttools
import numpy as np
from param_track import Parameters
from param_track.param_track_support import listify
from . import DATA_PATH
from os.path import join
import yaml


class Track(Parameters):
    header = ['observer', 'project_name', 'project_id', 'ants', 'focus', 'time_per_int_sec', 'backend', 'focus', 'attenuation', 'coord']
    short = ['freq', 'source', 'x', 'y', 'start', 'end', 'obs_time_sec']
    some_dtype_lists = {'freq': float, 'lo': str, 'attenuation': int}
    use_def = {'y': "use with ods",
               'n': "use without ods",
               's': "skip (don't use)"}

    def __init__(self, **kwargs):
        self.track_field_structure = yaml.safe_load(open(join(DATA_PATH, 'track_parameters.yaml')))
        self.fields = list(self.track_field_structure['fields'].keys())
        self.ods_mapping = self.track_field_structure['ods_mapping']
        super().__init__(ptnote='Track parameters', ptinit=self.fields, pttype=False, ptverbose=False)
        self.ptinit(['iobs', 'iref', 'istart', 'istop'])
        self.ptinit(['ra', 'dec', 'az', 'el', 'dist', 'utc'])
        self.update(**kwargs)

    def view(self, fields_to_show=None):
        self.ptshow(vals_only=True, include_par=fields_to_show)

    def update(self, **kwargs):
        """
        Update parameters, applying listify to some known dtypes.
        If changing time, must include 2 of 3 (and only 2) of [start, end, obs_time_sec] to keep them consistent.

        """
        for key in set(kwargs.keys()).intersection({'ra', 'dec', 'az', 'el', 'dist'}):
            if isinstance(kwargs[key], str):
                kwargs[key] = float(kwargs[key])
            else:
                kwargs[key] = kwargs[key]
        if 'utc' in kwargs:
            kwargs['utc'] = ttools.interpret_date(kwargs['utc'])
        for key in set(kwargs.keys()).intersection(self.some_dtype_lists.keys()):
            kwargs[key] = listify(kwargs[key], dtype=self.some_dtype_lists[key])
        timekeys = set(kwargs.keys()).intersection({'start', 'end', 'stop', 'obs_time_sec'})
        checked_timekeys = timekeys.copy()
        for key in timekeys:
            if kwargs[key] is None or not kwargs[key]:
                checked_timekeys.remove(key)
            elif key == 'stop':
                kwargs['end'] = kwargs.pop('stop')
        if checked_timekeys == {'start', 'end'} or checked_timekeys == {'start', 'stop'}:
            kwargs['obs_time_sec'] = int((kwargs['end'] - kwargs['start']).to_value('sec'))
        elif checked_timekeys == {'start', 'obs_time_sec'}:
            kwargs['end'] = ttools.t_delta(kwargs['start'], kwargs['obs_time_sec'], 's')
        elif checked_timekeys == {'end', 'obs_time_sec'}:
            kwargs['start'] = ttools.t_delta(kwargs['end'], -1.0 * kwargs['obs_time_sec'], 's')
        elif len(checked_timekeys) > 0:
            raise ValueError("When updating time parameters, must include 2 of 3 (and only 2) of 'start', 'end', 'obs_time_sec' to keep them consistent.")
        
        self._pt_set(**kwargs)

    def calc_properties(self):
        self.duration = self.utc[-1] - self.utc[0]
        self.imax = np.argmax(self.el)
        dt = np.diff(self.utc.mjd) * 24 * 3600
        self.daz = self.az.diff().to_value('deg')
        wrap = np.where(abs(self.daz) > 180.0)
        dwrap = wrap[0] - 1
        self.daz[wrap] = self.daz[dwrap]
        self.azdot = self.daz / dt
        self.azdot = np.insert(self.azdot, 0, self.azdot[0])
        self.azdot = abs(self.azdot)
        self.eldot = self.el.diff().to_value('deg') / dt
        self.eldot = np.insert(self.eldot, 0, self.eldot[0])