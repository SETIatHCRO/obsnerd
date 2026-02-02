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
        self._parse_track_field_structure()
        super().__init__(ptnote='Track parameters', ptinit=self.fields, pttype=False, ptverbose=False)
        self.update(**kwargs)

    def _parse_track_field_structure(self):
        self.track_field_structure = yaml.safe_load(open(join(DATA_PATH, 'track_parameters.yaml')))
        self.fields = list(self.track_field_structure['fields'].keys())
        self.ods_mapping = self.track_field_structure['ods_mapping']
        self.field_types = {}
        for key, value in self.track_field_structure['fields'].items():
            for tval in value.get('type', []):
                self.field_types.setdefault(tval, set())
                self.field_types[tval].add(key)

    def view(self, fields_to_show=None):
        self.ptshow(vals_only=True, include_par=fields_to_show)

    def update(self, **kwargs):
        """
        Update parameters, applying listify to some known dtypes.
        If changing time, must include 2 of 3 (and only 2) of [start, end, obs_time_sec] to keep them consistent.

        """
        newargs = {}
        for key, value in kwargs.items():
            dtype_info = self.track_field_structure[key]['type']
            if key in self.field_types['list']:
                if len(dtype_info) == 2:
                    ind = 0 if self.track_field_structure[key]['type'][1] == 'list' else 1
                    newargs[key] = listify(value, dtype=dtype_info[ind])
                else:
                    newargs[key] = listify(value)
            elif len(dtype_info) == 1:
                if dtype_info[0] in ['Time', 'TimeDelta']:
                    newargs[key] = ttools.interpret_date(value)
                else:
                    newargs[key] = eval(dtype_info[0])(value)
            else:
                newargs[key] = value

        timekeys = set(newargs.keys()).intersection({'start', 'end', 'stop', 'obs_time_sec'})
        checked_timekeys = timekeys.copy()
        for key in timekeys:
            if newargs[key] is None or not newargs[key]:
                checked_timekeys.remove(key)
            elif key == 'stop':
                newargs['end'] = newargs.pop('stop')
        if checked_timekeys == {'start', 'end'} or checked_timekeys == {'start', 'stop'}:
            newargs['obs_time_sec'] = int((newargs['end'] - newargs['start']))
        elif checked_timekeys == {'start', 'obs_time_sec'}:
            newargs['end'] = ttools.t_delta(newargs['start'], newargs['obs_time_sec'], 's')
        elif checked_timekeys == {'end', 'obs_time_sec'}:
            newargs['start'] = ttools.t_delta(newargs['end'], -1.0 * newargs['obs_time_sec'], 's')
        elif len(checked_timekeys) > 0:
            raise ValueError("When updating time parameters, must include 2 of 3 (and only 2) of 'start', 'end', 'obs_time_sec' to keep them consistent.")
        
        self._pt_set(**newargs)

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