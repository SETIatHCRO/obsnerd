###DELETED 2/18/2025 --> became on_track...

from odsutils import ods_timetools as ttools
from param_track import Parameters
from param_track.param_track_support import listify

class Record(Parameters):
    fields = [
        'observer', 'project_id', 'ants', 'freq', 'lo', 'attenuation', 'focus', 'backend',
        'source', 'x', 'y', 'coord',
        'start', 'end', 'obs_time', 'time_per_int'
    ]
    header = ['observer', 'project_id', 'ants', 'focus', 'time_per_int', 'backend', 'focus', 'attenuation', 'coord']
    short = ['freq', 'source', 'x', 'y', 'start', 'end', 'obs_time']
    some_dtype_lists = {'freq': float, 'lo': str, 'attenuation': int}

    def __init__(self, **kwargs):
        super().__init__(ptnote='Ono Record', ptverbose=False, pttype=False, ptinit=self.fields, **kwargs)

    def update(self, **kwargs):
        """
        Update parameters, applying listify to some known dtypes.
        If changing time, must include 2 of 3 (and only 2) of [start, end, obs_time] to keep them consistent.

        """
        dtypekeys = set(kwargs.keys()).intersection(self.some_dtype_lists.keys())
        for key in dtypekeys:
            kwargs[key] = listify(kwargs[key], dtype=self.some_dtype_lists[key])
        timekeys = set(kwargs.keys()).intersection({'start', 'end', 'obs_time'})
        if timekeys == {'start', 'end'}:
            kwargs['obs_time'] = int((ttools.t_delta(kwargs['end'], kwargs['start'])).to_value('sec'))
        elif timekeys == {'start', 'obs_time'}:
            kwargs['end'] = ttools.t_delta(kwargs['start'], kwargs['obs_time'], 's')
        elif timekeys == {'end', 'obs_time'}:
            kwargs['start'] = ttools.t_delta(kwargs['end'], -1.0 * kwargs['obs_time'], 's')
        elif len(timekeys) > 0:
            raise ValueError("When updating time parameters, must include 2 of 3 (and only 2) of 'start', 'end', 'obs_time' to keep them consistent.")
        
        self._pt_set(**kwargs)

    def view(self, fields_to_show=None):
        self.ptshow(vals_only=True, include_par=fields_to_show)
