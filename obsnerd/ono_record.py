from odsutils import ods_timetools as ttools
from odsutils import ods_tools as tools

class Record:
    fields = [
        'observer', 'project_name', 'project_id', 'ants', 'freq', 'lo', 'attenuation', 'focus', 'backend',
        'source_name', 'x', 'y', 'coord',
        'start', 'end', 'obs_time_sec', 'time_per_int_sec'
    ]
    header = ['observer', 'project_name', 'project_id', 'ants', 'focus', 'time_per_int_sec', 'backend', 'focus', 'attenuation', 'coord']
    short = ['freq', 'source_name', 'x', 'y', 'start', 'end', 'obs_time_sec']

    def __init__(self, **kwargs):
        self.update(**kwargs)

    def __repr__(self, use='fields'):
        return self.view(fields_to_show=getattr(self, use), bracket=['<', '>'], sep='  -- ')

    def __str__(self):
        return self.view(fields_to_show=self.fields, bracket=['', ''], sep='\n')

    def update(self, **kwargs):
        for key, val in kwargs.items():
            if key in self.fields:
                setattr(self, key, val)

    def _listify(self, key, dtype=None):
        try:
            setattr(self, key, tools.listify(getattr(self, key), dtype=dtype))
        except AttributeError:
            pass

    def _get_or_None(self, key):
        try:
            return getattr(self, key)
        except AttributeError:
            return None

    def proc(self):
        t = {}
        for key in ['start', 'end', 'obs_time_sec']:
            t[key] = self._get_or_None(key)
        if t['obs_time_sec'] is None or t['obs_time_sec'] == '-':
            self.obs_time_sec = int((self.end - self.start).to_value('sec'))
        elif t['end'] is None or t['end'] == '-':
            self.end = ttools.t_delta(self.start, self.obs_time_sec, 's')
        elif t['start'] is None or t['start'] == '-':
            self.start = ttools.t_delta(self.end, -1.0 * self.obs_time_sec, 's')
        for key, dtype in {'freq': None, 'lo': str, 'attenuation':int}.items():
            self._listify(key, dtype)

    def view(self, fields_to_show=None, bracket=['', ''], sep='\n'):
        if fields_to_show is None:
            fields_to_show = self.fields
        s = bracket[0]
        for fld in fields_to_show:
            try:
                val = getattr(self, fld)
                if isinstance(val, list):
                    val = ', '.join([str(x) for x in val])
                s += f'{fld}: {val}{sep}'
            except AttributeError:
                continue
        s = s.strip(sep)
        s += bracket[1]
        return s
    