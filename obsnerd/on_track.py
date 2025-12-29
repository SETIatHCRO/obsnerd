from odsutils import ods_timetools as ttools
from odsutils import ods_tools as tools
import numpy as np
from argparse import Namespace
import json



def read_obsinfo(filename):
    """
    Read an obsinfo file.

    Parameters
    ----------
    filename : str
        Filename of the obsinfo file to be read

    """
    if filename is None:
        return Namespace(filename=None, sources={})
    obsinfo = Namespace(filename=filename, sources={})
    with open(filename, 'r') as fp:
        data = json.load(fp)

    for key, value in data.items():
        if key != 'Sources':
            setattr(obsinfo, key.lower(), value)
    for key, value in data['Sources'].items():
        obsinfo.sources[key] = Track(source=key)
        obsinfo.sources[key].set_track(**value)
        for extra in ['off_time', 'off_angle']:
            try:
                setattr(obsinfo.sources[key], extra, value[extra])
            except KeyError:
                pass
    return obsinfo


class Track:
    fields = [
        'observer', 'project_name', 'project_id', 'ants', 'freq', 'lo', 'attenuation', 'focus', 'backend',
        'source', 'x', 'y', 'coord',
        'start', 'end', 'obs_time_sec', 'time_per_int_sec'
    ]
    header = ['observer', 'project_name', 'project_id', 'ants', 'focus', 'time_per_int_sec', 'backend', 'focus', 'attenuation', 'coord']
    short = ['freq', 'source', 'x', 'y', 'start', 'end', 'obs_time_sec']
    use_def = {'y': "use with ods",
               'n': "use without ods",
               's': "skip (don't use)"}
    ods_mapping = {
        'source': 'src_id',
        'x': 'src_ra_j2000_deg',
        'y': 'src_dec_j2000_deg',
        'start': 'src_start_utc',
        'end': 'src_end_utc',
        'time_per_int_sec': 'corr_integ_time_sec',
        'freq0': 'obs_freq_lo_mhz',
        'freq1': 'obs_freq_hi_mhz'}

    def __init__(self, **kwargs):
        self.update(**kwargs)
        self.iobs = None

    def __repr__(self, fprnt='fields'):
        return self.view(fields_to_show=getattr(self, fprnt), bracket=['<', '>'], sep='  -- ')

    def __str__(self):
        return self.view(fields_to_show=self.fields, bracket=['', ''], sep='\n')

    def to_dict(self):
        rec = {}
        for fld in self.fields:
            rec[fld] = self.listify(fld)

    def update(self, **kwargs):
        for key, val in kwargs.items():
            if key in self.fields:
                setattr(self, key, val)

    def listify(self, key, dtype=None):
        try:
            setattr(self, key, tools.listify(getattr(self, key), dtype=dtype))
        except AttributeError:
            pass

    def proc(self):
        t = {}
        for key in ['start', 'end', 'obs_time_sec']:
            t[key] = getattr(self, key, None)
        if t['obs_time_sec'] is None or t['obs_time_sec'] == '-':
            self.obs_time_sec = int((self.end - self.start).to_value('sec'))
        elif t['end'] is None or t['end'] == '-':
            self.end = ttools.t_delta(self.start, self.obs_time_sec, 's')
        elif t['start'] is None or t['start'] == '-':
            self.start = ttools.t_delta(self.end, -1.0 * self.obs_time_sec, 's')
        for key, dtype in {'freq': None, 'lo': str, 'attenuation':int}.items():
            self.listify(key, dtype)

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
    
    def index_tracker(self, **kwargs):
        """
        If the track is part of a bigger set, then this will keep track of the index in the bigger set.

        """
        for par, val in kwargs.items():
            if par not in ['iref', 'istart', 'istop']:
                continue
            setattr(self, par, val)

    def set_track(self, **kwargs):
        for par, val in kwargs.items():
            if par in ['ra', 'dec', 'az', 'el', 'dist']:
                if isinstance(val, str):
                    this_val = float(val)
                else:
                    this_val = val
            elif par == 'utc':
                this_val = ttools.interpret_date(val)
            else:
                this_val = val
            if this_val is not None:  # Probably not necessary anymore
                setattr(self, par, this_val)

    def set_par(self, **kwargs):
        for par, val in kwargs.items():
            setattr(self, par, val)

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