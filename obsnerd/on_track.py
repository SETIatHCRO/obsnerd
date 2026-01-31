from odsutils import ods_timetools as ttools
import numpy as np
from argparse import Namespace
from param_track import Parameters
from param_track.param_track_support import listify
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


class Track(Parameters):
    fields = [
        'observer', 'project_name', 'project_id', 'ants', 'freq', 'lo', 'attenuation', 'focus', 'backend',
        'source', 'x', 'y', 'coord',
        'start', 'end', 'obs_time_sec', 'time_per_int_sec'
    ]

    header = ['observer', 'project_name', 'project_id', 'ants', 'focus', 'time_per_int_sec', 'backend', 'focus', 'attenuation', 'coord']
    short = ['freq', 'source', 'x', 'y', 'start', 'end', 'obs_time_sec']
    some_dtype_lists = {'freq': float, 'lo': str, 'attenuation': int}
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
        super().__init__(ptnote='Track parameters', ptinit=self.fields, pttype=False, ptverbose=False, **kwargs)
        self.ptadd(iobs=None)

    def view(self, fields_to_show=None):
        self.ptshow(vals_only=True, include_par=fields_to_show)

    def update(self, **kwargs):
        """
        Update parameters, applying listify to some known dtypes.
        If changing time, must include 2 of 3 (and only 2) of [start, end, obs_time] to keep them consistent.

        """
        dtypekeys = set(kwargs.keys()).intersection(self.some_dtype_lists.keys())
        for key in dtypekeys:
            kwargs[key] = listify(kwargs[key], dtype=self.some_dtype_lists[key])
        timekeys = set(kwargs.keys()).intersection({'start', 'end', 'obs_time_sec'})
        if timekeys == {'start', 'end'}:
            kwargs['obs_time_sec'] = int((ttools.t_delta(kwargs['end'], kwargs['start'])).to_value('sec'))
        elif timekeys == {'start', 'obs_time_sec'}:
            kwargs['end'] = ttools.t_delta(kwargs['start'], kwargs['obs_time_sec'], 's')
        elif timekeys == {'end', 'obs_time_sec'}:
            kwargs['start'] = ttools.t_delta(kwargs['end'], -1.0 * kwargs['obs_time_sec'], 's')
        elif len(timekeys) > 0:
            raise ValueError("When updating time parameters, must include 2 of 3 (and only 2) of 'start', 'end', 'obs_time' to keep them consistent.")
        
        self._pt_set(**kwargs)

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