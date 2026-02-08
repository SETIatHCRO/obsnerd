from .on_track import Track
from . import on_sys
import astropy.units as u
from os.path import isfile
from param_track import param_track_timetools as ttools
import json


FREQ_UNIT_DEFAULT = 'MHz'
TUNING_BANDWIDTH_DEFAULT = 690.0

class Obsinfo:
    def __init__(self, filename=None, contents_only=False):
        if filename is not None:
            self.read(filename, contents_only=contents_only)

    def read(self, filename, contents_only=False):
        """
        Read an obsinfo file.

        Parameters
        ----------
        filename : str
            Filename of the obsinfo file to be read

        """
        self.filename = filename
        self.sources = {}
        self.ants = None

        if filename is None:
            return
        if not (filename.endswith('json') or filename.endswith('yaml') or filename.endswith('yml')):
            print(f"File must be a json or yaml file - you provided {filename}")
            return
        if not isfile(filename):
            raise FileNotFoundError(f"{filename} not found.")

        with open(filename, 'r') as fp:
            if filename.endswith('json'):
                self.contents = json.load(fp)
            elif filename.endswith('yaml') or filename.endswith('yml'):
                import yaml
                self.contents = yaml.safe_load(fp)

        if contents_only:
            return

        for key, value in self.contents.items():
            if key != 'Sources':
                setattr(self, key.lower(), value)
        self.proc_freqs()
        self.proc_filters()
        self.proc_sources()

    def proc_sources(self):
        """ Process sources from configuration. """
        if 'Sources' not in self.contents:
            return
        units = getattr(self, 'units', None)
        for key, value in self.contents['Sources'].items():
            self.sources[key] = Track(source=key, units=units)
            value = {k: ttools.interpret_date(v) if k in ['start', 'stop'] else v for k, v in value.items()}
            self.sources[key].update(**value)
            for extra in ['off_time', 'off_angle']:
                try:
                    setattr(self.sources[key], extra, value[extra])
                except KeyError:
                    pass

    def proc_freqs(self):
        """ Process frequency list from configuration. """
        self.tunings = self.contents.get('Tunings', None)
        self.freq_unit = u.Unit(self.contents.get('Freq_unit', FREQ_UNIT_DEFAULT))
        self.tuning_bandwidth = self.contents.get('Tuning_Bandwidth', TUNING_BANDWIDTH_DEFAULT) * self.freq_unit
        self.lo, self.freqs = [], []
        self.filters = {}
        if self.tunings is not None:
            for lo in sorted(self.tunings.keys()):
                self.tunings[lo] = self.tunings[lo] * self.freq_unit
                self.lo.append(lo)
                self.freqs.append(self.tunings[lo])
                self.filters[lo] = {}

    def proc_filters(self):
        """ Process filter settings from configuration. """
        self.filter_list = self.contents.get('Filters', [])
        for filt in self.filter_list:
            for i in range(len(filt['band'])):
                filt['band'][i] = filt['band'][i] * self.freq_unit
            for lo, freq in self.tunings.items():
                if filt['band'][0] >= freq - self.tuning_bandwidth / 2.0 and filt['band'][1] <= freq + self.tuning_bandwidth / 2.0:
                    self.filters[lo][filt['color']] = filt['band']
    
    def write_track_plan_to_obsinfo(self, obsinfofn, config_file, tracks):
        from astropy.coordinates import angular_separation
        import numpy as np
        self.read(config_file, contents_only=True)
        obsinfo = {'filename': obsinfofn, 'dir_data': 'data', "Sources": {}}
        for k, v in self.contents.items():
            obsinfo[k] = v
        timed_tracks = {}
        for satname in tracks:
            for i, track in enumerate(tracks[satname]):
                if track.iobs is not None:
                    timed_tracks[track.tstart.datetime] = {'satname': satname, 'index': i}
        for track_info in sorted(timed_tracks.keys()):
            satname = timed_tracks[track_info]['satname']
            ind = timed_tracks[track_info]['index']
            track = tracks[satname][ind]
            obsinfo['Sources'][track.source] = {'obsid': on_sys.make_obsid(track.source, track.time[track.istart].mjd)}
            obsinfo['Sources'][track.source]['ra'] = track.ra[track.iobs].to_value('deg')
            obsinfo['Sources'][track.source]['dec'] = track.dec[track.iobs].to_value('deg')
            obsinfo['Sources'][track.source]['time'] = track.time[track.iobs].datetime.isoformat(timespec='seconds')
            obsinfo['Sources'][track.source]['az'] = track.az[track.iobs].to_value('deg')
            obsinfo['Sources'][track.source]['el'] = track.el[track.iobs].to_value('deg')
            obsinfo['Sources'][track.source]['start'] = track.tstart.datetime.isoformat(timespec='seconds')
            obsinfo['Sources'][track.source]['stop'] = track.tstop.datetime.isoformat(timespec='seconds')
            obsinfo['Sources'][track.source]['off_time'] = []
            obsinfo['Sources'][track.source]['off_angle'] = []
            obsinfo['Sources'][track.source]['ods'] = True if track.use == 'y' else False
            for i in range(track.istart, track.istop+1):
                toff = (track.time[i] - track.tobs).to_value('sec')
                aoff = angular_separation(track.ra[i], track.dec[i], track.ra[track.iobs], track.dec[track.iobs]) * np.sign(toff)
                obsinfo['Sources'][track.source]['off_time'].append(np.round(toff, 1))
                obsinfo['Sources'][track.source]['off_angle'].append(np.round(aoff.to_value('deg'), 2))
        try:
            with open(obsinfofn, 'r') as fp:
                existing_obsinfo = json.load(fp)
            print(f"Updating {obsinfofn}")
        except FileNotFoundError:
            existing_obsinfo = {}
            print(f"Writing {obsinfofn}")
        existing_obsinfo.update(obsinfo)
        with open(obsinfofn, 'w') as fp:
            json.dump(existing_obsinfo, fp, indent=4)