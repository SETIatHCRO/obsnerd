from astropy.time import Time
from astropy.coordinates import EarthLocation, AltAz, SkyCoord
import numpy as np
import astropy.units as u
from argparse import Namespace
from datetime import timedelta
from copy import copy
import json


class Observatory:
    def __init__(self, name, location, **kwargs):
        self.name = name
        self.location = location
        for key, val in kwargs.items():
            setattr(self, key, val)

# From Wael
latitude = "40:49:02.75"
longitude = "-121:28:14.65"
altitude = 1019.222
ATA = Observatory('ATA',
                  # location = EarthLocation(lat = 40.814871 * u.deg, lon = -121.469001 * u.deg, height = 1019.0 * u.m),  # From SpaceX
                  location = EarthLocation(lat = 40.81743055556 * u.deg, lon = -121.5292638889 * u.deg, height = 1019.222 * u.m),  # From Wael
                  slew_speed = 1.5,  # deg/sec
                  startup = 15.0,  # number of seconds to actually start taking dat after started, sec
                  latency = 30.0  # time to acquire etc before/after slew, sec
)

class ObservationData:
    def __init__(self, name, **kwargs):
        self.name = name
        self.times = []
        self.ra = []
        self.dec = []
        self.az = []
        self.el = []
        for key, val in kwargs.items():
            setattr(self, key, val)


class Base:
    def __init__(self, observatory=ATA, obsmeta_file='obsmeta.json'):
        """
        Parameters
        ----------
        observatory : class Observatory
            Observatory information
        obsmeta_file : str
            Name of file containing meta data about the obsfile json files.

        """
        self.observatory = observatory
        self.obsmeta_file = obsmeta_file
        try:
            with open(obsmeta_file, 'r') as fp:
                self.obsmeta_data = json.load(fp)
        except FileNotFoundError:
            self.obsmeta_data = {}
        print(f"Handle observation ephemerides, tools etc for {observatory.name}")

    def read_obsinfo(self, obs):
        """
        This reads the observation meta json files generated.
        Makes an 'obsinfo' Namespace with 3 attributes: 'filename', 'array' and 'obsid'

        Parameter
        ---------
        obsid : str or None
            If a .json extension, reads that file.  If not, it checks that obsid as found in the obsid_file
        
        """
        self.obsinfo = Namespace(array=ObservationData(name=[]), obsid={}, dir_data='.')
        import json
        if obs.endswith('.json'):
            self.obsinfo.filename = obs  # Assume that this is a full obsinfo filepath/name
        else:
            from os.path import join
            dir2use = '.'
            for key, val in self.obsmeta_data.items():
                if key == 'dir_obsinfo':
                    dir2use = val
                elif obs in val:
                    self.obsinfo.filename = join(dir2use, key)
                    break
            else:
                self.obsinfo.filename = None
                print(f"No obs file found for {obs}")
                return
        
        parameters = {'name'}
        self.obsinfo.obsinfo_keys = ['obsinfo_keys']
        print(f"Reading {self.obsinfo.filename}")
        with open(self.obsinfo.filename, 'r') as fp:
            json_input = json.load(fp)
            for key, val in json_input.items():  # Get non-'Sources' info
                self.obsinfo.obsinfo_keys.append(key)
                if key != 'Sources':
                    setattr(self.obsinfo, key.lower(), val)
            for src, data in json_input['Sources'].items():
                self.obsinfo.array.name.append(src)
                self.obsinfo.obsid[src] = ObservationData(name=src)
                for key, value in data.items():
                    parameters.add(key)
                    if key[0] == 't': # Is Time format
                        try:
                            value = Time(value)
                        except ValueError:
                            raise ValueError(f"Parameters starting with 't' must be astropy.time.Time format:  {key}: {value}")
                    try:
                        getattr(self.obsinfo.array, key).append(value)
                    except AttributeError:
                        setattr(self.obsinfo.array, key, [value])
                    setattr(self.obsinfo.obsid[src], key, value)
        print(f"\t{', '.join(self.obsinfo.array.name)}")
        print(f"\t{', '.join(list(parameters))}")
        for par in [x for x in parameters if x[0]=='t']:
            value = Time(copy(getattr(self.obsinfo.array, par)))
            setattr(self.obsinfo.array, par, value)

    def write_sorted_obs_file(self, obslen=6.0, tz=0.0, fn='observe.dat', fmt=2):
        """
        obslen : float
            Length of observation in minutes
        tz : float
            Hours from UTC (-8.0 for PST)
        fn : str
            Name of output file
        fmt : int
            Just a way to keep track of arbitrary output format

        """
        from copy import copy
        sorter = {}
    
        for src, data in self.obsid.items():
            sorter[data.tref.datetime] = ObservationData(name=src, ra=data.ra, dec=data.dec, az=data.az, el=data.el)
        with open(fn, 'w') as fp:
            if fmt == 1:
                print("#Source,pause[s],obs[s],start,ref,ra,dec,az,el", file=fp)
            elif fmt == 2:
                print("#Source, Start")
            elif fmt == 3:
                print("#Source,ra,dec")
            for i, key in enumerate(sorted(sorter.keys())):
                if i:
                    daz = sorter[key].az - prev_az
                    acquire = abs(daz) / self.observatory.slew_speed + self.observatory.latency
                    delay = (key - prev_time).total_seconds() - 60.0 * obslen
                    pause = delay - acquire
                else:
                    acquire = 0.0
                    delay = 0.0
                    pause = delay - acquire
                if pause < 0.0:
                    print(f"Not enough to acquire {sorter[key].name}")
                prev_az = sorter[key].az
                prev_time = copy(key)
                startat = key - timedelta(hours=-1.0*tz, minutes=obslen/2.0, seconds=self.observatory.startup)
                if fmt == 1:
                    print(f"{sorter[key].name},{pause:.1f},{obslen*60.0:.1f},{startat.isoformat()},{key.isoformat()},", end='', file=fp)
                    print(f"{sorter[key].ra},{sorter[key].dec},{sorter[key].az},{sorter[key].el}", file=fp)
                elif fmt == 2:
                    print(f'["{sorter[key].name}"], ["{startat.isoformat()}"]', file=fp)
                elif fmt == 3:
                    ra = 360.0 + sorter[key].ra if sorter[key].ra < 0.0 else sorter[key].ra
                    ra = ra / 15.0
                    print(f"{sorter[key].name}: [{ra:.4f},{sorter[key].dec}]", file=fp)

    def find_sats_gen(self, efn, start_seconds=200, tlefile='tle/starlink.tle', duration=6.0):
        self.read_obsinfo(efn)
        with open('find_sats.sh', 'w') as fp:
            for src in self.obsinfo.obsid:
                start = self.obsinfo.obsid[src].tref.datetime - timedelta(seconds=start_seconds)
                print(f'find_sats.py -t "{start.isoformat()}" --tle_file {tlefile} -d {duration} --output_file --sat2write {src[1:-1]}', file=fp)

    def boresight_update_view(self, efn):
        import matplotlib.pyplot as plt
        self.read_obsinfo(efn)
        for src in self.obsinfo.obsid:
            plt.plot(self.obsinfo.obsid[src].off_times, self.obsinfo.obsid[src].off_boresight)

    def update_obsinfo_boresight(self, fn):
        """
        You need to have run SOPP and written satellite output files.

        """
        from astropy.coordinates import angular_separation
        self.read_obsinfo(fn)
        import json
        with open(fn, 'r') as fp:
            obsinfo_file_dict = json.load(fp)
        for src in self.obsinfo.obsid:
            satfile = f"{src[:-1]}.txt"  # strip off the last letter
            try:
                fp = open(satfile, 'r')
                print(f"Reading {satfile}")
            except FileNotFoundError:
                print(F"{satfile} not found")
                continue
            self.obsinfo.obsid[src].sopp = Namespace(dt=[], az=[], el=[], dist=[])
            ctr = 0
            center_found = False
            min_dt_off = 1E6
            for line in fp:
                data = [float(x) for x in line.strip().split(',')[1:]]
                t = Time(line.split(',')[0])
                dt = (t - self.obsinfo.obsid[src].tref).to('second').value
                if abs(dt) < abs(min_dt_off):
                    min_dt_off = dt
                if abs(dt) < 1E-6:
                    center_found = True
                    self.obsinfo.obsid[src].sopp.d_az = data[0] - self.obsinfo.obsid[src].az
                    self.obsinfo.obsid[src].sopp.d_el = data[1] - self.obsinfo.obsid[src].el
                    if abs(self.obsinfo.obsid[src].sopp.d_az) > 0.1 or abs(self.obsinfo.obsid[src].sopp.d_az) > 0.1:
                        print(f"{src} has large excursions:  {self.obsinfo.obsid[src].sopp.d_az}, {self.obsinfo.obsid[src].sopp.d_el}")
                    az0 = self.obsinfo.obsid[src].az * u.deg
                    el0 = self.obsinfo.obsid[src].el * u.deg
                if abs(dt) < 1E-6 or not ctr % 10:
                    self.obsinfo.obsid[src].sopp.dt.append(dt)
                    self.obsinfo.obsid[src].sopp.az.append(data[0])
                    self.obsinfo.obsid[src].sopp.el.append(data[1])
                    self.obsinfo.obsid[src].sopp.dist.append(data[2])
                ctr += 1
            if center_found:
                for i in range(len(self.obsinfo.obsid[src].sopp.az)):  # "Fix" for offset
                    self.obsinfo.obsid[src].sopp.az[i] -= self.obsinfo.obsid[src].sopp.d_az
                    self.obsinfo.obsid[src].sopp.el[i] -= self.obsinfo.obsid[src].sopp.d_el
                self.obsinfo.obsid[src].sopp.angsep = []
                for i in range(len(self.obsinfo.obsid[src].sopp.az)):  # "Fix" for offset
                    angsep = float(angular_separation(az0, el0, self.obsinfo.obsid[src].sopp.az[i]*u.deg, self.obsinfo.obsid[src].sopp.el[i]*u.deg).to(u.deg).value)
                    self.obsinfo.obsid[src].sopp.angsep.append(angsep * np.sign(self.obsinfo.obsid[src].sopp.dt[i]))
                obsinfo_file_dict['Sources'][src]['off_times'] = np.round(self.obsinfo.obsid[src].sopp.dt, 1).tolist()
                obsinfo_file_dict['Sources'][src]['off_boresight'] = np.round(self.obsinfo.obsid[src].sopp.angsep, 2).tolist()
                obsinfo_file_dict['Sources'][src]['off_distance'] = np.round(self.obsinfo.obsid[src].sopp.dist, 1).tolist()
                obsinfo_file_dict['Sources'][src]['off_az_offset'] = self.obsinfo.obsid[src].sopp.d_az
                obsinfo_file_dict['Sources'][src]['off_el_offset'] = self.obsinfo.obsid[src].sopp.d_el
                obsinfo_file_dict['Sources'][src]['off_min_dt'] = min_dt_off
            else:
                print(f"Center was not found - min dt offset was {min_dt_off}")
        self.write_obsinfo(fn, obsinfo_file_dict)

    def get_azel(self, name=None):
        if name is None:
            name = list(self.sats.keys())
        for this_name in name:
            aa = AltAz(location=self.observatory.location, obstime=self.sats[this_name].times)
            coord = SkyCoord(self.sats[this_name].ra * u.deg, self.sats[this_name].dec * u.deg)
            altazsky = coord.transform_to(aa)
            self.sats[this_name].az = altazsky.az.value
            self.sats[this_name].el = altazsky.alt.value

    def write_obsinfo(self, fn=None, package='obsinfo'):
        if isinstance(package, str):
            print("Package up self.obsinfo into dict...not done yet")
            return
        if fn is None:
            mjd = float(Time.now().jd) - 2400000.5
            fn = f"obsinfo_{mjd:.4f}.json"
        import json
        print(f"Writing {fn}")
        with open(fn, 'w') as fp:
            json.dump(package, fp, indent=2)