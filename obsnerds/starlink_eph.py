from astropy.time import Time
from astropy.coordinates import EarthLocation, AltAz, SkyCoord
import matplotlib.pyplot as plt
import numpy as np
import astropy.units as u
from argparse import Namespace
from datetime import datetime, timedelta
from copy import copy


#ATA = EarthLocation(lat = 40.8171 * u.deg, lon = -121.469 * u.deg, height = 1050.0 * u.m)
ATA = EarthLocation(lat = 40.814871 * u.deg, lon = -121.469001 * u.deg, height = 1019.0 * u.m)  # From SpaceX
DEGSEC = 1.5
LATENCY = 30.0


class Satellite:
    def __init__(self, satno, **kwargs):
        self.satno = satno
        self.times = []
        self.ra = []
        self.dec = []
        self.az = []
        self.el = []
        self.provided = Namespace(times=[], az=[], el=[])
        for key, val in kwargs.items():
            setattr(self, key, val)


class Eph:
    def __init__(self):
        print("Handle Starlink observation ephemerides etc")

    def read_SpaceX(self, fn, ftype='c'):
        """
        Since the files from SpaceX vary so much, pull the readers out but produce common self.sats dictionary.

        """
        from . import starlink_input
        self.SpaceX = fn
        if ftype == 'a':
            self.sats = starlink_input.reada(fn)
        if ftype == 'b':
            self.sats = starlink_input.readb(fn)
        if ftype == 'c':
            self.sats = starlink_input.readc(fn)
        self.get_azel()


    def read_feph(self, fn, show_plot=False):
        """
        This reads the "feph" json files written out below in self.filter.  Makes a 'feph' Namespace with
        2 attributes: 'array' and 'sources'

        """
        import json
        self.feph = Namespace(array=Satellite(satno=[]), sources={})
        parameters = set()
        print(f"Reading {fn}")
        with open(fn, 'r') as fp:
            feph_json_input = json.load(fp)
            for src, data in feph_json_input['Sources'].items():
                self.feph.array.satno.append(src)
                self.feph.sources[src] = Satellite(satno=src)
                for key, value in data.items():
                    parameters.add(key)
                    if key[0] == 't': # Is Time format
                        try:
                            value = Time(value)
                        except ValueError:
                            raise ValueError(f"Parameters starting with 't' must be astropy.time.Time format:  {key}: {value}")
                    try:
                        getattr(self.feph.array, key).append(value)
                    except AttributeError:
                        setattr(self.feph.array, key, [value])
                    setattr(self.feph.sources[src], key, value)
        parameters.add('satno')
        print(f"\t{', '.join(list(parameters))}")
        for par in parameters:
            if par[0] == 't':
                value = Time(copy(getattr(self.feph.array, par)))
                setattr(self.feph.array, par, value)
                
        if show_plot:
            plt.plot(self.feph.array.tref.datetime, self.feph.array.az, 'o')
            plt.plot(self.feph.array.tref.datetime, self.feph.array.el, 'x')

    def write_sorted_obs_file(self, obslen=6.0, startup=15.0, tz=0.0, fn='observe.dat', fmt=2):
        """
        obslen : float
            Length of observation in minutes
        startup : float
            Number of seconds it takes to start taking data
        tz : float
            Hours from UTC (-8.0 for PST)

        """
        from copy import copy
        sorter = {}
    
        for src, data in self.sources.items():
            sorter[data.tref.datetime] = Satellite(satno=src, ra=data.ra, dec=data.dec, az=data.az, el=data.el)
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
                    acquire = abs(daz) / DEGSEC + LATENCY
                    delay = (key - prev_time).total_seconds() - 60.0 * obslen
                    pause = delay - acquire
                else:
                    acquire = 0.0
                    delay = 0.0
                    pause = delay - acquire
                if pause < 0.0:
                    print(f"Not enough to acquire {sorter[key].satno}")
                prev_az = sorter[key].az
                prev_time = copy(key)
                startat = key - timedelta(hours=-1.0*tz, minutes=obslen/2.0, seconds=startup)
                if fmt == 1:
                    print(f"{sorter[key].satno},{pause:.1f},{obslen*60.0:.1f},{startat.isoformat()},{key.isoformat()},", end='', file=fp)
                    print(f"{sorter[key].ra},{sorter[key].dec},{sorter[key].az},{sorter[key].el}", file=fp)
                elif fmt == 2:
                    print(f'["{sorter[key].satno}"], ["{startat.isoformat()}"]', file=fp)
                elif fmt == 3:
                    ra = 360.0 + sorter[key].ra if sorter[key].ra < 0.0 else sorter[key].ra
                    ra = ra / 15.0
                    print(f"{sorter[key].satno}: [{ra:.4f},{sorter[key].dec}]", file=fp)

    def update_feph_boresight(self, fn):
        """
        This assumes you've read in the reada/readb stuff

        """
        import json
        self.readfeph(fn)
        with open(fn, 'r') as fp:
            feph_input_file = json.load(fp)
        for src in self.sources:
            satno = int(src[1:-1])  # strip off the first and last letters
            self.angular_sep(satno, self.sources[src].ra, self.sources[src].dec, ra_unit='deg')
            suptimes = [np.round(x, 1) for x in (self.sats[satno].times - self.sources[src].tref).to(u.second).value]
            feph_input_file['Sources'][src]['suptimes'] = suptimes
            off_boresight = [np.round(x, 3) for x in self.sats[satno].angsep.value]
            feph_input_file['Sources'][src]['off_boresight'] = off_boresight
        self.write_feph(fn, feph_input_file)

    def angular_sep(self, satno, ra, dec, ra_unit='hourangle'):
        from astropy.coordinates import angular_separation
        ra = ra * getattr(u, ra_unit)
        dec = dec * u.deg
        self.sats[satno].angsep = []
        for i in range(len(self.sats[satno].ra)):
            self.sats[satno].angsep.append(angular_separation(ra, dec,
                                           self.sats[satno].ra[i] * u.deg, self.sats[satno].dec[i] * u.deg).to(u.deg).value)
        self.sats[satno].angsep = np.array(self.sats[satno].angsep) * u.deg

    def angular_speed(self, satno, using='provided'):
        from astropy.coordinates import angular_separation
        self.sats[satno].angvel = [0.0]  # Just to keep the number same as times
        if using[0].lower() == 'p':
            print('Using provided azel')
            for i in range(1, len(self.sats[satno].provided.az)):
                az0 = self.sats[satno].provided.az[i-1] * u.deg
                el0 = self.sats[satno].provided.el[i-1] * u.deg
                az1 = self.sats[satno].provided.az[i] * u.deg
                el1 = self.sats[satno].provided.el[i] * u.deg
                dt = self.sats[satno].times[i] - self.sats[satno].times[i-1]
                self.sats[satno].angvel.append(
                    (angular_separation(az0, el0, az1, el1) / dt).to(u.deg / u.second).value)
        self.sats[satno].angvel = np.array(self.sats[satno].angvel) * u.deg / u.second  

    def get_azel(self, satno=None, observatory=ATA):
        if satno is None:
            satno = list(self.sats.keys())
        for this_satno in satno:
            aa = AltAz(location=observatory, obstime=self.sats[this_satno].times)
            coord = SkyCoord(self.sats[this_satno].ra * u.deg, self.sats[this_satno].dec * u.deg)
            altazsky = coord.transform_to(aa)
            self.sats[this_satno].az = altazsky.az.value
            self.sats[this_satno].el = altazsky.alt.value

    def filter(self, satno=None, ytime=['2024-11-06T23:00:00', '2024-11-15T01:00:00'], yel=[30, 70]):
        """
        use None, None, None to write out everything.

        """
        import string
        tags = string.ascii_lowercase
        if satno is None:
            satno = list(self.sats.keys())
        if ytime is not None:
            ytime[0] = Time(ytime[0], format='isot')
            ytime[1] = Time(ytime[1], format='isot')
        fephjson = {'Sources': {}}
        ctr = 0
        for this_sat in satno:
            for i in range(len(self.sats[this_sat].times)):
                use_it = True
                if ytime is not None:
                    this_time = self.sats[this_sat].times[i]
                    if this_time < ytime[0] or this_time > ytime[1]:
                        use_it = False
                if use_it and yel is not None:
                    this_el = self.sats[this_sat].el[i]
                    if this_el < yel[0] or this_el > yel[1]:
                        use_it = False
                source_name = f"S{this_sat}{tags[i]}"
                if use_it:
                    fephjson['Sources'][source_name] = {
                        'sat': this_sat,
                        'az': self.sats[this_sat].az[i],
                        'el': self.sats[this_sat].el[i],
                        'ra': self.sats[this_sat].ra[i],
                        'dec': self.sats[this_sat].dec[i],
                        'tref': self.sats[this_sat].times[i].datetime.isoformat(),
                        'norad': self.sats[this_sat].norad[i],
                        'bf_distance': self.sats[this_sat].bf_distance[i]
                    }
                    #print(f"S{this_sat}{tags[i]},{self.sats[this_sat].ra[i] / 15.0:.6f},{self.sats[this_sat].dec[i]:.6f}", file=fp)
                    print(f"S{this_sat}{tags[i]},{self.sats[this_sat].az[i]:.6f},{self.sats[this_sat].el[i]:.6f}")
                    #print(f"S{this_sat}{tags[i]},{self.sats[this_sat].times[i].datetime},{self.sats[this_sat].ra[i] / 15.0:.6f},{self.sats[this_sat].dec[i]:.6f},", file=fp, end='')
                    #print(f"{self.sats[this_sat].az[i]:.6f},{self.sats[this_sat].el[i]:.6f}", file=fp)
                    ctr += 1
        print(f"{ctr} entries")
        now = datetime.now()
        fn = f"feph_{now.strftime('%Y%m%dT%H%M')}.json"
        self.write_feph(fn=fn, feph_dict=fephjson)

    def write_feph(self, fn, feph_dict):
        import json
        print(f"Writing {fn}")
        with open(fn, 'w') as fp:
            json.dump(feph_dict, fp, indent=2)