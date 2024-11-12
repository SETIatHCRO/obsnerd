from astropy.time import Time, TimeDelta
from astropy.coordinates import EarthLocation, AltAz, SkyCoord
import matplotlib.pyplot as plt
import numpy as np
import astropy.units as u
from argparse import Namespace
from datetime import datetime, timedelta


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
    def __init__(self, fn=None):
        if fn is not None:
            print(f"Starlink ephemerides:  {fn}")
            self.fn = fn

    def readfeph(self, fn, show_plot=False):
        """
        This reads the "feph" json files written out below in self.filter.
        Manually added the Tx times from the file later from SpaceX

        """
        import json
        self.feph = Satellite([], tref=[], t0=[], t1=[])
        self.sources = {}
        with open(fn, 'r') as fp:
            fephjson = json.load(fp)
            for src, data in fephjson.items():
                data['tref'] = Time(data['tref'])
                self.feph.satno.append(src)
                for dval in ['tref', 'ra', 'dec', 'az', 'el']:
                    getattr(self.feph, dval).append(data[dval])
                data['t0'] = Time(data['t0']) if 't0' in data else 0
                data['t1'] = Time(data['t1']) if 't1' in data else 0
                self.sources[src] = Satellite(src)
                for dval in data:
                    setattr(self.sources[src], dval, data[dval])
        self.feph.tref = Time(self.feph.tref)
        if show_plot:
            plt.plot(self.feph.tref.datetime, self.feph.az, 'o')
            plt.plot(self.feph.tref.datetime, self.feph.el, 'x')

    def write_sorted_obs_file(self, obslen=6.0, startup=15.0, tz=0.0, fn='observe.dat'):
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
        with open('observe.dat', 'w') as fp:
            print("#Source,pause[s],obs[s],start,ref,ra,dec,az,el", file=fp)
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
                print(f"{sorter[key].satno},{pause:.1f},{obslen*60.0:.1f},{startat.isoformat()},{key.isoformat()},", end='', file=fp)
                print(f"{sorter[key].ra},{sorter[key].dec},{sorter[key].az},{sorter[key].el}", file=fp)

    def reada(self):
        """
        Read the version supplied for first set in October 2024

        """
        self.sats = {}
        with open(self.fn, 'r') as fp:
            self.T0 = Time(int(fp.readline().strip()), scale='utc', format='gps')
            print(f"T0 = {self.T0}")
            for line in fp:
                if line[0] != '(':
                    this_sat = int(line.strip())
                    self.sats[this_sat] = Satellite(satno=line.strip())
                elif line[0] == '(':
                    for _E in line.strip().split('('):
                        if len(_E):
                            #print(_E.strip(',').strip(')').split(','))
                            t, ra, dec = [float(x) for x in _E.strip(',').strip(')').split(',')]
                            self.sats[this_sat].times.append(self.T0 + TimeDelta(t, format='sec'))
                            #self.sats[this_sat].times.append(t)
                            if ra < 0.0:
                                ra += 360.0
                            self.sats[this_sat].ra.append(ra)
                            self.sats[this_sat].dec.append(dec)
        for this_sat in self.sats:
            self.sats[this_sat].ra = np.array(self.sats[this_sat].ra)
            self.sats[this_sat].dec = np.array(self.sats[this_sat].dec)

    def readb(self):
        """
        Read the version supplied for second set in November 2024.

        """
        self.sats = {}
        with open(self.fn, 'r') as fp:
            self.T0 = Time(int(fp.readline().strip()), scale='utc', format='gps')
            print(f"T0 = {self.T0}")
            for line in fp:
                if not line[0].isspace():
                    this_sat = int(line.strip())
                    self.sats[this_sat] = Satellite(satno=line.strip())
                else:
                    t, tstr, rastr, ra, decstr, dec, azstr, az, elstr, el = line.split()
                    self.sats[this_sat].times.append(self.T0 + TimeDelta(float(t), format='sec'))
                    self.sats[this_sat].provided.times.append(tstr)
                    ra, dec, az, el = float(ra.strip(',')), float(dec.strip(',')), float(az.strip(',')), float(el.strip(','))
                    if ra < 0.0:
                        ra += 360.0
                    if az < 0.0:
                        az += 360.0
                    self.sats[this_sat].ra.append(ra)
                    self.sats[this_sat].dec.append(dec)
                    self.sats[this_sat].provided.az.append(az)
                    self.sats[this_sat].provided.el.append(el)
        for this_sat in self.sats:
            self.sats[this_sat].ra = np.array(self.sats[this_sat].ra)
            self.sats[this_sat].dec = np.array(self.sats[this_sat].dec)
            self.sats[this_sat].provided.az = np.array(self.sats[this_sat].provided.az)
            self.sats[this_sat].provided.el = np.array(self.sats[this_sat].provided.el)
            self.sats[this_sat].times = Time(self.sats[this_sat].times)
        
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

    def plot(self, ptype='radec'):
        for sat, data in self.sats.items():
            if ptype[0].lower() == 'r':
                plt.plot(data.ra, data.dec)
                plt.xlabel('RA')
                plt.ylabel('Dec')
            elif ptype[0].lower() == 'a':
                plt.plot(data.az, data.el)
                plt.xlabel('Az')
                plt.ylabel('El')
            elif ptype[0].lower() == 'p':
                plt.plot(data.provided.az, data.provided.el)
                plt.xlabel('Az')
                plt.ylabel('El')
            elif ptype[0].lower() == 'd':
                plt.plot(data.az - data.provided.az, data.el - data.provided.el, '.')
                plt.xlabel('dAz')
                plt.ylabel('dEl')    

    def get_azel(self, satno='all', observatory=ATA):
        if satno == 'all':
            satno = list(self.sats.keys())
        for this_satno in satno:
            aa = AltAz(location=observatory, obstime=self.sats[this_satno].times)
            coord = SkyCoord(self.sats[this_satno].ra * u.deg, self.sats[this_satno].dec * u.deg)
            altazsky = coord.transform_to(aa)
            self.sats[this_satno].az = altazsky.az.value
            self.sats[this_satno].el = altazsky.alt.value

    def filter(self, satno='all', ytime=['2024-11-06T23:00:00', '2024-11-07T01:00:00'], yel=[40, 50]):
        import json
        now = datetime.now()
        fn = f"feph_{now.strftime('%Y%m%dT%H%M')}.json"
        import string
        tags = string.ascii_lowercase
        if satno == 'all':
            satno = list(self.sats.keys())
        if ytime is not None:
            ytime[0] = Time(ytime[0], format='isot')
            ytime[1] = Time(ytime[1], format='isot')
        fephjson = {}
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
                    fephjson[source_name] = {
                        'az': self.sats[this_sat].az[i],
                        'el': self.sats[this_sat].el[i],
                        'ra': self.sats[this_sat].ra[i],
                        'dec': self.sats[this_sat].dec[i],
                        'tref': self.sats[this_sat].times[i].datetime.isoformat()
                    }
                    #print(f"S{this_sat}{tags[i]},{self.sats[this_sat].ra[i] / 15.0:.6f},{self.sats[this_sat].dec[i]:.6f}", file=fp)
                    print(f"S{this_sat}{tags[i]},{self.sats[this_sat].az[i]:.6f},{self.sats[this_sat].el[i]:.6f}")
                    #print(f"S{this_sat}{tags[i]},{self.sats[this_sat].times[i].datetime},{self.sats[this_sat].ra[i] / 15.0:.6f},{self.sats[this_sat].dec[i]:.6f},", file=fp, end='')
                    #print(f"{self.sats[this_sat].az[i]:.6f},{self.sats[this_sat].el[i]:.6f}", file=fp)
                    ctr += 1
        with open(fn, 'w') as fp:
            json.dump(fephjson, fp, indent=2)
        print(f"Wrote {ctr} to {fn}")