from astropy.time import Time, TimeDelta
from astropy.coordinates import EarthLocation, AltAz, SkyCoord
import matplotlib.pyplot as plt
import numpy as np
import astropy.units as u
from argparse import Namespace
from datetime import datetime, timedelta


#ATA = EarthLocation(lat = 40.8171 * u.deg, lon = -121.469 * u.deg, height = 1050.0 * u.m)
ATA = EarthLocation(lat = 40.814871 * u.deg, lon = -121.469001 * u.deg, height = 1019.0 * u.m)
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

    def readfeph(self, fn='feph_241106.out', show_plot=False):
        """
        This reads the "feph" files written out below in self.filter.
        Manually added the Tx times from the file later from SpaceX
        """
        self.feph = Satellite([], t0=[], t1=[])
        self.sources = {}
        sorter = {}
        with open(fn, 'r') as fp:
            for line in fp:
                offset = 0
                data = line.split(',')
                if len(data) == 6:
                    sn, ti, ra, dec, az, el = data
                    t0, t0n, t1n = None, None, None
                elif len(data) == 8:
                    sn, ti, ra, dec, az, el, t0, t1 = data
                elif len(data) == 9:
                    sn, ti, ra, dec, az, el, t0, t1, offset = data
                    #offset = int(offset)
                self.feph.satno.append(sn)
                self.feph.times.append(ti)
                ra, dec, az, el = float(ra), float(dec), float(az), float(el)
                self.feph.ra.append(ra)
                self.feph.dec.append(dec)
                self.feph.az.append(az)
                self.feph.el.append(el)
                sorter[Time(ti).datetime] = [sn, ra, dec, az, el]
                if t0 is not None:
                    self.feph.t0.append(t0)
                    self.feph.t1.append(t1)
                    t0n, t1n = Time(t0), Time(t1)
                self.sources[sn] = Satellite(sn, ra=ra, dec=dec, az=az, el=el, times=Time(ti), t0=t0n, t1=t1n, offset=offset)
        self.feph.times = Time(self.feph.times)
        if len(self.feph.t0):
            self.feph.t0 = Time(self.feph.t0)
            self.feph.t1 = Time(self.feph.t1)
        if show_plot:
            plt.plot(self.feph.times.datetime, self.feph.az, 'o')
            plt.plot(self.feph.times.datetime, self.feph.el, 'x')
        self.sorted = Satellite([])
        self.sorted.times = sorted(sorter.keys())
        self.sorted.slew = []
        self.sorted.delay = []
        with open('observe.dat', 'w') as fp:
            for i, key in enumerate(sorted(self.sorted.times)):
                self.sorted.satno.append(sorter[key][0])
                self.sorted.ra.append(sorter[key][1])
                self.sorted.dec.append(sorter[key][2])
                self.sorted.az.append(sorter[key][3])
                self.sorted.el.append(sorter[key][4])
                if not i:
                    self.sorted.slew.append(0.0)
                    self.sorted.delay.append(0.0)
                else:
                    self.sorted.slew.append(abs(self.sorted.az[i] - self.sorted.az[i-1]) / DEGSEC + LATENCY)
                    self.sorted.delay.append((self.sorted.times[i] - self.sorted.times[i-1]).total_seconds())
                print(f"{self.sorted.satno[i]},{self.sorted.times[i] - timedelta(hours=8.0, minutes=1.5, seconds=15)}", file=fp)
                #print(f"{self.sorted.satno[i][1:-1]},{self.sorted.times[i]}", file=fp)

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

    def angular_speed(self, satno, using='provided'):
        from astropy.coordinates import angular_separation
        self.sats[satno].angvel = []
        if using[0].lower() == 'p':
            print('Using provided azel')
            for i in range(1, len(self.sats[satno].provided.az)):
                az0 = self.sats[satno].provided.az[i-1] * u.deg
                el0 = self.sats[satno].provided.el[i-1] * u.deg
                az1 = self.sats[satno].provided.az[i] * u.deg
                el1 = self.sats[satno].provided.el[i] * u.deg
                dt = (self.sats[satno].times[i] - self.sats[satno].times[i-1]).to(u.second).value
                print(dt)
                self.sats[satno].angvel.append(angular_separation(az0, el0, az1, el1).to(u.deg).value / dt)

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
        fn = 'feph.out'
        import string
        tags = string.ascii_lowercase
        if satno == 'all':
            satno = list(self.sats.keys())
        if ytime is not None:
            ytime[0] = Time(ytime[0], format='isot')
            ytime[1] = Time(ytime[1], format='isot')
        ctr = 0
        with open(fn, 'w') as fp:
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
                    if use_it:
                        # print(f"S{this_sat}{tags[i]},{self.sats[this_sat].ra[i] / 15.0:.6f},{self.sats[this_sat].dec[i]:.6f}", file=fp)
                        print(f"S{this_sat}{tags[i]},{self.sats[this_sat].az[i]:.6f},{self.sats[this_sat].el[i]:.6f}")
                        print(f"S{this_sat}{tags[i]},{self.sats[this_sat].times[i].datetime},{self.sats[this_sat].ra[i] / 15.0:.6f},{self.sats[this_sat].dec[i]:.6f},", file=fp, end='')
                        print(f"{self.sats[this_sat].az[i]:.6f},{self.sats[this_sat].el[i]:.6f}", file=fp)
                        ctr += 1
        print(f"Wrote {ctr} to {fn}")