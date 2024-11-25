from astropy.time import Time, TimeDelta
from obsnerds.starlink_eph import ObsData
import numpy as np

"""
This has the various input stuff from the SpaceX provided satellite info.
"""

def reada(fn):
    """
    Read the version supplied for first set in October 2024

    """
    sats_inp = {}
    with open(fn, 'r') as fp:
        T0 = Time(int(fp.readline().strip()), scale='utc', format='gps')
        print(f"T0 = {T0}")
        for line in fp:
            if line[0] != '(':
                this_sat = int(line.strip())
                sats_inp[this_sat] = ObsData(name=line.strip())
            elif line[0] == '(':
                for _E in line.strip().split('('):
                    if len(_E):
                        #print(_E.strip(',').strip(')').split(','))
                        t, ra, dec = [float(x) for x in _E.strip(',').strip(')').split(',')]
                        sats_inp[this_sat].times.append(T0 + TimeDelta(t, format='sec'))
                        #sats_inp[this_sat].times.append(t)
                        if ra < 0.0:
                            ra += 360.0
                        sats_inp[this_sat].ra.append(ra)
                        sats_inp[this_sat].dec.append(dec)
    for this_sat in sats_inp:
        sats_inp[this_sat].ra = np.array(sats_inp[this_sat].ra)
        sats_inp[this_sat].dec = np.array(sats_inp[this_sat].dec)
    return sats_inp


def readb(fn):
    """
    Read the version supplied for second set in November 2024.

    """
    sats_inp = {}
    with open(fn, 'r') as fp:
        T0 = Time(int(fp.readline().strip()), scale='utc', format='gps')
        print(f"T0 = {T0}")
        for line in fp:
            if not line[0].isspace():
                this_sat = int(line.strip())
                sats_inp[this_sat] = ObsData(name=line.strip())
            else:
                t, tstr, rastr, ra, decstr, dec, azstr, az, elstr, el = line.split()
                sats_inp[this_sat].times.append(T0 + TimeDelta(float(t), format='sec'))
                sats_inp[this_sat].provided.times.append(tstr)
                ra, dec, az, el = float(ra.strip(',')), float(dec.strip(',')), float(az.strip(',')), float(el.strip(','))
                if ra < 0.0:
                    ra += 360.0
                if az < 0.0:
                    az += 360.0
                sats_inp[this_sat].ra.append(ra)
                sats_inp[this_sat].dec.append(dec)
                sats_inp[this_sat].provided.az.append(az)
                sats_inp[this_sat].provided.el.append(el)
    for this_sat in sats_inp:
        sats_inp[this_sat].ra = np.array(sats_inp[this_sat].ra)
        sats_inp[this_sat].dec = np.array(sats_inp[this_sat].dec)
        sats_inp[this_sat].provided.az = np.array(sats_inp[this_sat].provided.az)
        sats_inp[this_sat].provided.el = np.array(sats_inp[this_sat].provided.el)
        sats_inp[this_sat].times = Time(sats_inp[this_sat].times)
    return sats_inp


def readc(fn):
    """
    Read in the later November version that is already "filtered"

    """
    sats_inp = {}
    with open(fn, 'r') as fp:
        for line in fp:
            if '#ObsData' in line:
                continue
            data = line.split(',')
            this_sat = int(data[0])
            sats_inp[this_sat] = ObsData(name=this_sat, norad=int(data[1]))
            da = data[4].split()[0].split('/')
            day = f"20{da[2]}-{da[0]}-{da[1]}"
            sats_inp[this_sat].times = [Time(f"{day}T{data[8]}")]
            sats_inp[this_sat].ra = [float(data[9])]
            sats_inp[this_sat].dec = [float(data[10])]
            sats_inp[this_sat].provided.az = [float(data[11])]
            sats_inp[this_sat].provided.el = [float(data[12])]
            sats_inp[this_sat].norad = [int(data[1])]
            sats_inp[this_sat].bf_distance = [int(data[-1])]
    return sats_inp

class Input:
    def __init__(self, tools):
        """
        Parameter
        ---------
        tools : class (nominally obsid_base.Base)
        """
        self.tools = tools

    def read_SpaceX(self, fn, ftype='c'):
        """
        Since the files from SpaceX vary so much, pull the readers out but produce common self.sats dictionary.

        """
        from . import starlink_input
        self.SpaceX = fn
        self.sats = getattr(starlink_input, f"read{ftype}")(fn)
        # if ftype == 'a':
        #     self.sats = starlink_input.reada(fn)
        # if ftype == 'b':
        #     self.sats = starlink_input.readb(fn)
        # if ftype == 'c':
        #     self.sats = starlink_input.readc(fn)
        self.tools.get_azel()

    def filter(self, name=None, ytime=['2024-11-06T23:00:00', '2024-11-15T01:00:00'], yel=[30, 70]):
        """
        use None, None, None to write out everything.

        """
        import string
        tags = string.ascii_lowercase
        if name == 'dump':
            name, ytime, yel = None, None, None
        if name is None:
            name = list(self.sats.keys())
        if ytime is not None:
            ytime[0] = Time(ytime[0], format='isot')
            ytime[1] = Time(ytime[1], format='isot')
        fephjson = {'Sources': {}}
        ctr = 0
        for this_sat in name:
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
                obsid_name = f"S{this_sat}{tags[i]}"
                if use_it:
                    fephjson['Sources'][obsid_name] = {
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
        self.write_feph(fn=None, feph_dict=fephjson)