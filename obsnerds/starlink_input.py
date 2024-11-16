from astropy.time import Time, TimeDelta
from obsnerds.starlink_eph import Satellite
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
                sats_inp[this_sat] = Satellite(satno=line.strip())
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
                sats_inp[this_sat] = Satellite(satno=line.strip())
            else:
                t, tstr, rastr, ra, decstr, dec, azstr, az, elstr, el = line.split()
                sats_inp[this_sat].times.append(self.T0 + TimeDelta(float(t), format='sec'))
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
            if '#Satellite' in line:
                continue
            data = line.split(',')
            this_sat = int(data[0])
            sats_inp[this_sat] = Satellite(satno=this_sat, norad=int(data[1]))
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