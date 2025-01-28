from astropy.time import Time, TimeDelta
from obsnerd.onv_base import ObservationData as ObsData
import numpy as np
from copy import copy

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

def readd(fn):
    """
    Read in the later December file

    """
    sats_inp = {}
    with open(fn, 'r') as fp:
        header = fp.readline().split(',')
        header = [x.strip() for x in header]
        print(', '.join(header))
        for line in fp:
            data = line.split(',')
            this_sat = int(data[1])
            sats_inp.setdefault(this_sat, ObsData(name=this_sat, center={}))
            sats_inp[this_sat].times.append(Time(data[3].strip()))
            t_ra = float(data[4])
            this_ra = t_ra if t_ra > 0 else 360.0 + t_ra
            sats_inp[this_sat].ra.append(this_ra)
            sats_inp[this_sat].dec.append(float(data[5]))
            sats_inp[this_sat].az.append(float(data[6]))
            sats_inp[this_sat].el.append(float(data[7]))
            if int(data[8]):
                sats_inp[this_sat].center['time'] = Time(data[3].strip())
                sats_inp[this_sat].center['ra'] = float(this_ra)
                sats_inp[this_sat].center['dec'] = float(data[5])
                sats_inp[this_sat].center['az'] = float(data[6])
                sats_inp[this_sat].center['el'] = float(data[7])
    return sats_inp

def write(sats, src_tag, fn='sources.json'):
    import json
    source_list = {}
    fpsl = open("source_list.txt", 'w')
    print(f"src_id src_ra_j2000_deg src_dec_j2000_deg src_start_utc src_end_utc", file=fpsl)
    for this_sat in sats:
        src_name = f"S{this_sat}_{src_tag}"
        if 'ra' in sats[this_sat].center:
            print(f"Found {this_sat} at {sats[this_sat].center['time']}")
            source_list[src_name] = {
                "ra": sats[this_sat].center['ra'],
                "dec": sats[this_sat].center['dec'],
                "time": sats[this_sat].center['time'].datetime.isoformat(timespec='seconds'),
                "tuning": {'a': 1980.0, 'b': 5500.0}
                }
            start_t = (sats[this_sat].center['time'] - TimeDelta(240, format='sec')).datetime.isoformat(timespec='seconds')
            end_t = (sats[this_sat].center['time'] + TimeDelta(240, format='sec')).datetime.isoformat(timespec='seconds')
            print(f"{src_name}  {sats[this_sat].center['ra']:.5f}  {sats[this_sat].center['dec']} {start_t} {end_t}", file=fpsl)
    fpsl.close()
    with open(fn, 'w') as fp:
        json.dump(source_list, fp, indent=4)


class Input:
    def __init__(self, tools=False):
        """
        Parameter
        ---------
        tools : class (nominally obs_base.Base)
        """
        self.tools = tools

    def read_SpaceX(self, fn, ftype='d', getazel=False):
        """
        Since the files from SpaceX vary so much, pull the readers out but produce (somewhat) common self.sats dictionary.

        """
        self.SpaceX = fn
        from . import starlink_io
        self.sats = getattr(starlink_io, f"read{ftype}")(fn)
        self.sort_in_place('times', ['ra', 'dec', 'az', 'el'])
        if self.tools:
            self.tools.get_azel()

    def sort_in_place(self, sortkey, otherkeys):
        allkeys = [sortkey] + otherkeys
        for this_sat in self.sats:
            sorter = {}
            for i, tt in enumerate(getattr(self.sats[this_sat], sortkey)):
                sorter[tt] = i
            this_map = [sorter[tt] for tt in sorted(sorter.keys())]
            for thiskey in allkeys:
                sorted_list = []
                for sorted_order in this_map:
                    next_entry = copy(getattr(self.sats[this_sat], thiskey)[sorted_order])
                    sorted_list.append(next_entry)
                setattr(self.sats[this_sat], thiskey, sorted_list)

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
        obsinfojson = {'Sources': {}}
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
                    obsinfojson['Sources'][obsid_name] = {
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
        self.write_obsinfo(fn=None, obsinfo_dict=obsinfojson)