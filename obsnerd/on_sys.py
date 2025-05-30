from odsutils.ods_tools import listify

ALL_CNODES = ['C0352', 'C0544', 'C0736', 'C0928', 'C1120', 'C1312', 'C1504']
ALL_LOS = ['A', 'B']
AXIS_OPTIONS = {'b': 'boresight', 'd': 'datetime', 's': 'seconds'}
FILTER_AXIS = {'time': 0, 'freq': 1}
OBS_START_DELAY = 10  # sec
FEED_FOCUS_SLEEP = 20
STARTUP_LATENCY = 15
LATENCY = 30
SLEW_SPEED = 1.5  # deg/sec

ANT_LISTS = {'old_feeds': ['']}


def parse_ods_notes(entry):
    """
    Parses notes and returns a dictionary with keys 'source', 'cnode', 'lo', 'obsid', 'mjd'.
    """
    if 'notes' not in entry:
        return {}
    data = {'ods': 'True'}  # default value to use ODS entry
    for item in entry['notes'].split(';'):
        if ':' in item:
            key, val = item.split(':')
            data[key.strip()] = val.strip()
    return data

def make_lo(lo):
    if lo.upper() not in ALL_LOS:
        raise ValueError(f"LO must be one of {ALL_LOS}")
    return lo.upper()

def make_cnode(cns):
    """
    Takes a list and makes a list of valid cnodes -> C####

    """
    cns = listify(cns, {'all': ALL_CNODES})
    try:
        return [f"C{int(x):04d}" for x in cns]
    except ValueError:
        return cns


def make_obsid(source, mjd):
    return f"{source}_{float(mjd):.5f}"


def split_obsid(obsid):
    if obsid is None:
        return None, None
    mjd = obsid.split('_')[-1]
    source = obsid[:(len(obsid) - len(mjd) - 1)]
    try:
        mjd = float(mjd)
        if mjd < 60000:
            mjd = None
    except ValueError:
        mjd = None
    return source, mjd


def make_obsrec(source, mjd, lo, cnode):
    cnode = make_cnode(cnode)
    obsid = make_obsid(source, mjd)
    return f"{obsid}_{lo.upper()}_{cnode[0]}"


def split_obsrec(obsrec):
    if obsrec.endswith('.npz'):
        obsrec = obsrec[:-4]
    data = obsrec.split('_')
    obsid = '_'.join(data[:-2])
    lo = data[-2]
    cnode = data[-1]
    source, mjd = split_obsid(obsid)
    return {'obsrec': obsrec, 'source': source, 'mjd': mjd, 'obsid': obsid, 'lo': lo, 'cnode': cnode}


def parse_uvh5_filename(fullfn):
    pieces = split_obsrec('-_-_-_-')
    pieces.update({'filename': fullfn})
    data = fullfn.split('/')
    fn = data[-1]
    if fn.startswith('uvh5_'):
        sdat = fn.split('_')
        mjd = float(f"{sdat[1]}.{sdat[2]}")
        src = '_'.join([x for x in sdat[4:-1]])
        pieces.update({'mjd': mjd, 'source': src, 'obsid': make_obsid(src, mjd)})
        if len(data) > 1 and data[-2].startswith('Lo'):
            lo = data[-2][2]
            cnode = data[-2].split('.')[-1]
            pieces.update({'lo': lo, 'cnode': cnode, 'obsrec': make_obsrec(src, mjd, lo, cnode)})
    return pieces