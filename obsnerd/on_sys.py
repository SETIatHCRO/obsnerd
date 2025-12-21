from odsutils.ods_tools import listify
from os import path as op
from glob import glob


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

def make_lo(lo):
    if lo is None:
        return ''
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


class InputMetadata:
    """
    This holds metadata for a Look instance supplied input.

    It does not track other attributes of the Look class.

    From the README.md file:
        source - a unique source name.
        obsid - a unique observation identifier: <SOURCE>_<MJD{.5f}>
        obsrec - a unique observation/hardware identifier: <OBSID>_<LO>_<CNODE>
        obsrec file - the filename holding the obsrec information (generally <OBSREC>.npz note)
        experiment - a session looking at sources (typically within an MJD day or two)
        obsinfo - a file containing information on obsids for a given experiment: obsinfo_<MJD>.json

    Parameters
    ----------
    oinput : str
        oinput to parse

    """
    parameters = ['input', 'source', 'mjd', 'obsid', 'obsrec', 'obsinfo', 'input_type', 'is_file', 'lo', 'cnode']
    def __init__(self, oinput):
        for p in self.parameters:
            setattr(self, p, None)
        self.input = oinput
        self.parse()

    def __str__(self):
        s = 'Metadata:\n'
        for p in self.parameters:
            s += f"  {p}: {getattr(self, p)}\n"
        return s

    def __repr__(self):
        s = '<Metadata: '
        for p in self.parameters:
            s += f"{p}={getattr(self, p)}, "
        s = s[:-2] + '>'
        return s

    def parse(self):
        """
        This parses the input to set obsid, source, mjd, obsrec_files.

        Parameters
        ----------
        oinput : str
            oinput to use to find obsid etc...

        Attributes
        ----------
        obsid : str
            ObsID found from input
        source : str
            Source name found from input
        mjd : float
            MJD found from input
        obsrec_files : list of str
            List of obsrec files found from input
        """
        if self.input is None:
            self.input_type = None
            return
        try:
            self.mjd = float(self.input)
            self.input_type = 'mjd'
            return
        except (ValueError, TypeError):
            pass
        if self.input.endswith('.uvh5'):
            self.is_file = True
            self.input_type = 'uvh5'
            for k, v in  parse_uvh5_filename(self.input).items():
                setattr(self, k, v)
        elif self.input.endswith('.npz'):
            self.is_file = True
            self.input_type = 'npz'
            for k, v in  split_obsrec(self.input).items():
                setattr(self, k, v)
        elif self.input.startswith('obsinfo') or self.input.endswith('.json'):
            self.is_file = True
            self.input_type = 'obsinfo'
            if not self.input.endswith('.json'):
                self.obsinfo = f"{self.input}.json"
            else:
                self.obsinfo = self.input
            self.mjd =  self.obsinfo[:-5].split('_')[-1]
        else:  # obsid or source or obsrec
            data = self.input.split('_')
            if len(data) == 4:
                self.input_type = 'obsrec'
                for k, v in  split_obsrec(self.input).items():
                    setattr(self, k, v)
            elif len(data) == 2:
                self.input_type = 'obsid'
                self.source, self.mjd =  split_obsid(self.input)
                self.obsid = self.input
            else:
                self.input_type = 'source'
                self.source = self.input
        try:
            self.mjd = float(self.mjd)
        except (ValueError, TypeError):
            pass


def get_obsinfo_filename_from_oinput(oinput):
    """
    Get the obsinfo filename from an input.
    This function tries to determine the obsinfo filename based on the input provided.

    Parameters
    ----------
    oinput : str
        Either an obsid, mjd or an obsinfo filename.

    """
    meta = InputMetadata(oinput)
    if meta.input_type == 'obsinfo':
        if op.exists(meta.obsinfo):
            return meta.obsinfo
        return None
    if meta.mjd is None:
        return None
    for p in range(len(oinput)-5, 0, -1):
        mjd_option = make_obsinfo_filename(meta.mjd, decimal_places=p)
        if op.exists(mjd_option):
            return mjd_option
    return None


def get_obsid_from_source(source, data_dir='.'):
    """
    Get the obsid from a source.
 
    Parameters
    ----------
    source : str
        Source name.
    data_dir : str
        Directory where the data is stored.

    """
    for npzfnfp in glob(f'{data_dir}/*.npz'):
        npzfn = op.basename(npzfnfp)
        if source in npzfn:
            return  split_obsrec(npzfn)['obsid'] 
    return None


def make_obsid(source, mjd, decimal_places=5):
    smjd = make_mjd_for_filename(mjd, decimal_places=decimal_places)
    return f"{source}_{smjd}"


def make_obsinfo_filename(mjd, decimal_places=1):
    smjd = make_mjd_for_filename(mjd, decimal_places=decimal_places)
    return f"obsinfo_{smjd}.json"


def make_mjd_for_filename(mjd, decimal_places=5):
    from numpy import floor
    mult = 10 ** decimal_places
    mjd = int(floor(mjd * mult)) / mult
    return f"{float(mjd):.{decimal_places}f}"


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