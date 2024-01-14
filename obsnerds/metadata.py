from datetime import datetime
import yaml
from . import onutil


ONLOG_FILENAME = 'onlog.log'
META_FILENAME = 'metadata.yaml'


# Log functions
def onlog(notes):
    """
    Add notes to the log.

    Parameter
    ---------
    notes : str or list
        entries to add
    """
    if isinstance(notes, str):
        notes = [notes]
    ts = datetime.now()
    with open(ONLOG_FILENAME, 'a') as fp:
        for note in notes:
            print(f"{ts} -- {note}", file=fp)


def get_latest_value(param, parse=False):
    """
    Return the latest entry for given param in line.

    Parameters
    ----------
    param : str
        string to search for
    parse : str or False
        if 'timestamp' uses the timestamp
        if str will split on that string and return last index
    """
    metadata = {}
    indentry = 0 if parse == 'timestamp' else -1
    with open(ONLOG_FILENAME, 'r') as fp:
        for line in fp:
            if param in line:
                data = [x.strip() for x in line.split('--')]
                metadata[data[0]] = data[indentry]
    ts = sorted(metadata)
    if not len(ts):
        return None
    val = metadata[ts[-1]]
    if parse:
        val = val.split(parse)[-1].strip()
    return val


# Metadata functions
def get_meta():
    with open(META_FILENAME, 'r') as fp:
        meta = yaml.safe_load(fp)
    for key, val in meta.items():
        tval = onutil.make_datetime(date=val)
        if isinstance(tval, datetime):
            meta.update({key: tval})
    return meta


def start(samp_rate, decimation, nfft):
    move = get_latest_value('move to', parse=':')
    data = {
        'tstart': datetime.now().isoformat(),
        'fcen': float(get_latest_value('fcen', parse=':')),
        'bw': samp_rate / 1E6,
        'decimation': decimation,
        'nfft': nfft,
        'tle': get_latest_value('TLEs', parse='timestamp'),
        'source': get_latest_value('source', parse=':'),
        'expected': get_latest_value('expected', parse=' '),
        'move': move,
        'move_data': get_latest_value(move, parse=':')
    }
    add_value(initialize=True, **data)
    onlog(['tstart', f"bw: {samp_rate}"])


def stop():
    add_datetimestamp('tstop')
    onlog('tstop')


def add_value(initialize=False, **kwargs):
    if initialize:
        meta = kwargs
    else:
        meta = get_meta()
        meta.update(kwargs)
    with open(META_FILENAME, 'w') as fp:
        yaml.safe_dump(meta, fp)


def add_datetimestamp(kw):
    add_value(**{kw: datetime.now().isoformat()})
