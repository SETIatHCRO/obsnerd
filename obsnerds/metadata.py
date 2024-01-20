from datetime import datetime, timezone, timedelta
import yaml
from . import onutil


ONLOG_FILENAME = 'onlog.log'
META_FILENAME = 'metadata.yaml'
UTC = timezone(timedelta(0), 'UTC')
PST = timezone(timedelta(hours=-8), 'PST')
PDT = timezone(timedelta(hours=-7), 'PDT')


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
    ts = datetime.now().astimezone(UTC).isoformat()
    with open(ONLOG_FILENAME, 'a') as fp:
        for note in notes:
            print(f"{ts} -- {note}", file=fp)


def read_onlog(params=['tstart', 'source:', 'expected:', 'azel:', 'fcen:', 'bw:', 'tstop']):
    logdata = {}
    for param in params:
        logdata[param.strip(':')] = {}
    with open(ONLOG_FILENAME, 'r') as fp:
        for line in fp:
            for param in params:
                if param in line:
                    data = [x.strip() for x in line.split('--')]
                    if param in ['tstart', 'tstop']:
                        logdata[param][data[0]] = data[0]
                    else:
                        logdata[param.strip(':')][data[0]] = ' -- '.join(data[1:])
    return logdata

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
    logdata = read_onlog()
    ts = sorted(logdata[param])
    if not len(ts):
        return None
    val = logdata[param][ts[-1]]
    if parse:
        val = val.split(parse)[-1].strip()
    return val


def get_summary():
    from copy import copy
    logdata = read_onlog()
    tstart = list(logdata['tstart'].keys())
    source = list(logdata['source'].keys())
    azel = list(logdata['azel'].keys())
    expected = list(logdata['expected'].keys())
    fcen = list(logdata['fcen'].keys())
    bw = list(logdata['bw'].keys())
    tstop = list(logdata['tstop'].keys())
    tsall = set(tstart + source + azel + expected + fcen + bw + tstop)
    tsall = sorted(list(tsall))
    table_data = []
    row = ['' for x in range(8)]
    previous_fcen = ''
    for this_ts in tsall:
        if this_ts in tstart:
            row[0] = logdata['tstart'][this_ts]
        if this_ts in source:
            row[2] = logdata['source'][this_ts].split(':')[-1]
        if this_ts in expected:
            row[3] = logdata['expected'][this_ts][9:].strip()
        if this_ts in azel:
            payload = logdata['azel'][this_ts].split(':')[-1].split(',')
            row[4] = payload[0]
            row[5] = payload[1]
        if this_ts in fcen:
            row[6] = logdata['fcen'][this_ts].split(':')[-1]
            previous_fcen = copy(row[6])
        if this_ts in bw:
            row[7] = float(logdata['bw'][this_ts].split(':')[-1]) / 1E6
        if this_ts in tstop:
            row[1] = logdata['tstop'][this_ts]
            if not len(row[6]):
                row[6] = previous_fcen
            table_data.append(row)
            row = ['' for x in range(8)]
    from tabulate import tabulate
    print(tabulate(table_data, headers=['start', 'stop', 'source', 'expected', 'az', 'el', 'fcen', 'bw']))
    return table_data
    

# Metadata functions
def get_meta():
    with open(META_FILENAME, 'r') as fp:
        meta = yaml.safe_load(fp)
    for key in ['tstart', 'tstop', 'tle', 'expected']:
        if key in meta:
            meta.update({key: onutil.make_datetime(**{key: meta[key]})})
    return meta


def start(samp_rate, decimation, nfft):
    move = get_latest_value('move to', parse=':')
    data = {
        'tstart': datetime.now().astimezone(UTC).isoformat(),
        'fcen': float(get_latest_value('fcen', parse=':')),
        'bw': samp_rate / 1E6,
        'decimation': decimation,
        'nfft': nfft,
        'tle': onutil.make_datetime(date=get_latest_value('TLEs', parse='timestamp'), tz=0.0),
        'source': get_latest_value('source', parse=':'),
        'expected': onutil.make_datetime(date=get_latest_value('expected', parse=' '), tz=0.0),
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
    add_value(**{kw: datetime.now().astimezone(UTC).isoformat()})
