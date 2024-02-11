from datetime import datetime, timezone, timedelta
import yaml
from . import onutil


ONLOG_FILENAME = 'onlog.log'
META_FILENAME = 'metadata.yaml'
UTC = timezone(timedelta(0), 'UTC')
PST = timezone(timedelta(hours=-8), 'PST')
PDT = timezone(timedelta(hours=-7), 'PDT')
LOGFILE_DELIMITER = '--'
LOG_ENTRIES_TO_GET = ['tstart', 'source:', 'expected:', 'azel:', 'fcen:', 'bw:', 'session start:',
                      'end:', 'traj:', 'track:', 'Writing', 'move to:', 'TLEs', 'tstop']


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

class Onlog:
    def __init__(self, entries=LOG_ENTRIES_TO_GET, delimiter=LOGFILE_DELIMITER, auto_read=False):
        self.file = ONLOG_FILENAME
        self.delimiter = delimiter
        self.entries = entries
        self.data = {'other': {}}
        self.latest = {}
        self.keys = [p.strip(':') for p in self.entries]
        for key in self.keys:
            self.data[key] = {}
            self.latest[key] = ''
        if auto_read:
            self.read()

    def read(self):
        self.all_timestamps = set()
        with open(self.file, 'r') as fp:
            for line in fp:
                par_not_found = True
                linedata = [x.strip() for x in line.split(self.delimiter)]
                timestamp, payload = linedata[0], self.delimiter.join(linedata[1:])
                self.all_timestamps.add(timestamp)
                for key, entry in zip(self.keys, self.entries):
                    if entry in linedata[1]:
                        par_not_found = False
                        if key in ['tstart', 'tstop', 'TLEs']:
                            self.data[key][timestamp] = timestamp
                            self.latest[key] = timestamp
                        elif key == 'track':
                            self.data[key].setdefault(timestamp, [])
                            self.data[key][timestamp].append(payload)
                            self.latest[key] = payload
                        else:
                            self.data[key][timestamp] = payload
                            self.latest[key] = payload
                if par_not_found:
                    self.data['other'][timestamp] = payload
        self.all_timestamps = sorted(list(self.all_timestamps))

    def get_latest_value(self, key, parse=False):
        val = self.latest[key]
        if parse:
            val = val.split(parse)[-1].strip()
        return val


def get_latest_value(key, parse=False):
    """
    Return the latest entry for given entry in line.

    Parameters
    ----------
    key : str
        key to use
    parse : str or False
        if 'timestamp' uses the timestamp
        if str will split on that string and return last index
    """
    log = Onlog(auto_read=True)
    return log.get_latest_value(key=key, parse=parse)


def get_summary():
    from copy import copy
    log = Onlog(auto_read=True)
    headers=['obs', 'start', 'stop', 'source', 'expected', 'az', 'el', 'fcen', 'bw']
    I = {}
    for i, hdr in enumerate(headers):
        I[hdr] = i
    table_data = []
    row = ['' for x in range(len(headers))]
    previous = {'fcen': '', 'obs': '', 'az': '', 'el': ''}
    for this_ts in log.all_timestamps:
        if this_ts in log.data['session start']:
            row[I['obs']] = (log.data['session start'][this_ts].split(log.delimiter)[0][15:]).strip()
            previous['obs'] = copy(row[I['obs']])
        if this_ts in log.data['tstart']:
            row[I['start']] = log.data['tstart'][this_ts]
        if this_ts in log.data['source']:
            row[I['source']] = log.data['source'][this_ts].split(':')[-1]
        if this_ts in log.data['expected']:
            row[I['expected']] = log.data['expected'][this_ts][9:].strip()
        if this_ts in log.data['azel']:
            payload = log.data['azel'][this_ts].split(':')[-1].split(',')
            row[I['az']] = payload[0]
            row[I['el']] = payload[1]
            previous['az'] = copy(row[I['az']])
            previous['el'] = copy(row[I['el']])
        if this_ts in log.data['fcen']:
            row[I['fcen']] = log.data['fcen'][this_ts].split(':')[-1]
            previous['fcen'] = copy(row[I['fcen']])
        if this_ts in log.data['bw']:
            row[I['bw']] = float(log.data['bw'][this_ts].split(':')[-1]) / 1E6
        if this_ts in log.data['tstop']:
            row[I['stop']] = log.data['tstop'][this_ts]
            for par in ['obs', 'fcen', 'az', 'el']:
                if not len(row[I[par]]):
                    row[I[par]] = previous[par]
            table_data.append(row)
            row = ['' for x in range(len(headers))]
    from tabulate import tabulate
    print(tabulate(table_data, headers=headers))
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
    log = Onlog(auto_read=True)
    move = log.get_latest_value('move to', parse=':')
    print(f"Move type is {move}")
    data = {
        'tstart': datetime.now().astimezone(UTC).isoformat(),
        'fcen': float(log.get_latest_value('fcen', parse=':')),
        'bw': samp_rate / 1E6,
        'decimation': decimation,
        'nfft': nfft,
        'tle': onutil.make_datetime(date=log.get_latest_value('TLEs', parse='timestamp'), tz=0.0),
        'source': log.get_latest_value('source', parse=':'),
        'expected': onutil.make_datetime(date=log.get_latest_value('expected', parse=' '), tz=0.0),
        'move': move,
        'move_data': log.get_latest_value(move, parse=':')
    }
    if move == 'traj':
        import yaml
        from .trajectory_engine import TRACK_LOG_FILENAME
        with open(TRACK_LOG_FILENAME, 'r') as fp:
            move_data = yaml.safe_load(fp)
        data['track'] = yaml.safe_dump(move_data)
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
