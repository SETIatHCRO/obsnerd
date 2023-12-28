from datetime import datetime
"""
This writes out datetime as UTC
"""

ONLOG_FILENAME = 'onlog.log'
META_FILENAME = 'metadata.yaml'


def onlog(notes):
    if isinstance(notes, str):
        notes = [notes]
    ts = datetime.now()
    with open(ONLOG_FILENAME, 'a') as fp:
        for note in notes:
            print(f"{ts} -- {note}", file=fp)


def get_latest(param):
    metadata = {}
    with open(ONLOG_FILENAME, 'r') as fp:
        for line in fp:
            if param in line:
                data = [x.strip() for x in line.split('--')]
                metadata[data[0]] = data[1]
    ts = sorted(metadata)
    return metadata[ts[-1]]


def start(samp_rate, decimation, nfft):
    with open(META_FILENAME, 'w') as fp:
        print(f"tstart: {datetime.now()}", file=fp)
        print(get_latest('fcen'), file=fp)
        print(f"bw: {samp_rate / 1E6}", file=fp)
        print(f"decimation: {decimation}", file=fp)
        print(f"nfft: {nfft}", file=fp)
    onlog(['tstart', f"bw: {samp_rate}"])


def stop():
    with open(META_FILENAME, 'a') as fp:
        print(f"tstop: {datetime.now()}", file=fp)
    onlog('tstop')