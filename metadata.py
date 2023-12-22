from datetime import datetime

def start(samp_rate, fn='metadata.yaml'):
    with open(fn, 'a') as fp:
        print(f"tstart: {datetime.now()}", file=fp)
        print(f"bw: {samp_rate / 1E6}", file=fp)

def stop():
    with open(fn, 'a') as fp:
        print(f"tstop: {datetime.now()}", file=fp)
