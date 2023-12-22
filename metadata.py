from datetime import datetime


def start(samp_rate, fn='metadata.yaml'):
    with open(fn, 'a') as fp:
        print(f"tstart: {datetime.now()}", file=fp)
        print(f"bw: {samp_rate / 1E6}", file=fp)
    with open('onlog.log', 'a') as fp:
        print(f"{datetime.now()} -- tstart", file=fp)
        print(f"{datetime.now()} -- bw: {samp_rate}", file=fp)


def stop(fn='metadata.yaml'):
    with open(fn, 'a') as fp:
        print(f"tstop: {datetime.now()}", file=fp)
    with open('onlog.log', 'a') as fp:
        print(f"{datetime.now()} -- tstop", file=fp)