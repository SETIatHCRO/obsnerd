#! /usr/bin/env python
import argparse
from obsnerd import obs_sources


ap = argparse.ArgumentParser()
ap.add_argument('sources', help='Target list or json file with inputs')
ap.add_argument('integrations', help="Integration times in seconds for each source", nargs='?', default=None)
ap.add_argument('starts', help="Start times of observations (seconds or datetime.isoformat)", nargs='?', default=0)
ap.add_argument('-a', '--freq_loA', help='Freq LO A in MHz', default=1500, type=float)
ap.add_argument('-b', '--freq_loB', help='Freq LO B in MHz', default=6000, type=float)
ap.add_argument('-t', '--integration_times', help="Integration times for sources [sec]", default='300')
ap.add_argument('--ant_list', help="List of antennas or group list", default='rfsoc_active')
ap.add_argument('--known_bad', help='List of known bad antennas', default='')
ap.add_argument('--fiddle', help="Extra fiddle time in seconds", default=5)
ap.add_argument('--focus', help="Focus parameter", choices=['a', 'b', 'max'], default='max')
ap.add_argument('--data_record', help="Data recording to use", choices=['gnuradio', 'hpguppi'], default='hpguppi')
args = ap.parse_args()

if args.sources.endswith('.json'):
    import json
    with open(args.sources, 'r') as fp:
        vals = json.load(fp)
    args.sources = None
    argv = vars(args)
    argv.update(vals)
    args = argparse.Namespace(**argv)
if args.sources is None or args.integrations is None:
    raise ValueError("Must specify sources and integration times.")

freqs = {'a': args.freq_loA,
         'b': args.freq_loB}

observer = obs_sources.Observer(args.sources, args.integrations, start_times=args.starts, freqs=freqs,
                                ant_list=args.ant_list, known_bad=args.known_bad, focus_on=args.focus,
                                obs_fiddle=args.fiddle, data_record=args.data_record)
observer.setup_session()
observer.step_obs()