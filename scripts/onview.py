#! /usr/bin/env python
import argparse
from obsnerds import onview_engine as oe

ap = argparse.ArgumentParser()
ap.add_argument('fn', help="Name of hdf5 datafile")
ap.add_argument('-o', '--output_type', help='wf, series, spectra [wf]', choices=['wf', 'series', 'spectra', 'showdata'], default='wf')
ap.add_argument('-x', '--xticks', help="Number of xticks in waterfall [10]", type=int, default=10)
ap.add_argument('-y', '--yticks', help="Number of yticks to use in waterfall [4]", type=int, default=4)
ap.add_argument('-c', '--colorbar', help="Flag to hide colorbar", action='store_false')
ap.add_argument('-b', '--beamfit', help='Flag to fit for the beam', action='store_true')
ap.add_argument('-f', '--freq', help='Frequency [range in csv] to use.', default=None)
ap.add_argument('-t', '--time', help='Time [range in csv] to use.', default=None)
ap.add_argument('-l', '--log', help="Flag to take log10 of data", action='store_true')
ap.add_argument('-d', '--dB', help="Flag to convert to dB", action='store_true')
ap.add_argument('-P', '--total_power', help="Show total power (for series)", action='store_true')
ap.add_argument('-n', '--norm', help="Norm to apply for total power (Total, [/bin], /full)", choices=['Total', '/bin', '/full'], default='/bin')
ap.add_argument('--tz', help="Timezone offset to UTC in hours [-8.0]", type=float, default=0.0)

args = ap.parse_args()
obs = oe.Data(args.fn, args.tz)
getattr(obs, args.output_type)(**vars(args))
oe.plt.show()
