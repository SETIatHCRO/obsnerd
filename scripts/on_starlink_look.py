import matplotlib.pyplot as plt
import argparse
ap = argparse.ArgumentParser()
from obsnerds import starlink_look

if __name__ == '__main__':
    ap.add_argument('filename', help="Name of file to dump.")
    ap.add_argument('-a', '--ants', help='csv list of ants', default=None)
    ap.add_argument('-p', '--pols', help="List of polarizations to use for some options...", default='xx,yy,yx,xy')
    ap.add_argument('-w', '--waterfall', help="Flag to generate all of the waterfalls.", action='store_true')
    ap.add_argument('-s', '--save', help="Save the generated plots", action='store_true')
    ap.add_argument('-t', '--time_axis', help="Type of 'time' axis for dashboard", default='boresight')
    ap.add_argument('-f', '--feph_file', help="Name of feph file", default=None)
    ap.add_argument('--dash', help="Generate the dashboard", action='store_true')
    ap.add_argument('--dump_autos', help="Flag to dump the autos", action='store_true')
    args = ap.parse_args()
    if args.ants is not None:
        args.ants = args.ants.split(',')
    args.pols = args.pols.split(',')
    file_type = args.filename.split('.')[-1]
    sl = starlink_look.Look(eph=args.feph_file)
    if file_type == 'uvh5':
        sl.read_uvh5(args.filename)
    elif file_type == 'npz':
        sl.read_npz(args.filename)
    if args.dump_autos:
        sl.dump_autos(ants=args.ants, pols=args.pols)
    if args.waterfall:
        for pol in args.pols:
            sl.all_wf(pol=pol, save=args.save)
        if not args.save:
            plt.show()
    if args.dash:
        sl.dashboard(ant=args.ants[0], pol=args.pols[0], save=args.save, time_axis=args.time_axis)
        if not args.save:
            plt.show()