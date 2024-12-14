#! /usr/bin/env python
from obsnerds import obs_dump


if __name__ == '__main__':
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument('filename', help="Name of UVH5 file to dump.")
    ap.add_argument('--ants', help="Ant list to drop or'all'", default='all')
    ap.add_argument('--pols', help="Pols to use", default='all')

    args = ap.parse_args()
    look = obs_dump.Dump(args.filename)
    look.dump_autos(args.ants, args.pols)
