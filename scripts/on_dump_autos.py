#! /usr/bin/env python
from obsnerd import onv_dump


if __name__ == '__main__':
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument('filename', help="Name of UVH5 file to dump.")
    ap.add_argument('--ants', help="Ant list to drop or'all'", default='all')
    ap.add_argument('--pols', help="Pols to use", default='all')

    args = ap.parse_args()
    look = onv_dump.Dump(args.filename)
    look.dump_autos(args.ants, args.pols)
