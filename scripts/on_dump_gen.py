#! /usr/bin/env python
from os.path import join
from os import walk, listdir
import argparse


def generate(base_path, query=False, script_filename='dump_autos.sh'):
    if query:
        print(f"Available observation dates in {base_path}:")
        for x in listdir(base_path):
            print(f"\t{x}")
        return
    print(f"Retrieving from {base_path}")
    files = {}
    for basedir, _, filelist in walk(base_path):
        if base_path in basedir and '/Lo' in basedir:
            lolo, cnode = basedir.split('/')[-1].split('.')
            lo = lolo[2:]
            for fn in filelist:
                if fn.startswith('uvh5_'):
                    _, mjd1, mjd2, _, src, _ = fn.split('_')
                    mjd = float(f"{mjd1}.{mjd2}")
                    obsid = f"{src}_{mjd:.4f}"
                    obsrec = f"{obsid}_{lo}_{cnode}"
                    dfn = join(basedir, fn)
                    files[obsrec] = [dfn, lo, cnode]
                
    with open(script_filename, 'w') as fp:
        for obsrec, data in files.items():
            print(f"python on_dump_autos.py {data[0]} --lo {data[1]} --cnode {data[2]}", file=fp)
            print(f"Adding {obsrec}")

if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('date', help='Date-name of directory or ? to list options')
    ap.add_argument('--base_dir', help='Base directory', default='/mnt/primary/ata/projects/p054/')
    ap.add_argument('--script', help="Name of script file", default='dump_autos.sh')
    args = ap.parse_args()
else:
    args = argparse.Namespace(date='2024-11-19-03:22:01', base_dir='/mnt/primary/ata/projects/p054/', script='dump_autos.sh')

if args.date == '?':
    base_path = args.base_dir
    query = True
else:
    base_path = join(args.base_dir, args.date)
    query = False
generate(base_path=base_path, query=query, script_filename=args.script)