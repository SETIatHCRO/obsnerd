#! /usr/bin/env python
from os.path import join
from os import walk
import argparse


def generate(base_path):
    print(f"Retrieving from {base_path}")
    files = {}
    for dinfo in walk(base_path):
        if base_path in dinfo[0] and '/Lo' in dinfo[0]:
            data = dinfo[0].split('/')
            lolo, cnode = data[-1].split('.')
            lo = lolo[2:]
            src = data[-2].split('_')[-2]
            obsid = f"{src}_{lo}_{cnode}"
            for fn in dinfo[2]:
                if fn.startswith('uvh5_'):
                    dfn = join(dinfo[0], fn)
                    files[obsid] = [dfn, lo, cnode]
                
    with open('dump_auto.sh', 'w') as fp:
        for src, data in files.items():
            print(f"python on_dump_autos.py {data[0]} --lo {data[1]} --cnode {data[2]}", file=fp)

if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('date', help='Date-name of directory')
    ap.add_argument('--base_dir', help='Base directory', default='/mnt/primary/ata/projects/p054/')
    args = ap.parse_args()
else:
    args = argparse.Namespace(date='2024-11-19-03:22:01', base_dir='/mnt/primary/ata/projects/p054/')

base_path = join(args.base_dir, args.date)
generate(base_path=base_path)