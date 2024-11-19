#! /usr/bin/env python
from os.path import join
from os import walk

files = {}
lsmethod = 2
base_dir = '/mnt/primary/ata/projects/p054/'
base_date = '2024-11-19-03:22:01'
base_path = join(base_dir, base_date)
print(f"Retrieving from {base_path}")

obs_ctr = {}
for dinfo in walk(base_path):
    if base_path in dinfo[0] and '/Lo' in dinfo[0]:
        data = dinfo[0].split('/')
        lolo, cnode = data[-1].split('.')
        lo = lolo[2:]
        src = data[-2].split('_')[-2]
        obs_ctr.setdefault(src, 1)
        obs_ctr[src] += 1
        obsid = f"{src}-{obs_ctr[src]}_{lo}_{cnode}"
        for fn in dinfo[2]:
            if fn.startswith('uvh5_'):
                dfn = join(dinfo[0], fn)
                files[obsid] = dfn
with open('dump_auto.sh', 'w') as fp:
    for src, this_file in files.items():
        print(f"python on_dump_autos.py {this_file} -o {src}.npz", file=fp)
