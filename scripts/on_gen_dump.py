#! /usr/bin/env python
import argparse
from obsnerds import obs_look


ap = argparse.ArgumentParser()
ap.add_argument('date', help='Date-name of directory or ? to list options')
ap.add_argument('--base_dir', help='Base directory', default='/mnt/primary/ata/projects/p054/')
ap.add_argument('--script', help="Name of script file", default='dump_autos.sh')
ap.add_argument('--ants', help="Ants to use - csv list of 'all'", default='2a,2b,2h,4e')
ap.add_argument('--pols', help="Pols to use", default='xx,xy')
ap.add_argument('--cnodes', help="CNODES to use.", default='all')
ap.add_argument('--los', help="LOs to use", default='all')
args = ap.parse_args()

#on_dump_autos.py /mnt/primary/ata/projects/p054/2024-12-12-22:22:01/uvh5_60656_83668_34122680_S11105_1212_0001/LoB.C0928/uvh5_60656_83668_34122680_S11105_1212_0001.uvh5 --lo B --cnode C0928 --ants 2a,2b,2h,4e
obs_look.Look(None, cnoee)
obs_look.gen_dump_script(date_path=args.date, base_path=args.base_dir, script_filename=args.script, ants=args.ants, pols=args.pols)