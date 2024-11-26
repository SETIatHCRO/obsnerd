#! /usr/bin/env python
import argparse
from obsnerds import obs_look


ap = argparse.ArgumentParser()
ap.add_argument('date', help='Date-name of directory or ? to list options')
ap.add_argument('--base_dir', help='Base directory', default='/mnt/primary/ata/projects/p054/')
ap.add_argument('--script', help="Name of script file", default='dump_autos.sh')
args = ap.parse_args()


obs_look.gen_dump_script(date_path=args.date, base_path=args.base_dir, script_filename=args.script)