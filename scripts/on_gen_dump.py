#! /usr/bin/env python
import argparse
from obsnerd import onv_dump


ap = argparse.ArgumentParser()
ap.add_argument('date', help='Date-name of directory or ? to list options')
ap.add_argument('--base-dir', dest='base_dir', help='Base directory', default='/mnt/primary/ata/projects/p054/')
ap.add_argument('--script', help="Name of script file", default='dump_autos.sh')
ap.add_argument('--ants', help="Ants to use - csv list or 'all'", default='2a,1e,2d,2e,4e')
ap.add_argument('--pols', help="Pols to use", default='xx,xy,yy')
ap.add_argument('--cnodes', help="CNODES to use.", default='all')
ap.add_argument('--los', help="LOs to use", default='all')
ap.add_argument('--download-script', dest='download_script_filename',
                help="Name of download script file", default='download_files.sh')
ap.add_argument('--copy-script', dest='copy_script_filename', 
                help="Name of copy script file", default=None)  # 'copy_files.sh'
args = ap.parse_args()

if args.date == 'check':
    needed = onv_dump.cull_tracking_file()
    print(f"Needed files")
    for cline in needed:
        filen = cline.split()[1].strip('"').split('/')[-1]
        print(f"\t{filen}")
else:
    onv_dump.gen_uvh5_dump_script(args.date, base_path=args.base_dir,
                                  ants=args.ants, pols=args.pols,
                                  LOs=args.los, CNODEs=args.cnodes,
                                  dump_script_filename=args.script,
                                  download_script_filename=args.download_script_filename,
                                  copy_script_filename=args.copy_script_filename)