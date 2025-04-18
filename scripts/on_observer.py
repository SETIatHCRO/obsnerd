#! /usr/bin/env python
import argparse
from obsnerd import ono_observer

"""
This reads in a source file and step through observations
"""

ap = argparse.ArgumentParser()
ap.add_argument('source', help='Target list or json file with inputs', nargs='?', default="/opt/mnt/share/ods_project/ods_rados.json")
ap.add_argument('-t', '--source_type', help="Type of source above", choices=['ods'], default='ods')
ap.add_argument('-o', '--observer', help="Name of observer", default="Was a loser and didn't give an observer")
ap.add_argument('-n', '--project_name', help="Name of project", default="Was a loser and didn't give a project name")
ap.add_argument('-p', '--project_id', help="Project ID for observations", default="p054")
ap.add_argument('-a', '--ant_list', help="List of antennas or group list", default='rfsoc_active')
ap.add_argument('--embargo', help='List of known antennas to not include', default='1k')
ap.add_argument('--focus', help="Focus parameter", choices=['a', 'b', 'max'], default='max')
ap.add_argument('--go', help="Flag to go 'live' for now for testing", action='store_true')
ap.add_argument('-c', '--calendar', help="Flag to add to calendar", action='store_true')
# ap.add_argument('--data_record', help="Data recording to use", choices=['gnuradio', 'hpguppi'], default='hpguppi')
args = ap.parse_args()

observer = ono_observer.Observer(observer=args.observer, project_name=args.project_name, project_id=args.project_id, ants=args.ant_list, embargo=args.embargo)
if args.source_type == 'ods':
    observer.get_ods(args.source, defaults='defaults.json')
    observer.get_obs_from_ods(add_to_calendar=args.calendar)
    observer.update_ods("https://www.seti.org/sites/default/files/HCRO/ods.json", "/opt/mnt/share/ods_upload/ods.json")
observer.observe(args.go)