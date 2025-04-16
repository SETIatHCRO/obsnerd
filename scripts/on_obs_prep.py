#! /usr/bin/env python
from obsnerd import ono_observer
import argparse

ap = argparse.ArgumentParser()
ap.add_argument('--add-to-calendar', dest='add_to_calendar', help="Add the observations to the calendar", action='store_true')
ap.add_argument('--ods2use', help="ODS file to use for observations", default='/opt/mnt/share/ods_rados/ods_rados.json')
ap.add_argument('--ods-upload', dest='ods_upload', help="ODS file to upload for TBA", default="/opt/mnt/share/ods_upload/ods.json")
ap.add_argument('--ods-active', dest='ods_active', help="Active ODS location", default="https://www.seti.org/sites/default/files/HCRO/ods.json")
args = ap.parse_args()

# Initialize the observer
observer = ono_observer.Observer()
observer.observe_prep(add_to_calendar=args.add_to_calendar,
                     ods2use=args.ods2use,
                     ods_upload=args.ods_upload,
                     ods_active=args.ods_active)
