#! /usr/bin/env python
from obsnerd import ono_observer
import argparse

ap = argparse.ArgumentParser()
ap.add_argument('--add-to-calendar', dest='add_to_calendar', help="Add the observations to the calendar", action='store_true')
ap.add_argument('--input-obs', dest='input_obs', help="Input obsinfo file to use for observation preparation", default='obsinfo_rados.json')
ap.add_argument('--ods-rados', dest='ods_rados', help="ODS file to use for ODS publishing", default='/opt/mnt/share/ods_project/ods_rados.json')
ap.add_argument('--ods-upload', dest='ods_upload', help="ODS file to upload for TBA", default="/opt/mnt/share/ods_upload/ods.json")
# Note the inverted logic...
ap.add_argument('--skip-source-database', dest='update_source_database', help="[SKIP] Update the source database with the ODS entries", action='store_false')
ap.add_argument('--skip-ods-assembly', dest='ods_assembly', help="[SKIP] Assembling and publishing the ODS file", action='store_false')
args = ap.parse_args()

# Initialize the observer
observer = ono_observer.Observer(obsfile=args.input_obs, project_name='SatSpot', project_id='p054')
observer.observe_prep(add_to_calendar=args.add_to_calendar,
                      ods_rados=args.ods_rados,
                      ods_upload=args.ods_upload,
                      update_source_database=args.update_source_database,
                      ods_assembly=args.ods_assembly)
