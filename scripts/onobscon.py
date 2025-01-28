#!/usr/bin/env python
from obsnerd import ono_engine as oe
import argparse

ap = argparse.ArgumentParser()
ap.add_argument('cmd', help="Action [start, freq,  move, end, note, source, summary]", choices=['start', 'freq', 'move', 'end', 'note', 'source', 'summary'])
ap.add_argument('-i', '--initials', help="Initials of user (cmd=start)", default=None)
ap.add_argument('-p', '--project_id', help="Project ID designator (cmd=start;backend, default=p054)", default='p054')
ap.add_argument('-g', '--group_ants', help="csv-list of antennas to reserve in group (cmd=start;end, default=rfsoc_active)", default='rfsoc_active')
ap.add_argument('-u', '--use_ants', help="csv-list of antennas to use or None (cmd=move;freq;end, default=<group_ants>)", default=None)
ap.add_argument('-f', '--frequency', help="Frequency of observation [GHz] (cmd=freq)", default=None)
ap.add_argument('-a', '--attenuation', help="Attenuation [dB] (cmd=freq, default=20)", type=int, default=20)
ap.add_argument('--lo', help="LO to set (cmd=freq, default=A)", choices=['A', 'B', 'C', 'D'], default='A')
ap.add_argument('-b', '--back_end', help="Type of backend (cmd=backend, default=xgpu)", choices=['xgpu'], default='xgpu')
ap.add_argument('-l', '--location', help="Location for coord_type (cmd=move)", default=None)
ap.add_argument('-c', '--coord_type', help="For <move> coordinate type: azel/radec/source/traj [azel] (cmd=move, default=azel)",
                choices=['azel', 'radec', 'source', 'traj'], default='azel')
ap.add_argument('-n', '--notation', help="Note to add (cmd=note)", default=None)
ap.add_argument('-s', '--source', help="Source name for log (cmd=source)", default=None)
ap.add_argument('-d', '--datetime', help="Reference datetime for source (cmd=source)", default=None)
args = ap.parse_args()

session = oe.CommandHandler()
getattr(session, args.cmd)(**vars(args))