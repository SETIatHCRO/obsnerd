#!/usr/bin/env python
try:
    from ATATools import ata_control
except ImportError:
    ata_control = None
import atexit
import metadata
import onutil


class CommandHandler:
    def __init__(self, payload=None, arg2=None, coord_type=None, group_ants=['1a'], use_ants=['1a']):
        self.payload = payload
        self.arg2 = arg2
        self.coord_type = coord_type
        self.group_ants = group_ants
        self.use_ants = use_ants
    
    def start(self, initials=None):
        self.initials = self.payload if initials is None else initials
        if self.initials is None:
            print("Please include your name or initials.")
        else:
            ata_control.move_ant_group(self.group_ants, 'none', 'atagr')
            metadata.onlog(f"session start: {self.initials} -- {', '.join(self.use_ants)} / {', '.join(self.group_ants)}")

    def end(self):
        atexit.register(ata_control.move_ant_group, self.group_ants, 'atagr', 'none')
        atexit.register(ata_control.park_antennas, self.use_ants)
        metadata.onlog(f"end: {', '.join(self.use_ants)} / {', '.join(self.group_ants)}")

    def freq(self, freq=None, att=20):
        self.freq = float(self.payload) if freq is None else freq
        self.att = att
        metadata.onlog(f"fcen: {self.freq}")
        ata_control.set_freq(self.freq, self.use_ants, lo='d')
        ata_control.autotune(self.use_ants)
        ata_control.rf_switch_thread(self.use_ants)
        ata_control.set_atten_thread([[f'{ant}x', f'{ant}y'] for ant in self.use_ants], [[self.att, self.att] for ant in self.use_ants])

    def move(self, location=None, coord_type=None):
        self.location = self.payload if location is None else location
        self.coord_type = self.coord_type if coord_type is None else coord_type

        if ',' in self.location:
            x, y = [float(_v) for _v in self.location.split(',')]

        metadata.onlog(f'move to: {self.coord_type}')
        if self.coord_type == 'azel':
            ata_control.set_az_el(self.use_ants, x, y)
            metadata.onlog(f"azel: {x},{y}")
        elif self.coord_type == 'radec':
            source = ata_control.track_source(self.use_ants, radec=[x, y])
            metadata.onlog(f"radec: {x},{y}")
        elif self.coord_type == 'source':
            source = ata_control.track_source(self.use_ants, source=self.location)
            metadata.onlog(f"source: {self.location}")
        elif self.coord_type == 'traj':
            ephem = ata_control.upload_ephemeris(self.location)
            ata_control.track_ephemeris(ephem, self.use_ants, wait=True)
            metadata.onlog(f"traj: {self.location}")
            try:
                with open('track.log', 'r') as fp:
                    for line in fp:
                        if len(line) > 2:
                            metadata.onlog(line.strip())
            except FileNotFoundError:
                pass
    
    def note(self, notation=None):
        self.notation = self.payload if notation is None else notation
        if self.notation is None:
            print("Need to include a note.")
        else:
            metadata.onlog(f"note: {self.notation}")

    def source(self, name=None, datestamp=None):
        self.name = self.payload if name is None else name
        self.datestamp = self.arg2 if datestamp is None else datestamp
        if self.datestamp is None:
            from datetime import datetime
            self.datestamp = datetime.now()
            print(f"Using current as date: {self.datestamp}")
        else:
            self.datestamp = onutil.make_datetime(date=self.datestamp)
        metadata.onlog([f'source: {self.name}', f'expected: {self.datestamp.isoformat()}'])


if __name__ == '__main__':
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument('cmd', help="Action [start, freq,  move, end, note, source]", choices=['start', 'freq', 'move', 'end', 'note', 'source'])
    ap.add_argument('payload', help="Argument for command.", nargs='?', default=None)
    ap.add_argument('arg2', help="If source, associated datetime.", nargs='?', default=None)
    ap.add_argument('-c', '--coord-type', dest='coord_type', help="For <move> coordinate type: azel/radec/source/traj [azel]",
                    choices=['azel', 'radec', 'source', 'traj'], default='azel')
    args = ap.parse_args()

    session = CommandHandler(payload=args.payload, arg2=args.arg2, coord_type=args.coord_type)
    getattr(session, args.cmd)()