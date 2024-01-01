#!/usr/bin/python3
from ATATools import ata_control
import atexit
import metadata


class CommandHandler:
    def __init__(self, payload=None, coord_type=None, group_ants=['1a', '1f', '5c'], use_ants=['1a']):
        self.payload = payload
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

        if self.coord_type == 'azel':
            ata_control.set_az_el(self.use_ants, x, y)
            metadata.onlog(f"az: {x}, el: {y}")
        elif self.coord_type == 'radec':
            source = ata_control.track_source(self.use_ants, radec=[x, y])
            metadata.onlog(f"ra: {x}, dec: {y}")
        elif self.coord_type == 'source':
            source = ata_control.track_source(self.use_ants, source=self.location)
            metadata.onlog(f"source {self.location}")
        elif self.coord_type == 'traj':
            ephem = ata_control.upload_ephemeris(self.location)
            ata_control.track_ephemeris(ephem, self.use_ants, wait=True)
            metadata.onlog(f"trajectory {self.location}")

    def note(self, notation=None):
        self.notation = self.payload if notation is None else notation
        if self.notation is None:
            print("Need to include a note.")
        else:
            metadata.onlog(f"note: {self.notation}")


# if args.cmd == 'start':
#     if args.payload is None:
#         print("Please include your name or initials.")
#     else:
#         ata_control.move_ant_group(ants, 'none', 'atagr')
#         metadata.onlog(f"session start: {args.payload}")
# elif args.cmd == 'end':
#     atexit.register(ata_control.move_ant_group, ants, 'atagr', 'none')
#     atexit.register(ata_control.park_antennas, ants)
#     metadata.onlog(f"end: {', '.join(ants)}")
# elif args.cmd == 'freq':
#     freq = defaults.freq if args.payload is None else float(args.payload)
#     metadata.onlog(f"fcen: {freq}")
#     ata_control.set_freq(freq, [ants[0]], lo='d')
#     ata_control.autotune([ants[0]])
#     att = 20  # Attenuation in dB
#     ata_control.rf_switch_thread(ants)
#     ata_control.set_atten_thread([[f'{ant}x', f'{ant}y'] for ant in ants], [[att, att] for ant in ants])
# elif args.cmd == 'move':
#     if args.payload is None:
#         x, y = defaults.x, defaults.y
#     elif ',' in args.payload:
#         x, y = [float(_v) for _v in args.payload.split(',')]
#     else:
#         x = args.payload

#     if args.coord_type == 'azel':
#         ata_control.set_az_el([ants[0]], x, y)
#         metadata.onlog(f"az: {x}, el: {y}")
#     elif args.coord_type == 'radec':
#         source = ata_control.track_source([ants[0]], radec=[x, y])
#         metadata.onlog(f"ra: {x}, dec: {y}")
#     elif args.coord_type == 'source':
#         source = ata_control.track_source([ants[0]], source=x)
#         metadata.onlog(f"source {x}")
#     elif args.coord_type == 'traj':
#         print("MAKE THIS ONE WORK")
#         metadata.onlog(f"trajectory {x}")
# elif args.cmd == 'note':
#     if args.payload is None:
#         print("Need to include a note.")
#     else:
#         metadata.onlog(f"note: {args.payload}")


if __name__ == '__main__':
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument('cmd', help="Action [start, freq,  move, end, note]", choices=['start', 'freq', 'move', 'end', 'note'])
    ap.add_argument('payload', help="Argument for command.", nargs='?', default=None)
    ap.add_argument('-c', '--coord-type', dest='coord_type', help="Coordinate type: azel/radec/source/traj [azel]",
                    choices=['azel', 'radec', 'source', 'traj'], default='azel')
    args = ap.parse_args()

    session = CommandHandler(payload=args.payload, coord_type=args.coord_type)
    getattr(session, args.cmd)()