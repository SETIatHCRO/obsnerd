#!/usr/bin/python3
from ATATools import ata_control
import argparse
import atexit
import metadata


class CommandHandler:
    def __init__(self, cmd, payload, coord_type, ants=['1a', '1f', '5c']):
        self.cmd = cmd
        self.payload = payload
        self.coord_type = coord_type
        self.ants = ants
        self.use_ants = [ants[0]]  # Hard-coded for now
    
    def start(self, initials=None):
        initials = self.payload if initials is None else initials
        if initials is None:
            print("Please include your name or initials.")
        else:
            ata_control.move_ant_group(self.ants, 'none', 'atagr')
            metadata.onlog(f"session start: {self.payload}")

    def end(self):
        atexit.register(ata_control.move_ant_group, self.ants, 'atagr', 'none')
        atexit.register(ata_control.park_antennas, self.ants)
        metadata.onlog(f"end: {', '.join(self.ants)}")

    def freq(self, freq=None):
        freq = float(self.payload) if freq is None else freq
        metadata.onlog(f"fcen: {freq}")
        ata_control.set_freq(freq, self.use_ants, lo='d')
        ata_control.autotune(self.use_ants)
        att = 20  # Attenuation in dB
        ata_control.rf_switch_thread(self.ants)
        ata_control.set_atten_thread([[f'{ant}x', f'{ant}y'] for ant in self.ants], [[att, att] for ant in self.ants])

    def move(self, location=None, coord_type=None):
        location = self.payload if location is None else location
        coord_type = self.coord_type if coord_type is None else coord_type

        if ',' in location:
            x, y = [float(_v) for _v in location.split(',')]
        else:
            x = location

        if coord_type == 'azel':
            ata_control.set_az_el(self.use_ants, x, y)
            metadata.onlog(f"az: {x}, el: {y}")
        elif coord_type == 'radec':
            source = ata_control.track_source(self.use_ants, radec=[x, y])
            metadata.onlog(f"ra: {x}, dec: {y}")
        elif coord_type == 'source':
            source = ata_control.track_source(self.use_ants, source=x)
            metadata.onlog(f"source {x}")
        elif coord_type == 'traj':
            print("MAKE THIS ONE WORK")
            # metadata.onlog(f"trajectory {x}")

    def note(self, notation=None):
        notation = self.payload if notation is None else notation
        if self.payload is None:
            print("Need to include a note.")
        else:
            metadata.onlog(f"note: {notation}")


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
    ap = argparse.ArgumentParser()
    ap.add_argument('cmd', help="Action [start, freq,  move, end, note]", choices=['start', 'freq', 'move', 'end', 'note'])
    ap.add_argument('payload', help="Argument for command.", nargs='?', default=None)
    ap.add_argument('-c', '--coord-type', dest='coord_type', help="Coordinate type: azel/radec/source/traj [azel]",
                    choices=['azel', 'radec', 'source', 'traj'], default='azel')
    args = ap.parse_args()

    session = CommandHandler(cmd=args.cmd, payload=args.payload, coord_type=args.coord_type)
    getattr(session, args.cmd)()