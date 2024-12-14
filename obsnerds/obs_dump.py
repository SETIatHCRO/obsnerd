import numpy as np
from copy import copy
from . import obs_sys as OS
from . import obs_look


def gen_uvh5_dump_script(date_path, base_path='/mnt/primary/ata/projects/p054/',
                         ants='all', pols='xx,xy,yy,yx',
                         LOs='all', CNODEs='all', script_filename='dump_autos.sh'):
    from os import walk, listdir, path
    if date_path == '?':
        print(f"Available observation dates in {base_path}:")
        for x in listdir(base_path):
            print(f"\t{x}")
        return
    LOs = OS.listify(LOs, {'all': OS.ALL_LOS})
    CNODEs = OS.make_cnode(CNODEs)

    dbase_path = path.join(base_path, date_path)
    print(f"Retrieving from {dbase_path}")
    files = {}
    for basedir, _, filelist in walk(dbase_path):
        if base_path in basedir and '/Lo' in basedir:
            for fn in filelist:
                dfn = path.join(basedir, fn)
                X = OS.parse_uvh5_filename(dfn)
                if X['lo'] in LOs and X['cnode'] in CNODEs:
                    files[X['obsrec']] = copy(X)
                
    with open(script_filename, 'w') as fp:
        for obsrec, data in files.items():
            print(f"on_dump_autos.py {data['filename']} --ants {ants} --pols {pols}]", file=fp)
            print(f"Adding {obsrec}")


class Dump:
    def __init__(self, obsinput=None):
        """

        Parameters
        ----------
        obsinput : str
            File to use

        """
        self.look = obs_look.Look(obsinput)

    def dump_autos(self, ants='all', pols='all'):
        ants = OS.listify(ants, {'all': self.look.ant_names})
        pols = OS.listify(pols, {'all': ['xx', 'xy', 'yy', 'yx']})
        outdata = {'ants': ants, 'freqs': self.look.freqs, 'pols': pols, 'source': self.look.source, 'uvh5': self.look.fn, 'freq_unit': self.look.freq_unit}
        print(f"Dumping autos in {self.look.fn} for {','.join(ants)} {','.join(pols)}", end=' ... ')
        for ant in ants:
            for pol in pols:
                self.look.get_bl(ant, pol=pol)
                outdata[f"{ant}{pol}"] = self.look.data
        outdata['times'] = self.look.times.jd  # This assumes that all times in the UVH5 file are the same...
        obsrec_file = f"{self.look.uvh5_pieces['obsrec']}.npz"
        print(f"writing {obsrec_file}")
        np.savez(obsrec_file, **outdata)