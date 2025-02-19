import numpy as np
from copy import copy
from . import on_sys
from . import onv_look, on_track
from odsutils import ods_tools as tools


def gen_uvh5_dump_script(date_path, base_path='/mnt/primary/ata/projects/p054/',
                         ants='all', pols='xx,xy,yy,yx',
                         LOs='all', CNODEs='all', script_filename='dump_autos.sh'):
    from os import walk, listdir, path
    if date_path == '?':
        print(f"Available observation dates in {base_path}:")
        for x in listdir(base_path):
            print(f"\t{x}")
        return
    LOs = tools.listify(LOs, {'all': on_sys.ALL_LOS})
    CNODEs = on_sys.make_cnode(CNODEs)

    dbase_path = path.join(base_path, date_path)
    print(f"Retrieving from {dbase_path}")
    files = {}
    fp = open('download_files.sh', 'w')
    for basedir, _, filelist in walk(dbase_path):
        if base_path in basedir and '/Lo' in basedir:
            for fn in filelist:
                dfn = path.join(basedir, fn)
                X = on_sys.parse_uvh5_filename(dfn)
                if X['lo'] in LOs and X['cnode'] in CNODEs:
                    files[X['obsrec']] = copy(X)
                    fp.write(f'scp "sonata@obs-node1.hcro.org:./rfsoc_obs_scripts/p054/{X["obsrec"]}.npz" .\n')
    fp.close()
    with open(script_filename, 'w') as fp:
        for obsrec, data in files.items():
            print(f"on_dump_autos.py {data['filename']} --ants {ants} --pols {pols}", file=fp)
            print(f"Adding {obsrec}")


def redump(obsid, ants, pols='xx,yy,xy,yx'):
    from glob import glob
    import shutil
    keys2use = ['times', 'freqs', 'source', 'uvh5', 'freq_unit']
    usepols = set()
    useants = set()
    for fn in glob(f'{obsid}*.npz'):
        print(f"Redumping {fn}")
        shutil.copy(fn, f"{fn}.arc")
        outdata = {}
        data = np.load(fn)
        for key in keys2use:
            outdata[key] = data[key]
        for ant in ants.split(','):
            for pol in pols.split(','):
                if f"{ant}{pol}" in data:
                    usepols.add(pol)
                    useants.add(ant)
                    outdata[f"{ant}{pol}"] = data[f"{ant}{pol}"]
        outdata['pols'] = list(usepols)
        outdata['ants'] = list(useants)
        np.savez(f"{fn}", **outdata)


def cull_tracking_file(filename='download_files.sh'):
    from glob import glob
    npz = glob('*.npz')
    needed_files = []
    with open(filename, 'r') as fp:
        for line in fp:
            line_present = False
            for fn in npz:
                if fn in line:
                    line_present = True
                    break
            if not line_present:
                needed_files.append(line.strip())
    return needed_files

class Dump:
    def __init__(self, obsinput=None, lo='A', cnodes='all'):
        """

        Parameters
        ----------
        obsinput : str
            File to use generally an obsid
        lo : str, list, 'all'
        cnodes : str, list, 'all'

        """
        self.obsinput = obsinput
        self.lo = lo
        self.cnodes = cnodes

    def dump_autos(self, ants='all', pols='all'):
        """
        self.obsinput should be an obsid

        """
        self.look = onv_look.Look(self.obsinput, lo=self.lo, cnode=self.cnodes)
        ants = tools.listify(ants, {'all': self.look.ant_names})
        pols = tools.listify(pols, {'all': ['xx', 'xy', 'yy', 'yx']})
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

    def dump_jupyter(self, ants='2b', pols='xx,xy,yy'):
        """
        self.obsinput should be an obsinfo file

        """
        filters = {}
        filters['on'] = onv_look.Filter(ftype='time', unit='degrees', lo=-5, hi=5, norm=True, color='r')
        filters['off'] = onv_look.Filter(ftype='time', unit='degrees', lo=-5, hi=5, norm=True, color='k', invert=True)
        filters['adjacent_feature'] = onv_look.Filter(color='r', ftype='freq', unit='MHz', lo=1975.0, hi=1985.0)
        filters['dtz'] = onv_look.Filter(color='r', ftype='freq', unit='MHz', lo=1990.0, hi=1995.0, norm=True)
        filters['low'] = onv_look.Filter(color='r', ftype='freq', unit='MHz', lo=1910.0, hi=1915.0, norm=True)

        for i, obsid in enumerate(base.obsinfo.array.name):
            print(f"Reading {obsid}")
            look = onv_look.Look(obsid, self.lo, cnode=self.cnodes)
            look.get_time_axes()
            if not i:
                ants = tools.listify(ants, {'all': look.ant_names})
                pols = tools.listify(pols, {'all': ['xx', 'xy', 'yy', 'yx']})
                outdata = {'ants': ants, 'freqs': look.freqs, 'pols': pols}
            else:
                if abs(outdata['freqs'][0] - look.freqs[0]) > 1.0:
                    print(f"Skipping {obsid}:  Frequencies don't match: {outdata['freqs'][0]} vs {look.freqs[0]}")
                    continue
            outdata[obsid] = {'tref': look.obs.obsinfo.obsid[obsid].tref.datetime.isoformat(timespec='seconds'),
                              'boresight': look.taxes['boresight']['values'], 'seconds': look.taxes['seconds']['values']}
            for ant in ants:
                for pol in pols:
                    look.get_bl(ant, pol=pol)
                    outdata[obsid][f"{ant}{pol}"] = np.abs(look.data)
                    for key in ['on', 'off']:
                        filters[key].apply(look.taxes['boresight']['values'], look.data)
                        outdata[obsid][key] = np.abs(filters[key].power)
                    for key in ['low', 'adjacent_feature', 'dtz']:
                        filters[key].apply(look.freqs, look.data)
                        outdata[obsid][key] = np.abs(filters[key].power)
                
        fnout = f"{self.obsinput.split('.')[0]}.npz"
        print(f"Writing {fnout}")
        np.savez(fnout, **outdata, allow_pickle=True)

        # 
        # print(f"Dumping autos in {self.look.fn} for {','.join(ants)} {','.join(pols)}", end=' ... ')
        # for ant in ants:
        #     for pol in pols:
        #         self.look.get_bl(ant, pol=pol)
        #         outdata[f"{ant}{pol}"] = self.look.data
        # outdata['times'] = self.look.times.jd  # This assumes that all times in the UVH5 file are the same...
        # obsrec_file = f"{self.look.uvh5_pieces['obsrec']}.npz"
        # print(f"writing {obsrec_file}")
        # np.savez(obsrec_file, **outdata)