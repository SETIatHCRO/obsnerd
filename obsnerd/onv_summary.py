from . import onv_logs
from glob import glob
import json
from tabulate import tabulate


class Summary:
    def __init__(self, ods, obsinfo, log, directory='data'):
        """
        Parameters
        ----------
        obsinput : str
            File to use generally an obsid
        lo : str, list, 'all'
        cnodes : str, list, 'all'

        """
        self.obsinfo_file = obsinfo
        self.obsinfo = json.load(open(obsinfo, 'r'))
        self.ods_file = ods
        self.ods = json.load(open(ods, 'r'))
        self.obs_files = glob(f"{directory}/*.npz")
        self.log = onv_logs.TBALog(log)

    def process(self):
        """
        Process the data and create a summary.
        """
        oidict = {}
        for sat, obsinfo in self.obsinfo['Sources'].items():
            if sat.startswith('STARLINK'):
                this_sat = int(sat.split('K')[1].split('D')[0])
                oidict[this_sat] = {'utc': obsinfo['utc'], 'ods': obsinfo['ods']}
        oddict = {}
        for sat in self.ods['ods_data']:
            if sat['src_id'].startswith('STARLINK'):
                this_sat = int(sat['src_id'].split('K')[1].split('D')[0])
                oddict[this_sat] = {'utc': sat['src_start_utc'], 'ods': sat['notes'].split(';')[1].split(':')[1].strip()}
        obdict = {}
        for sat in self.obs_files:
            if 'STARLINK' in sat:
                this_sat = int(sat.split('K')[1].split('D')[0])
                obdict.setdefault(this_sat, {})
                # obdict[this_sat].setdefault('cnode', [])
                # obdict[this_sat]['cnode'].append(int(sat.split('_')[-1].split('.')[0][1:]))
                obdict[this_sat].setdefault('cnode', 0)
                obdict[this_sat]['cnode'] += 1
        self.summary = {}
        allsats = set(oidict.keys()).union(set(oddict.keys())).union(set(obdict.keys()))
        for sat in allsats:
            self.summary[sat] = {'obsinfo': {}, 'ods': {}, 'obs': {}, 'inner': [], 'outer': []}
            if sat in oidict:
                self.summary[sat]['obsinfo'] = oidict[sat]
            else:
                self.summary[sat]['obsinfo'] = {'utc': '', 'ods': ''}
            if sat in oddict:
                self.summary[sat]['ods'] = oddict[sat]
            else:
                self.summary[sat]['ods'] = {'utc': '', 'ods': ''}
            if sat in obdict:
                self.summary[sat]['obs'] = obdict[sat]
            else:
                self.summary[sat]['obs'] = {'cnode': 0}
        for scope in ['inner', 'outer']:
            for t, tsat in zip(self.log.times[scope], self.log.sats[scope]):
                this_sat = int(tsat)
                if this_sat in allsats:
                    self.summary[this_sat][scope].append(t.datetime.isoformat())

    def show(self):
        """
        Show the summary.
        """
        header = ['Satellite', 'OBSINFO', 'ODS', 'OBS', 'Inner', 'Outer']
        table_data = []
        for sat, data in self.summary.items():
            try:
                a = data['obsinfo']['utc']
            except KeyError:
                a = ''
            try:
                b = data['ods']['utc']
            except KeyError:
                b = ''
            try:
                c = data['obs']['cnode']
            except KeyError:
                c = 0
            try:
                d = ','.join(data['inner'])
            except KeyError:
                d = ''
            try:
                e = ','.join(data['outer'])
            except KeyError:
                e = ''
            row = [sat, a, b, c, d, e]
            table_data.append(row)
        print(tabulate(table_data, headers=header, tablefmt='grid'))
