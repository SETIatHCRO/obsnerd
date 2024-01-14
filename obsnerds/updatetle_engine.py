import requests
from html.parser import HTMLParser
from os import path
from obsnerd import metadata


class DataParser(HTMLParser):
    def handle_data(self, data):
        self.description = data


def updatetle(base_path, base_url):
    master_file = requests.get(base_url)
    master = master_file.text.splitlines()
    tlefiles = {}
    base_path = path.expanduser(base_path)

    ignore = ['debris', 'cesium']

    print('Reading Celestrak master file')
    for line in master:
        data = line.split('"')
        for word in data:
            if '.txt' in word and word.startswith('/NORAD/elements/'):
                hpp = DataParser()
                hpp.feed(line)
                this_file = word.split('/')[-1]
                tlefiles[this_file] = hpp.description

    with open(path.join(base_path, 'master.dat'), 'w') as master:
        for lll in tlefiles:
            useThis = True
            for ig in ignore:
                if ig in tlefiles[lll].lower():
                    useThis = False
            if useThis:
                a = lll.split('.')
                outfile = a[0]+'.tle'
                print(f'Reading {lll}:  {tlefiles[lll]}')
                sat = requests.get(path.join(base_url, lll))
                if sat.status_code == 404:
                    print(f"Normal request not working for {lll}", end=' -- ')
                    new_url = base_url + f"gp.php?GROUP={lll.split('.txt')[0]}&FORMAT=tle"
                    print(f"Trying {new_url}", end=' -- ')
                    sat = requests.get(new_url)
                    print(sat.status_code)

                if sat.status_code == 200:
                    sat_text = sat.text.splitlines()
                    with open(path.join(base_path, outfile), 'w') as fp:
                        for line in sat_text:
                            print(line, file=fp)
                print("{}:  {}".format(outfile, tlefiles[lll]), file=master)

def update_log():
    metadata.onlog('Updating TLEs')
