import requests
from html.parser import HTMLParser
from os import path
import logging
from odsutils import logger_setup
from . import LOG_FILENAME


logger = logging.getLogger(__name__)
logger.setLevel('DEBUG')  # Set to lowest
logger_setup.Logger(logger, conlog='INFO', filelog='INFO', log_filename=LOG_FILENAME, path='.')

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
        #if 'Active Satellites' in line:
        if '.txt' in line:
            print(line)
        data = line.split('"')
        for word in data:
            if '.txt' in word and word.startswith('/NORAD/elements/'):
                hpp = DataParser()
                hpp.feed(line)
                this_file = word.split('/')[-1]
                tlefiles[this_file] = hpp.description

    print(f"Updating from {base_url}")
    with open(path.join(base_path, 'master.dat'), 'w') as master:
        for lll in tlefiles:
            useThis = True
            for ig in ignore:
                if ig in tlefiles[lll].lower():
                    useThis = False
            if useThis:
                a = lll.split('.')
                outfile = a[0]+'.tle'
                reading = f'{tlefiles[lll]:30s}  {lll}'
                sat = requests.get(path.join(base_url, lll))
                if sat.status_code == 404:
                    new_url = base_url + f"gp.php?GROUP={lll.split('.txt')[0]}&FORMAT=tle"
                    reading = f"{tlefiles[lll]:30s}  {new_url.split('?')[1]}"
                    sat = requests.get(new_url)

                if sat.status_code != 200:
                    print(f"Invalid code for {tlefiles[lll]}:  {sat.status_code}")
                else:
                    print(reading)
                    sat_text = sat.text.splitlines()
                    with open(path.join(base_path, outfile), 'w') as fp:
                        for line in sat_text:
                            print(line, file=fp)
                print("{}:  {}".format(outfile, tlefiles[lll]), file=master)

def update_log():
    logger.info('Updating TLEs')
