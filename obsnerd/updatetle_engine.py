import requests
from bs4 import BeautifulSoup
from os import path
import logging
from odsutils import logger_setup
from . import LOG_FILENAME


logger = logging.getLogger(__name__)
logger.setLevel('DEBUG')  # Set to lowest
logger_setup.Logger(logger, conlog='INFO', filelog='INFO', log_filename=LOG_FILENAME, path='.')


def make_tle_filename(tle_name):
    """
    Create a TLE filename from the TLE name.
    """
    tle_name = tle_name.strip()
    tle_name = tle_name.replace(' ', '_').replace('/', '_')
    tle_name = tle_name.replace('(', '').replace(')', '')
    tle_name = tle_name.replace("'", '').replace('"', '')
    tle_name = tle_name.replace(',', '').replace('.', '')
    tle_name = tle_name.replace('?', '').replace('!', '')
    tle_name = tle_name.replace('&', 'and')
    return f"{tle_name}.tle"


def updatetle(group='*', base_path='./tle', base_url='https://celestrak.org/NORAD/elements/'):
    if group == '*':
        group = ''
    master_file = requests.get(base_url)
    soup = BeautifulSoup(master_file.text, 'html.parser')
    for td in soup.find_all('td'):
        this_href = td.find('a')
        try:
            ttype = this_href.get('title')
        except AttributeError:
            continue
        if "debris" in td.text.lower() or "cesium" in td.text.lower():
            continue
        if ttype == 'TLE Data' and group in td.text:
            actual_href = this_href.get('href')
            groupname = actual_href.split('=')[1].split('&')[0]
            tlefilename = path.join(base_path, make_tle_filename(groupname))
            print(f"{td.text} - {tlefilename}:  {actual_href}")
            tle_url = path.join(base_url, actual_href)
            try:
                tle_file = requests.get(tle_url, timeout=8)
            except Exception as e:
                print(e)
                continue
            with open(tlefilename, 'w') as f:
                f.write(tle_file.text)

 
def update_log():
    logger.info('Updating TLEs')
