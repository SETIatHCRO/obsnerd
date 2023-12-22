import requests
from html.parser import HTMLParser
from os import path
import argparse

ap = argparse.ArgumentParser()
ap.add_argument('--base-url', dest='base_url', help="Base url for tles",
                default='http://www.celestrak.com/NORAD/elements/')
ap.add_argument('--base-path', dest='base_path', help="Base path for tles",
                default='./tle')
args = ap.parse_args()


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
                print('Reading %s:  %s' % (lll, tlefiles[lll]))
                sat = requests.get(path.join(base_url, lll)).text.splitlines()
                with open(path.join(base_path, outfile), 'w') as fp:
                    for line in sat:
                        print(line, file=fp)
                print("{}:  {}".format(outfile, tlefiles[lll]), file=master)


if __name__ == '__main__':
    updatetle(base_path=args.base_path, base_url=args.base_url)

