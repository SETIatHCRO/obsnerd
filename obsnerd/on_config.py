from os.path import isfile, join
import yaml
from . import DATA_PATH
import astropy.units as u


"""
This takes the old ant_file.dat, filter.json, and LO tunings and makes a config class
out of it. Dec 2025
"""

FREQ_UNIT_DEFAULT = 'MHz'
TUNING_BANDWIDTH_DEFAULT = 690.0

class Config:
    def __init__(self, config_file=None):
        """ Initialize configuration by loading from a YAML file.

        Args:
            config_file (str): Path to configuration file. If None, loads default config.yaml from data directory.
        """

        self.get_config(config_file)

    def get_config(self, config_file):
        if config_file is None:
            config_file = join(DATA_PATH, 'config.yaml')
        if not isfile(config_file):
            raise FileNotFoundError(f"Configuration file {config_file} not found.")
        self.config_file = config_file
        self.contents = yaml.safe_load(open(config_file, 'r'))
        self.proc_ants()
        self.proc_freqs()
        self.proc_filters()

    def proc_ants(self):
        """ Process antenna list from configuration. """
        self.ants = self.contents.get('Ants', None)

    def proc_freqs(self):
        """ Process frequency list from configuration. """
        self.tunings = self.contents.get('Tunings', None)
        self.freq_unit = u.Unit(self.contents.get('Freq_unit', FREQ_UNIT_DEFAULT))
        self.tuning_bandwidth = self.contents.get('Tuning_Bandwidth', TUNING_BANDWIDTH_DEFAULT) * self.freq_unit
        self.lo, self.freqs = [], []
        self.filters = {}
        if self.tunings is not None:
            for lo in sorted(self.tunings.keys()):
                self.tunings[lo] = self.tunings[lo] * self.freq_unit
                self.lo.append(lo)
                self.freqs.append(self.tunings[lo])
                self.filters[lo] = {}

    def proc_filters(self):
        """ Process filter settings from configuration. """
        self.filter_list = self.contents.get('Filters', [])
        for filt in self.filter_list:
            for i in range(len(filt['band'])):
                filt['band'][i] = filt['band'][i] * self.freq_unit
            for lo, freq in self.tunings.items():
                if filt['band'][0] >= freq - self.tuning_bandwidth / 2.0 and filt['band'][1] <= freq + self.tuning_bandwidth / 2.0:
                    self.filters[lo][filt['color']] = filt['band']