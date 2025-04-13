from importlib.metadata import version
__version__ = version('obsnerd')
from os.path import join

DATA_PATH = join(__path__[0], 'data')
LOG_FORMATS = {'conlog_format': "{asctime} - {levelname} - {module} - {message}",
               'filelog_format': "{asctime} - {levelname} - {module} - {message}"}
LOG_FILENAME = 'onlog'