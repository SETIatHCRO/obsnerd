from importlib.metadata import version
__version__ = version('obsnerd')
from logging import config
from os.path import join

DATA_PATH = join(__path__[0], 'data')
LOG_FORMATS = {'conlog_format': "{asctime} - {levelname} - {module} - {message}",
               'filelog_format': "{asctime} - {levelname} - {module} - {message}"}
LOG_FILENAME = 'onlog'


def get_config(config_file=None ):
    """ Load configuration file. If no file is specified, load default config.yaml
    from data directory.

    Args:
        config_file (str): Path to configuration file.
    Returns:
        dict: Configuration dictionary.
    Raises:
        FileNotFoundError: If the configuration file does not exist.
    """
    from os.path import isfile
    import yaml

    if config_file is None:
        config_file = join(DATA_PATH, 'config.yaml')
    if not isfile(config_file):
        raise FileNotFoundError(f"Configuration file {config_file} not found.")
    return yaml.safe_load(open(config_file, 'r'))
