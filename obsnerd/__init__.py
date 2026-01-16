from importlib.metadata import version
__version__ = version('obsnerd')
from os.path import join
import math


DATA_PATH = join(__path__[0], 'data')
LOG_FORMATS = {'conlog_format': "{asctime} - {levelname} - {module} - {message}",
               'filelog_format': "{asctime} - {levelname} - {module} - {message}"}
LOG_FILENAME = 'onlog'


def approx_equal(a, b, rel_tol=1e-9, abs_tol=0.0):
    """
    Return True if a and b are approximately equal, using a mix of relative and absolute tolerance.

    - rel_tol controls comparison relative to magnitude (good for big numbers)
    - abs_tol controls a fixed tolerance floor (good near zero)
    """
    # Fast path for exact equality (also handles +0.0 == -0.0)
    if a == b:
        return True

    # Treat NaNs as never close; infinities only close if they are exactly equal (handled above)
    if math.isnan(a) or math.isnan(b):
        return False
    if math.isinf(a) or math.isinf(b):
        return False

    diff = abs(a - b)
    return diff <= max(rel_tol * max(abs(a), abs(b)), abs_tol)
