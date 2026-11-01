from ._low_index import *

# Version 1.1 was used for the cython implemention
# Version 1.2 (and later) is used for the C++/pybind11 implementation
# Version 1.3 fixed performance regressions between 1.1 and 1.2.
# 1.3+seeded.covers.2 is a downstream seeded-cover extension of 1.3.

__version__ = '1.3+seeded.covers.2'

def version():
    return __version__
