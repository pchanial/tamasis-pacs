# Copyrights 2010-2011 Pierre Chanial
# All rights reserved
#
try:
    import pydebug
    pydebug.listen()
    del pydebug
except:
    pass
from operators import *
from .core import *
from .madcap import *
from .pacs import *

del datatypes, mappers, numpyutils, observations, processing, quantity, tamasisfortran
del core, madcap, pacs

__all__ = [ f for f in dir() if f[0] != '_' and f not in ('mpiutils', 'stringutils', 'tmf', 'var', 'wcsutils')]

