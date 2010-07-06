import os
__verbose__ = False
tamasis_dir = os.path.dirname(__file__) + '/../' if os.path.dirname(__file__) != '' else '../'
del os

__version_info__ = (1, 0, 4)
__version__ = '.'.join((str(i) for i in __version_info__))

__all__ = [ 'tamasis_dir', '__verbose__', '__version__', '__version_info__' ]

