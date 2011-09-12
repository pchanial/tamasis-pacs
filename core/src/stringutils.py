from __future__ import division

__all__ = []

def strenum(choices, last='or'):
    """
    Enumerates elements in a list as strings.

    Parameters
    ----------
    choices : list of string
        list of elements to be enumerated
    last : string
        last separator

    Examples
    --------
    >>> strenum(['blue', 'red', 'yellow'], 'and')
    "'blue', 'red' and 'yellow'"
    """
    choices = [ "'" + choice + "'" for choice in choices ]
    return ', '.join(choices[0:-1]) + ' ' + last + ' ' + choices[-1]
