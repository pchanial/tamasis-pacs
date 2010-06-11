import numpy
import re

__all__ = ['DerivedUnits', 'Quantity', 'UnitError']

_re_unit = re.compile(r' *([/*])? *([a-zA-Z_]+|\?+)(\^-?[0-9]+(\.[0-9]*)?)? *')

class UnitError(Exception): pass

def _extract_unit(string):
    """
    Convert the input string into a unit as a dictionary
    """

    if string is None or isinstance(string, dict):
        return string

    if not isinstance(string, str):
        raise TypeError('Invalid unit type. Valid types are string or dictionary.')

    string = string.strip()
    result = None
    start = 0
    while start < len(string):
        match = _re_unit.match(string, start)
        if match is None:
            raise ValueError("Unit '"+string[start:]+"' cannot be understood.")
        op = match.group(1)
        u = match.group(2)
        exp = match.group(3)
        if exp is None:
            exp = 1.
        else:
            exp = float(exp[1:])
        if op == '/':
            exp = -exp
        result = _multiply_unit_inplace(result, u, exp)
        start = start + len(match.group(0))
    return result


def _multiply_unit_inplace(unit, key, val):
    """
    Multiply a unit as a dictionary by a non-composite one
    Unlike _divide_unit, _multiply_unit and _power_unit,
    the operation is done in place, to speed up main
    caller _extract_unit.
    """
    if unit is None:
        return { key : val }
    if unit.has_key(key):
        if unit[key] == -val:
            del unit[key]
        else:
            unit[key] += val
    else:
        unit[key] = val
    return unit

def _normalize_quantity(value):
    """
    Recursively convert a Quantity by looking up the DerivedUnits table
    """
    if value._unit is None:
        return value
    result = Quantity(value, '')
    for key, val in value._unit.iteritems():
        if (key,val) in DerivedUnits.table:
            result *= _normalize_quantity(DerivedUnits.table[(key,val)])
        elif (key,-val) in DerivedUnits.table:
            result /= _normalize_quantity(DerivedUnits.table[(key,-val)])
        elif key in DerivedUnits.table:
            if val == 1.:
                result *= _normalize_quantity(DerivedUnits.table[key])
            elif val == -1.:
                result /= _normalize_quantity(DerivedUnits.table[key])
            else:
                result *= _normalize_quantity(DerivedUnits.table[key]**val)
        else:
            if result._unit is None:
                result._unit = { key: val }
            elif key in result._unit:
                result._unit[key] += val
            else:
                result._unit[key] = val
    return result    

def _power_unit(unit, power):
    """
    Raise to power a unit as a dictionary
    """
    if unit is None or power == 0:
        return None
    result = unit.copy()
    for key in unit:
        result[key] *= power
    return result

def _multiply_unit(unit1, unit2):
    """
    Multiplication of units as dictionary. 
    Unlike _power_unit, _divide_unit, the operation is done in place
    of unit1, to speed up main caller _extract_unit.
    """
    if unit2 is None:
        return unit1
    if unit1 is None:
        return unit2
    unit = unit1.copy()
    for key, val in unit2.iteritems():
        if unit.has_key(key):
            if unit[key] == -val:
                del unit[key]
            else:
                unit[key] += val
        else:
            unit[key] = val
    return unit

def _divide_unit(unit1, unit2):
    """
    Division of units as dictionary
    """
    if unit2 is None:
        return unit1
    if unit1 is None:
        return _power_unit(unit2, -1)
    unit = unit1.copy()
    for key, val in unit2.iteritems():
        if unit.has_key(key):
            if unit[key] == val:
                del unit[key]
            else:
                unit[key] -= val
        else:
            unit[key] = -val
    return unit

def _get_units(units):
    return map(lambda q: getattr(q, '_unit', None), units)

def _strunit(unit):
    """
    Convert a unit as dictionary into a string
    """
    if unit is None:
        return ''
    result = ''
    has_pos = False

    for key, val in unit.iteritems():
        if val >= 0:
            has_pos = True
            break

    for key, val in sorted(unit.iteritems()):
        if val < 0 and has_pos:
            continue
        result += ' ' + key
        if val != 1:
            ival = int(val)
            if abs(val-ival) <= 1e-7 * abs(val):
                val = ival
            result += '^' + str(val)

    if not has_pos:
        return result[1:]

    for key, val in sorted(unit.iteritems()):
        if val >= 0:
            continue
        val = -val
        result += ' / ' + key
        if val != 1:
            ival = int(val)
            if abs(val-ival) <= 1e-7 * abs(val):
                val = ival
            result += '^' + str(val)

    return result[1:]


class Quantity(numpy.ndarray):
    """
    Represent a quantity, i.e. a scalar or array associated with a unit.

    Examples
    --------

    A Quantity can be converted to a different unit:
    >>> a = Quantity(18, 'm')
    >>> a.unit = 'km'
    >>> a
    Quantity(0.018, 'km')

    A more useful conversion:
    >>> sky = Quantity(4*numpy.pi, 'sr')
    >>> print sky.inunit('deg^2')
    41252.9612494 deg^2
    
    Quantities can be compared:
    >>> Quantity(0.018, 'km') > Quantity(10., 'm')
    True
    >>> minimum(Quantity(1, 'm'), Quantity(0.1, 'km'))
    Quantity(1.0, 'm')

    Quantities can be operated on:
    >>> time = Quantity(0.2, 's')
    >>> a / time
    Quantity(0.09, 'km / s')

    Units do not have to be standard and ? can be used as a 
    non-standard one:
    >>> value = Quantity(1., '?/detector')
    >>> value *= Quantity(100, 'detector')
    >>> value
    Quantity(100.0, '?')

    A derived unit can be converted to a SI one:
    >>> print Quantity(2, 'Jy').SI
    2e-26 kg / s^2

    It is also possible to add a new derived unit:
    >>> DerivedUnits.new('kloug', Quantity(3, 'kg'))
    >>> print Quantity(2, 'kloug').SI
    6.0 kg

    In fact, a non-standard unit is not handled differently than one of
    the 7 SI units, as long as it is not in the DerivedUnits table:
    >>> DerivedUnits.new('krouf', Quantity(0.5, 'broug^2'))
    >>> print Quantity(1, 'krouf').SI
    0.5 broug^2
    """

    __slots__ = ('_unit',)

    def __new__(cls, data, unit=None, dtype=numpy.float64, copy=True, order='C', subok=False, ndmin=0):

        # get a new Quantity instance (or a subclass if subok is True)
        result = numpy.array(data, dtype, copy=copy, order=order, subok=True, ndmin=ndmin)
        if not subok and result.__class__ is not cls or not issubclass(result.__class__, cls):
            result = result.view(cls)

        # copy unit attribute
        if unit is not None:
            result._unit = _extract_unit(unit)
            
        return result

    def mean(self, axis=None, dtype=None, out=None):
        result = super(self.__class__, self).mean(axis, dtype, out)
        if numpy.rank(result) == 0:
            result = Quantity(result, self._unit)
        return result
    mean.__doc__ = numpy.ndarray.mean.__doc__

    def sum(self, axis=None, dtype=None, out=None):
        result = super(self.__class__, self).sum(axis, dtype, out)
        if numpy.rank(result) == 0:
            result = Quantity(result, self._unit)
        return result
    sum.__doc__ = numpy.ndarray.sum.__doc__

    def std(self, axis=None, dtype=None, out=None, ddof=0):
        result = super(self.__class__, self).std(axis, dtype, out, ddof)
        if numpy.rank(result) == 0:
            result._unit = self._unit
        return result
    std.__doc__ = numpy.ndarray.std.__doc__

    def var(self, axis=None, dtype=None, out=None, ddof=0):
        result = super(self.__class__, self).var(axis, dtype, out, ddof)
        if numpy.rank(result) == 0:
            result._unit = _power_unit(self._unit, 2)
        return result
    var.__doc__ = numpy.ndarray.var.__doc__

    def __array_finalize__(self, obj):
        if obj is None: return
        self._unit = getattr(obj, '_unit', None)

    @property
    def __array_priority__(self):
        return 1. if self._unit is None else 1.5

    def __array_prepare__(self, array, context=None):
        """
        Homogenise ufunc's argument units and cast array to the class of the
        argument of highest __array_priority__
        """

        if self is not array:
            array = array.view(type(self))
        if context is None or self._unit is None:
            return array
       
        ufunc = context[0]
        if ufunc in (numpy.add, numpy.subtract, numpy.maximum, numpy.minimum,
                     numpy.greater, numpy.greater_equal, numpy.less,
                     numpy.less_equal, numpy.equal, numpy.not_equal):
            for arg in context[1]:
                u = getattr(arg, '_unit', None)
                if u is None or u == self._unit:
                    continue

                print "Warning: applying function '" + str(ufunc) + "' to Quant\
ities of different units may have changed operands to common unit '" + \
                    _strunit(self._unit) + "'."
                arg.unit = self._unit

        return array

    def __array_wrap__(self, array, context=None):
        """
        Set unit of the result of a ufunc.

        Since different quantities can share their _unit attribute, a
        a change in unit of the ufunc result must be preceded by a copy
        of the argument unit.
        Not all ufuncs are currently handled. For an exhaustive list of
        ufuncs, see http://docs.scipy.org/doc/numpy/reference/ufuncs.html
        """

        if context is None:
            return array

        ufunc = context[0]
        args = context[1]

        if ufunc in (numpy.add, numpy.subtract, numpy.maximum, numpy.minimum):
            if self is not array:
                # self has highest __array_priority__
                array._unit = self._unit
            else:
                # inplace operation
                if array._unit is None:
                    array._unit = getattr(args[1], '_unit', None)

        elif ufunc is numpy.reciprocal:
            array._unit = _power_unit(args[0]._unit, -1)

        elif ufunc is numpy.sqrt:
            array._unit = _power_unit(args[0]._unit, 0.5)

        elif ufunc in [numpy.square, numpy.var]:
            array._unit = _power_unit(args[0]._unit, 2)

        elif ufunc is numpy.power:
            array._unit = _power_unit(args[0]._unit, args[1])

        elif ufunc in [numpy.multiply, numpy.vdot]:
            units = _get_units(args)
            array._unit = _multiply_unit(units[0], units[1])

        elif ufunc is numpy.divide:
            units = _get_units(args)
            array._unit = _divide_unit(units[0], units[1])

        elif ufunc in (numpy.greater, numpy.greater_equal, numpy.less,
                       numpy.less_equal, numpy.equal, numpy.not_equal,
                       numpy.iscomplex, numpy.isfinite, numpy.isinf, 
                       numpy.isnan, numpy.isreal):
            if numpy.rank(array) == 0:
                return bool(array)
            array = array.view(numpy.ndarray)

        elif ufunc in (numpy.arccos, numpy.arccosh, numpy.arcsin, numpy.arcsinh,
                       numpy.arctan, numpy.arctanh, numpy.arctan2, numpy.cosh,
                       numpy.sinh, numpy.tanh, numpy.cos, numpy.sin, numpy.tan,
                       numpy.log, numpy.log2, numpy.log10, numpy.exp, 
                       numpy.exp2):
            array._unit = None

        elif ufunc not in (numpy.abs, numpy.negative ):
            print 'Quantity: unhandled ufunc ', ufunc, 'with', len(args), 'args'
            array._unit = None

        else:
            array._unit = self._unit
            
        return array

    def __getitem__(self, key):
        item = super(Quantity, self).__getitem__(key)
        if isinstance(item, Quantity):
            return item
        return Quantity(item, unit=self._unit, copy=False)

    @property
    def unit(self):
        """
        Return the Quantity unit as a string
        
        Example
        -------
        >>> Quantity(32., 'm/s').unit
        'm / s'
        """
        return _strunit(self._unit)

    @unit.setter
    def unit(self, unit):
        """
        Convert a Quantity into a new unit

        Parameters
        ----------
        unit : string or None
             A string representing the unit into which the Quantity
             should be converted

        Example
        -------
        >>> q = Quantity(1., 'km')
        >>> q.unit = 'm'
        >>> print q
        1000.0 m
        """
        newunit = _extract_unit(unit)
        if self._unit is None or newunit is None:
            self._unit = newunit
            return

        q1 = Quantity(1., self._unit).SI
        q2 = Quantity(1., newunit).SI
        if q1._unit != q2._unit:
            raise UnitError("Units '"+self.unit+"' and '" + _strunit(newunit)+"' are incompatible.")
        numpy.ndarray.__imul__(self.view(numpy.ndarray), numpy.asscalar(q1) / numpy.asscalar(q2))
        self._unit = newunit

    def inunit(self, unit):
        """
        Return a copy of the Quantity in the specified unit

        Parameters
        ----------
        unit : string or None
             A string representing the unit into which the Quantity
             should be converted

        Returns
        -------
        res : Quantity
            A copy of the Quantity in the new unit

        Example
        -------
        >>> print Quantity(1., 'km').inunit('m')
        1000.0 m
        """
        q = Quantity(self, subok=True)
        q.unit = unit
        return q

    @property
    def SI(self):
        """
        Return a copy of the Quantity in SI unit

        Example
        -------
        >>> print Quantity(1., 'km').SI
        1000.0 m
        """
        return _normalize_quantity(self)

    def __repr__(self):
        return self.__class__.__name__+'(' + str(numpy.array(self,copy=False)) + ", '" + _strunit(self._unit) + "')"

    def __str__(self):
        result = super(Quantity,self).__str__()
        if self._unit is None:
            return result
        return result + ' ' + _strunit(self._unit)


class DerivedUnits(object):
    """
    This class holds a table of derived units.

    This table is used to compared Quantitiies of different units and to
    convert a Quantity into the International System.
    """

    table = {
        'arcmin' : Quantity(1./60., 'deg'), 
        'arcsec' : Quantity(1./3600., 'deg'),
        'AU'     : Quantity(149597870700., 'm'),
        'C'      : Quantity(1., 'A s'),
        'deg'    : Quantity(numpy.pi/180., 'rad'),
        'Hz'     : Quantity(1., 's^-1'),
        'J'      : Quantity(1., 'kg m^2 / s^2'),
        'N'      : Quantity(1., 'kg m / s^2'),
        'Ohm'    : Quantity(1., 'V/A'),
        'pc'     : Quantity(30.857e15, 'm'),
        ('rad',2): Quantity(1., 'sr'),
        'V'      : Quantity(1., 'kg m^2 / A / s^3'),
        'W'      : Quantity(1., 'kg m^2 / s^3'),
        'km'     : Quantity(1000., 'm'),
        'uJy'    : Quantity(1.e-6, 'Jy'),
        'mJy'    : Quantity(1.e-3, 'Jy'),
        'Jy'     : Quantity(1.e-26, 'W/Hz/m^2'),
        'MJy'    : Quantity(1.e6, 'Jy'),
        }

    @staticmethod
    def new(name, value):
        """
        Adds a new derived unit to the current list
        
        Parameters
        ----------
        name : string
             The name of the new derived unit to be added to the table

        value : Quantity
             The Quantity corresponding to the new derived unit

        Example
        -------
        >>> DerivedUnits.new('mile', Quantity(1.609344, 'km'))
        >>> Quantity(1, 'm').inunit('mile')
        Quantity(0.000621371192237, 'mile')
        """
        if not isinstance(value, Quantity):
            raise TypeError('2nd argument should be a Quantity.')
        DerivedUnits.table[name] = value.SI

    @staticmethod
    def delete(name):
        """
        Delete a derived unit from the table.

        Parameters
        ----------
        name : string
             The name of the unit to be deleted

        Example
        -------
        >>> DerivedUnits.delete('mynewunit')
        """
        del DerivedUnits.table[name]

def _grab_doc(doc, func):
    doc = doc.replace(func + '(shape, ', func + '(shape, unit=None, ')
    doc = doc.replace('\n    dtype : ','\n    unit : string\n        Unit of the new array, e.g. ``W``. Default is None for scalar.\n    dtype : ''')
    doc = doc.replace('np.'+func, 'unit.'+func)
    doc = doc.replace('\n    out : ndarray', '\n    out : Quantity')
    return doc

def empty(shape, unit=None, dtype=None, order=None):
    return Quantity(numpy.empty(shape, dtype, order), unit, copy=False)
empty.__doc__ = _grab_doc(numpy.empty.__doc__, 'empty')

def ones(shape, unit=None, dtype=None, order=None):
    return Quantity(numpy.ones(shape, dtype, order), unit, copy=False)
empty.__doc__ = _grab_doc(numpy.ones.__doc__, 'ones')

def zeros(shape, unit=None, dtype=None, order=None):
    return Quantity(numpy.zeros(shape, dtype, order), unit, copy=False)
empty.__doc__ = _grab_doc(numpy.zeros.__doc__, 'zeros')
