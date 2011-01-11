import numpy
import re
from . import var

__all__ = ['Quantity', 'UnitError', 'units']

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
    if key in unit:
        if unit[key] == -val:
            del unit[key]
        else:
            unit[key] += val
    else:
        unit[key] = val
    return unit

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
    for key, val in unit2.items():
        if key in unit:
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
    for key, val in unit2.items():
        if key in unit:
            if unit[key] == val:
                del unit[key]
            else:
                unit[key] -= val
        else:
            unit[key] = -val
    return unit

def _get_units(units):
    return [getattr(q, '_unit', None) for q in units]

def _strunit(unit):
    """
    Convert a unit as dictionary into a string
    """
    if unit is None:
        return ''
    result = ''
    has_pos = False

    for key, val in unit.items():
        if val >= 0:
            has_pos = True
            break

    for key, val in sorted(unit.items()):
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

    for key, val in sorted(unit.items()):
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
    if dtype is not specified, the quantity is upcasted to float64 if its dtype
    is not real nor complex

    Examples
    --------

    A Quantity can be converted to a different unit:
    >>> a = Quantity(18, 'm')
    >>> a.unit = 'km'
    >>> a
    Quantity(0.018, 'km')

    A more useful conversion:
    >>> sky = Quantity(4*numpy.pi, 'sr')
    >>> print(sky.tounit('deg^2'))
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
    >>> print(Quantity(2, 'Jy').SI)
    2e-26 kg / s^2

    It is also possible to add a new derived unit:
    >>> unit['kloug'] = Quantity(3, 'kg')
    >>> print(Quantity(2, 'kloug').SI)
    6.0 kg

    In fact, a non-standard unit is not handled differently than one of
    the 7 SI units, as long as it is not in the unit table:
    >>> unit['krouf'] = Quantity(0.5, 'broug^2')
    >>> print(Quantity(1, 'krouf').SI)
    0.5 broug^2
    """

    __slots__ = ('_unit','derived_units')

    def __new__(cls, data, unit=None, derived_units=None, dtype=None, copy=True, order='C', subok=False, ndmin=0):

        if dtype is None:
            dtype = var.get_default_dtype(data)

        # get a new Quantity instance (or a subclass if subok is True)
        result = numpy.array(data, dtype, copy=copy, order=order, subok=True, ndmin=ndmin)
        if not subok and type(result) is not cls or not isinstance(result, cls):
            result = result.view(cls)

        # set derived_units attribute
        if derived_units is not None:
            for key in derived_units.keys():
                if not isinstance(derived_units[key], Quantity) and not type(derived_units[key]) is tuple:
                    raise UnitError("The user derived unit '" + key + "' is not a Quantity.")
            result.derived_units.update(derived_units)

        # copy unit attribute
        if unit is not None:            
            result._unit = _extract_unit(unit)

        return result

    def min(self, axis=None, out=None):
        return _wrap_func(numpy.min, self, self._unit, axis, out)
    min.__doc__ = numpy.ndarray.min.__doc__

    def max(self, axis=None, out=None):
        return _wrap_func(numpy.max, self, self._unit, axis, out)
    max.__doc__ = numpy.ndarray.max.__doc__

    def sum(self, axis=None, dtype=None, out=None):
        return _wrap_func(numpy.sum, self, self._unit, axis, dtype, out)
    sum.__doc__ = numpy.ndarray.sum.__doc__

    def mean(self, axis=None, dtype=None, out=None):
        return _wrap_func(numpy.mean, self, self._unit, axis, dtype, out)
    mean.__doc__ = numpy.ndarray.mean.__doc__

    def std(self, axis=None, dtype=None, out=None, ddof=0):
        return _wrap_func(numpy.std, self, self._unit, axis, dtype, out)
    std.__doc__ = numpy.ndarray.std.__doc__

    def var(self, axis=None, dtype=None, out=None, ddof=0):
        return _wrap_func(numpy.var, self, _power_unit(self._unit,2), axis, dtype, out)
    var.__doc__ = numpy.ndarray.var.__doc__

    def __array_finalize__(self, obj):
        # for some numpy methods (append): the result doesn't go through __new__ and
        # obj is None in __array_finalize__
        #if obj is None: return
        self._unit = getattr(obj, '_unit', None)
        self.derived_units = getattr(obj, 'derived_units', {})

    @property
    def __array_priority__(self):
        return 1. if self._unit is None else 1.5

    def __array_prepare__(self, array, context=None):
        """
        Homogenise ufunc's argument units and cast array to the class of the
        argument of highest __array_priority__
        """

        if not isinstance(array, type(self)):
            array = array.view(type(self))

        #XXX derived_units should be obtained from all arguments
        array.derived_units = self.derived_units
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

                print("Warning: applying function '" + str(ufunc) + "' to Quant\
ities of different units may have changed operands to common unit '" + \
                    _strunit(self._unit) + "'.")
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

        if numpy.__version__ < '1.4.1' and self is not array:
            array = array.view(type(self))

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
            array = array.magnitude

        elif ufunc in (numpy.arccos, numpy.arccosh, numpy.arcsin, numpy.arcsinh,
                       numpy.arctan, numpy.arctanh, numpy.arctan2, numpy.cosh,
                       numpy.sinh, numpy.tanh, numpy.cos, numpy.sin, numpy.tan,
                       numpy.log, numpy.log2, numpy.log10, numpy.exp, 
                       numpy.exp2):
            array._unit = None

        elif ufunc in (numpy.bitwise_or, numpy.invert, numpy.logical_and,
                       numpy.logical_not, numpy.logical_or, numpy.logical_xor):
            array._unit = None

        elif ufunc not in (numpy.abs, numpy.negative ):
            print('Quantity: unhandled ufunc ', ufunc, 'with', len(args), 'args')
            array._unit = None

        else:
            array._unit = self._unit
            
        return array

    def __getitem__(self, key):
        item = super(Quantity, self).__getitem__(key)
        if isinstance(item, Quantity):
            return item
        return Quantity(item, self._unit, self.derived_units, copy=False)

    @property
    def magnitude(self):
        """
        Return the magnitude of the quantity
        """
        return self.view(numpy.ndarray)

    @magnitude.setter
    def magnitude(self, value):
        """
        Set the magnitude of the quantity
        """
        if numpy.rank(self) == 0:
            self.shape = (1,)
            try:
                self.view(numpy.ndarray)[:] = value
            except ValueError as error:
                raise error
            finally:
                self.shape = ()
        else:
            self.view(numpy.ndarray)[:] = value

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
        >>> print(q)
        1000.0 m
        """
        newunit = _extract_unit(unit)
        if self._unit is None or newunit is None:
            self._unit = newunit
            return

        q1 = Quantity(1., self._unit, self.derived_units).SI
        q2 = Quantity(1., newunit, self.derived_units).SI
        if q1._unit != q2._unit:
            raise UnitError("Units '" + self.unit + "' and '" + \
                            _strunit(newunit) + "' are incompatible.")
        factor = q1.magnitude / q2.magnitude
        if numpy.rank(self) == 0:
            self.magnitude = self.magnitude * factor
        else:
            self.magnitude.T[:] *= factor.T
        self._unit = newunit

    def tounit(self, unit):
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
        >>> print(Quantity(1., 'km').tounit('m'))
        1000.0 m
        """
        q = self.copy()
        q.unit = unit
        return q

    @property
    def SI(self):
        """
        Return a copy of the Quantity in SI unit

        Example
        -------
        >>> print(Quantity(1., 'km').SI)
        1000.0 m
        """
        if self._unit is None:
            return self
        result = self.copy()
        result.unit = ''
        factor = Quantity(1., '')
        for key, val in self._unit.items():

            # check if the unit is a local derived unit
            newfactor = _check_du(self, key, val, self.derived_units)

            # check if the unit is a global derived unit
            if newfactor is None:
                newfactor = _check_du(self, key, val, units_derived)
                
            if newfactor is None:
                # the unit is not derived, check if the unit is SI
                if key in units_SI:
                    if factor._unit is None:
                        factor._unit = { key : val }
                    elif key in factor._unit:
                        factor._unit[key] += val
                        if factor._unit[key] == 0:
                            del factor._unit[key]
                    else:
                        factor._unit[key] = val
                    continue
            
                # the unit is not derived and not SI, we add it to the dictionary
                if factor._unit is None:
                    factor._unit = { key : val }
                elif key in factor._unit:
                    factor._unit[key] += val
                    if factor._unit[key] == 0:
                        del factor._unit[key]
                else:
                    factor._unit[key] = val
                continue

            # factor may be broadcast
            factor = factor * newfactor

            # if the unit factor is scalar, nothing more to do
            if newfactor.size == 1:
                continue
            
            # broadcast the result if it is a scalar
            if numpy.rank(self) == 0:
                newfactor[:] = 1.
                result = result * newfactor
                continue

            # make sure dimensions are compatible
            if self.shape[0] != 1 and self.shape[0] != newfactor.size:
                raise ValueError("The derived unit '" + key + "' has a number of elements '" + str(newfactor.size) + "' which is incompatible with that along the first axis '" + str(shape[0]) + "'.")

            # broadcast the result if necessary
            if self.shape[0] == 1:
                shape = numpy.ones(numpy.rank(self), dtype=int)
                shape[0] = newfactor.size
                newfactor = newfactor.reshape(shape)
                newfactor[:] = 1.
                result = result * newfactor

        if numpy.rank(result) == 0:
            result = result * factor
        else:
            result.T[:] *= factor.magnitude.T
        result._unit = factor._unit
        return result

    def __repr__(self):
        return type(self).__name__+'(' + str(numpy.asarray(self)) + ", '" + _strunit(self._unit) + "')"

    def __str__(self):
        result = str(numpy.asarray(self))
        if self._unit is None:
            return result
        return result + ' ' + _strunit(self._unit)


def _get_du(input, key, d):
    if type(d[key]) is tuple:
        return d[key][0](input) * Quantity(d[key][1], derived_units=d, copy=False)
    return Quantity(d[key], derived_units=d, copy=False)

def _check_du(input, key, val, derived_units):
    if len(derived_units) == 0:
        return None
    if (key,val) in derived_units:
        return _get_du(input, (key,val), derived_units).SI
    if (key,-val) in derived_units:
        return (1. / _get_du(input, (key,-val), derived_units)).SI
    if key not in derived_units:
        return None
    if val == 1.:
        return _get_du(input, key, derived_units).SI
    if val == -1.:
        return (1. / _get_du(input, key, derived_units)).SI
    return (_get_du(input, key, derived_units)**val).SI

units_SI = {
    'A'   : Quantity(1., 'A'),
    'cd'  : Quantity(1., 'cd'),
    'K'   : Quantity(1., 'K'),
    'kg'  : Quantity(1., 'kg'),
    'm'   : Quantity(1., 'm'),
    'mol' : Quantity(1., 'mol'),
    'rad' : Quantity(1., 'rad'),
    's'   : Quantity(1., 's'),
    'sr'  : Quantity(1., 'sr'),
}

units_derived = {
    
    # angle
    'arcmin' : Quantity(1./60., 'deg'), 
    'arcsec' : Quantity(1./3600., 'deg'),
    'deg'    : Quantity(numpy.pi/180., 'rad'),

    # flux_densities
    'uJy'    : Quantity(1.e-6, 'Jy'),
    'mJy'    : Quantity(1.e-3, 'Jy'),
    'Jy'     : Quantity(1.e-26, 'W/Hz/m^2'),
    'MJy'    : Quantity(1.e6, 'Jy'),

    # force
    'N'      : Quantity(1., 'kg m / s^2'),
    
    # frequency
    'Hz'     : Quantity(1., 's^-1'),

    # energy
    'J'      : Quantity(1., 'kg m^2 / s^2'),

    # length
    'AU'     : Quantity(149597870700., 'm'),
    'km'     : Quantity(1000., 'm'),
    'pc'     : Quantity(30.857e15, 'm'),

    # power
    'W'      : Quantity(1., 'kg m^2 / s^3'),

    # pressure
    'atm'    : Quantity(101325., 'Pa'),
    'bar'    : Quantity(1e5, 'Pa'),
    'mmHg'   : Quantity(1./760, 'atm'),
    'cmHg'   : Quantity(10, 'mmHg'),
    'Pa'     : Quantity(1., 'kg / m / s^2'),
    'Torr'   : Quantity(1., 'mmHg'),
    
    # resistivity
    'Ohm'    : Quantity(1., 'V/A'),
    
    # solid angle
    ('rad',2): Quantity(1., 'sr'),
    
    # misc
    'C'      : Quantity(1., 'A s'),
    'V'      : Quantity(1., 'kg m^2 / A / s^3'),
}

class Unit(dict):
    def __init__(self):
        
        for k, v in units_SI.items():
            self[k] = v
            setattr(self, k, v)
        for k, v in units_derived.items():
            self[k] = v
            if isinstance(k, str):
                setattr(self, k, v)

units = Unit()

def _wrap_func(func, array, unit, *args):
    result = func(array.magnitude, *args)
    result = Quantity(result, unit, array.derived_units, copy=False)
    return result

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
ones.__doc__ = _grab_doc(numpy.ones.__doc__, 'ones')

def zeros(shape, unit=None, dtype=None, order=None):
    return Quantity(numpy.zeros(shape, dtype, order), unit, copy=False)
zeros.__doc__ = _grab_doc(numpy.zeros.__doc__, 'zeros')
