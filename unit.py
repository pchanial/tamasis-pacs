import numpy
import re
re_unit = re.compile(r' *([/*])? *([a-zA-Z]+|\?)(\^-?[0-9]+(\.[0-9]*)?)? *')

class UnitError(Exception): pass

def extract_unit(string):

    if string is None:
        return {}

    if not isinstance(string, str) and not isinstance(string, dict):
        raise TypeError('Invalid unit type. Valid types are string or a dictionary.')

    if isinstance(string, dict):
        return string

    string = string.strip()
    result = {}
    start = 0
    while start < len(string):
        match = re_unit.match(string, start)
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
        multiply_unit(result, {u:exp})
        start = start + len(match.group(0))
    return result

def normalize_quantity(value):
    if len(value._unit) == 0:
        return value
    result = Quantity(value, '')
    for key, val in value._unit.iteritems():
        if (key,val) in DerivedUnit._units:
            result *= normalize_quantity(DerivedUnit._units[(key,val)])
        elif (key,-val) in DerivedUnit._units:
            result /= normalize_quantity(DerivedUnit._units[(key,-val)])
        elif key in DerivedUnit._units:
            if val == 1.:
                result *= normalize_quantity(DerivedUnit._units[key])
            elif val == -1.:
                result /= normalize_quantity(DerivedUnit._units[key])
            else:
                result *= normalize_quantity(DerivedUnit._units[key]**val)
        else:
            if key in result._unit:
                result._unit[key] += val
            else:
                result._unit[key] = val
    return result

def add_unit(unit1, unit2):
    if len(unit1) == 0:
        return unit2
    if len(unit2) == 0:
        return unit1
    if unit1 == unit2:
        return unit1
    

    
def power_unit(unit, power):
    if len(unit) == 0 or power == 0:
        return {}
    result = unit.copy()
    for key in unit:
        result[key] *= power
    return result

def multiply_unit(unit1, unit2):
    if len(unit2) == 0:
        return
    for key, val in unit2.iteritems():
        if unit1.has_key(key):
            if unit1[key] == -val:
                del unit1[key]
            else:
                unit1[key] += val
        else:
            unit1[key] = val

def divide_unit(unit1, unit2):
    if len(unit2) == 0:
        return
    for key, val in unit2.iteritems():
        if unit1.has_key(key):
            if unit1[key] == val:
                del unit1[key]
            else:
                unit1[key] -= val
        else:
            unit1[key] = -val

        
class Quantity(numpy.ndarray):

    def __new__(cls, data, unit=None, dtype=numpy.float64, order='C', copy=True):
        from copy import copy as cp
        result = numpy.array(data, dtype=dtype, order=order, copy=copy, subok=True)
        if unit is not None:
            unit = extract_unit(unit)
        if not isinstance(result, Quantity):
            result = result.view(cls)
            if unit is not None:
                result._unit = unit
        else:
            result._unit = data._unit.copy() if copy else data._unit
            if unit is not None:
                result.unit = unit
        return result

    def __array_finalize__(self, obj):
        if obj is None: return
        self._unit = getattr(obj, '_unit', {})

    def __array_wrap__(self, obj, context=None):
        result = numpy.ndarray.__array_wrap__(self, obj, context).view(type(self))

        if context is None:
            return result
        
        ufunc = context[0]
        if ufunc in [numpy.add, numpy.subtract]:
            pass
        elif ufunc is numpy.reciprocal:
            result._unit = power_unit(context[1][0]._unit, -1)
        elif ufunc is numpy.sqrt:
            result._unit = power_unit(context[1][0]._unit, 0.5)
        elif ufunc in [numpy.square, numpy.var]:
            result._unit = power_unit(context[1][0]._unit, 2)
        elif ufunc is numpy.power:
            result._unit = power_unit(context[1][0]._unit, context[1][1])
        elif ufunc in [numpy.multiply, numpy.vdot]:
            result._unit = getattr(context[1][0], '_unit', {}).copy()
            multiply_unit(result._unit, getattr(context[1][1], '_unit', {}))
        elif ufunc is numpy.divide:
            result._unit = getattr(context[1][0], '_unit', {}).copy()
            divide_unit(result._unit, getattr(context[1][1], '_unit', {}))
        elif ufunc in (numpy.arccos, numpy.arccosh, numpy.arcsin, numpy.arcsinh, numpy.arctan, numpy.arctanh, numpy.arctan2, numpy.cosh, numpy.sinh, numpy.tanh, numpy.cos, numpy.sin, numpy.tan, numpy.log, numpy.log2, numpy.log10, numpy.exp, numpy.exp2):
            result._unit = {}
        else:
            result._unit = context[1][0]._unit
            
        return result

    def __add__(self, other):
        if not isinstance(other, Quantity) or len(other._unit) == 0 or self._unit == other._unit:
            return super(Quantity, self).__add__(other)

        if len(self._unit) == 0:
            result = super(Quantity, self).__add__(other)
            result._unit = other._unit
            return result
        
        q1 = Quantity(1., self._unit).unit_si
        q2 = Quantity(1., other._unit).unit_si
        if q1._unit != q2._unit:
            raise UnitError("Units '"+self.unit+"' and '"+other.unit+"' are incompatible.")

        return super(Quantity, self).__add__(other.view(numpy.ndarray) * (numpy.asscalar(q2) / numpy.asscalar(q1)))

    def __radd__(self, other):
        return super(Quantity, self).__add__(other)

    def __sub__(self, other):
        if not isinstance(other, Quantity) or len(other._unit) == 0 or self._unit == other._unit:
            return super(Quantity, self).__sub__(other)

        if len(self._unit) == 0:
            result = super(Quantity, self).__sub__(other)
            result._unit = other._unit
            return result
        
        q1 = Quantity(1., self._unit).unit_si
        q2 = Quantity(1., other._unit).unit_si
        if q1._unit != q2._unit:
            raise UnitError("Units '"+str(self.unit)+"' and '"+str(other.unit)+"' are incompatible.")

        return super(Quantity, self).__sub__(other.view(numpy.ndarray) * (numpy.asscalar(q2) / numpy.asscalar(q1)))

    def __rsub__(self, other):
        return -super(Quantity, self).__sub__(other)

    def __getitem__(self, key):
        item = super(Quantity, self).__getitem__(key)
        if isinstance(item, Quantity):
            return item
        return Quantity(item, unit=self._unit, copy=False)

    @property
    def unit(self):
        return self._strunit(self._unit)

    @unit.setter
    def unit(self, unit):
        
        newunit = extract_unit(unit)
        if len(self._unit) == 0 or len(newunit) == 0:
            self._unit = newunit
            return

        q1 = Quantity(1., self._unit).unit_si
        q2 = Quantity(1., newunit).unit_si
        if q1._unit != q2._unit:
            raise UnitError("Units '"+self.unit+"' and '"+self._strunit(newunit)+"' are incompatible.")
        numpy.ndarray.__imul__(self.view(numpy.ndarray), numpy.asscalar(q1) / numpy.asscalar(q2))
        self._unit = newunit

    @property
    def unit_si(self):
        return normalize_quantity(Quantity(1., self._unit))

    def __repr__(self):
        return self.__class__.__name__+'(' + str(numpy.array(self,copy=False)) + ", '" + self._strunit(self._unit) + "')"

    def __str__(self):
        result = super(Quantity,self).__str__()
        if len(self._unit) == 0:
            return result
        return result + ' ' + self._strunit(self._unit)

    @staticmethod
    def _strunit(unit):
        if len(unit) == 0:
            return ''
        result = ''
        has_pos = False

        if unit.has_key('?'):
            result +=' ?'

        for key, val in unit.iteritems():
            if val >= 0 or key == '?':
                has_pos = True
                break

        for key, val in sorted(unit.iteritems()):
            if val < 0 and has_pos or key == '?':
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
            if val >= 0 or key == '?':
                continue
            val = -val
            result += ' / ' + key
            if val != 1:
                ival = int(val)
                if abs(val-ival) <= 1e-7 * abs(val):
                    val = ival
                result += '^' + str(val)

        return result[1:]


class DerivedUnit(object):

    _units = {
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
    def new(name, value, unit=None):
        if unit is None and not isinstance(value, Quantity):
            raise ValueError('')
        if unit is None:
            DerivedUnit._units[name] = value
        else:
            DerivedUnit._units[name] = Quantity(value, unit)

    @staticmethod
    def delete(name):
        del DerivedUnit._units[name]
