import numpy

from tamasis.core import *
from tamasis.numpyutils import _my_isscalar

__all__ = ['FlatField', 'MaskPolicy', 'Pointing']

class FlatField(object):
    def __init__(self, optical, detector):
        self.optical = Map(optical, origin='upper')
        self.detector = Map(detector, origin='upper')

    @property
    def total(self):
        return Map(self.optical * self.detector, origin='upper')


#-------------------------------------------------------------------------------


class MaskPolicy(object):
    KEEP   = 0
    MASK   = 1
    REMOVE = 2
    def __init__(self, flags, values, description=None):
        self._description = description
        if _my_isscalar(flags):
            if isinstance(flags, str):
                flags = flags.split(',')
            else:
                flags = (flags,)
        if _my_isscalar(values):
            if isinstance(values, str):
                values = values.split(',')
            else:
                values = (values,)
        if len(flags) != len(values):
            raise ValueError('The number of policy flags is different from the number of policies.')

        self._policy = []
        for flag, value in zip(flags, values):
            if flag[0] == '_':
                raise ValueError('A policy flag should not start with an underscore.')
            value = value.strip().lower()
            if value not in ('keep', 'mask', 'remove'):
                raise KeyError("Invalid policy "+flag+"='" + value + "'. Valid policies are 'keep', 'mask' or 'remove'.")
            self._policy.append({flag:value})
            setattr(self, flag, value)
        self._policy = tuple(self._policy)

    def __array__(self, dtype=int):
        conversion = { 'keep':self.KEEP, 'mask':self.MASK, 'remove':self.REMOVE }
        return numpy.array([conversion[policy.values()[0]] for policy in self._policy], dtype=dtype)

    def __str__(self):
        str = self._description + ': ' if self._description is not None else ''
        str_ = []
        for policy in self._policy:
            str_.append(policy.values()[0] + " '" + policy.keys()[0] + "'")
        str += ', '.join(str_)
        return str


#-------------------------------------------------------------------------------


class Pointing(numpy.ndarray):
    INSCAN     = 1
    TURNAROUND = 2
    OTHER      = 3
    def __new__(cls, time, ra, dec, pa, info=None, policy=None, nsamples=None):
        if nsamples is None:
            nsamples = (time.size,)
        else:
            if _my_isscalar(nsamples):
                nsamples = (nsamples,)
            else:
                nsamples = tuple(nsamples)
        nsamples_tot = numpy.sum(nsamples)

        if info is None:
            info = Pointing.INSCAN

        if policy is None:
            policy = MaskPolicy.KEEP

        time   = numpy.asarray(time)
        ra     = numpy.asarray(ra)
        dec    = numpy.asarray(dec)
        pa     = numpy.asarray(pa)
        info   = numpy.asarray(info)
        policy = numpy.asarray(policy)
        if numpy.any(map(lambda x: x.size not in (1, nsamples_tot), [ra,dec,pa,info,policy])):
            raise ValueError('The pointing inputs do not have the same size.')

        ftype = config.get_default_dtype_float()
        dtype = numpy.dtype([('time', ftype), ('ra', ftype), ('dec', ftype), ('pa', ftype), ('info', numpy.int64), ('policy', numpy.int64)])
        result = numpy.recarray(nsamples_tot, dtype=dtype)
        result.time = time
        result.ra   = ra
        result.dec  = dec
        result.pa   = pa
        result.info = info
        result.policy = policy
        result.nsamples = nsamples
        return result

    @property
    def velocity(self):
        dra  = Quantity(numpy.diff(self.ra), 'deg')
        ddec = Quantity(numpy.diff(self.dec), 'deg')
        dtime = Quantity(numpy.diff(self.time), 's')
        vel = numpy.sqrt((dra*numpy.cos(self.dec[0:-1].inunit('rad')))**2 + ddec**2) / dtime
        vel.unit = 'arcsec/s'
        u = vel._unit
        vel = numpy.append(vel, vel[-1])
        # BUG: append eats the unit...
        vel._unit = u
        return vel


#-------------------------------------------------------------------------------


