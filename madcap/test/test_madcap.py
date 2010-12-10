import numpy
import pyfits
import os
import tamasis

from tamasis import *

tamasis.__verbose__ = False
path = os.path.abspath(os.path.dirname(__file__)) + '/data/madmap1/'
obs = MadMap1Observation(path+'todSpirePsw_be', path+'invnttSpirePsw_be', 
                         path+'madmapSpirePsw.fits[coverage]', 'big_endian',
                         135, missing_value=numpy.nan)
obs.instrument.name = 'SPIRE/PSW'

tod = obs.get_tod()
tod.unit = 'Jy/beam'
invNtt = InvNtt(len(tod.nsamples)*(1024,), obs.get_filter_uncorrelated())
fft = FftHalfComplex(len(tod.nsamples)*(1024,))
padding = Padding(left=invNtt.ncorrelations, right=1024-numpy.array(tod.nsamples)-invNtt.ncorrelations)
projection = Projection(obs)
packing = Unpacking(obs.info.mapmask, field=numpy.nan).T

map_naive = mapper_naive(tod, projection*packing)
map_ref = pyfits.fitsopen(path+'naivemapSpirePsw.fits')['image'].data
if any_neq(map_naive,map_ref): print('FAILED: mapper_naive madcap 1')

map_naive_1d = mapper_naive(tod, projection)
map_naive_2d = packing.transpose(map_naive_1d)
if any_neq(map_naive_2d,map_ref): print('FAILED: mapper_naive madcap 2')

packing = Unpacking(obs.info.mapmask).T

#M = 1 / map_naive.coverage
#M[~numpy.isfinite(M)] = numpy.max(M[numpy.isfinite(M)])
#map_rlsw1 = mapper_rls(tod, projection*packing, padding.T * fft.T * invNtt * fft * padding, hyper=0, tol=1.e-5, M=M, solver=cg)
#print('Elapsed time: ' + str(map_rlsw1.header['time']))

M = packing(1/map_naive.coverage)
if numpy.any(~numpy.isfinite(M)):
    raise ValueError()

def callback(x):
    pass

map_rlsw2_packed = mapper_rls(tod, projection, padding.T * fft.T * invNtt * fft * padding, hyper=0, tol=1.e-7, M=M, callback=callback)
print 'Elapsed time:', map_rlsw2_packed.header['TIME']

map_rlsw2 = packing.T(map_rlsw2_packed)
