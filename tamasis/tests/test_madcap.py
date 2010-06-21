import numpy
import pyfits

from tamasis import *

observation = MadMap1Observation(tamasis_dir+'tests/madmap1/todSpirePsw_be', tamasis_dir+'tests/madmap1/invnttSpirePsw_be', 
                                 tamasis_dir+'tests/madmap1/madmapSpirePsw.fits[coverage]', 'big_endian', 135, missing_value=numpy.nan)

tod = observation.get_tod()
tod.unit = 'Jy/beam'
invNtt= InvNtt((observation.ndetectors, len(tod.nsamples)*(1024,)), tamasis_dir+'tests/madmap1/invnttSpirePsw_be', convert='big_endian')
fft = Fft(len(tod.nsamples)*(1024,))
padding = Padding(left=invNtt.ncorrelations, right=1024-numpy.array(tod.nsamples)-invNtt.ncorrelations)
projection = Projection(observation)
packing = Unpacking(observation.mapmask, field=numpy.nan).T

map_naive = mapper_naive(tod, projection*packing)
map_ref = pyfits.fitsopen(tamasis_dir+'tests/madmap1/naivemapSpirePsw.fits')['image'].data
if any_neq(map_naive,map_ref,15): print 'FAILED: mapper_naive madcap 1'

map_naive_1d = mapper_naive(tod, projection)
map_naive_2d = packing.transpose(map_naive_1d)
if any_neq(map_naive_2d,map_ref,15): print 'FAILED: mapper_naive madcap 2'

packing = Unpacking(observation.mapmask).T

M = 1 / map_naive.coverage
M[numpy.isfinite(M) == False] = numpy.max(M[numpy.isfinite(M)])
map_rlsw1 = mapper_rls(tod, projection*packing, padding.T * fft.T * invNtt * fft * padding, hyper=0, tol=1.e-5, M=M)
print map_rlsw1.header['time']

M = packing(1/map_naive.coverage)
M[numpy.isfinite(M) == False] = numpy.max(M[numpy.isfinite(M)])
M = M.reshape(numpy.product(M.shape))
map_rlsw2 = packing.T(mapper_rls(tod, projection, padding.T * fft.T * invNtt * fft * padding, hyper=0, tol=1.e-5, M=M))

print 'OK.'
