import numpy
from   tamasis import *
import pyfits
from   scipy.sparse import dia_matrix

observation = MadMap1Observation(tamasis_dir+'tests/madmap1/todSpirePsw_be', tamasis_dir+'tests/madmap1/invnttSpirePsw_be', 
                                 tamasis_dir+'tests/madmap1/madmapSpirePsw.fits[coverage]', 'big_endian', 135, missing_value=numpy.nan)

tod = observation.get_tod()
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

M = packing(1/map_naive.coverage)
M[numpy.where(numpy.isfinite(M) == False)] = numpy.max(M[numpy.where(numpy.isfinite(M))])
M = M.reshape(numpy.product(M.shape))
M  = dia_matrix((M, 0), shape=2*M.shape)
#M = None
map_rlsw1 = packing.T(mapper_rls(tod, projection, padding.T * fft.T * invNtt * fft * padding, hyper=0, maxiter=2000, M=M))
dfsdfsdf
M = 1/map_naive.coverage
M[numpy.where(numpy.isfinite(M) == False)] = numpy.max(M[numpy.where(numpy.isfinite(M))])
M = M.reshape(numpy.product(M.shape))
M  = dia_matrix((M, 0), shape=2*M.shape)

map_rlsw2 = mapper_rls(tod, projection*packing, padding.T * fft.T * invNtt * fft * padding, hyper=0, maxiter=2000, M=M)


# we solve M x = y0
# M^T M x = M^T y0
# M^T N^-1 M x = (M^T F^T s^T) (s F M) y0 = (M^T F^T s^T) s F y0
#model = sqrtInvNtt * fft * padding * projection * packing
#rhs = model.transpose(sqrtInvNtt.direct(fft.direct(tod)))
#map_filtered = mapper_rls(tod, model, rhs=rhs)
#if numpy.any(numpy.isnan(map_filtered)): print 'FAILED rls'

print 'OK.'
