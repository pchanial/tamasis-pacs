from   tamasis import *
import pyfits

observation = MadMap1Observation(tamasis_dir+'tests/madmap1/todSpirePsw_be', tamasis_dir+'tests/madmap1/invnttSpirePsw_be', 
                                 tamasis_dir+'tests/madmap1/madmapSpirePsw.fits[coverage]', 'big_endian', 135, missing_value=numpy.nan)

tod = observation.get_tod()
sqrtInvNtt = SqrtInvNtt(observation, tamasis_dir+'tests/madmap1/invnttSpirePsw_be', convert='big_endian')

fft = Fft(tod.nsamples)
#padding = Padding(left = sqrtInvNtt.ncorrelations, right = 1024 - tod.shape[-1] - sqrtInvNtt.ncorrelations)
projection = Projection(observation)
packing = Packing(observation.mapmask)

model = projection * packing

map_naive = mapper_naive(tod, model)
map_ref = pyfits.fitsopen(tamasis_dir+'tests/madmap1/naivemapSpirePsw.fits')['image'].data
if any_neq(map_naive,map_ref,15): print 'FAILED: mapper_naive madcap 1'

map_naive_1d = mapper_naive(tod, projection)
map_naive_2d = packing.transpose(map_naive_1d)
if any_neq(map_naive_2d,map_ref,15): print 'FAILED: mapper_naive madcap 2'

# we solve M x = y0
# M^T M x = M^T y0
# M^T N^-1 M x = (M^T F^T s^T) (s F M) y0 = (M^T F^T s^T) s F y0
#model = sqrtInvNtt * fft * padding * projection * packing
#rhs = model.transpose(sqrtInvNtt.direct(fft.direct(tod)))
#map_filtered = mapper_rls(tod, model, rhs=rhs)
#if numpy.any(numpy.isnan(map_filtered)): print 'FAILED rls'

print 'OK.'
