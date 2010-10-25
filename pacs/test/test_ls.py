from tamasis import *
import numpy

obs = PacsObservation(tamasis_dir+'pacs/test/data/frames_blue.fits', fine_sampling_factor=1)
tod = obs.get_tod()

telescope    = Identity('Telescope PSF')
projection   = Projection(obs, oversampling=False, npixels_per_sample=6)
multiplexing = CompressionAverage(obs.instrument.fine_sampling_factor, 'Multiplexing')
crosstalk    = Identity('Crosstalk')
compression  = CompressionAverage(obs)
masking_tod  = Masking(tod.mask)

model = masking_tod * crosstalk * multiplexing * projection * telescope
print model

# naive map
map_naive = mapper_naive(tod, model)
map_mask = map_naive.coverage == 0

# iterative map, restricting oneself to observed map pixels
unpacking = Unpacking(map_mask)
M = unpacking.T(1./map_naive.coverage)
map_iter1 = unpacking(mapper_ls(tod, model * unpacking, tol=1.e-4, M=M))

# iterative map, taking all map pixels
unpacking = Masking(map_mask)
M = 1./map_naive.coverage
M[map_mask] = numpy.max(M[map_mask == False])
M0 = unpacking.transpose(1./map_naive.coverage)
map_iter2 = mapper_ls(tod, model * unpacking, tol=1.e-4, maxiter=200, M=M)
print map_iter2.header['time']
