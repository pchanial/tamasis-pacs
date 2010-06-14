from tamasis import *
from copy import copy
import numpy
import scipy
from scipy.sparse import dia_matrix
from scipy.sparse.linalg import LinearOperator, cgs

pacs = PacsObservation(tamasis_dir+'tests/frames_blue.fits')
tod = pacs.get_tod()

telescope    = Identity('Telescope PSF')
projection   = Projection(pacs, oversampling=False, npixels_per_sample=6)
multiplexing = CompressionAverage(pacs.fine_sampling_factor, 'Multiplexing')
crosstalk    = Identity('Crosstalk')
compression  = CompressionAverage(pacs.compression_factor)
masking_tod  = Masking(tod.mask)

model = masking_tod * crosstalk * multiplexing * projection * telescope
print model

# naive map
map_naive = mapper_naive(tod, model)
map_mask = map_naive.coverage == 0


# iterative map, restricting oneself to observed map pixels
unpacking = Unpacking(map_mask)
M = unpacking.T(1./map_naive.coverage)
shape = 2*(numpy.sum(map_mask == False),)
M0  = dia_matrix((unpacking.transpose(1./map_naive.coverage), 0), shape=shape)
map_iter1 = unpacking(mapper_ls(tod, model * unpacking, tol=1.e-4, M=M0))


# iterative map, taking all map pixels
unpacking = Masking(map_mask)
M = 1./map_naive.coverage
M[numpy.where(map_mask)] = numpy.max(M[numpy.where(map_mask == False)])
M0 = unpacking.transpose(1./map_naive.coverage)
M0 = M0.reshape(map_naive.size)
M0 = dia_matrix((M0, 0), shape=2*(map_naive.size,))
map_iter2 = mapper_ls(tod, model * unpacking, tol=1.e-4, maxiter=200, M=M0)
print map_iter2.header['time']
