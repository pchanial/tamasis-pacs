from tamasis import *
from copy import copy
import numpy
import scipy
from scipy.sparse import dia_matrix
from scipy.sparse.linalg import LinearOperator, cgs

pacs = PacsObservation(tamasis_dir+'tests/frames_blue.fits')
tod = pacs.get_tod()

telescope    = Identity('Telescope PSF')
projection   = Projection(pacs, finer_sampling=False, npixels_per_sample=6)
multiplexing = CompressionAverage(pacs.fine_sampling_factor, 'Multiplexing')
crosstalk    = Identity('Crosstalk')
compression  = CompressionAverage(pacs.compression_factor)
masking_tod  = Masking(tod.mask)

model = masking_tod * crosstalk * multiplexing * projection * telescope
print model

# naive map
map_naive, weights = mapper_naive(tod, model, weights=True)
map_mask = weights == 0
map_naive.mask = map_mask
backmap = model.transpose(tod)
backmap.mask = map_mask


# iterative map, restricting oneself to observed map pixels
unpacking = Unpacking(map_mask)

shape = 2*(numpy.sum(map_mask == False),)
matvec = LeastSquareMatvec(model, unpacking)
operator = LinearOperator(matvec=matvec, dtype=numpy.float64, shape=shape)
b  = unpacking.transpose(backmap)
x0 = unpacking.transpose(map_naive)
M  = dia_matrix((unpacking.transpose(1./weights), 0), shape=shape)
solution, nit = cgs(operator, b, x0=x0, M=M, tol=1.e-4, maxiter=20, callback=PcgCallback())
map_iter1 = unpacking.direct(Map(solution))


# iterative map, taking all map pixels
unpacking = Masking(map_mask) * Reshaping(numpy.product(map_naive.shape), map_naive.shape)

shape = 2*(map_naive.size,)
matvec = LeastSquareMatvec(model, unpacking)
operator = LinearOperator(matvec=matvec, dtype=numpy.float64, shape=shape)
b  = unpacking.transpose(backmap)
x0 = unpacking.transpose(map_naive)
M  = dia_matrix((unpacking.transpose(1./weights), 0), shape=shape)
solution, nit = cgs(operator, b, x0=x0, M=M, tol=1.e-4, maxiter=200, callback=PcgCallback())
map_iter2 = unpacking.direct(Map(solution))
