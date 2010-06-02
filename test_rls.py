from tamasis import *
from copy import copy
import numpy
import scipy
from scipy.sparse import dia_matrix
from scipy.sparse.linalg import LinearOperator, cgs

pacs = PacsObservation(filename=tamasis_dir+'tests/frames_blue.fits',
                       fine_sampling_factor=1,
                       keep_bad_detectors=False)

tod = pacs.get_tod()

telescope    = Identity('Telescope PSF')
projection   = Projection(pacs, resolution=3.2, finer_sampling=False, npixels_per_sample=6)
multiplexing = CompressionAverage(pacs.fine_sampling_factor, 'Multiplexing')
crosstalk    = Identity('Crosstalk')
compression  = CompressionAverage(pacs.compression_factor)
masking_tod  = Masking(tod.mask)

model = masking_tod * crosstalk * multiplexing * projection * telescope
print model

# naive map * masking_map
map_naive = mapper_naive(tod, model)
weights = map_naive.coverage
map_mask = weights == 0
map_naive.mask = map_mask
backmap = model.transpose(tod)
backmap.mask = map_mask

# iterative map, taking all map pixels
reshaping = Reshaping(numpy.product(map_naive.shape), map_naive.shape)

shape = 2*(map_naive.size,)
matvec = RegularizedLeastSquareMatvec(model, reshaping, hyper=1e1)
operator = LinearOperator(matvec=matvec, dtype=numpy.float64, shape=shape)
b  = reshaping.transpose(backmap)
x0 = reshaping.transpose(map_naive)
x0[numpy.isnan(x0)] = 0.
#m = reshaping.transpose(1./weights)
#m[numpy.isfinite(m) == False] = 1.
#m = numpy.maximum(m,1.)
#M = dia_matrix((m, 0), shape=shape)
solution, nit = cgs(operator, b, x0=x0, tol=1.e-4, maxiter=200, callback=PcgCallback())
map_iter = Map(reshaping.direct(solution))

map_rls = mapper_rls(tod, model)

dX = DiscreteDifference(axis=1)
dY = DiscreteDifference(axis=0)
C = model.T * model + 1.e1 * (dX.T * dX + dY.T * dY)
operator = C.aslinearoperator()
b = operator.unpacking.T(backmap)
x0 = operator.unpacking.T(map_naive)
x0[numpy.where(numpy.isfinite(x0) == False)] = 0
solution, nit = cgs(operator, b, x0=x0, tol=1.e-4, maxiter=200, callback=PcgCallback())
map_iter2 = Map(operator.unpacking(solution))
