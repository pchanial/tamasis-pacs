from tamasis import *
from copy import copy
import numpy
import scipy
from scipy.sparse import dia_matrix
from scipy.sparse.linalg import LinearOperator, cgs

pacs = PacsObservation(filename='tests/frames_blue.fits',
                       resolution=3.2,
                       fine_sampling_factor=1,
                       npixels_per_sample=6, 
                       keep_bad_detectors=False)

tod = pacs.get_tod()

telescope    = Identity('Telescope PSF')
projection   = Projection(pacs, finer_sampling=False)
multiplexing = CompressionAverage(pacs.fine_sampling_factor, 'Multiplexing')
crosstalk    = Identity('Crosstalk')
compression  = CompressionAverage(pacs.compression_factor)
masking      = Masking(tod.mask)

model = masking * crosstalk * multiplexing * projection * telescope
print model

# naive map
backmap = model.transpose(tod)
tod[:] = 1
weights = model.transpose(tod)
map_naive = Map(backmap / weights, mask=weights <=1)

# iterative map
backmap.mask = map_naive.mask
shape = 2*(map_naive.count(),)
matvec = RLSMatvec(1e-3, model, map_naive.mask)
operator = LinearOperator(matvec=matvec, dtype=numpy.float64, shape=shape)
b  = backmap.compressed()
x0 = map_naive.compressed()
M  = dia_matrix(((1./weights)[weights > 1], 0), shape=shape)
solution, nit = cgs(operator, b, x0=x0, M=M, tol=1.e-4, maxiter=20, callback=PcgCallback())
map_iter = map_naive.copy()
map_iter[map_naive.mask == False] = solution

#ra0  = 20.
#dec0 = 0.1
#time = numpy.arange(0.,100., 1./40)
#simulation = PacsSimulation(inputmap           = 
#                            time               = time \
#                            ra                 = numpy.linspace(ra0, ra0+0.1, nsamples)   \
#                            dec                = numpy.linspace(dec0, dec0+0.1, nsamples) \
#                            pa                 = numpy.zeros(nsamples) \
#                            chop               = numpy.zeros(nsamples) \
#                            array              = 'blue'        \
#                            npixelsPerSample   = 9             \
#                            observingMode      = 'transparent' \
#                            fineSamplingFactor = 1             \
#                            compressionFactor  = 1             \
#                            keepBadDetectors   = True)

