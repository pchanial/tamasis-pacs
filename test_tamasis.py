from tamasis import Identity, CompressionAverage, PacsProjectionSharpEdges, PacsMultiplexing, Masking, PacsObservation, PacsSimulation, Map, Tod, hcss_photproject, LeastSquareMatvec, RLSMatvec,  PcgCallback
from copy import copy
import numpy
import scipy
from scipy.sparse import dia_matrix
from scipy.sparse.linalg import LinearOperator, cgs

datadir = '/home/pchanial/work/pacs/data/transparent/'
pacs = PacsObservation(filename=datadir+'NGC6946/1342184520_blue',
                       first=20000,
                       last=86000,
                       resolution=3.,
                       fine_sampling_factor=1,
                       npixels_per_sample=6)

tod = pacs.get_tod()

telescope    = Identity('Telescope PSF')
projection   = PacsProjectionSharpEdges(pacs)
multiplexing = CompressionAverage(pacs.fine_sampling_factor, 'Multiplexing')
crosstalk    = Identity('Crosstalk')
compression  = CompressionAverage(pacs.compression_factor)
masking      = Masking(tod.mask)

model = crosstalk * multiplexing * projection * telescope
print model

tod = pacs.get_tod()
backmap = copy(model.transpose(tod))
backmap.mask=numpy.ma.nomask
tod[:] = 1.
weights = copy(model.transpose(tod))
map0 = Map(backmap / weights, mask=weights <= 1)
backmap.mask = map0.mask

shape = 2*(map0.count(),)
matvec = RLSMatvec(1e-3, model, map0.mask)
operator = LinearOperator(matvec=matvec, dtype=numpy.float64, shape=shape)
b  = backmap.compressed()
x0 = map0.compressed()
M  = dia_matrix(((1./weights)[weights > 1], 0), shape=shape)
solution, nit = cgs(operator, b, x0=x0, M=M, tol=1.e-4, maxiter=20, callback=PcgCallback())
mymap = copy(map0)
mymap[mymap.mask == False] = solution

#map = hcss_photproject(pacs)

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

