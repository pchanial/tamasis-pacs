import numpy
import os
import scipy
import tamasis

from scipy.sparse import dia_matrix
from scipy.sparse.linalg import LinearOperator, cgs
from tamasis import *

class TestFailure(Exception): pass

tamasis.var.verbose = False
data_dir = os.path.dirname(__file__) + '/data/'
obs = PacsObservation(filename=data_dir+'frames_blue.fits', fine_sampling_factor=1)

tod = obs.get_tod()

telescope    = Identity('Telescope PSF')
projection   = Projection(obs, resolution=3.2, oversampling=False, npixels_per_sample=6)
multiplexing = CompressionAverage(obs.instrument.fine_sampling_factor, 'Multiplexing')
crosstalk    = Identity('Crosstalk')
compression  = CompressionAverage(obs)
masking_tod  = Masking(tod.mask)

model = masking_tod * crosstalk * multiplexing * projection * telescope
print(model)

# naive map * masking_map
map_naive = mapper_naive(tod, model)
weights = map_naive.coverage
map_mask = weights == 0


# iterative map, taking all map pixels
map_iter = mapper_rls(tod, model, hyper=1., tol=1.e-4)
print 'Elapsed time: '+str(map_iter.header['TIME'])
if map_iter.header['NITER'] > 63:
    raise TestFailure()
