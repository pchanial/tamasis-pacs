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
projection   = Projection(pacs, resolution=3.2, oversampling=False, npixels_per_sample=6)
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
map_iter = mapper_rls(tod, model, hyper=1e1, tol=1.e-4)
