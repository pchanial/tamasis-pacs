import numpy as np
import os
import tamasis
from tamasis import *

class TestFailure(Exception): pass

tamasis.var.verbose = False
profile = None#'test_ls.png'
data_dir = os.path.dirname(__file__) + '/data/'
obs = PacsObservation(data_dir+'frames_blue.fits', fine_sampling_factor=1)
tod = obs.get_tod()

telescope    = Identity(description='Telescope PSF')
projection   = Projection(obs, oversampling=False, npixels_per_sample=6)
multiplexing = CompressionAverage(obs.instrument.fine_sampling_factor,
                                  description='Multiplexing')
crosstalk    = Identity(description='Crosstalk')
compression  = CompressionAverage(obs.slice.compression_factor)
masking_tod  = Masking(tod.mask)

model = masking_tod * crosstalk * multiplexing * projection * telescope
#model = masking_tod * crosstalk * projection * telescope
print(model)

# naive map
map_naive = mapper_naive(tod, model)
map_mask = map_naive.coverage == 0

# iterative map, restricting oneself to observed map pixels
unpacking = Unpacking(map_mask)
old_settings = np.seterr(divide='ignore')
M = unpacking.T(1./map_naive.coverage)
np.seterr(**old_settings)
#map_iter1 = mapper_ls(tod, model * unpacking, tol=1.e-4, M=M)
#if map_iter1.header['NITER'] > 11:
#    raise TestFailure()

# iterative map, taking all map pixels
unpacking = Masking(map_mask)
old_settings = np.seterr(divide='ignore')
M = 1./map_naive.coverage
np.seterr(**old_settings)
M[map_mask] = np.max(M[map_mask == False])
old_settings = np.seterr(divide='ignore')
M0 = unpacking.T(1./map_naive.coverage)
np.seterr(**old_settings)
class Callback():
    def __init__(self):
        self.niterations = 0
    def __call__(self, x):
        self.niterations += 1
map_iter2 = mapper_ls(tod, model * unpacking, tol=1.e-4, maxiter=200, M=M, callback=Callback(), profile=profile)
if profile is None:
    print(map_iter2.header['time'])
    if map_iter2.header['NITER'] > 11:
        raise TestFailure()
