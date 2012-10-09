import numpy as np
import pyoperators
import os
import scipy
from tamasis import (PacsObservation, CompressionAverageOperator,
                     DiagonalOperator, IdentityOperator, MaskOperator,
                     ProjectionOperator, UnpackOperator, mapper_ls,
                     mapper_naive)

pyoperators.memory.verbose = False
profile = None#'test_ls.png'
solver = scipy.sparse.linalg.bicgstab
tol = 1.e-6 if profile else 1.e-4
maxiter = 10
data_dir = os.path.dirname(__file__) + '/data/'
obs = PacsObservation(data_dir + 'frames_blue.fits', fine_sampling_factor=1,
                      reject_bad_line=False)
tod = obs.get_tod()

telescope   = IdentityOperator()
projection  = ProjectionOperator(obs, downsampling=True, npixels_per_sample=6)
compression = CompressionAverageOperator(obs.slice.compression_factor)
masking_tod = MaskOperator(tod.mask)
masking_map = MaskOperator(projection.get_mask())

model = masking_tod * projection * telescope * masking_map

# naive map
map_naive = mapper_naive(tod, model)

# iterative map, restricting oneself to observed map pixels
unpacking = UnpackOperator(projection.get_mask())
old_settings = np.seterr(divide='ignore')
M = DiagonalOperator(unpacking.T(1./map_naive.coverage))
np.seterr(**old_settings)
#map_iter1 = mapper_ls(tod, model * unpacking, tol=1.e-4, M=M)
#if map_iter1.header['NITER'] > 11:
#    raise TestFailure()

# iterative map, taking all map pixels
class Callback():
    def __init__(self):
        self.niterations = 0
    def __call__(self, x):
        self.niterations += 1
callback = Callback()
#callback=None

map_iter2 = mapper_ls(tod, model,
                      tol=tol,
                      maxiter=maxiter,
                      M=DiagonalOperator(masking_map(1./map_naive.coverage)),
                      callback=callback,
                      solver=solver,
                      profile=profile)
if profile is None:
    print 'Elapsed time:', map_iter2.header['TIME']

def test():
    assert map_iter2.header['NITER'] <= 10
