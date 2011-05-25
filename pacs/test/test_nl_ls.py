import numpy as np
import os
import scipy
import tamasis
from tamasis import *

class TestFailure(Exception): pass

tamasis.var.verbose = False
data_dir = '~/work/tamasis/tamasis-dev/pacs/test/data/'
obs = PacsObservation(data_dir+'frames_blue.fits', fine_sampling_factor=1)
tod = obs.get_tod(flatfielding=False)

projection   = Projection(obs, oversampling=False, npixels_per_sample=6)
masking_tod  = Masking(tod.mask)
masking_map  = Masking(projection.mask)

model = masking_tod * projection * masking_map

# naive map
map_naive = mapper_naive(tod, model)

# iterative map, taking all map pixels
class Callback():
    def __init__(self):
        self.niterations = 0
    def __call__(self, x):
        self.niterations += 1

invntt = Diagonal(1/obs.get_detector_stddev(100)**2)
norm = tamasis.linalg.norm2_ellipsoid(invntt)
precond = 1./map_naive.coverage
precond[precond > 1] = 0
map_nl = mapper_nl(tod, model,
                   tol=1.e-6,
                   maxiter=1000,
                   norms=[norm],
                   M=Diagonal(precond),
                   callback=Callback())

print('Elapsed time:', map_nl.header['TIME'])
if map_nl.header['NITER'] > 21:
    raise TestFailure()
