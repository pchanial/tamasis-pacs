import os
import scipy
import tamasis

from scipy.sparse.linalg import cgs
from tamasis import *
from tamasis.linalg import norm2, norm2_ellipsoid

class TestFailure(Exception): pass

tamasis.var.verbose = False
profile = None#'test_rls.png'
data_dir = os.path.dirname(__file__) + '/data/'
obs = PacsObservation(filename=data_dir+'frames_blue.fits',
                      fine_sampling_factor=1)
obs.pointing.chop[:] = 0
tod = obs.get_tod(flatfielding=False)

projection   = Projection(obs, resolution=3.2, oversampling=False,
                          npixels_per_sample=6)
masking_tod  = Masking(tod.mask)

model = masking_tod * projection

# iterative map, taking all map pixels
class Callback():
    def __init__(self):
        self.niterations = 0
    def __call__(self, x):
        self.niterations += 1

invntt = Diagonal(1/obs.get_detector_stddev(100)**2)
invntt = Identity()
map_nl = mapper_nl(tod, model, hypers=2*[1.],
                   norms=[norm2_ellipsoid(invntt)] + 2*[norm2],
                   tol=1.e-4, maxiter=1000,
                   callback=None if tamasis.var.verbose else Callback(),
                   )

if profile is None:
    print 'Elapsed time: ' + str(map_nl.header['TIME']) + \
        ' after ' + str(map_nl.header['NITER']) + ' iterations.'
#    if map_nl.header['NITER'] > 48:
#        raise TestFailure()
