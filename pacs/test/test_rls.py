import os
import tamasis
import pyoperators

from scipy.sparse.linalg import cgs
from tamasis import *
from tamasis.numpyutils import any_neq

class TestFailure(Exception): pass

pyoperators.memory.verbose = False
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

map_rls = mapper_rls(tod, model, hyper=1., tol=1.e-4, profile=profile,
                     callback=None if tamasis.var.verbose else Callback(),
                     maxiter=10 if profile else 1000, solver=cgs)

if profile is None:
    print 'Elapsed time: ' + str(map_rls.header['TIME'])
    if map_rls.header['NITER'] > 48:
        raise TestFailure()

ref = Map(data_dir + 'frames_blue_map_rls_cgs_tol1e-6.fits')
ref.derived_units = map_rls.derived_units
cov = ref.coverage > 80
if any_neq(ref[cov], map_rls[cov], 1.e-1): raise TestFailure()
cov = ref.coverage > 125
if any_neq(ref[cov], map_rls[cov], 1.e-2): raise TestFailure()
