import os
import scipy
import tamasis

from mpi4py import MPI
from scipy.sparse.linalg import cgs
from tamasis import *

class TestFailure(Exception): pass

tamasis.var.verbose = True
tamasis.var.comm_tod = MPI.COMM_WORLD
tamasis.var.comm_map = MPI.COMM_SELF

profile = None#'test_rls.png'
data_dir = os.path.dirname(__file__) + '/data/'
rank = MPI.COMM_WORLD.Get_rank()
size = MPI.COMM_WORLD.Get_size()
obs = PacsObservation(filename=data_dir + 'frames_blue.fits',
                      fine_sampling_factor=1)
obs.pointing.chop[:] = 0
tod = obs.get_tod(flatfielding=False)

map_ref = Map(data_dir + 'frames_blue_map_rls_cgs_tol1e-6.fits')

masking_tod  = Masking(tod.mask)
multiplexing = CompressionAverage(obs.instrument.fine_sampling_factor,
                                  description='Multiplexing')
projection   = Projection(obs, resolution=3.2, oversampling=False,
                          npixels_per_sample=6, header=map_ref.header)
telescope    = Identity(description='Telescope PSF')

model = masking_tod * multiplexing * projection * telescope
print(model)

map_naive = mapper_naive(tod, model)
map_naive_ref = Map(data_dir+'../../../core/test/data/frames_blue_map_naive.fits')
if any_neq(map_naive, map_naive_ref, 1.e-8): raise TestFailure()

# iterative map, including all map pixels
class Callback():
    def __init__(self):
        self.niterations = 0
    def __call__(self, x):
        self.niterations += 1

map_iter = mapper_rls(tod, model, hyper=1., tol=1.e-6, profile=profile,
                      callback=None if tamasis.var.verbose else Callback(),
                      maxiter=10 if profile else 1000, solver=cgs)

if profile is None:
    print 'Elapsed time: ' + str(map_iter.header['TIME'])
    if map_iter.header['NITER'] > 121:
        raise TestFailure()

if any_neq(map_ref, map_iter, 1.e-8): raise TestFailure()
