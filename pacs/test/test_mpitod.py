import os
import tamasis

from scipy.sparse.linalg import cgs
from tamasis import (PacsObservation, Map, MaskOperator, ProjectionOperator,
                     MPI, mapper_naive, mapper_rls)
from tamasis.numpyutils import assert_all_eq

tamasis.var.verbose = True
comm_tod = MPI.COMM_WORLD
comm_map = MPI.COMM_SELF
rank = MPI.COMM_WORLD.rank
size = MPI.COMM_WORLD.size

profile = None#'test_rls.png'
data_dir = os.path.dirname(__file__) + '/data/'
obs = PacsObservation(filename=data_dir + 'frames_blue.fits',
                      reject_bad_line=False)
obs.pointing.chop[:] = 0
tod = obs.get_tod(flatfielding=False)

map_ref = Map(data_dir + 'frames_blue_map_rls_cgs_tol1e-6.fits', comm=comm_map)

masking_tod = MaskOperator(tod.mask)
projection  = ProjectionOperator(obs, resolution=3.2, downsampling=True,
                                 npixels_per_sample=6, header=map_ref.header,
                                 commin=comm_map, commout=comm_tod)
model = masking_tod * projection

map_naive = mapper_naive(tod, model)
map_naive_ref = Map(data_dir + '../../../core/test/data/frames_blue_map_naive.f'
                    'its', comm=comm_map)
def test_naive():
    assert_all_eq(map_naive, map_naive_ref.magnitude, 1.e-8)

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

def test_rls():
    assert map_iter.header['NITER'] < 121
    assert_all_eq(map_ref, map_iter.magnitude, 1.e-8)
