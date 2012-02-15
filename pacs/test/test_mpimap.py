import numpy as np
import os
import tamasis
from pyoperators.utils.mpi import distribute_slice
from tamasis import PacsObservation, Map, DistributionGlobalOperator, DistributionLocalOperator, MaskOperator, PackOperator, ProjectionOperator, mapper_naive, MPI
from tamasis.numpyutils import assert_all_eq

from nose.plugins.skip import SkipTest
raise SkipTest

data_dir = os.path.dirname(__file__) + '/data/'
rank = MPI.COMM_WORLD.Get_rank()
size = MPI.COMM_WORLD.Get_size()
tamasis.var.comm_tod = MPI.COMM_WORLD
tamasis.var.comm_map = MPI.COMM_WORLD
tamasis.var.verbose = True

map_ref = Map(data_dir + '../../../core/test/data/frames_blue_map_naive.fits')
obs = PacsObservation(data_dir + 'frames_blue.fits')
tod = obs.get_tod(flatfielding=False)

proj = ProjectionOperator(obs, downsampling=True, packed=True,
                          header=map_ref.header, npixels_per_sample=6).blocks[0]
masking = MaskOperator(tod.mask)
packing = PackOperator(proj.mask, attrin = {'header':proj.header})
distrib = DistributionGlobalOperator(proj.mask.shape, share=True)
model = masking * proj * packing * distrib

m = mapper_naive(tod, model)
assert_all_eq(m.magnitude, map_ref.magnitude, 1e-11)
assert_all_eq(m.coverage, map_ref.coverage, 1e-11)

proj_local = DistributionLocalOperator(proj.mask)
proj_global = DistributionGlobalOperator(proj.mask.shape)
model = masking * proj * proj_local
mlocal = model.T(np.ones(tod.shape))
mglobal = proj_global.T(mlocal)
assert_all_eq(mglobal, map_ref.coverage, 1.e-11)

mlocal2 = proj_global(mglobal)
assert_all_eq(mlocal2, mlocal)

mlocal[:] = rank + 1
mlocalpacked = proj_local(mlocal)
mlocalunpacked = packing.T(mlocalpacked)

pb = False
for irank in range(size):
    s = distribute_slice(proj.mask.shape[0], rank=irank)
    if np.any((mlocalunpacked[s] != 0) & (mlocalunpacked[s] != irank+1)):
        pb = True
        print 'Problem in rank ' + str(rank) + ' with local ' + str(irank)
assert pb is False
