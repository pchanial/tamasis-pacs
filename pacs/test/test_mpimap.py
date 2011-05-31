import numpy as np
import os
import tamasis
from mpi4py import MPI
from tamasis import *

class TestFailure(Exception): pass

data_dir = os.path.dirname(__file__) + '/data/'
rank = MPI.COMM_WORLD.Get_rank()
size = MPI.COMM_WORLD.Get_size()
tamasis.var.comm_tod = MPI.COMM_WORLD
tamasis.var.comm_map = MPI.COMM_WORLD
tamasis.var.verbose = True

map_ref = Map(data_dir + '../../../core/test/data/frames_blue_map_naive.fits')
obs = PacsObservation(data_dir + 'frames_blue.fits')
tod = obs.get_tod(flatfielding=False)

proj = Projection(obs, oversampling=False, packed=True, header=map_ref.header,
                  npixels_per_sample=6).blocks[0]
masking = Masking(tod.mask)
packing = Packing(proj.mask, attrin = {'header':proj.header})
distrib = DistributionGlobal(proj.mask.shape, share=True)
model = masking * proj * packing * distrib

m = mapper_naive(tod, model)
if any_neq(m.magnitude, map_ref.magnitude, 1e-11): raise TestFailure()
if any_neq(m.coverage, map_ref.coverage, 1e-11): raise TestFailure()

proj_local = DistributionLocal(proj.mask)
proj_global = DistributionGlobal(proj.mask.shape)
model = masking * proj * proj_local
mlocal = model.T(np.ones(tod.shape))
mglobal = proj_global.T(mlocal)
if any_neq(mglobal, map_ref.coverage, 1.e-11): raise TestFailure()

mlocal2 = proj_global(mglobal)
if any_neq(mlocal2, mlocal): raise TestFailure()

mlocal[:] = rank + 1
mlocalpacked = proj_local(mlocal)
mlocalunpacked = packing.T(mlocalpacked)

pb = False
for irank in range(size):
    s = tamasis.mpiutils.split_work(proj.mask.shape[0], rank=irank)
    if np.any((mlocalunpacked[s] != 0) & (mlocalunpacked[s] != irank+1)):
        pb = True
        print 'Problem in rank ' + str(rank) + ' with local ' + str(irank)
if pb: raise TestFailure()
