import numpy as np
import os
import tamasis
from tamasis import *

class TestFailure(Exception): pass

data_dir = os.path.dirname(__file__) + '/data/'
rank = tamasis.var.mpi_comm.Get_rank()
size = tamasis.var.mpi_comm.Get_size()

map_ref = Map(data_dir + 'frames_blue_map_naive.fits')
obs = PacsObservation(data_dir + 'frames_blue.fits')
tod = obs.get_tod(flatfielding=False)

proj = Projection(obs, oversampling=False, packed=True, header=map_ref.header, npixels_per_sample=6)
masking = Masking(tod.mask)
packing = Packing(proj.mask, attrin = {'header':proj.header})
model = masking * proj * packing

m = mapper_naive(tod, model)
if any_neq(m.magnitude, map_ref.magnitude, 1e-11): raise TestFailure()
if any_neq(m.coverage, map_ref.coverage, 1e-11): raise TestFailure()

reduce = AllReduceLocal(proj.mask)
gather = AllGather(proj.mask.shape)

o = masking(np.ones(tod.shape))
mlocal = reduce(proj.T(o)) # input is packed, output is local
mglobal = gather(mlocal)
if any_neq(mglobal, map_ref.coverage, 1.e-11): raise TestFailure()

mlocal2 = gather.T(mglobal)
if any_neq(mlocal2, mlocal): raise TestFailure()

mlocal[:] = rank + 1
mlocalpacked = reduce.T(mlocal)
mlocalunpacked = packing.T(mlocalpacked)

pb = False
for irank in range(size):
    s = tamasis.mpiutils.split_work(tamasis.var.mpi_comm, proj.mask.shape[0], rank=irank)
    if np.any((mlocalunpacked[s] != 0) & (mlocalunpacked[s] != irank+1)):
        pb = True
        print 'Problem in rank ' + str(rank) + ' with local ' + str(irank)
if pb: raise TestFailure()
