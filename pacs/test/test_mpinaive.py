import numpy as np
import os
import tamasis

from pyoperators import MPIDistributionGlobalOperator
from pyoperators.utils.mpi import MPI, distribute_slice
from pysimulators import Map
from tamasis import (PacsInstrument, PacsObservation, UnpackOperator,
                     mapper_naive)
from tamasis.utils import assert_all_eq


data_dir = os.path.dirname(__file__) + '/data/'
rank = MPI.COMM_WORLD.rank
size = MPI.COMM_WORLD.size
tamasis.var.verbose = True

PacsInstrument.info.CALFILE_BADP = tamasis.var.path + '/pacs/PCalPhotometer_Ba'\
                                   'dPixelMask_FM_v5.fits'
PacsInstrument.info.CALFILE_RESP = tamasis.var.path + '/pacs/PCalPhotometer_Re'\
                                   'sponsivity_FM_v5.fits'

map_ref_local = Map(data_dir + '../../../core/test/data/frames_blue_map_naive.fits')
mask_ref_local = map_ref_local.coverage == 0
header_ref_local = map_ref_local.header

map_ref_global = Map(data_dir + '../../../core/test/data/frames_blue_map_naive.fits', comm=MPI.COMM_SELF)
mask_ref_global = map_ref_global.coverage == 0
header_ref_global = map_ref_global.header

def check_map_global(m):
    assert_all_eq(m.magnitude, map_ref_global.magnitude, 1e-10)
    assert_all_eq(m.coverage, map_ref_global.coverage, 1e-10)
    assert_all_eq(m.tounit('Jy/arcsec^2'), map_ref_global.tounit('Jy/arcsec^2'),
                  1e-10)

def check_map_local(m):
    assert_all_eq(m.magnitude, map_ref_local.magnitude, 1e-10)
    assert_all_eq(m.coverage, map_ref_local.coverage, 1e-10)
    assert_all_eq(m.tounit('Jy/arcsec^2'), map_ref_local.tounit('Jy/arcsec^2'),
                  1e-10)

obs1 = PacsObservation(data_dir + 'frames_blue.fits', reject_bad_line=False)
obs1.pointing.chop = 0
tod1 = obs1.get_tod()
obs2 = PacsObservation([data_dir + 'frames_blue.fits[1:176]',
                        data_dir + 'frames_blue.fits[177:360]'],
                       reject_bad_line=False)
obs2.pointing.chop = 0
tod2 = obs2.get_tod()

tolocal = MPIDistributionGlobalOperator(map_ref_global.shape,
                                        attrin={'header':header_ref_global})

# non-distributed map, distributed TOD
def test1():
    comm_map = MPI.COMM_SELF
    for obs, tod in ((obs1, tod1), (obs2, tod2)):
        proj = obs.get_projection_operator(downsampling=True,
                   npixels_per_sample=6, header=map_ref_global.header,
                   commin=comm_map)
        proj.apply_mask(tod.mask)
        m = mapper_naive(tod, proj)
        yield check_map_global, m
        assert_all_eq(proj.get_mask(), mask_ref_global)

# non-distributed map, packed projection, distributed TOD
def test2():
    comm_map = MPI.COMM_SELF
    for obs, tod in ((obs1, tod1), (obs2, tod2)):
        proj = obs.get_projection_operator(downsampling=True,
                   npixels_per_sample=6, header=map_ref_global.header,
                   packed=True, commin=comm_map)
        proj.apply_mask(tod.mask)
        m = mapper_naive(tod, proj)
        yield check_map_global, m

# distributed map, distributed TOD
def test3():
    comm_map = MPI.COMM_WORLD
    for obs, tod in ((obs1, tod1), (obs2, tod2)):
        proj = obs.get_projection_operator(downsampling=True,
                   npixels_per_sample=6, header=map_ref_global.header,
                   commin=comm_map)
        proj.apply_mask(tod.mask)
        m = mapper_naive(tod, proj)
        yield check_map_local, m

        if size == 1:
            continue

        mlocal = m
        mlocal[:] = rank + 1
        mlocalpacked = proj.operands[1](mlocal)
        mglobal = UnpackOperator(proj.operands[1].mask)(mlocalpacked)
        def func(s, irank):
            assert np.all((mglobal[s] == 0) | (mglobal[s] == irank+1))
        for irank in range(size):
            s = distribute_slice(map_ref_global.shape[0], rank=irank)
            yield func, s, irank

# same map for all processors and distributed as local maps, distributed TOD
def test4():
    for obs, tod in ((obs1, tod1), (obs2, tod2)):
        proj = obs.get_projection_operator(downsampling=True,
                   npixels_per_sample=6, header=map_ref_local.header,
                   commin=MPI.COMM_WORLD)
        proj.apply_mask(tod.mask)
        model = proj * tolocal
        m = mapper_naive(tod, model)
        yield check_map_global, m
