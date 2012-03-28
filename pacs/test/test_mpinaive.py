import numpy as np
import os
import tamasis
from pyoperators.utils.mpi import distribute_slice
from pyoperators import DistributionGlobalOperator
from tamasis import PacsObservation, Map, ProjectionOperator, UnpackOperator, mapper_naive, MPI
from tamasis.numpyutils import assert_all_eq


data_dir = os.path.dirname(__file__) + '/data/'
rank = MPI.COMM_WORLD.rank
size = MPI.COMM_WORLD.size
tamasis.var.verbose = True

map_ref_local = Map(data_dir + '../../../core/test/data/frames_blue_map_naive.fits')
mask_ref_local = map_ref_local.coverage == 0
header_ref_local = map_ref_local.header

map_ref_global = Map(data_dir + '../../../core/test/data/frames_blue_map_naive.fits', comm=MPI.COMM_SELF)
mask_ref_global = map_ref_global.coverage == 0
header_ref_global = map_ref_global.header

def check_map_global(m, mask):
    assert_all_eq(m.magnitude, map_ref_global.magnitude, 1e-10)
    assert_all_eq(m.coverage, map_ref_global.coverage, 1e-10)
    assert_all_eq(mask, mask_ref_global)
    assert_all_eq(m.tounit('Jy/arcsec^2'), map_ref_global.tounit('Jy/arcsec^2'),
                  1e-10)

def check_map_local(m, mask):
    assert_all_eq(m.magnitude, map_ref_local.magnitude, 1e-10)
    assert_all_eq(m.coverage, map_ref_local.coverage, 1e-10)
    assert_all_eq(mask, mask_ref_local)
    assert_all_eq(m.tounit('Jy/arcsec^2'), map_ref_local.tounit('Jy/arcsec^2'),
                  1e-10)

def mask_pmatrix(pmatrix, tod):
    if not isinstance(pmatrix, (list, tuple)):
        pmatrix = (pmatrix,)
    dest = 0
    for p in pmatrix:
        n = p.shape[-2]
        p.value.T[...] *= 1 - tod[:,dest:dest+n].mask.T
        dest += n

obs1 = PacsObservation(data_dir + 'frames_blue.fits', reject_bad_line=False)
obs1.pointing.chop = 0
tod1 = obs1.get_tod(flatfielding=False, subtraction_mean=False)
obs2 = PacsObservation([data_dir + 'frames_blue.fits[1:176]',
                        data_dir + 'frames_blue.fits[177:360]'],
                       reject_bad_line=False)
obs2.pointing.chop = 0
tod2 = obs2.get_tod(flatfielding=False, subtraction_mean=False)

tolocal = DistributionGlobalOperator(map_ref_global.shape,
                                     attrin={'header':header_ref_global})

# non-distributed map, distributed TOD
def test1():
    comm_tod = MPI.COMM_WORLD
    comm_map = MPI.COMM_SELF
    for obs, tod in ((obs1, tod1), (obs2, tod2)):
        pmatrix = obs.get_pointing_matrix(header=map_ref_global.header,
                      downsampling=True, npixels_per_sample=6)
        mask_pmatrix(pmatrix, tod)
        proj = ProjectionOperator(pmatrix, commin=comm_map, commout=comm_tod)
        m = mapper_naive(tod, proj)
        mask = proj.get_mask()
        yield check_map_global, m, mask

# non-distributed map, packed projection, distributed TOD
def test2():
    comm_tod = MPI.COMM_WORLD
    comm_map = MPI.COMM_SELF
    for obs, tod in ((obs1, tod1), (obs2, tod2)):
        pmatrix = obs.get_pointing_matrix(header=map_ref_global.header,
                      downsampling=True, npixels_per_sample=6)
        mask_pmatrix(pmatrix, tod)
        proj = ProjectionOperator(pmatrix, commin=comm_map, commout=comm_tod,
                                  packed=True)
        m = mapper_naive(tod, proj)
        mask = proj.get_mask()
        yield check_map_global, m, mask

# distributed map, distributed TOD
def test3():
    comm_tod = MPI.COMM_WORLD
    comm_map = MPI.COMM_WORLD
    for obs, tod in ((obs1, tod1), (obs2, tod2)):
        pmatrix = obs.get_pointing_matrix(header=map_ref_local.header,
                      downsampling=True, npixels_per_sample=6)
        mask_pmatrix(pmatrix, tod)
        proj = ProjectionOperator(pmatrix, commin=comm_map, commout=comm_tod)
        m = mapper_naive(tod, proj)
        mask = proj.get_mask()
        yield check_map_local, m, mask

        if size <= 1:
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
    comm_tod = MPI.COMM_WORLD
    comm_map = MPI.COMM_SELF
    for obs, tod in ((obs1, tod1), (obs2, tod2)):
        pmatrix = obs.get_pointing_matrix(header=header_ref_local,
                      downsampling=True, npixels_per_sample=6)
        mask_pmatrix(pmatrix, tod)
        proj = ProjectionOperator(pmatrix, commin=MPI.COMM_WORLD,
                                  commout=comm_tod)
        model = proj * tolocal
        m = mapper_naive(tod, model)
        yield check_map_global, m, tolocal.T((proj.get_mask()))