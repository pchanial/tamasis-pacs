import os
import tamasis

from pyoperators import (DiagonalOperator, DistributionGlobalOperator,
                         MaskOperator)
from pyoperators.utils.mpi import MPI
from tamasis import PacsObservation, mapper_ls, mapper_naive
from tamasis.utils import assert_all_eq


#solver = scipy.sparse.linalg.bicgstab
solver = tamasis.solvers.cg
tol = 1.e-2
mtol = 1.e-10
maxiter = 100

rank = MPI.COMM_WORLD.rank
tamasis.var.verbose = True
profile = None#'test_ls.png'
data_dir = os.path.dirname(__file__) + '/data/'

# reference map (no communication)
comm_tod = MPI.COMM_SELF
comm_map = MPI.COMM_SELF
obs_ref = PacsObservation(data_dir + 'frames_blue.fits', comm=comm_tod)
obs_ref.pointing.chop = 0
tod_ref = obs_ref.get_tod()
model_ref = MaskOperator(tod_ref.mask) * \
            obs_ref.get_projection_operator(downsampling=True,
                                            npixels_per_sample=6,
                                            commin=comm_map)
map_naive_ref = mapper_naive(tod_ref, model_ref, unit='Jy/arcsec^2')
map_ref_global = mapper_ls(tod_ref, model_ref, tol=tol, maxiter=maxiter,
                           solver=solver, M=DiagonalOperator(
                           1/map_naive_ref.coverage))
cov_ref_global = map_ref_global.coverage
mask_ref_global = map_ref_global.coverage == 0
header_ref_global = map_ref_global.header
tolocal = DistributionGlobalOperator(map_naive_ref.shape,
                                     attrin={'header':header_ref_global})
map_ref_local = tolocal(map_ref_global)
cov_ref_local = tolocal(map_ref_global.coverage)
mask_ref_local = tolocal(mask_ref_global)

def check_map_global(m):
    assert_all_eq(m.magnitude, map_ref_global.magnitude, mtol)
    assert_all_eq(m.coverage, cov_ref_global, mtol)

def check_map_local(m):
    assert_all_eq(m.magnitude, map_ref_local.magnitude, mtol)
    assert_all_eq(m.coverage, cov_ref_local, mtol)


obs1 = PacsObservation(data_dir + 'frames_blue.fits')
obs1.pointing.chop = 0
tod1 = obs1.get_tod()
obs2 = PacsObservation([data_dir + 'frames_blue.fits[1:176]',
                        data_dir + 'frames_blue.fits[177:360]'])
obs2.pointing.chop = 0
tod2 = obs2.get_tod()

# non-distributed map, distributed TOD
def test1():
    comm_map = MPI.COMM_SELF
    def func(obs, tod):
        masking = MaskOperator(tod.mask)
        proj = obs.get_projection_operator(downsampling=True,
                                           npixels_per_sample=6,
                                           header=header_ref_global,
                                           commin=comm_map)
        model = masking * proj
        m = mapper_ls(tod, model, tol=tol, maxiter=maxiter, solver=solver,
                      M=DiagonalOperator(1/cov_ref_global))
        check_map_global(m)
    for obs, tod in ((obs1,tod1), (obs2,tod2)):
        yield func, obs, tod

# non-distributed map, packed projection, distributed TOD
def test2():
    comm_map = MPI.COMM_SELF
    for obs, tod in ((obs1, tod1), (obs2, tod2)):
        masking = MaskOperator(tod.mask)
        proj = obs.get_projection_operator(downsampling=True,
                                           npixels_per_sample=6,
                                           header=header_ref_global,
                                           commin=comm_map,
                                           packed=True)
        model = masking * proj
        m = mapper_ls(tod, model, tol=tol, maxiter=maxiter, solver=solver,
                      M=DiagonalOperator(1/cov_ref_global))
        yield check_map_global, m

# distributed map, distributed TOD
def test3():
    comm_map = MPI.COMM_WORLD
    for obs, tod in ((obs1,tod1), (obs2,tod2)):
        masking = MaskOperator(tod.mask) 
        proj = obs.get_projection_operator(downsampling=True,
                                           npixels_per_sample=6,
                                           header=header_ref_global,
                                           commin=comm_map)
        model = masking * proj
        m = mapper_ls(tod, model, tol=tol, maxiter=maxiter, solver=solver,
                      M=DiagonalOperator(1/cov_ref_local))
        yield check_map_local, m

# same map for all processes but also distributed as local maps, distributed TOD
def test4():
    comm_map = MPI.COMM_SELF
    for obs, tod in ((obs1,tod1), (obs2,tod2)):
        masking = MaskOperator(tod.mask)
        proj = obs.get_projection_operator(downsampling=True,
                                           npixels_per_sample=6,
                                           header=header_ref_global,
                                           commin=MPI.COMM_WORLD)
        model = masking * proj * tolocal
        m = mapper_ls(tod, model, tol=tol, maxiter=maxiter, solver=solver,
                      M=DiagonalOperator(1/cov_ref_global))
        yield check_map_global, m
