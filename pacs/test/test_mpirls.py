import os
import tamasis
from tamasis import (MPI, PacsObservation, DistributionGlobalOperator,
                     MaskOperator, ProjectionOperator, create_fitsheader,
                     mapper_rls)
from tamasis.numpyutils import assert_all_eq


#solver = scipy.sparse.linalg.bicgstab
solver = tamasis.solvers.cg
tol = 1.e-2
mtol = 1.e-10
maxiter = 100
hyper = 1

rank = MPI.COMM_WORLD.rank
tamasis.var.verbose = True
profile = None#'test_rls.png'
data_dir = os.path.dirname(__file__) + '/data/'

# reference map (no communication)
header_ref_global = create_fitsheader((97,108), cdelt=3.2/3600,
                                      crval=(245.998427916727,61.5147650744551))
comm_tod = MPI.COMM_SELF
comm_map = MPI.COMM_SELF
obs_ref = PacsObservation(data_dir + 'frames_blue.fits', comm_tod=comm_tod)
obs_ref.pointing.chop = 0
tod_ref = obs_ref.get_tod(flatfielding=False, subtraction_mean=False)
model_ref = MaskOperator(tod_ref.mask) * \
            ProjectionOperator(obs_ref, downsampling=True, npixels_per_sample=6,
                               commout=comm_tod, commin=comm_map,
                               header=header_ref_global)
map_ref_global = mapper_rls(tod_ref, model_ref, tol=tol, maxiter=maxiter,
                            solver=solver, hyper=hyper)
header_ref_global = map_ref_global.header
cov_ref_global = map_ref_global.coverage
mask_ref_global = map_ref_global.coverage == 0
tolocal = DistributionGlobalOperator(map_ref_global.shape,
                                     attrin={'header':header_ref_global})
map_ref_local = tolocal(map_ref_global)
cov_ref_local = tolocal(map_ref_global.coverage)
mask_ref_local = tolocal(mask_ref_global)

if rank == 0:
    map_ref_global.save('map_ref_global_v4.fits', comm=MPI.COMM_SELF)

def check_map_global(m):
    assert_all_eq(m.magnitude, map_ref_global.magnitude, mtol)
    assert_all_eq(m.coverage, cov_ref_global, mtol)

def check_map_local(m):
    assert_all_eq(m.magnitude, map_ref_local.magnitude, mtol)
    assert_all_eq(m.coverage, cov_ref_local, mtol)


obs1 = PacsObservation(data_dir + 'frames_blue.fits')
obs1.pointing.chop = 0
tod1 = obs1.get_tod(flatfielding=False, subtraction_mean=False)
obs2 = PacsObservation([data_dir + 'frames_blue.fits[1:176]',
                        data_dir + 'frames_blue.fits[177:360]'])
obs2.pointing.chop = 0
tod2 = obs2.get_tod(flatfielding=False, subtraction_mean=False)

# non-distributed map, distributed TOD
def test1():
    comm_tod = MPI.COMM_WORLD
    comm_map = MPI.COMM_SELF
    for obs, tod in ((obs1,tod1), (obs2,tod2)):
        masking = MaskOperator(tod.mask)
        proj = ProjectionOperator(obs, downsampling=True, npixels_per_sample=6,
                                  header=header_ref_global, commin=comm_map,
                                  commout=comm_tod)
        model = masking * proj
        m = mapper_rls(tod, model, tol=tol, maxiter=maxiter, solver=solver,
                       hyper=hyper)
        yield check_map_global, m

# non-distributed map, packed projection, distributed TOD
def test2():
    comm_tod = MPI.COMM_WORLD
    comm_map = MPI.COMM_SELF
    for obs, tod in ((obs1, tod1), (obs2, tod2)):
        masking = MaskOperator(tod.mask)
        proj = ProjectionOperator(obs, downsampling=True, npixels_per_sample=6,
                                  header=header_ref_global, commin=comm_map,
                                  commout=comm_tod, packed=True)
        model = masking * proj
        m = mapper_rls(tod, model, tol=tol, maxiter=maxiter, solver=solver,
                       hyper=hyper)
        yield check_map_global, m

# distributed map, distributed TOD
def test3():
    comm_tod = MPI.COMM_WORLD
    comm_map = MPI.COMM_WORLD
    for obs, tod in ((obs1,tod1), (obs2,tod2)):
        masking = MaskOperator(tod.mask)
        proj = ProjectionOperator(obs, downsampling=True, npixels_per_sample=6,
                                  header=header_ref_global, commin=comm_map,
                                  commout=comm_tod)
        model = masking * proj
        m = mapper_rls(tod, model, tol=tol, maxiter=maxiter, solver=solver,
                       hyper=hyper)
        yield check_map_local, m

# same map for all processes but also distributed as local maps, distributed TOD
def test4():
    comm_tod = MPI.COMM_WORLD
    comm_map = MPI.COMM_SELF
    for obs, tod in ((obs1,tod1), (obs2,tod2)):
        masking = MaskOperator(tod.mask)
        proj = ProjectionOperator(obs, downsampling=True, npixels_per_sample=6,
                                  header=header_ref_global,
                                  commin=MPI.COMM_WORLD, commout=comm_tod)
        model = masking * proj * tolocal
        m = mapper_rls(tod, model, tol=tol, maxiter=maxiter, solver=solver,
                       hyper=hyper)
        yield check_map_global, m
