import os
import tamasis
from tamasis import PacsObservation, DistributionGlobalOperator, MaskOperator, ProjectionOperator, create_fitsheader, mapper_naive, mapper_rls
from tamasis import MPI
from tamasis.numpyutils import assert_all_eq

#solver = scipy.sparse.linalg.bicgstab
solver = tamasis.solvers.cg
tol = 1.e-3
mtol = 1.e-8
maxiter = 100
hyper = 1.

rank = MPI.COMM_WORLD.Get_rank()
tamasis.var.verbose = True
profile = None#'test_rls.png'
data_dir = os.path.dirname(__file__) + '/data/'

header_ref = create_fitsheader((97,108), cdelt=3.2/3600,
                               crval=(245.998427916727,61.5147650744551))

# no communication for the reference map
tamasis.var.comm_tod = MPI.COMM_SELF
tamasis.var.comm_map = MPI.COMM_SELF
obs_ref = PacsObservation(data_dir + 'frames_blue.fits')
tod_ref = obs_ref.get_tod(flatfielding=False)
model_ref = MaskOperator(tod_ref.mask) * \
            ProjectionOperator(obs_ref, downsampling=True, npixels_per_sample=6,
                               header=header_ref)
map_naive_ref = mapper_naive(tod_ref, model_ref, unit='Jy/arcsec^2')
map_rls_ref = mapper_rls(tod_ref, model_ref, tol=tol, maxiter=maxiter,
                         solver=solver, hyper=hyper)

# LS map, not sliced
tamasis.var.comm_tod = MPI.COMM_WORLD
tamasis.var.comm_map = MPI.COMM_SELF
obs = PacsObservation(data_dir + 'frames_blue.fits')
tod = obs.get_tod(flatfielding=False)
masking = MaskOperator(tod.mask)
proj = ProjectionOperator(obs, downsampling=True, npixels_per_sample=6,
                          header=header_ref)
model = masking * proj

map_naive_g1 = mapper_naive(tod, model, unit='Jy/arcsec^2')
map_rls_g1 = mapper_rls(tod, model, tol=tol, maxiter=maxiter, solver=solver,
                        hyper=hyper)

# LS map, sliced
tamasis.var.comm_tod = MPI.COMM_WORLD
tamasis.var.comm_map = MPI.COMM_WORLD
proj = ProjectionOperator(obs, downsampling=True, npixels_per_sample=6,
                          header=header_ref)
model = masking * proj

map_naive_local = mapper_naive(tod, model, unit='Jy/arcsec^2')
map_rls_local = mapper_rls(tod, model, tol=tol, maxiter=maxiter, solver=solver,
                           hyper=hyper)
proj_global = DistributionGlobalOperator(proj.mask.shape)
map_naive_g2 = proj_global.T(map_naive_local)
map_rls_g2 = proj_global.T(map_rls_local)

# LS map, same for all processor and distributed as local maps
tamasis.var.comm_tod = MPI.COMM_WORLD
tamasis.var.comm_map = MPI.COMM_SELF
model = masking * proj * proj_global
map_naive_g3 = mapper_naive(tod, model, unit='Jy/arcsec^2')
map_rls_g3 = mapper_rls(tod, model, tol=tol, maxiter=maxiter, solver=solver,
                        comm_map=MPI.COMM_SELF, hyper=hyper)

assert_all_eq(map_naive_ref, map_naive_g1.magnitude, mtol)
assert_all_eq(map_naive_ref, map_naive_g2.magnitude, mtol)
assert_all_eq(map_naive_ref, map_naive_g3.magnitude, mtol)

assert_all_eq(map_naive_ref.coverage, map_naive_g1.coverage)
assert_all_eq(map_naive_ref.coverage, map_naive_g2.coverage)
assert_all_eq(map_naive_ref.coverage, map_naive_g3.coverage.magnitude)

assert_all_eq(map_rls_ref, map_rls_g1.magnitude, mtol)
assert_all_eq(map_rls_ref, map_rls_g2.magnitude, mtol)
assert_all_eq(map_rls_ref, map_rls_g3.magnitude, mtol)
