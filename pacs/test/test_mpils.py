import os
import tamasis
from tamasis import *
from mpi4py import MPI


#solver = scipy.sparse.linalg.bicgstab
solver = tamasis.solvers.cg
tol = 1.e-3
mtol = 1.e-8
maxiter = 100

class TestFailure(Exception): pass

rank = MPI.COMM_WORLD.Get_rank()
tamasis.var.verbose = True
profile = None#'test_ls.png'
data_dir = os.path.dirname(__file__) + '/data/'

# no communication for the reference map
tamasis.var.comm_tod = MPI.COMM_SELF
tamasis.var.comm_map = MPI.COMM_SELF
obs_ref = PacsObservation(data_dir + 'frames_blue.fits')
tod_ref = obs_ref.get_tod(flatfielding=False)
model_ref = Masking(tod_ref.mask) * \
            Projection(obs_ref, downsampling=True, npixels_per_sample=6)
map_naive_ref = mapper_naive(tod_ref, model_ref, unit='Jy/arcsec^2')
map_ls_ref = mapper_ls(tod_ref, model_ref, tol=tol, maxiter=maxiter,
                       solver=solver, M=DiagonalOperator(1/map_naive_ref.coverage))
header_ref = map_naive_ref.header

# LS map, not sliced
tamasis.var.comm_tod = MPI.COMM_WORLD
tamasis.var.comm_map = MPI.COMM_SELF
obs = PacsObservation(data_dir + 'frames_blue.fits')
tod = obs.get_tod(flatfielding=False)
masking = Masking(tod.mask)
proj = Projection(obs, downsampling=True, npixels_per_sample=6,
                  header=header_ref)
#proj_global = DistributionGlobal(proj_tod.mask.shape, share=True,
#                                 comm=MPI.COMM_WORLD)
model = masking * proj

map_naive_g1 = mapper_naive(tod, model, unit='Jy/arcsec^2')
map_ls_g1 = mapper_ls(tod, model, tol=tol, maxiter=maxiter, solver=solver,
                      M=DiagonalOperator(1/map_naive_ref.coverage))

# LS map, sliced
tamasis.var.comm_tod = MPI.COMM_WORLD
tamasis.var.comm_map = MPI.COMM_WORLD
proj = Projection(obs, downsampling=True, npixels_per_sample=6,
                  header=header_ref)
model = masking * proj

map_naive_local = mapper_naive(tod, model, unit='Jy/arcsec^2')
map_ls_local = mapper_ls(tod, model, tol=tol, maxiter=maxiter, solver=solver,
                         M=DiagonalOperator(1/map_naive_local.coverage))
proj_global = DistributionGlobal(proj.mask.shape)
map_naive_g2 = proj_global.T(map_naive_local)
map_ls_g2 = proj_global.T(map_ls_local)

# LS map, same for all processor and distributed as local maps
tamasis.var.comm_tod = MPI.COMM_WORLD
tamasis.var.comm_map = MPI.COMM_SELF
model = masking * proj * proj_global
map_naive_g3 = mapper_naive(tod, model, unit='Jy/arcsec^2')
map_ls_g3 = mapper_ls(tod, model, tol=tol, maxiter=maxiter, solver=solver,
                      M=DiagonalOperator(1/map_naive_ref.coverage),
                      comm_map=MPI.COMM_SELF)

if any_neq(map_naive_ref, map_naive_g1.magnitude, mtol): raise TestFailure()
if any_neq(map_naive_ref, map_naive_g2.magnitude, mtol): raise TestFailure()
if any_neq(map_naive_ref, map_naive_g3.magnitude, mtol): raise TestFailure()

if any_neq(map_naive_ref.coverage, map_naive_g1.coverage): raise TestFailure()
if any_neq(map_naive_ref.coverage, map_naive_g2.coverage): raise TestFailure()
if any_neq(map_naive_ref.coverage, map_naive_g3.coverage.magnitude): raise TestFailure()

if any_neq(map_ls_ref, map_ls_g1.magnitude, mtol): raise TestFailure()
if any_neq(map_ls_ref, map_ls_g2.magnitude, mtol): raise TestFailure()
if any_neq(map_ls_ref, map_ls_g3.magnitude, mtol): raise TestFailure()
