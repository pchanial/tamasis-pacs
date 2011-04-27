import numpy as np
import tamasis
from kapteyn import wcs
from mpi4py import MPI
from tamasis import any_neq, DistributionGlobal, utils as u, mpiutils as mu, \
                    wcsutils as wu

rank = MPI.COMM_WORLD.Get_rank()
size = MPI.COMM_WORLD.Get_size()

class TestFailure(Exception):
    def __init__(self, *args):
        print('Test failure: rank ' + str(rank) + '\n' + ',\n'.join([repr(a) for a in args]))
        raise Exception()

vec = np.ones(rank+1)

# sum, dot
if mu.sum(vec) != size*(size+1)/2: raise TestFailure()
if mu.dot(vec,-vec) != -size*(size+1)/2: raise TestFailure()

# gradient
n = np.random.random_integers(10)
x = np.random.random_sample(n) + 1
dsum = mu.gradient(mu.sum)
if any_neq(dsum(x), 1., 1.e-6): raise TestFailure()
def _prod(x, comm=None):
    if comm is None: comm = MPI.COMM_WORLD
    output = np.product(x)
    return comm.allreduce(output, op=MPI.PROD)
if any_neq(_prod([2,2]), 4**size): raise TestFailure()
dprod = mu.gradient(_prod)
if any_neq(dprod(x), _prod(x)/x, 1.e-6): raise TestFailure()

# norms
a = 4 * (np.random.random_sample(10) - 0.5)
aglobal = np.empty(10*size)
MPI.COMM_WORLD.Allgather([a, MPI.DOUBLE], [aglobal, MPI.DOUBLE])

norm_l1_ref = np.sum(np.abs(aglobal))
norm_l2_ref = np.sqrt(np.sum(aglobal**2))
norm_lp_ref = (np.sum(np.abs(aglobal)**-0.11))**(-1/0.11)
norm_linf_ref = np.max(np.abs(aglobal))
norm2_ref = np.sum(aglobal**2)
normp_ref = np.sum(np.abs(aglobal)**2.8)

tamasis.var.verbose = True
if any_neq(norm_l1_ref, mu.norm_l1(a)): raise TestFailure()
if any_neq(norm_l2_ref, mu.norm_l2(a)): raise TestFailure()
if any_neq(norm_lp_ref, mu.norm_lp(-0.11)(a)): raise TestFailure()
if any_neq(norm_linf_ref, mu.norm_linf(a)): raise TestFailure()
if any_neq(norm2_ref, mu.norm2(a)): raise TestFailure()
if any_neq(normp_ref, mu.normp(2.8)(a)): raise TestFailure()

# norm gradients

dnorm_l1_ref = mu.gradient(mu.norm_l1)(a)
dnorm_l2_ref = mu.gradient(mu.norm_l2)(a)
dnorm_lp_ref = mu.gradient(mu.norm_lp(-0.11))(a)
dnorm2_ref   = mu.gradient(mu.norm2)(a)
dnormp_ref   = mu.gradient(mu.normp(-0.11))(a)
dnorm_huber_ref = mu.gradient(mu.norm_huber(1.2))(a)

tol = 1.e-3
if any_neq(dnorm_l1_ref, mu.dnorm_l1(a), tol): raise TestFailure()
if any_neq(dnorm_l2_ref, mu.dnorm_l2(a), tol): raise TestFailure()
if any_neq(dnorm_lp_ref, mu.dnorm_lp(-0.11)(a), tol): raise TestFailure()
if any_neq(dnorm2_ref, mu.dnorm2(a), tol): raise TestFailure()
if any_neq(dnormp_ref, mu.dnormp(-0.11)(a), tol): raise TestFailure()
if any_neq(dnorm_huber_ref, mu.dnorm_huber(1.2)(a), tol): raise TestFailure()
tamasis.var.verbose = False

# diff, diffT, diffTdiff
shape = np.array((4,8,4))
for r in range(1,shape.size+1):
    for axis in range(r):
        a=np.random.random_integers(1, 10, size=shape[0:r]).astype(float)
        MPI.COMM_WORLD.Bcast([a, MPI.DOUBLE], root=0)
        distribution = DistributionGlobal(shape[0:r], comm=MPI.COMM_WORLD)
        for f in (u.diff, u.diffT, u.diffTdiff):
            ref = a.copy()
            f(ref, axis=axis, comm=MPI.COMM_SELF)
            reflocal = distribution(ref)
            local = distribution(a)
            f(local, axis=axis)
            if any_neq(local, reflocal):
                print 'ERROR: ', r, axis, f
                print 'ERROR: input', a
                print 'ERROR: output', ref
                print 'ERROR: ', rank, 'reflocal', reflocal
                print 'ERROR: ', rank, 'local:', local
                raise TestFailure()

