import numpy as np
import tamasis

from tamasis import linalg as la, MPI
from tamasis.utils import assert_all_eq

from nose.plugins.skip import SkipTest
raise SkipTest

rank = MPI.COMM_WORLD.Get_rank()
size = MPI.COMM_WORLD.Get_size()

class TestFailure(Exception):
    def __init__(self, *args):
        print('Test failure: rank ' + str(rank) + '\n' + ',\n'.join([repr(a) for a in args]))
        raise Exception()

vec = np.ones(rank+1)

# sum, dot
if la.sum(vec) != size*(size+1)/2: raise TestFailure()
if la.dot(vec,-vec) != -size*(size+1)/2: raise TestFailure()

# gradient
n = np.random.random_integers(10)
x = np.random.random_sample(n) + 1
dsum = la.gradient(la.sum)
if any_neq(dsum(x), 1., 1.e-6): raise TestFailure()
def _prod(x, comm=None):
    if comm is None: comm = MPI.COMM_WORLD
    output = np.product(x)
    return comm.allreduce(output, op=MPI.PROD)
if any_neq(_prod([2,2]), 4**size): raise TestFailure()
dprod = la.gradient(_prod)
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
if any_neq(norm_l1_ref, la.norm_l1(a)): raise TestFailure()
if any_neq(norm_l2_ref, la.norm_l2(a)): raise TestFailure()
if any_neq(norm_lp_ref, la.norm_lp(-0.11)(a)): raise TestFailure()
if any_neq(norm_linf_ref, la.norm_linf(a)): raise TestFailure()
if any_neq(norm2_ref, la.norm2(a)): raise TestFailure()
if any_neq(normp_ref, la.normp(2.8)(a)): raise TestFailure()

# norm gradients

dnorm_l1_ref = la.gradient(la.norm_l1)(a)
dnorm_l2_ref = la.gradient(la.norm_l2)(a)
dnorm_lp_ref = la.gradient(la.norm_lp(-0.11))(a)
dnorm2_ref   = la.gradient(la.norm2)(a)
dnormp_ref   = la.gradient(la.normp(-0.11))(a)
dnorm_huber_ref = la.gradient(la.norm_huber(1.2))(a)

tol = 1.e-3
if any_neq(dnorm_l1_ref, la.norm_l1.D(a), tol): raise TestFailure()
if any_neq(dnorm_l2_ref, la.norm_l2.D(a), tol): raise TestFailure()
if any_neq(dnorm_lp_ref, la.norm_lp(-0.11).D(a), tol): raise TestFailure()
if any_neq(dnorm2_ref, la.norm2.D(a), tol): raise TestFailure()
if any_neq(dnormp_ref, la.normp(-0.11).D(a), tol): raise TestFailure()
if any_neq(dnorm_huber_ref, la.norm_huber(1.2).D(a), tol): raise TestFailure()

