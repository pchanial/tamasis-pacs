import numpy as np
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

# sum, norm2, dot
if mu.sum(vec) != size*(size+1)/2: raise TestFailure()
if mu.norm2(vec) != size*(size+1)/2: raise TestFailure()
if mu.dot(vec,-vec) != -size*(size+1)/2: raise TestFailure()

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

