import numpy as np
from numpy.testing import assert_equal
from tamasis import DistributionGlobal, utils as u, MPI
from tamasis.mpiutils import distribute_shape, distribute_slice

rank = MPI.COMM_WORLD.rank
size = MPI.COMM_WORLD.size

def test_mpiutils():
    if size > 1:
        return
    for n in range(10):
        for sz in range(1,7):
            work = np.zeros(n, int)
            for i in range(n):
                work[i] = i % sz
            a = np.zeros(sz, int)
            for r in range(sz):
                a[r] = sum(work==r)
            stop = tuple(np.cumsum(a))
            start = (0,) + stop[:-1]
            for r in range(sz):
                yield assert_equal, a[r], distribute_shape((n,), rank=r, size=sz)[0]
                s = slice(start[r], stop[r])
                yield assert_equal, s, distribute_slice(n, rank=r, size=sz)

# diff, diffT, diffTdiff
def test_diff():
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
                yield assert_equal, local, reflocal, 'r={0}, axis={1}, f={2}, rank={3}, size={4}'.format(r, axis, f.__name__, rank, size)

