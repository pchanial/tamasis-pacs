import numpy as np
from numpy.testing import assert_equal
from tamasis import DistributionGlobalOperator, utils as u, MPI

rank = MPI.COMM_WORLD.rank
size = MPI.COMM_WORLD.size

# diff, diffT, diffTdiff
def test_diff():
    shape = np.array((4,8,4))
    for r in range(1,shape.size+1):
        for axis in range(r):
            a=np.random.random_integers(1, 10, size=shape[0:r]).astype(float)
            MPI.COMM_WORLD.Bcast([a, MPI.DOUBLE], root=0)
            distribution = DistributionGlobalOperator(shape[0:r],
                                                      commout=MPI.COMM_WORLD)
            for f in (u.diff, u.diffT, u.diffTdiff):
                ref = a.copy()
                f(ref, ref, axis=axis, comm=MPI.COMM_SELF)
                reflocal = distribution(ref)
                local = distribution(a)
                f(local, local, axis=axis)
                yield assert_equal, local, reflocal, 'r={0}, axis={1}, f={2}, rank={3}, size={4}'.format(r, axis, f.__name__, rank, size)

