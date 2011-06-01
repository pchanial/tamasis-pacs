import numpy as np
import tamasis
from tamasis import *

class TestFailure(Exception): pass
tamasis.var.verbose = True

rank = MPI.COMM_WORLD.Get_rank()
size = MPI.COMM_WORLD.Get_size()

# test DistributedArray
shapes = [ (10,), (5,4), (3,4,2) ]
dtypes = (int, float, complex)
for shape in shapes:
    for dtype in dtypes:
        for order in ('c', 'f'):
            a = np.random.random_sample(shape).astype(dtype)
            a = MPI.COMM_WORLD.bcast(a)
            objs = [
                FitsArray(a, dtype=dtype, order=order),
                Map(a, dtype=dtype, order=order, error=a, coverage=a),
                Tod(a, dtype=dtype, order=order, mask=a.astype(np.bool))
                ]
            for obj in objs:
                if isinstance(obj, Tod) and obj.ndim == 1:
                    continue
                local = obj.tolocal()
                globa = local.toglobal()
                if any_neq(obj, globa):
                    raise TestFailure()
