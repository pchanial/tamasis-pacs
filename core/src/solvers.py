# Copyrights 2010-2011 Pierre Chanial

import numpy as np
import tamasisfortran as tmf

from mpi4py import MPI
from scipy.sparse.linalg import aslinearoperator

from . import var
from .mpiutils import norm2, dot
from .acquisitionmodels import Identity, asacquisitionmodel

__all__ = []

def cg(A, b, x0, tol=1.e-5, maxiter=300, callback=None, M=None, comm=None):
    """OpenMPI/MPI hybrid conjugate gradient solver"""

    if comm is None:
        comm = MPI.COMM_WORLD

    rank = comm.Get_rank()

    A = asacquisitionmodel(A)
    if M is None:
        M = Identity()
    M = asacquisitionmodel(M)

    maxRelError = tol**2

    n = b.size
    x = np.empty(n)
#   d = np.empty(n)
#   q = np.empty(n)
    r = np.empty(n)
#   s = np.empty(n)
    xfinal = np.zeros(n)

    if x0 is None:
        x[:] = 0
    else:
        x[:] = x0.ravel()

    norm = norm2(b, comm)
    if norm == 0:
        return xfinal
    norm = 1./norm

    r[:] = b
    r -= A.matvec(x)
    epsilon = norm * norm2(r, comm)
    minEpsilon = epsilon

    d = M.matvec(r)
    delta0 = dot(r, d, comm)
    deltaNew = delta0

    for i in range(maxiter):
        if epsilon <= maxRelError:
            break
        
        q = d.copy()
        q = A.matvec(q, True, True, True)

        alpha = deltaNew / dot(d, q, comm)
        x += alpha * d
        r -= alpha * q
        epsilon = norm * norm2(r, comm)

        qnorm = np.sqrt(norm2(q, comm))
        if rank == 0:
            print("ccPCG:  Iteration %s, epsilon = %s, J(x) = %s." % (str(i+1),
                  str(np.sqrt(epsilon)), str(qnorm)))

        if epsilon < minEpsilon:
            xfinal[:] = x
            minEpsilon = epsilon

        s = r.copy()
        s = M.matvec(s, True, True, True)

        deltaOld = deltaNew

        deltaNew = dot(r, s, comm)
        beta = deltaNew / deltaOld
        d *= beta
        d += s
    
    if rank == 0:
        if minEpsilon > maxRelError:
            print("ccPCG terminated without reaching specified tolerance")
            print("after %s iterations.\n" % str(i+1))
        else:
            print("ccPCG terminated after reaching specified tolerance\n")
            print("in %s iterations.\n" % str(i+1))

    minEpsilon  = np.sqrt(minEpsilon)
    maxRelError = np.sqrt(maxRelError)
    print("Relative error is %s, " % str(minEpsilon))
    print("requested relative error is %s.\n" % str(maxRelError))

    if callback is not None:
        callback.niterations = i+1
        callback.residual = minEpsilon
        callback.criterion = np.sqrt(norm2(q-b, comm))

    return xfinal, 0
