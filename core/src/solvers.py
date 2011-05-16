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
    """OpenMPI/MPI hybrid conjugate gradient solver with preconditioning."""

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

    if rank == 0:
        print('Iteration\tResiduals')

    for i in xrange(maxiter):
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
            print("%s\t%s" % (str(i+1).rjust(9), repr(np.sqrt(epsilon))))

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

def nlcg(criterion, n, linesearch, descent_method='pr', x0=None, M=None,
         tol=1.e-6, maxiter=500, callback=None, comm=None):
    """Non-linear conjugate gradient with preconditioning.

    Parameters
    ----------
    criterion : function
        cost function to be minimised
    n : integer
        size of criterion input vector
    linesearch : function
        line search function, minimising the criterion along a descent
    
    x0 : vector
        initial guess for the solution
    descent_method : string
        method for the descent vector update. Available methods are 
            - 'fs' (Fletcher-Reeves)
            - 'pr' (Polak-Ribiere)
            - 'hs' (Hestenes-Stiefel)
    callback : function
        callback function, called after each iteration
    tol : float
        Relative tolerance to achieve before terminating.
    maxiter : integer
        Maximum number of iterations.  Iteration will stop after maxiter
        steps even if the specified tolerance has not been achieved.
    M : {sparse matrix, dense matrix, LinearOperator}
        Preconditioner for A.  The preconditioner should approximate the
        inverse of A.  Effective preconditioning dramatically improves the
        rate of convergence, which implies that fewer iterations are needed
        to reach a given error tolerance.
    callback : function
        User-supplied function to call after each iteration.  It is called
        as callback(x), where x is the current solution vector.
    comm : MPI communicator

    Returns
    --------
    x : solution

    """

    # first guess
    if x0 is None:
        x = np.zeros(n)
    else:
        x = x0.flatten()

    # check preconditioner
    if M is not None:
        if not hasattr(M, 'matvec'):
            raise TypeError('The preconditioner is not a linear operator.')

    # initialise gradient and its conjugate
    r = np.empty(n)
    s = np.empty(n)
    d = np.empty(n)

    if comm is None:
        comm = MPI.COMM_WORLD

    rank = comm.Get_rank()

    # tolerance
    Js = criterion(x, gradient=r)
    J = sum(Js)
    r[:] = -r
    if M is not None:
        s[:] = M.matvec(r)
    else:
        s[:] = r
    d[:] = s
    delta_new = dot(r, d, comm)
    delta_0 = delta_new
    Jnorm = J
    resid = np.sqrt(delta_new/delta_0)

    iteration = 0

    if rank == 0:
        print('Iteration\tResiduals\tCriterion')

    while iteration < maxiter and resid > tol:
        iteration += 1

        # step
        a = linesearch(d, r)

        # update
        x += a * d

        # criterion and gradient at new position
        Js = criterion(x, gradient=r)
        J = sum(Js)
        r[:] = -r

        delta_old = delta_new
        if descent_method in ('pr', 'hs'):
            delta_mid = dot(r, s, comm)
        if M is not None:
            s[:] = M.matvec(r)
        else:
            s[:] = r
        delta_new = dot(r, s, comm)

        # descent direction
        if descent_method == 'fr':
            b = delta_new / delta_old
        elif descent_method == 'pr':
            b = (delta_new - delta_mid) / delta_old
        elif descent_method == 'hs':
            b = (delta_new - delta_mid) / (delta_mid - delta_old)
        else:
            raise ValueError("The descent update method must be 'fs' (Fletcher"\
                "-Reeves), 'pr' (Polak-Ribiere) or 'hs' (Hestenes-Stiefel).")

        if (iteration  % 50) == 1 or b <= 0:
            d[:] = s # reset to steepest descent
        else:
            d[:] = s + b * d

        resid = np.sqrt(delta_new/delta_0)

        if rank == 0:
            Jstr = repr(J)
            if len(Js) > 1:
                Jstr += ' ' + str(Js)
            print("%s\t%s\t%s" % (str(iteration).rjust(9), str(resid), Jstr))
        
        if callback is not None:
            callback(x)

    return x
