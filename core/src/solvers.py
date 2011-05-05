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
            print("Iteration %s \tepsilon = %s \tJ(x) = %s." % (str(i+1),
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

def nlcg(criterion, n, linesearch, tol=1e-6, x0=None, maxiter=300,
         callback=None):
    """Non-linear conjugate gradient
    
    Parameters
    ----------
    criterion : function
        cost function to be minimised
    n : integer
        size of criterion input vector
    linesearch : function
        line search function, minimising the criterion along a descent
    x0 : vector
        initial guess
    maxiter : integer
        maximum number of iterations.
    callback : function
        callback function, called after each iteration

    Returns
    --------
    x : solution

    """
    if callback is None:
        callback = CallbackFactory(verbose=True, criterion=True)

    # first guess
    if x0 is None:
        x = np.zeros(n)
    else:
        x = x0.flatten()

    # tolerance
    Js, g, ng = criterion(x, gradient=True)
    J = sum(Js)
    Jnorm = J
    resid = 2 * tol

    # maxiter
    if maxiter is None:
        maxiter = x.size
    iter_ = 0

    while iter_ < maxiter and resid > tol:
        iter_ += 1

        # descent direction
        if (iter_  % 10) == 1:
            d = - g
        else:
            b = ng / ng_old
            d = - g + b * d
        ng_old = ng

        # step
        a = linesearch(d, g)

        # update
        x += a * d

        # criterion
        J_old = J
        Js, g, ng = criterion(x, gradient=True)
        J = sum(Js)
        resid = (J_old - J) / Jnorm
        callback(x)

    # define output
    if resid > tol:
        info = resid
    else:
        info = 0

    return x#, info

# To create callback functions
class CallbackFactory():
    def __init__(self, verbose=False, criterion=False):
        self.iter_ = []
        self.resid = []
        if criterion:
            self.criterion = []
        else:
            self.criterion = False
        self.verbose = verbose
    def __call__(self, x):
        import inspect
        parent_locals = inspect.stack()[1][0].f_locals
        self.iter_.append(parent_locals['iter_'])
        self.resid.append(parent_locals['resid'])
        if self.criterion is not False:
            self.criterion.append(parent_locals['J'])
        if self.verbose:
            # print header at first iteartion
            if len(self.iter_) == 1:
                header = 'Iteration \t Residual'
                if self.criterion is not False:
                    header += '\t Criterion'
                    print(header)
            # print status
            report = "\t%i \t %e" % (self.iter_[-1], self.resid[-1])
            if self.criterion is not False:
                report += '\t %e' % (self.criterion[-1])
            print(report)
