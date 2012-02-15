# Copyrights 2010-2011 Pierre Chanial
from __future__ import division

import numpy as np

from pyoperators import IdentityOperator, asoperator
from pyoperators.utils.mpi import MPI

from .linalg import dot, norm2


__all__ = []

def cg(A, b, x0, tol=1.e-5, maxiter=300, callback=None, M=None, comm=None):
    """OpenMPI/MPI hybrid conjugate gradient solver with preconditioning."""

    if comm is None:
        comm = MPI.COMM_WORLD

    A = asoperator(A)
    if M is None:
        M = IdentityOperator()
    M = asoperator(M)

    maxRelError = tol**2

    n = b.size
    x = np.empty(n)
    d = np.empty(n)
    q = np.empty(n)
    r = np.empty(n)
    s = np.empty(n)
    xfinal = np.zeros(n)

    if x0 is None:
        x[:] = 0
    else:
        x[:] = x0.ravel()

    norm = norm2(b, comm=comm)
    if norm == 0:
        return xfinal, 0

    r[:] = b
    r -= A.matvec(x)
    epsilon = norm2(r, comm=comm) / norm
    minEpsilon = epsilon

    M.matvec(r, d)
    delta0 = dot(r, d, comm=comm)
    deltaNew = delta0

    for iter_ in xrange(maxiter):
        if epsilon <= maxRelError:
            break
        
        A.matvec(d, q)

        alpha = deltaNew / dot(d, q, comm=comm)
        x += alpha * d
        r -= alpha * q
        epsilon = norm2(r, comm=comm) / norm

        if callback is not None:
            resid = np.sqrt(epsilon)
            callback(x)

        if epsilon < minEpsilon:
            xfinal[:] = x
            minEpsilon = epsilon

        qnorm = np.sqrt(norm2(q, comm=comm))
        M.matvec(r, s)

        deltaOld = deltaNew

        deltaNew = dot(r, s, comm=comm)
        beta = deltaNew / deltaOld
        d *= beta
        d += s
    
    minEpsilon  = np.sqrt(minEpsilon)
    maxRelError = np.sqrt(maxRelError)

    return xfinal, int(minEpsilon > maxRelError)

def nlcg(objfunc, n, linesearch, descent_method='pr', x0=None, M=None,
         tol=1.e-6, maxiter=500, callback=None, comm=None):
    """Non-linear conjugate gradient with preconditioning.

    Parameters
    ----------
    objfunc : function
        objective function to be minimised
    n : integer
        input size of the objective function
    linesearch : function
        line search function, minimising the objective function along
        a direction
    x0 : vector
        initial guess for the solution
    descent_method : string
        method for the descent vector update. Available methods are 
            - 'fs' (Fletcher-Reeves)
            - 'pr' (Polak-Ribiere)
            - 'hs' (Hestenes-Stiefel)
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
    g = np.empty(n)
    s = np.empty(n)
    d = np.empty(n)

    if comm is None:
        comm = MPI.COMM_WORLD

    # tolerance
    f = np.array(0.)
    work = []
    objfunc(x, out=f, outwork=work)
    objfunc.D(x, out=g, inwork=work)
    if M is not None:
        s[:] = -M.matvec(g)
    else:
        s[:] = -g
    d[:] = s
    delta_new = dot(g, d, comm=comm)
    delta_0 = delta_new
    resid = np.sqrt(delta_new/delta_0)

    iter_ = 0

    while iter_ < maxiter and resid > tol:
        iter_ += 1

        # update position and gradient
        info = linesearch(objfunc, d, x, f, g)

        delta_old = delta_new
        if descent_method in ('pr', 'hs'):
            delta_mid = dot(g, s, comm=comm)
        if M is not None:
            s[:] = -M.matvec(g)
        else:
            s[:] = -g
        delta_new = dot(g, s, comm=comm)

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

        if (iter_  % 50) == 1 or b <= 0:
            d[:] = s # reset to steepest descent
        else:
            d[:] = s + b * d

        resid = np.sqrt(delta_new/delta_0)

        if callback is not None:
            callback(x)

    return x

class LineSearch(object):
    def __call__(self, d, x, f, g):
        """
        Parameters
        ----------
        d : array
            descent vector
        x : array
            origin vector (updated after call)
        f : float
            objective function value at origin (updated after call)
        g : array
            gradient at origin (updated after call)
        
        """
        raise NotImplementedError()

class QuadraticStep(LineSearch):

    def __init__(self, hypers, norms, operators, comms, unpacking=None, comm=MPI.COMM_SELF, ):
        self.hypers = hypers
        self.norms = norms
        self.operators = operators
        self.comms = comms
        self.unpacking=unpacking
        self.comm = comm
        
    def __call__(self, objfunc, d, x, f, g):
        if self.unpacking is not None:
            d = self.unpacking.matvec(d)
            g = self.unpacking.matvec(g)
        a = - 0.5 * dot(d, g, comm=self.comm) / sum([ h * n(M * d, comm=c) \
            for h, n, M, c in zip(self.hypers, self.norms, self.operators,
                                  self.comms)])
        x += a * d
        work = []
        fs = objfunc(x, out=f, outwork=work)
        objfunc.D(x, out=g, inwork=work)
        return { 'alpha' : a, 'terms' : fs }

class BacktrackingStep(LineSearch):
    def __init__(self, alpha=1., rho=0.1, beta=0.1, comm=MPI.COMM_SELF):
        self.alpha = alpha
        self.rho = rho
        self.beta = beta
        self.comm = comm

    def __call__(self, objfunc, d, x, f, g):
        x0 = x.copy()
        f0 = f.copy()
        d = d / dot(d, d)
        s0 = dot(g, d, comm=self.comm)
        a = self.alpha
        work = []

        # Armijo rule
        niterations = 0
        while True:
            niterations += 1
            x[:] = x0 + a * d
            fs = objfunc(x, out=f, outwork=work)
            if f <= f0 + self.rho * a * s0:
                break
            a *= self.beta

        objfunc.D(x, out=g, inwork=work)
        return { 'alpha' : a, 'terms' : fs, 'niterations' : niterations }
        
class StrongWolfePowellStep(LineSearch):

    def __init__(self, alpha=1, rho=0.1, sigma=0.4, alpha_min=0,
                 alpha_max=1, alpha_limit = 0.1, comm=MPI.COMM_SELF):
        self.comm = comm
        self.rho = rho
        self.sigma = sigma
        self.alpha = alpha
        self.alpha_min = alpha_min
        self.alpha_max = alpha_max
        self.alpha_limit = alpha_limit
        self.x0 = None

    def __call__(self, objfunc, d, x, f, g):

        a  = self.alpha
        a1 = self.alpha_min
        a2 = self.alpha_max
        d = d / dot(d, d)

        # allocate space to copy the origin
        if self.x0 is None or self.x0.size != x.size:
            self.x0 = np.empty(x.size)
        self.x0[:] = x

        f0 = f.copy()
        fp0 = dot(g, d, comm=self.comm)

        f1 = f0
        fp1 = fp0
        work = []
        niterations = 0
        while True:

            niterations += 1
            x[:] = self.x0 + a * d
            fs = objfunc(x, out=f, outwork=work)

            condition1 = f <= f0 + self.rho * a * fp0
            if not condition1:
                a, a2 = a1 + 0.5 * (a - a1) / (1 + (f1 - f)/((a - a1) * fp1)), a
                if abs(a - a1) < self.alpha_limit * abs(a2 - a1):
                    a = a1 + self.alpha_limit * abs(a2 - a1)
                continue

            objfunc.D(x, out=g, inwork=work)
            fp = dot(g, d, comm=self.comm)

            condition2 = abs(fp) <= self.sigma * abs(fp0)
            if not condition2:
                a, a1 = a + (a - a1) * fp / (fp1 - fp), a
                f1 = f.copy()
                fp1 = fp
                continue

            return { 'alpha' : a, 'terms' : fs, 'niterations' : niterations }
