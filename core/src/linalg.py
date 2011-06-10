import numpy as np
import tamasisfortran as tmf
from mpi4py import MPI
from scipy.sparse.linalg import aslinearoperator
from . import var

__all__ = []


def sum(x, comm=None):
    """Parallelised sum"""
    x = np.array(x, copy=False, order='c', dtype=var.FLOAT_DTYPE)
    if comm is None:
        comm = MPI.COMM_WORLD
    output = np.array(tmf.sum(x.T))
    comm.Allreduce(MPI.IN_PLACE, [output, MPI.DOUBLE], op=MPI.SUM)
    return float(output)


#-------------------------------------------------------------------------------


def dot(x1, x2, comm=None):
    """Parallelised dot product"""
    x1 = np.array(x1, copy=False, order='c', dtype=var.FLOAT_DTYPE)
    x2 = np.array(x2, copy=False, order='c', dtype=var.FLOAT_DTYPE)
    if comm is None:
        comm = MPI.COMM_WORLD
    output = np.array(tmf.dot(x1.T, x2.T))
    comm.Allreduce(MPI.IN_PLACE, [output, MPI.DOUBLE], op=MPI.SUM)
    return float(output)


#-------------------------------------------------------------------------------


class Function:
    def __init__(self, f, df=None):
        self.__call__ = f
        if df is not None:
            self.D = Function(df)
    def __call__(x, out=None, inwork=None, outwork=None, comm=None):
        raise NotImplementedError()
    def D(x, out=None, inwork=None, outwork=None, comm=None):
        raise NotImplementedError()
    def set_out(self, out, value):
        if out is not None:
            if isinstance(value, (list, tuple)):
                if out.size != value[0].size:
                    raise ValueError('Incompatible size.')
                out.flat = sum(value)
            else:
                if out.size != value.size:
                    raise ValueError('Incompatible size.')
                out.flat = value
        return value
    def set_outwork(self, outwork, value):
        if outwork is None:
            return value
        for i in xrange(len(outwork)):
            outwork.pop()
        outwork.extend(value)
        return outwork
        


#-------------------------------------------------------------------------------


def gradient(f):
    """Parallelised three-point estimation of the gradient"""
    def df(x, out=None, comm=None):
        x = np.array(x, copy=False)
        if out is None:
            out = np.empty(x.size)
        elif len(out.shape) == 0:
            out.shape = (1,)
        if comm is None:
            comm = MPI.COMM_WORLD
        dx = np.sqrt(np.finfo(float).eps) * x
        rank = comm.Get_rank()
        size = comm.Get_size()
        for r in range(size):
            n = comm.bcast(x.size, root=r)
            if r == rank:
                x_dx = x.copy()
                for i in xrange(n):
                    x_dx[i] = x[i] + dx[i]
                    try:
                        f2 = f(x_dx, comm=comm)
                    except TypeError:
                        f2 = f(x_dx)
                    x_dx[i] = x[i] - dx[i]
                    try:
                        f1 = f(x_dx, comm=comm)
                    except TypeError:
                        f1 = f(x_dx)
                    out[i] = (f2-f1) / (2*dx[i])
                    x_dx[i] = x[i]
            else:
                for i in xrange(n):
                    try:
                        f(x, comm=comm)
                        f(x, comm=comm)
                    except TypeError:
                        f(x)
                        f(x)
        out.shape = x.shape
        return out
    return df


#-------------------------------------------------------------------------------

def f(x, out=None, inwork=None, outwork=None, comm=None):
    """Parallelised L-1 norm"""
    x = np.ascontiguousarray(x)
    if out is None:
        out = np.empty(())
    if comm is None:
        comm = MPI.COMM_WORLD
    out.flat = tmf.norm_l1(x.T, comm.py2f())
    return float(out)

def df(x, out=None, inwork=None, outwork=None, comm=None):
    x = np.ascontiguousarray(x)
    if out is None:
        out = np.empty(x.shape)
    tmf.dnorm_l1(x.T, out.T)
    return out

norm_l1 = Function(f, df)

#-------------------------------------------------------------------------------


def f(x, out=None, inwork=None, outwork=None, comm=None):
    """Parallelised euclidian norm"""
    x = np.array(x, copy=False, order='c', dtype=var.FLOAT_DTYPE)
    if out is None:
        out = np.empty(())
    if comm is None:
        comm = MPI.COMM_WORLD
    out.flat = np.sqrt(tmf.norm2(x.T, comm.py2f()))
    return float(out)

def df(x, out=None, inwork=None, outwork=None, comm=None):
    x = np.array(x, copy=False, order='c', dtype=var.FLOAT_DTYPE)
    if out is None:
        out = np.empty(x.shape)
    if comm is None:
        comm = MPI.COMM_WORLD
    tmf.dnorm_l2(x.T, comm.py2f(), out.T)
    return out

norm_l2 = Function(f, df)


#-------------------------------------------------------------------------------


class norm_lp(Function):
    """Parallelised L-p norm"""
    def __init__(self, p):
        def f(x, out=None, inwork=None, outwork=None, comm=None):
            x = np.array(x, copy=False, order='c', dtype=var.FLOAT_DTYPE)
            if out is None:
                out = np.empty(())
            if comm is None:
                comm = MPI.COMM_WORLD
            out.flat = tmf.normp(x.T, p, comm.py2f())**(1./p)
            return float(out)
        def df(x, out=None, inwork=None, outwork=None, comm=None):
            x = np.array(x, copy=False, order='c', dtype=var.FLOAT_DTYPE)
            if out is None:
                out = np.empty(x.shape)
            if comm is None:
                comm = MPI.COMM_WORLD
            tmf.dnorm_lp(x.T, p, comm.py2f(), out.T)
            return out
        self.__call__ = f
        self.D = Function(df)


#-------------------------------------------------------------------------------


def f(x, out=None, inwork=None, outwork=None, comm=None):
    """Parallelised L-inf norm"""
    x = np.array(x, copy=False, order='c', dtype=var.FLOAT_DTYPE)
    if out is None:
        out = np.empty(())
    if comm is None:
        comm = MPI.COMM_WORLD
    out.flat = tmf.norm_linf(x.T, comm.py2f())
    return float(out)

norm_linf = Function(f)


#-------------------------------------------------------------------------------


class norm_huber(Function):
    """Parallelised Huber's norm"""
    def __init__(self, delta):
        if delta < 0:
            raise ValueError("Huber's norm delta must be positive.")
        if delta == 0:
            self.__call__ = norm2.__call__
            self.D = norm2.D
            return
        def f(x, out=None, inwork=None, outwork=None, comm=None):
            x = np.array(x, copy=False, order='c', dtype=var.FLOAT_DTYPE)
            if out is None:
                out = np.empty(())
            if comm is None:
                comm = MPI.COMM_WORLD
            out.flat = tmf.norm_huber(x.T, delta, comm.py2f())
            return float(out)
        def df(x, out=None, inwork=None, outwork=None, comm=None):
            x = np.array(x, copy=False, order='c', dtype=var.FLOAT_DTYPE)
            if out is None:
                out = np.empty(x.shape)
            tmf.dnorm_huber(x.T, delta, out.T)
            return out

        self.__call__ = f
        self.D = Function(df)
        

#-------------------------------------------------------------------------------


def f(x, out=None, inwork=None, outwork=None, comm=None):
    """Parallelised squared euclidian norm"""
    x = np.array(x, copy=False, order='c', dtype=var.FLOAT_DTYPE)
    if out is None:
        out = np.empty(())
    if comm is None:
        comm = MPI.COMM_WORLD
    out.flat = tmf.norm2(x.T, comm.py2f())
    return float(out)

def df(x, out=None, inwork=None, outwork=None, comm=None):
    x = np.array(x, copy=False, order='c', dtype=var.FLOAT_DTYPE)
    if out is None:
        out = np.empty(x.shape)
    tmf.dnorm2(x.T, out.T)
    return out

norm2 = Function(f, df)


#-------------------------------------------------------------------------------


class normp(Function):
    """Parallelised L-p norm"""
    def __init__(self, p):
        def f(x, out=None, inwork=None, outwork=None, comm=None):
            x = np.array(x, copy=False, order='c', dtype=var.FLOAT_DTYPE)
            if out is None:
                out = np.empty(())
            if comm is None:
                comm = MPI.COMM_WORLD
            out.flat = tmf.normp(x.T, p, comm.py2f())
            return float(out)
        def df(x, out=None, inwork=None, outwork=None, comm=None):
            x = np.array(x, copy=False, order='c', dtype=var.FLOAT_DTYPE)
            if out is None:
                out = np.empty(x.shape)
            tmf.dnormp(x.T, p, out.T)
            return out
        self.__call__ = f
        self.D = Function(df)


#-------------------------------------------------------------------------------


class norm2_ellipsoid(Function):
    """Ellipsoid norm

    Returns x^T A x, where A is a definite positive symmetric matrix
    """
    def __init__(self, A):
        A = aslinearoperator(A)
        def f(x, out=None, inwork=None, outwork=None, comm=None):
            if out is None:
                out = np.empty(())
            if comm is None:
                comm = MPI.COMM_WORLD
            
            if outwork is not None:
                if len(outwork) == 0:
                    outwork.append(A.matvec(x.ravel()))
                else:
                    outwork[0][:] = A.matvec(x.ravel())
                Ax = outwork[0]
            else:
                Ax = A.matvec(x.ravel())
            out.flat = dot(x, Ax, comm=comm)
            return float(out)
        def df(x, out=None, inwork=None, outwork=None, comm=None):
            if comm is None:
                comm = MPI.COMM_WORLD
            if inwork is not None:
                Ax = inwork[0]
            else:
                Ax = A.matvec(x.ravel())
            if out is not None:
                out[:] = 2 * Ax
                return out
            return 2 * Ax
        self.__call__ = f
        self.D = Function(df)

del f, df
