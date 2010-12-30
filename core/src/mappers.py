import numpy
import scipy
import time

from mpi4py import MPI
from .acquisitionmodels import asacquisitionmodel, Diagonal, DiscreteDifference, Identity, Masking, AllReduce
from .config import get_default_dtype_float
from .datatypes import Map, Tod, create_fitsheader
from .quantity import Quantity, UnitError


__all__ = [ 'mapper_naive', 'mapper_ls', 'mapper_rls' ]

def mapper_naive(tod, model, unit=None):
    """
    Returns a naive map, i.e.: model.transpose(tod) / model.transpose(1)

    The unit of the input Time Ordered Data (TOD) is not restricted, but because
    by construction, the unit of the output map is the same, it is advisable to 
    have the TOD in unit of surface brightness.
    """

    # make sure the input is a surface brightness
    if tod._unit is not None and 'detector' in tod._unit:
        tod = tod.tounit(tod.unit + ' detector / arcsec^2')
        copy = False
    elif tod._unit is not None and 'detector_reference' in tod._unit:
        tod = tod.tounit(tod.unit + ' detector_reference / arcsec^2')
        copy = False
    else:
        copy = True

    model = model * AllReduce().T
    if tod.mask is not None:
        model = Masking(tod.mask) * model

    mymap = model.T(tod)
    unity = Tod(tod, copy=copy)
    unity[:] = 1.
    map_weights = model.T(unity, reusein=True)
    map_weights.unit = ''
    old_settings = numpy.seterr(divide='ignore', invalid='ignore')
    mymap /= map_weights
    numpy.seterr(**old_settings)
    mymap.coverage = map_weights
   
    if unit is not None:
        mymap.unit = unit
    
    return mymap


#-------------------------------------------------------------------------------


def mapper_ls(tod, model, weight=None, unpacking=None, x0=None, tol=1.e-5, maxiter=300, M=None, solver=None, verbose=True, callback=None):
    return mapper_rls(tod, model, weight=weight, unpacking=unpacking, hyper=0, x0=x0, tol=tol, maxiter=maxiter, M=M, solver=solver, verbose=verbose, callback=callback)


#-------------------------------------------------------------------------------


def mapper_rls(tod, model, weight=None, unpacking=None, hyper=1.0, x0=None, tol=1.e-5, maxiter=300, M=None, solver=None, verbose=True, callback=None):

    if weight is None:
        weight = Identity('Weight')

    if solver is None:
        solver = scipy.sparse.linalg.cgs

    C = model.T * weight * model

    if hyper != 0:
        ntods = MPI.COMM_WORLD.allreduce(numpy.sum(tod.mask == 0), op=MPI.SUM)
        nmaps = C.aslinearoperator(unpacking=unpacking).shape[0]
        if MPI.COMM_WORLD.Get_rank() == 0:
            hyper = numpy.array(hyper * ntods / nmaps, dtype=get_default_dtype_float())
            dX = DiscreteDifference(axis=1)
            dY = DiscreteDifference(axis=0)
            C += hyper * ( dX.T * dX + dY.T * dY )

    C = AllReduce() * C
    C = C.aslinearoperator(unpacking=unpacking)

    if M is None:
        M = Identity('Preconditioner')
    elif isinstance(M, numpy.ndarray):
        M = M.copy()
        M[~numpy.isfinite(M)] = numpy.max(M[numpy.isfinite(M)])
        M = Diagonal(M, description='Preconditioner')
    else:
        M = asacquisitionmodel(M)
    M0 = M.aslinearoperator(C.shape, unpacking=unpacking)

    rhs = C.packing * AllReduce() * model.T * weight * tod
    if not numpy.all(numpy.isfinite(rhs)):
        raise ValueError('RHS contains not finite values.')
    if rhs.shape != C.shape[1]:
        raise ValueError("Incompatible size for RHS: '"+str(rhs.shape)+"' instead of '"+str(C.shape[1])+"'.")

    if x0 is not None:
        x0 = C.packing(x0)
        x0[numpy.isnan(x0)] = 0.
        if x0.shape != C.shape[0]:
            raise ValueError("Incompatible size for x0: '"+str(x0.shape)+"' instead of '"+str(C.shape[0])+"'.")

    class PcgCallback():
        def __init__(self):
            self.niterations = 0
            self.residual = 0.
        def __call__(self, x):
            import inspect
            parent_locals = inspect.stack()[1][0].f_locals
            self.niterations = parent_locals['iter_'] - 1
            self.residual = parent_locals['resid']
            if self.residual < tol:
                self.niterations += 1
            if verbose and MPI.COMM_WORLD.Get_rank() == 0: 
                print('Iteration ' + str(self.niterations) + ': ' + str(self.residual))

    if callback is None:
        callback = PcgCallback()
    
    time0 = time.time()
    solution, info = solver(C, rhs, x0=x0, tol=tol, maxiter=maxiter, callback=callback, M=M0)
    if info < 0:
        raise RuntimeError('Solver failure (code='+str(info)+' after '+str(callback.niterations)+' iterations).')

    if info > 0:
        callback.niterations += 1
        if callback.niterations != maxiter:
            print('Warning: mapper_rls: maxiter != niter...')

    output = Map(C.unpacking(solution))
    if output.header is None:
        output.header = create_fitsheader(solution)

    output.header.update('time', time.time() - time0)

    if hasattr(callback, 'niterations'):
        output.header.update('niter', callback.niterations)
    output.header.update('nitermax', maxiter)
    if hasattr(callback, 'residual'):
        output.header.update('residual', callback.residual)
    output.header.update('tol', tol)
    output.header.update('solver', solver.__name__)

    map_naive = mapper_naive(tod, model)
    output.coverage = map_naive.coverage

    return output
