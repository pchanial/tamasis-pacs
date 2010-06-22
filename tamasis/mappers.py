import numpy
import scipy
import time

from acquisitionmodels import asacquisitionmodel, Diagonal, DiscreteDifference, Identity
from datatypes import Map, Tod
from unit import Quantity, UnitError
from utils import create_fitsheader

__all__ = [ 'mapper_naive', 'mapper_ls', 'mapper_rls' ]

def mapper_naive(tod, model, unit=None):
    """
    Returns a naive map, i.e.: model.transpose(tod) / model.transpose(1)

    The unit of the input Time Ordered Data (TOD) is not restricted, but because
    by construction, the unit of the output map is the same, it is advisable to 
    have the TOD in unit of surface brightness.
    """
    mymap = Map(model.transpose(tod), copy=False)
    unity = tod.copy()
    unity[:] = 1.
    unity.unit = ''
    map_weights = model.transpose(unity, reusein=True)
    mymap /= map_weights
    mymap.coverage = map_weights
   
    if unit is None:
        return mymap

    newunit = Quantity(1., unit)
    newunit_si = newunit.SI._unit
    if 'pixel' in newunit_si and newunit_si['pixel'] == -1:
        h = mymap.header
        try:
            cd = numpy.array([ [h['cd1_1'],h['cd1_2']], [h['cd2_1'],h['cd2_2']] ])
        except:
            raise KeyError('The pixel size cannot be determined from the map header: The CD matrix is missing.')
        pixel_area = Quantity(abs(numpy.linalg.det(cd)), 'deg^2/pixel')
        mymap *= pixel_area
       
    mymap.unit = newunit._unit
       
    return mymap


#-------------------------------------------------------------------------------


def mapper_ls(tod, model, weight=None, unpacking=None, x0=None, tol=1.e-5, maxiter=300, M=None, solver=None, verbose=True):
    return mapper_rls(tod, model, weight=weight, unpacking=unpacking, hyper=0, x0=x0, tol=tol, maxiter=maxiter, M=M, solver=solver, verbose=verbose)


#-------------------------------------------------------------------------------


def mapper_rls(tod, model, weight=None, unpacking=None, hyper=1.0, x0=None, tol=1.e-5, maxiter=300, M=None, solver=None, verbose=True):

    if weight is None:
        weight = Identity('Weight')

    if unpacking is None:
        unpacking = Identity()

    if solver is None:
        solver = scipy.sparse.linalg.cgs

    C = model.T * weight * model

    if hyper != 0:
        dX = DiscreteDifference(axis=1)
        dY = DiscreteDifference(axis=0)
        C += hyper * ( dX.T * dX + dY.T * dY )

    C = (unpacking.T * C * unpacking).aslinearoperator()

    if M is None:
        M = Identity('Preconditioner')
    elif isinstance(M, numpy.ndarray):
        M = M.copy()
        M[~numpy.isfinite(M)] = numpy.max(M[numpy.isfinite(M)])
        M = Diagonal(M, description='Preconditioner')
    else:
        M = asacquisitionmodel(M)
    M0 = (unpacking.T * M * unpacking).aslinearoperator(C.shape)

    rhs = unpacking.T * model.T * weight * tod
    rhs = C.packing(rhs)
    if not numpy.all(numpy.isfinite(rhs)): raise ValueError('RHS contains not finite values.')

    if x0 is not None:
        x0 = C.packing * unpacking.T * x0
        x0 = C.packing(x0)
        x0[numpy.isnan(x0)] = 0.

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
            if verbose: 
                print 'Iteration ' + str(self.niterations) + ': ' + str(self.residual)
                

    callback = PcgCallback()
    
    time0 = time.time()
    solution, info = solver(C, rhs, x0=x0, tol=tol, maxiter=maxiter, callback=callback, M=M0)
    if info < 0:
        raise RuntimeError('Solver failure (code='+str(info)+' after '+str(callback.niterations)+' iterations).')

    if info > 0:
        callback.niterations += 1
        if callback.niterations != maxiter:
            print 'Warning: mapper_rls: maxiter != niter...'

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
