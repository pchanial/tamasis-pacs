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


def mapper_ls(tod, model, weight=None, tol=1.e-5, maxiter=300, M=None, unpacking=None, verbose=True, solver=None):
    return mapper_rls(tod, model, weight=weight, hyper=0, tol=tol, maxiter=maxiter, M=M, verbose=verbose, solver=solver, unpacking=unpacking)


#-------------------------------------------------------------------------------


def mapper_rls(tod, model, weight=None, hyper=1.0, tol=1.e-5, maxiter=300, M=None, verbose=True, solver=None, unpacking=None):

    if solver is None:
        solver = scipy.sparse.linalg.cgs

    if unpacking is None:
        unpacking = Identity()

    if weight is None:
        weight = Identity('Weight')

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
    map_naive = mapper_naive(tod, model)
    x0 = unpacking.T * map_naive

    rhs = C.packing(rhs)
    if not numpy.all(numpy.isfinite(rhs)): raise ValueError('RHS contains not finite values.')
    x0 = C.packing(x0)
    x0[numpy.isnan(x0)] = 0.
    x0[:] = 0

    class PcgCallback():
        def __init__(self):
            self.niterations = 0
            self.residuals = 0.
        def __call__(self, x):
            import inspect
            parent_locals = inspect.stack()[1][0].f_locals
            self.niterations = parent_locals['iter_'] - 1
            self.residuals = parent_locals['resid']
            if self.residuals < tol:
                self.niterations += 1
            if verbose: 
                print 'Iteration ' + str(self.niterations) + ': ' + str(self.residuals)
                

    callback = PcgCallback()
    
    time0 = time.time()
    solution, nit = solver(C, rhs, x0=x0, tol=tol, maxiter=maxiter, callback=callback, M=M0)
    if nit < 0:
        raise ValueError('Invalid input for the solver.')

    if nit > 0:
        callback.niterations += 1
        if callback.niterations != maxiter:
            print 'Warning: mapper_rls: maxiter != niter...'

    output = Map(C.unpacking(solution))

    if output.header is None:
        output.header = create_fitsheader(solution)

    if hasattr(callback, 'niterations'):
        output.header.update('niter', callback.niterations)
    output.header.update('nitermax', maxiter)
    if hasattr(callback, 'residuals'):
        output.header.update('residual', callback.residuals)
    output.header.update('tol', tol)
    output.header.update('solver', solver.__name__)

    output.header.update('time', time.time() - time0)
    output.coverage = map_naive.coverage

    return output
