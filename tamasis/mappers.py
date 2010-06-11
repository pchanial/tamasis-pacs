import numpy
import scipy
import time

from acquisitionmodels import Diagonal, DiscreteDifference, Identity
from datatypes import Map, Tod
from unit import Quantity, UnitError

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


def mapper_ls(tod, model, weight=None, tol=1.e-6, maxiter=300, M=None, unpacking=None, verbose=True, M0=None):
    return mapper_rls(tod, model, weight=weight, hyper=0, tol=tol, maxiter=maxiter, M=M, verbose=verbose, M0=M0)


#-------------------------------------------------------------------------------


def mapper_rls(tod, model, weight=None, hyper=1.0, tol=1.e-6, maxiter=300, M=None, verbose=True, M0=None):

    if weight is None:
        weight = Identity('Weight')

    C = model.T * weight * model

    if hyper != 0:
        dX = DiscreteDifference(axis=1)
        dY = DiscreteDifference(axis=0)
        C += hyper * ( dX.T * dX + dY.T * dY )

    if M is None:
        M = Identity('Preconditioner')
    elif isinstance(M, numpy.ndarray):
        M = Diagonal(M, description='Preconditioner')
    else:
        M = asacquisitionmodel(M)

    C   = M * C
    rhs = (M * model.T * weight)(tod)

    operator = C.aslinearoperator()
    rhs = operator.packing(rhs)
    if numpy.any(numpy.isnan(rhs)) or numpy.any(numpy.isinf(rhs)): raise ValueError('RHS contains not finite values.')

    map_naive = mapper_naive(tod, model)
    x0 = operator.packing(map_naive)
    x0[numpy.isnan(x0)] = 0.

    if M0 is True:
        M0 = operator.packing(1./numpy.maximum(0.01,map_naive.coverage))
        M0[numpy.isfinite(M0) == False] = 1./0.01
        M0 = scipy.sparse.dia_matrix((M0,0), shape = 2*(M0.size,))

    if verbose:
        def callback(x):
            import inspect
            parent_locals = inspect.stack()[1][0].f_locals
            iter_ = parent_locals['iter_']
            resid = parent_locals['resid']
            if resid > parent_locals['tol']:
                iter_ -= 1
                print 'Iteration ' + str(iter_) + ': ' + str(resid)

    else:
        callback = None
    
    time0 = time.time()
    solution, nit = scipy.sparse.linalg.cgs(operator, rhs, x0=x0, tol=tol, maxiter=maxiter, callback=callback, M=M0)
    output = Map(operator.unpacking(solution))
    output.time = time.time() - time0
    output.coverage = map_naive.coverage

    return output
