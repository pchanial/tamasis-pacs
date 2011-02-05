import numpy as np
import scipy
import time

from mpi4py import MPI
from . import var
from .acquisitionmodels import Diagonal, DiscreteDifference, Identity, Masking,\
     AllReduce, Reshaping
from .datatypes import Map, Tod, create_fitsheader
from .quantity import Quantity, UnitError


__all__ = [ 'mapper_naive', 'mapper_ls', 'mapper_rls' ]

def mapper_naive(tod, model, unit=None):
    """
    Returns a naive map, i.e.: model.T(tod) / model.T(1)

    The unit of the input Time Ordered Data (TOD) is not restricted, but because
    by construction, the unit of the output map is the same, it is advisable to 
    have the TOD in unit of surface brightness.
    """

    # make sure the input is a surface brightness
    if 'detector' in tod._unit:
        tod = tod.tounit(tod.unit + ' detector / arcsec^2')
        copy = False
    elif 'detector_reference' in tod._unit:
        tod = tod.tounit(tod.unit + ' detector_reference / arcsec^2')
        copy = False
    else:
        copy = True

    model = model * AllReduce().T
    if tod.mask is not None:
        model = Masking(tod.mask) * model

    # model.T expects a quantity / detector, we hide our units to 
    # prevent a unit validation exception
    tod_ = tod.view()
    tod_.unit = ''
    mymap = model.T(tod_)

    unity = Tod(tod_, copy=copy)
    unity[:] = 1.
    map_weights = model.T(unity, reusein=True)
    old_settings = np.seterr(divide='ignore', invalid='ignore')
    mymap /= map_weights
    mymap.unit = tod.unit

    np.seterr(**old_settings)
    mymap.coverage = map_weights
   
    if unit is not None:
        mymap.inunit(unit)
        return mymap

    # make sure that the map unit is compatible with the input unit of the model
    if all([u in mymap._unit and mymap._unit[u] == v \
            for u,v in model.unitin.items()]):
        return mymap

    for u, v in zip(['sr', 'rad', 'deg', 'arcmin', "'", 'arcsec', '"',
                     'pixel_reference'], [1,2,2,2,2,2,2,1]):
        if u in mymap._unit and mymap._unit[u] == -v:
            newunit = mymap.unit + ' ' + u + '^' + str(v) + ' ' + \
                      Quantity(1,model.unitin).unit
            try:
                mymap.inunit(newunit)
            except:
                pass
            return mymap
    
    print("Warning: cannot set naive map to a unit compatible with that of th" \
          "e model '" + Quantity(1,model.unitin).unit + "'.")

    return mymap


#-------------------------------------------------------------------------------


def mapper_ls(tod, model, weight=None, unpacking=None, x0=None, tol=1.e-5,
              maxiter=300, M=None, solver=None, verbose=True, callback=None):
    return mapper_rls(tod, model, weight=weight, unpacking=unpacking, hyper=0,
                      x0=x0, tol=tol, maxiter=maxiter, M=M, solver=solver,
                      verbose=verbose, callback=callback)


#-------------------------------------------------------------------------------


def mapper_rls(tod, model, weight=None, unpacking=None, hyper=1.0, x0=None,
               tol=1.e-5, maxiter=300, M=None, solver=None, verbose=True,
               callback=None):

    
    # make sure that the tod unit is compatible with the model's output unit
    if tod.unit == '':
        tod = tod.view()
        tod.unit = model.unitout
    else:
        if any([u not in tod._unit or tod._unit[u] != v
                for u,v in model.unitout.items()]):
            return UnitError("The input tod has a unit '" + tod.unit + \
                "' incompatible with that of the model '" + Quantity(1, 
                model.unitout).unit + "'.")

    if weight is None:
        weight = Identity(description='Weight')

    if unpacking is None:
        unpacking = Reshaping(shapein=model.shape[1], shapeout=model.shapein)

    if solver is None:
        solver = scipy.sparse.linalg.cgs

    C = model.T * weight * model

    if hyper != 0:
        ntods = var.mpi_comm.allreduce(np.sum(tod.mask == 0), op=MPI.SUM)
        nmaps = unpacking.shape[1]
        if var.mpi_comm.Get_rank() == 0:
            hyper = np.array(hyper * ntods / nmaps, dtype=var.FLOAT_DTYPE)
            dX = DiscreteDifference(axis=1)
            dY = DiscreteDifference(axis=0)
            C += hyper * ( dX.T * dX + dY.T * dY )

    C   = AllReduce() * unpacking.T * C * unpacking
    rhs = AllReduce() * unpacking.T * model.T * weight * tod

    if not np.all(np.isfinite(rhs)):
        raise ValueError('RHS contains not finite values.')
    if rhs.shape != (C.shape[1],):
        raise ValueError("Incompatible size for RHS: '" + str(rhs.shape) + \
                         "' instead of '" + str(C.shape[1]) + "'.")

    if M is not None:
        if isinstance(M, np.ndarray):
            M = M.copy()
            M[~np.isfinite(M)] = np.max(M[np.isfinite(M)])
            M = Diagonal(M, description='preconditioner')
        M = unpacking.T * M

    if x0 is not None:
        x0 = unpacking.T * x0
        x0[np.isnan(x0)] = 0.
        if x0.shape != (C.shape[1],):
            raise ValueError("Incompatible size for x0: '" + str(x0.shape) + \
                             "' instead of '" + str(C.shape[1]) + "'.")

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
            if verbose and var.mpi_comm.Get_rank() == 0: 
                print('Iteration ' + str(self.niterations) + ': ' + \
                      str(self.residual))

    if callback is None:
        callback = PcgCallback()

    time0 = time.time()
    solution, info = solver(C, rhs, x0=x0, tol=tol, maxiter=maxiter,
                            callback=callback, M=M)
    if info < 0:
        raise RuntimeError('Solver failure (code=' + str(info) + ' after ' + \
                           str(callback.niterations) + ' iterations).')

    if info > 0:
        callback.niterations += 1
        if callback.niterations != maxiter:
            print('Warning: mapper_rls: maxiter != niter...')

    output = Map(unpacking(solution, reusein=True, reuseout=True), copy=False)
    output.unit = tod.unit + ' ' + (1/Quantity(1, model.unitout)).unit + ' ' + \
                  Quantity(1, model.unitin).unit
    map_naive = mapper_naive(tod, model)

    output.header = map_naive.header
    output.coverage = map_naive.coverage

    output.header.update('time', time.time() - time0)

    if hasattr(callback, 'niterations'):
        output.header.update('niter', callback.niterations)
    output.header.update('maxiter', maxiter)
    if hasattr(callback, 'residual'):
        output.header.update('residual', callback.residual)
    output.header.update('tol', tol)
    output.header.update('solver', solver.__name__)

    return output
