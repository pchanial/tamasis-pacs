# Copyrights 2010-2011 Pierre Chanial
# All rights reserved
#
import cProfile
import numpy as np
import os
import scipy
import tamasisfortran as tmf
import time

from copy import copy
from mpi4py import MPI
from . import var
from .acquisitionmodels import Addition, DdTdd, Diagonal, DiscreteDifference, \
    Identity, Masking, Scalar, asacquisitionmodel
from .datatypes import Map, Tod, flatten_sliced_shape
from .linalg import dot, norm2, norm2_ellipsoid
from .quantity import Quantity, UnitError
from .solvers import cg, nlcg


def mapper_naive(tod, model, unit=None, local_mask=None):
    """
    Returns a naive map, i.e.: map = model.T(tod) / model.T(1)

    This equation is valid for a map and a Time Ordered Data (TOD) expressed as
    surface brightness. When the input Time Ordered Data (TOD) does not meet
    this requirement, this method performs a unit conversion if the input is
    a quantity per detector.

    Parameters
    ----------

    tod : Tod
        The input Time Ordered Data

    model : AcquisitionModel
        The instrument model such as tod = model(map)

    unit : string
        Output map unit. By default, the output map unit is chosen to be
        compatible with the model's input unit model.unitin (usually pixel^-1)
    """

    # make sure the input is a surface brightness
    if 'detector' in tod._unit:
        tod = tod.tounit(tod.unit + ' detector / arcsec^2')
        inplace = True
    elif 'detector_reference' in tod._unit:
        tod = tod.tounit(tod.unit + ' detector_reference / arcsec^2')
        inplace = True
    else:
        inplace = False

    mask = getattr(tod, 'mask', None)
    tod = Masking(mask)(tod, inplace=inplace)

    # model.T expects a quantity / detector, we hide our units to 
    # prevent a unit validation exception
    todunit = tod.unit
    tod.unit = ''

    # compute model.T(tod)/model.T(one)
    mymap = model.T(tod, True, True, False)
    tod[:] = 1
    map_weights = model.T(tod, True, True, False)
    old_settings = np.seterr(divide='ignore', invalid='ignore')
    tmf.divide_inplace(mymap.T, map_weights.T)
    mymap.unit = todunit
    np.seterr(**old_settings)
    mymap.coverage = map_weights
   
    if unit is not None:
        mymap.inunit(unit)
        return mymap

    # make sure that the map unit is compatible with the input unit of the model
    if len(mymap._unit) == 0:
        mymap._unit = model.unitin
        return mymap

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
                return mymap
            except:
                pass
    
    print("Warning: cannot set naive map to a unit compatible with that of th" \
          "e model '" + Quantity(1,model.unitin).unit + "'.")

    return mymap


#-------------------------------------------------------------------------------


def mapper_ls(tod, model, invntt=None, unpacking=None, x0=None, tol=1.e-5,
              maxiter=300, M=None, solver=None, verbose=True, callback=None,
              criterion=True, profile=None, comm_map=None):

    tod = _validate_tod(tod, model)

    if invntt is None:
        invntt = Identity()

    A = model.T * invntt * model
    b = (model.T * invntt)(tod, inplace=True)
    if M is not None:
        M = asacquisitionmodel(M, shapein=model.shapein, shapeout=model.shapein)

    return _solver(A, b, tod, model, invntt, hyper=0, x0=x0, tol=tol,
                   maxiter=maxiter, M=M, solver=solver, verbose=verbose,
                   callback=callback, criterion=criterion, profile=profile,
                   unpacking=unpacking, comm_map=comm_map)


#-------------------------------------------------------------------------------


def mapper_rls(tod, model, invntt=None, unpacking=None, hyper=1.0, x0=None,
               tol=1.e-5, maxiter=300, M=None, solver=None, verbose=True,
               callback=None, criterion=True, profile=None, comm_map=None):

    tod = _validate_tod(tod, model)

    if invntt is None:
        invntt = Identity()

    if comm_map is None:
        comm_map = var.comm_map

    A = model.T * invntt * model

    hyper = float(hyper)
    ntods = int(np.sum(~tod.mask)) if getattr(tod, 'mask', None) is not None \
            else tod.size
    ntods = var.comm_tod.allreduce(ntods, op=MPI.SUM)
    nmaps = A.shape[1] * comm_map.Get_size()

    npriors = len(model.shapein)
    priors = [ DiscreteDifference(axis=axis, shapein=model.shapein,
               comm=var.comm_map) for axis in range(npriors) ]
    prior = Addition([DdTdd(axis=axis, scalar=hyper * ntods / nmaps,
        comm=comm_map) for axis in range(npriors)])
    if hyper != 0 and (comm_map.Get_rank() == 0 or comm_map.Get_size() > 1):
# commenting this for now, matvec is not updated
#        A += prior
        A = A + prior

    b = (model.T * invntt)(tod)

    return _solver(A, b, tod, model, invntt, hyper=hyper, priors=priors, x0=x0,
                   tol=tol, maxiter=maxiter, M=M, solver=solver,
                   verbose=verbose, callback=callback, criterion=criterion,
                   profile=profile, unpacking=unpacking, comm_map=comm_map)


#-------------------------------------------------------------------------------


def mapper_nl(tod, model, unpacking=None, priors=[], hypers=[], norms=[],
              comms=[], x0=None, tol=1.e-6, maxiter=300, solver=None,
              descent_method='pr', M=None, verbose=True, callback=None,
              criterion=True, profile=None):

    if len(priors) == 0 and len(hypers) != 0:
        npriors = len(model.shapein)
        priors = [ DiscreteDifference(axis=axis, shapein=model.shapein,
                   comm=var.comm_map) for axis in range(npriors) ]
    else:
        npriors = len(priors)

    if np.isscalar(hypers):
        hypers = npriors * [hypers]
    elif npriors != len(hypers):
        raise ValueError('There should be one hyperparameter per prior.')

    if len(norms) == 0:
        norms = (npriors+1) * [norm2]

    if len(comms) == 0:
        comms = [var.comm_tod] + npriors * [var.comm_map]

    if isinstance(M, Diagonal):
        tmf.remove_nonfinite(M.diagonal.T)

    hypers = np.asarray(hypers, dtype=var.FLOAT_DTYPE)
    ntods = int(np.sum(~tod.mask)) if getattr(tod, 'mask', None) is not None \
            else tod.size
    ntods = var.comm_tod.allreduce(ntods, op=MPI.SUM)
    nmaps = model.shape[1] * var.comm_map.Get_size()
    hypers /= nmaps
    hypers = np.hstack([1./ntods, hypers])
    
    tod = _validate_tod(tod, model)
    y = tod.view(np.ndarray).ravel()

    if criterion:
        def criter(x, gradient=None):
            if unpacking is not None:
                x = unpacking.matvec(x)
            rs = [ model * x - y] + [p * x for p in priors ]
            Js = [ h * n(r,comm=c) for h, n, r, c in zip(hypers, norms, rs,
                   comms) ]
            if gradient is None:
                return Js
            gradient[:] = sum([h * M.T * n.D(r,comm=c) for h, M, n, r, c in \
                               zip(hypers, [model] + priors, norms, rs, comms)])
            return Js
    else:
        raise NotImplementedError()
        criter = None

    def quadratic_optimal_step(d, r):
        if unpacking is not None:
            d = unpacking.matvec(d)
            r = unpacking.matvec(r)
        a = 0.5 * dot(d, r, comm=var.comm_map) / sum([ h * n(M * d, comm=c) \
            for h, n, M, c in zip(hypers, norms, [model] + priors, comms)])
        return a

    if callback is None:
        callback = CgCallback(verbose=verbose, criterion=criter)
    
    time0 = time.time()

    solution = nlcg(criter, model.shape[1], M=M, descent_method=descent_method,
                    linesearch=quadratic_optimal_step, maxiter=maxiter, tol=tol,
                    callback=callback, comm=var.comm_map)

    time0 = time.time() - time0
    Js = criter(solution)

    if unpacking is not None:
        solution = unpacking(solution)
    output = Map(solution.reshape(model.shapein), copy=False)
    output.unit = tod.unit + ' ' + (1/Quantity(1, model.unitout)).unit + ' ' + \
                  Quantity(1, model.unitin).unit
    coverage = Map(model.T(np.ones(tod.shape), True, True, True), copy=False)
    output.coverage = coverage
    output.header = coverage.header
    output.header.update('likeliho', Js[0])
    output.header.update('criter', sum(Js))
    output.header.update('hyper', str(hypers))
    output.header.update('nsamples', ntods)
    output.header.update('npixels', model.shape[1])
    output.header.update('time', time0)
    if hasattr(callback, 'niterations'):
        output.header.update('niter', callback.niterations)
    output.header.update('maxiter', maxiter)
    if hasattr(callback, 'residual'):
        output.header.update('residual', callback.residual)
    output.header.update('tol', tol)
    output.header.update('solver', 'nlcg')
    return output


#-------------------------------------------------------------------------------


def _validate_tod(tod, model):
    # make sure that the tod unit is compatible with the model's output unit
    mask = getattr(tod, 'mask', None)
    tod = Masking(mask)(tod)

    if tod.unit == '':
        tod.unit = model.unitout
    else:
        if any([u not in tod._unit or tod._unit[u] != v
                for u,v in model.unitout.items()]):
            return UnitError("The input tod has a unit '" + tod.unit + \
                "' incompatible with that of the model '" + Quantity(1, 
                model.unitout).unit + "'.")
        
    return tod


#-------------------------------------------------------------------------------


def _solver(A, b, tod, model, invntt, priors=[], hyper=0, x0=None, tol=1.e-5,
            maxiter=300, M=None, solver=None, verbose=True, callback=None,
            criterion=True, profile=None, unpacking=None, comm_map=None):

    npriors = len(priors)

    if isinstance(M, Diagonal):
        tmf.remove_nonfinite(M.diagonal.T)

    if unpacking is not None:
        A = unpacking.T * A * unpacking
        b = unpacking.T(b)
        if x0 is not None:
            x = unpacking.T(b)
        if M is not None:
            M = unpacking.T * M * unpacking

    b = b.ravel()
    if not np.all(np.isfinite(b)):
        raise ValueError('RHS contains not finite values.')
    if b.size != A.shape[1]:
        raise ValueError("Incompatible size for RHS: '" + str(b.size) + \
                         "' instead of '" + str(A.shape[1]) + "'.")

    if comm_map is None:
        comm_map = var.comm_map

    ntods = int(np.sum(~tod.mask)) if getattr(tod, 'mask', None) is not None \
            else tod.size
    ntods = var.comm_tod.allreduce(ntods, op=MPI.SUM)

    if hyper != 0:
        hc = np.hstack([1., npriors * [hyper]]) / ntods
    else:
        hc = [ 1. / ntods ]
    norms = [norm2_ellipsoid(invntt)] + npriors * [norm2]
    comms = [var.comm_tod] + npriors * [comm_map]
    def criter(x):
        if unpacking is not None:
            x = unpacking * x
        rs = [model * x - tod.view(np.ndarray).ravel() ] + \
             [ p * x for p in priors]
        Js = [h * n(r,comm=c) for h, n, r, c in zip(hc, norms, rs, comms)]
        return Js

    if callback is None:
        if comm_map.Get_rank() == 0 and verbose:
            print('Iteration\tResiduals' + ('\tCriterion' if criterion else ''))
        callback = CgCallback(verbose=verbose, criterion=criter if criterion \
            else None)

    if solver is None:
        solver = cg

    if var.verbose or profile:
        print('')
        print('Model:')
        print(A)
        if M is not None:
            print('Preconditioner:')
            print(M)

    time0 = time.time()
    if profile is not None:

        def run():
            try:
                solution, info = solver(A, b, x0=x0, tol=tol, maxiter=maxiter,
                                        M=M, comm=comm_map)
            except TypeError:
                solution, info = solver(A, b, x0=x0, tol=tol, maxiter=maxiter,
                                        M=M)
            if info < 0:
                print('Solver failure: info='+str(info))

        cProfile.runctx('run()', globals(), locals(), profile+'.prof')
        print('Profile time: '+str(time.time()-time0))
        os.system('python -m gprof2dot -f pstats -o ' + profile +'.dot ' + \
                  profile + '.prof')
        os.system('dot -Tpng ' + profile + '.dot > ' + profile)
        os.system('rm -f ' + profile + '.prof' + ' ' + profile + '.dot')
        return None

    try:
        solution, info = solver(A, b, x0=x0, tol=tol, maxiter=maxiter, M=M,
                                callback=callback, comm=comm_map)
    except TypeError:
        solution, info = solver(A, b, x0=x0, tol=tol, maxiter=maxiter, M=M,
                                callback=callback)
        

    time0 = time.time() - time0
    if info < 0:
        raise RuntimeError('Solver failure (code=' + str(info) + ' after ' + \
                           str(callback.niterations) + ' iterations).')

    if info > 0:
        print('Warning: Solver reached maximum number of iterations without r' \
              'eaching specified tolerance.')

    Js = criter(solution)

    if unpacking is not None:
        solution = unpacking(solution)

    output = Map(solution.reshape(model.shapein), copy=False)

    output.unit = tod.unit + ' ' + (1/Quantity(1, model.unitout)).unit + ' ' + \
                  Quantity(1, model.unitin).unit
    coverage = Map(model.T(np.ones(tod.shape), True, True, True), copy=False)
    output.coverage = coverage
    output.header = coverage.header
    output.header.update('likeliho', Js[0])
    output.header.update('criter', sum(Js))
    output.header.update('hyper', hyper)
    output.header.update('nsamples', ntods)
    output.header.update('npixels', A.shape[1])
    output.header.update('time', time0)
    if hasattr(callback, 'niterations'):
        output.header.update('niter', callback.niterations)
    output.header.update('maxiter', maxiter)
    if hasattr(callback, 'residual'):
        output.header.update('residual', callback.residual)
    output.header.update('tol', tol)
    output.header.update('solver', solver.__name__)

    return output

class CgCallback():
    def __init__(self, verbose=True, criterion=None, comm=None):
        self.niterations = 0
        self.residual = 0.
        self.verbose = verbose
        self.criterion = criterion
        self.comm = comm or MPI.COMM_WORLD
    def __call__(self, x):
        import inspect
        parent_locals = inspect.stack()[1][0].f_locals
        self.niterations = parent_locals['iter_'] - 1
        self.residual = parent_locals['resid']
        if self.niterations == 0 and self.comm.Get_rank() == 0:
            print('Iteration\tResiduals' + ('\tCriterion' if self.criterion \
                  else ''))
        if self.verbose: 
            if self.criterion is not None:
                Js = self.criterion(x)
                Jstr = repr(sum(Js))
                if len(Js) > 1:
                    Jstr += ' ' + str(Js)
            else:
                Jstr = ''
            if self.comm.Get_rank() == 0:
                print("%s\t%s\t%s" % (str(self.niterations).rjust(9),
                                      str(self.residual), Jstr))

