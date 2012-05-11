# Copyrights 2010-2011 Pierre Chanial
# All rights reserved
#
from __future__ import division

import cProfile
import numpy as np
import os
import time

from pyoperators import DiagonalOperator, IdentityOperator, MaskOperator
from pyoperators.utils.mpi import MPI

from . import var
from .acquisitionmodels import DiscreteDifferenceOperator
from .datatypes import Map, Tod
from .linalg import Function, norm2, norm2_ellipsoid
from .processing import filter_nonfinite
from .solvers import cg, nlcg, QuadraticStep
from .wcsutils import create_fitsheader

__all__ = [ 'mapper_naive',
            'mapper_ls',
            'mapper_rls',
            'mapper_nl',
]

def mapper_naive(tod, model, unit=None):
    """
    Returns a naive map, i.e.: map = model.T(tod) / model.T(1)

    This equation is valid for a map and a Time Ordered Data (TOD) expressed as
    a surface brightness, so when the TOD does not meet this requirement and
    is a quantity per detector, a unit conversion is attempted.

    Parameters
    ----------

    tod : Tod
        The input Time Ordered Data

    model : Operator
        The instrument model such as tod = model(map)

    unit : string
        Output map unit. By default, the output map unit is chosen to be
        compatible with the model (usually pixel^-1)
    """

    # apply mask
    if hasattr(tod, 'mask') and tod.mask is not None:
        mask = MaskOperator(tod.mask)
    else:
        mask = IdentityOperator()
    tod = mask(tod)

    # get tod units
    if not hasattr(tod, '_unit') or len(tod._unit) == 0:
        attr = {'_unit' : {'?' : 1.}}
        model.propagate_attributes(None, attr)
        u = getattr(attr, '_unit', {})
        if 'detector' in u and u['detector'] == -1:
            u = {'detector' : -1.}
        elif u == {'?':1.}:
            u = {}
        elif len(u) > 1:
            raise ValueError('The timeline units are not known and cannot be in'
                             'ferred from the model.')
        tod_du = getattr(attr, '_derived_units', {})
        tod = Tod(tod.magnitude, unit=u, derived_units=tod_du, copy=False)
    else:
        attr = {'_unit' : {'?':1}}
        model.T.propagate_attributes(None, attr)
        u = attr['_unit']
        if 'detector' not in tod._unit and 'detector' in u and u['detector']==1:
            raise ValueError("The model is incompatible with input units '{0}'"\
                             .format(tod.unit))

    tod_unit = tod._unit
    tod_du = tod._derived_units

    # make sure the input is a surface brightness
    if 'detector' in tod._unit:
        tod.inunit(tod.unit + ' detector / arcsec^2')
    elif 'detector_reference' in tod._unit:
        tod.inunit(tod.unit + ' detector_reference / arcsec^2')

    # compute model.T(tod)/model.T(one)
    mymap = model.T(tod.magnitude)
    tod[...] = 1
    mask(tod, tod)
    map_weights = model.T(tod.magnitude)
    old_settings = np.seterr(divide='ignore', invalid='ignore')
    mymap /= map_weights
    mymap.unit = tod.unit
    np.seterr(**old_settings)
    mymap.coverage = Map(map_weights.magnitude, header=mymap.header, copy=False)
   
    if unit is not None:
        mymap.inunit(unit)
        return mymap

    # set map units according to model
    attr = {'_unit' : tod_unit, '_derived_units' : tod_du}
    model.T.propagate_attributes(None, attr)
    if '_derived_units' in attr:
        mymap.derived_units = attr['_derived_units']
    if '_unit' in attr:
        mymap.inunit(attr['_unit'])

    return mymap


#-------------------------------------------------------------------------------


def mapper_ls(y, H, invntt=None, unpacking=None, x0=None, tol=1.e-5,
              maxiter=300, M=None, solver=None, verbose=True, callback=None,
              criterion=True, profile=None):
    """
    Solve the linear equation
        y = H(x)
    where H is an Operator representing the  acquisition, and y is the observed
    time series.
    x is the least square solution as given by:
        x = argmin (Hx-y)^T N^-1 (Hx-y)
    or:
        x = (H^T N^-1 H)^-1 H^T N^-1 y

    """
    return mapper_rls(y, H, invntt=invntt, unpacking=unpacking, x0=x0,
                      tol=tol, maxiter=maxiter, M=M, solver=solver, hyper=0,
                      verbose=verbose, callback=callback, criterion=criterion,
                      profile=profile)


#-------------------------------------------------------------------------------


def mapper_rls(y, H, invntt=None, unpacking=None, hyper=1.0, x0=None,
               tol=1.e-5, maxiter=300, M=None, solver=None, verbose=True,
               callback=None, criterion=True, profile=None):
    """
    Solve the linear equation
        y = H(x)
    where H is an Operator representing the  acquisition, and y is the observed
    time series.
    x is the regularised least square solution as given by:
        x = argmin (Hx-y)^T N^-1 (Hx-y) + h ||D(x)||^2
    or:
        x = (H^T N^-1 H + h D1^T D1 + h D2^T D2)^-1 H^T N^-1 y

    """
    comm_map = H.commin or MPI.COMM_WORLD
    comm_tod = H.commout or comm_map

    tod = _validate_tod(y)

    ntods_ = int(np.sum(~tod.mask)) if getattr(tod, 'mask', None) is not None \
            else tod.size
    nmaps_ = unpacking.shape[1] if unpacking is not None else None
    if nmaps_ is None:
        nmaps_ = H.shape[1]
    if nmaps_ is None:
        raise ValueError('The model H has not an explicit input shape.')
    ntods = comm_tod.allreduce(ntods_)
    nmaps = comm_map.allreduce(nmaps_)

    # get A
    if invntt is None:
        invntt = IdentityOperator()

    A = H.T * invntt * H

    npriors = len(H.shapein)
    priors = [DiscreteDifferenceOperator(axis=axis, shapein=H.shapein,
              commin=comm_map) for axis in range(npriors)]
    if comm_map.rank == 0 or comm_map.size > 1:
        A += sum((hyper * ntods / nmaps) * p.T * p for p in priors)

    # get b
    b = (H.T * invntt)(tod)
    if not np.all(np.isfinite(b)):
        raise ValueError('RHS contains non-finite values.')
    if b.size != A.shape[1]:
        raise ValueError("Incompatible size for RHS: '" + str(b.size) +
                         "' instead of '" + str(A.shape[1]) + "'.")
    if np.min(b) == np.max(b) == 0:
        print('Warning: in equation Ax=b, b is zero.')

    # unpack input
    if unpacking is None:
        unpacking = IdentityOperator()
    A = unpacking.T * A * unpacking
    b = unpacking.T(b).ravel()
    if x0 is not None:
        x0 = unpacking.T(x0).ravel()
    if M is not None:
        if isinstance(M, DiagonalOperator):
            filter_nonfinite(M.data, out=M.data)
        M = unpacking.T * M * unpacking
    H_ = H * unpacking
    priors = [p * unpacking for p in priors]

    # criterion
    if hyper != 0:
        hc = np.hstack([1, npriors * [hyper]]) / ntods
    else:
        hc = [ 1 / ntods ]
    norms = [norm2_ellipsoid(invntt)] + npriors * [norm2]
    comms = [comm_tod] + npriors * [comm_map]

    def criter(x):
        rs = [H_*x - tod.view(np.ndarray).ravel()] + [p*x for p in priors]
        Js = [h * n(r,comm=c) for h, n, r, c in zip(hc, norms, rs, comms)]
        return Js

    if callback is None:
        callback = CgCallback(verbose=verbose, objfunc=criter if criterion
                              else None)

    if solver is None:
        solver = cg

    if (verbose or profile) and comm_map.rank == 0:
        print('')
        print('H.T * N^-1 * H:')
        print(repr(A))
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
        os.system('python -m gprof2dot -f pstats -o ' + profile +'.dot ' +
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
        raise RuntimeError('Solver failure (code=' + str(info) + ' after ' +
                           str(callback.niterations) + ' iterations).')

    if info > 0 and comm_map.rank == 0:
        print('Warning: Solver reached maximum number of iterations without rea'
              'ching specified tolerance.')

    Js = criter(solution)

    if isinstance(unpacking, IdentityOperator):
        solution = solution.reshape(H.shapein)
    else:
        solution = unpacking(solution)
 
    tod[...] = 1
    coverage = H.T(tod)
    unit = getattr(coverage, 'unit', None)
    derived_units = getattr(coverage, 'derived_units', None)
    coverage = coverage.view(Map)

    header = getattr(coverage, 'header', None)
    if header is None:
        header = create_fitsheader(fromdata=coverage)
    header.update('likeliho', Js[0])
    header.update('criter', sum(Js))
    header.update('hyper', hyper)
    header.update('nsamples', ntods)
    header.update('npixels', nmaps)
    header.update('time', time0)
    if hasattr(callback, 'niterations'):
        header.update('niter', callback.niterations)
    header.update('maxiter', maxiter)
    if hasattr(callback, 'residual'):
        header.update('residual', callback.residual)
    header.update('tol', tol)
    header.update('solver', solver.__name__)

    output = Map(solution,
                 header=header,
                 coverage=coverage,
                 unit=unit,
                 derived_units=derived_units,
                 copy=False)

    return output


#-------------------------------------------------------------------------------


def mapper_nl(tod, model, unpacking=None, priors=[], hypers=[], norms=[],
              comms=[], x0=None, tol=1.e-6, maxiter=300, solver=None,
              linesearch=None, descent_method='pr', M=None, verbose=True,
              callback=None, profile=None):

    comm_map = model.commin or MPI.COMM_WORLD
    comm_tod = model.commout or comm_map

    tod = _validate_tod(tod)

    if len(priors) == 0 and len(hypers) != 0:
        npriors = len(model.shapein)
        priors = [ DiscreteDifferenceOperator(axis=axis, shapein=model.shapein,
                   commin=comm_map) for axis in range(npriors) ]
    else:
        npriors = len(priors)

    if np.isscalar(hypers):
        hypers = npriors * [hypers]
    elif npriors != len(hypers):
        raise ValueError('There should be one hyperparameter per prior.')

    if len(norms) == 0:
        norms = (npriors+1) * [norm2]

    if len(comms) == 0:
        comms = [comm_tod] + npriors * [comm_map]

    if isinstance(M, DiagonalOperator):
        filter_nonfinite(M.data, out=M.data)

    hypers = np.asarray(hypers, dtype=var.FLOAT_DTYPE)
    ntods = int(np.sum(~tod.mask)) if getattr(tod, 'mask', None) is not None \
            else tod.size
    ntods = comm_tod.allreduce(ntods)
    nmaps = comm_map.allreduce(model.shape[1])
    hypers /= nmaps
    hypers = np.hstack([1/ntods, hypers])
    
    y = tod.view(np.ndarray).ravel()

    class ObjectiveFunction(Function):
        def __init__(self):
            pass
        def __call__(self, x, inwork=None, outwork=None, out=None):
            outwork = self.set_outwork(outwork, self.get_work(x))
            return self.set_out(out, [ h * n(r,comm=c) for h, n, r, c in \
                zip(hypers, norms, outwork, comms) ])
        def D(self, x, inwork=None, outwork=None, out=None):
            inwork = self.get_work(x, inwork)
            return self.set_out(out, sum([h * M.T * n.D(r,comm=c) for h, M, n,
                r, c in zip(hypers, [model] + priors, norms, inwork, comms)]))
        def get_work(self, x, work=None):
            if work is not None:
                return work
            if unpacking is not None:
                x = unpacking.matvec(x)
            return [ model * x - y] + [p * x for p in priors ]

    objfunc = ObjectiveFunction()

    if linesearch is None:
        linesearch = QuadraticStep(hypers, norms, [model] + priors, comms,
                                   unpacking, comm_map)

    if callback is None:
        callback = CgCallback(verbose=verbose)
    
    time0 = time.time()

    solution = nlcg(objfunc, model.shape[1], M=M, maxiter=maxiter, tol=tol,
                    descent_method=descent_method,  linesearch=linesearch,
                    callback=callback, comm=comm_map)

    time0 = time.time() - time0
    Js = objfunc(solution)

    solution = unpacking(solution)

    tod[...] = 1
    coverage = model.T(tod)
    header = getattr(coverage, 'header', None)
    if header is None:
        header = create_fitsheader(fromdata=coverage)
    unit = getattr(coverage, 'unit', None)
    derived_units = getattr(coverage, 'derived_units', None)
    coverage = Map(coverage.magnitude, header=header, copy=False)

    header.update('likeliho', Js[0])
    header.update('criter', sum(Js))
    header.update('hyper', str(hypers))
    header.update('nsamples', ntods)
    header.update('npixels', nmaps)
    header.update('time', time0)
    if hasattr(callback, 'niterations'):
        header.update('niter', callback.niterations)
    header.update('maxiter', maxiter)
    if hasattr(callback, 'residual'):
        header.update('residual', callback.residual)
    header.update('tol', tol)
    header.update('solver', 'nlcg')

    output = Map(solution.reshape(model.shapein),
                 header=header,
                 coverage=coverage,
                 unit=unit,
                 derived_units=derived_units,
                 copy=False)

    return output


#-------------------------------------------------------------------------------


def _validate_tod(tod):
    # make sure that the tod is masked
    if hasattr(tod, 'mask') and tod.mask is not None:
        return MaskOperator(tod.mask)(tod)
    return tod.copy()


#-------------------------------------------------------------------------------


class CgCallback():
    def __init__(self, verbose=True, objfunc=None, comm=None):
        self.niterations = 0
        self.residual = 0.
        self.verbose = verbose
        self.objfunc = objfunc
        self.comm = comm or MPI.COMM_WORLD
    def __call__(self, x):
        import inspect
        plocals = inspect.stack()[1][0].f_locals
        self.niterations = plocals['iter_'] - 1
        self.residual = plocals['resid']
        if self.niterations == 0 and self.comm.Get_rank() == 0:
            docriterion = self.objfunc is not None or 'info' in plocals and \
                isinstance(plocals['info'], dict) and 'terms' in plocals['info']
            if MPI.COMM_WORLD.Get_rank() == 0 and self.verbose:
                print('Iteration\tResiduals' + ('\tCriterion' if docriterion \
                      else ''))
        if self.verbose: 
            if 'info' in plocals and isinstance(plocals['info'], dict):
                info = plocals['info']
                Js = info['terms'] if 'terms' in info else None
                niterls = info['niterations'] if 'niterations' in info else 0
            else:
                Js = None
                niterls = 0

            if Js is None and self.objfunc is not None:
                Js = self.objfunc(x)

            if Js is None:
                Jstr = ''
            else:
                Jstr = repr(sum(Js))
                if len(Js) > 1:
                    Jstr += ' ' + str(Js)

            if self.comm.Get_rank() == 0:
                if niterls < 1:
                    print("%s\t%s\t%s" % (str(self.niterations).rjust(9),
                                          str(self.residual), Jstr))
                else:
                    print("%s (%s)\t%s\t%s" % (str(self.niterations).rjust(9),
                         str(niterls), str(self.residual), Jstr))

