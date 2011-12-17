import numpy as np
import pyfits
import os
import tamasis
from tamasis import *

class TestFailure(Exception): pass

tamasis.var.verbose = True
path = os.path.abspath(os.path.dirname(__file__)) + '/data/madmap1/'
obs = MadMap1Observation(path+'todSpirePsw_be', path+'invnttSpirePsw_be', 
                         path+'madmapSpirePsw.fits[coverage]', 'big_endian',
                         135, missing_value=np.nan)
obs.instrument.name = 'SPIRE/PSW'

tod = obs.get_tod(unit='Jy/beam')
projection = Projection(obs)
packing = Unpacking(obs.info.mapmask, field=np.nan).T

map_naive = mapper_naive(tod, projection)
map_naive = mapper_naive(tod, projection*packing)
map_ref = pyfits.open(path+'naivemapSpirePsw.fits')['image'].data
if any_neq(map_naive,map_ref): raise TestFailure('mapper_naive madcap 1')

map_naive_1d = mapper_naive(tod, projection)
map_naive_2d = packing.T(map_naive_1d)
if any_neq(map_naive_2d, map_ref): raise TestFailure('mapper_naive madcap 2')

packing = Unpacking(obs.info.mapmask).T

#M = 1 / map_naive.coverage
#M[~np.isfinite(M)] = np.max(M[np.isfinite(M)])
#map_rlsw1 = mapper_rls(tod, projection*packing, padding.T * fft.T * invNtt * fft * padding, hyper=0, tol=1.e-5, M=M, solver=cg)
#print('Elapsed time: ' + str(map_rlsw1.header['time']))

M = Diagonal(packing(1/map_naive.coverage))
if np.any(~np.isfinite(M.data)):
    raise TestFailure()

def callback(x):
    pass
#callback=None
map_lsw2_packed = mapper_ls(tod, projection, invntt=InvNtt(obs), tol=1.e-7, M=M, callback=callback, criterion=False)
print 'Elapsed time:', map_lsw2_packed.header['TIME']

map_lsw2 = packing.T(map_lsw2_packed)
