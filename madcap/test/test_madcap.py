import numpy as np
import pyfits
import os
import tamasis
from tamasis import (MadMap1Observation, DiagonalOperator, InvNtt, Packing,
                     Projection, mapper_naive, mapper_ls, Tod)
from tamasis.numpyutils import all_eq

class TestFailure(Exception): pass

tamasis.var.verbose = False
profile=None #'test_madcap.png'
path = os.path.abspath(os.path.dirname(__file__)) + '/data/madmap1/'
obs = MadMap1Observation(path+'todSpirePsw_be', path+'invnttSpirePsw_be', 
                         path+'madmapSpirePsw.fits[coverage]', 'big_endian',
                         135, missing_value=np.nan)
obs.instrument.name = 'SPIRE/PSW'

tod = obs.get_tod(unit='Jy/beam')
projection = Projection(obs)
packing = Packing(obs.info.mapmask)

model = projection*packing
map_naive = mapper_naive(tod, projection*packing)
map_ref = pyfits.open(path+'naivemapSpirePsw.fits')['image'].data

def test_madcap1():
    assert all_eq(map_naive, map_ref)

map_naive_1d = mapper_naive(tod, projection)
map_naive_2d = packing.T(map_naive_1d)
map_naive_2d = packing.T(projection.T(tod)/projection.T(Tod.ones(tod.shape)))
map_naive_2d[obs.info.mapmask] = np.nan

def test_madcap2():
    assert all_eq(map_naive_2d, map_ref)

M = DiagonalOperator(packing(1/map_naive.coverage))
if np.any(~np.isfinite(M.data)):
    raise TestFailure()

class Callback:
    def __init__(self):
        self.niterations = 0
    def __call__(self, x):
        self.niterations += 1

invntt = InvNtt(obs)

def test_madcap3():
    map_lsw2_packed = mapper_ls(tod, projection, invntt=invntt, tol=1.e-7, M=M, callback=Callback(), criterion=False, profile=profile)
    if profile:
        return
    print 'Elapsed time:', map_lsw2_packed.header['TIME']
    assert map_lsw2_packed.header['NITER'] < 50

if __name__ == '__main__':
    test_madcap3()
