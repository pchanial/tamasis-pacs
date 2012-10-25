import os
import tamasis
import pyoperators

from pyoperators import DiagonalOperator, MaskOperator
from tamasis import PacsObservation, mapper_nl
from tamasis.linalg import norm2, norm2_ellipsoid

class TestFailure(Exception): pass

pyoperators.memory.verbose = False
tamasis.var.verbose = False
profile = None#'test_rls.png'
data_dir = os.path.dirname(__file__) + '/data/'
obs = PacsObservation(filename=data_dir+'frames_blue.fits',
                      fine_sampling_factor=1)
obs.pointing.chop[:] = 0
tod = obs.get_tod()

projection = obs.get_projection_operator(resolution=3.2, downsampling=True,
                                         npixels_per_sample=6)
masking_tod = MaskOperator(tod.mask)

model = masking_tod * projection

# iterative map, taking all map pixels
class Callback():
    def __init__(self):
        self.niterations = 0
    def __call__(self, x):
        self.niterations += 1

invntt = DiagonalOperator(1/obs.get_detector_stddev(100)**2,
                          broadcast='rightward')

def test():
    map_nl = mapper_nl(tod, model, hypers=2*[1.],
                       norms=[norm2_ellipsoid(invntt)] + 2*[norm2],
                       tol=1.e-4, maxiter=1000,
                       callback=None if tamasis.var.verbose else Callback(),
                       )

    print 'Elapsed time: ' + str(map_nl.header['TIME']) + ' after ' + \
        str(map_nl.header['NITER']) + ' iterations.'
    if map_nl.header['NITER'] > 150:
        raise TestFailure()

if __name__ == '__main__':
    test()
