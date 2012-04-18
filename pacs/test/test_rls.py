import os
import tamasis
import pyoperators

from scipy.sparse.linalg import cgs
from tamasis import *
from tamasis.numpyutils import all_eq

pyoperators.memory.verbose = False
tamasis.var.verbose = False
profile = None#'test_rls.png'
data_dir = os.path.dirname(__file__) + '/data/'

PacsInstrument.info.CALFILE_BADP = tamasis.var.path + '/pacs/PCalPhotometer_Ba'\
                                   'dPixelMask_FM_v5.fits'
PacsInstrument.info.CALFILE_RESP = tamasis.var.path + '/pacs/PCalPhotometer_Re'\
                                   'sponsivity_FM_v5.fits'

obs = PacsObservation(filename=data_dir+'frames_blue.fits',
                      fine_sampling_factor=1, reject_bad_line=False)
obs.pointing.chop[:] = 0
tod = obs.get_tod(subtraction_mean=True)

projection  = ProjectionOperator(obs, resolution=3.2, downsampling=True,
                                 npixels_per_sample=6)
masking_tod = MaskOperator(tod.mask)
model = masking_tod * projection

# iterative map, taking all map pixels
class Callback():
    def __init__(self):
        self.niterations = 0
    def __call__(self, x):
        self.niterations += 1

map_rls = mapper_rls(tod, model, hyper=1., tol=1.e-4, profile=profile,
                     callback=None if tamasis.var.verbose else Callback(),
                     solver=cgs)

if profile is None:
    print 'Elapsed time: ' + str(map_rls.header['TIME'])

def test():
    assert map_rls.header['NITER'] < 49
    ref = Map(data_dir + 'frames_blue_map_rls_cgs_tol1e-6.fits')
    ref.derived_units = map_rls.derived_units
    cov = ref.coverage > 80
    assert all_eq(ref[cov], map_rls[cov], 1.e-1)
    cov = ref.coverage > 125
    assert all_eq(ref[cov], map_rls[cov], 1.e-2)
