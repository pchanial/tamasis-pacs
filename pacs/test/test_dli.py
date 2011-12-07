import os
import tamasis
import pyoperators

from tamasis import *


class TestFailure(Exception): pass

pyoperators.memory.verbose=False
tamasis.var.verbose = True
data_dir = os.path.dirname(__file__) + '/data/'
obs = PacsObservation(filename=data_dir+'frames_blue.fits')
obs.pointing.chop[:] = 0
tod = obs.get_tod(flatfielding=False)

projection  = Projection(obs, resolution=3.2, downsampling=True,
                         npixels_per_sample=6)
masking_tod = Masking(tod.mask)
model = masking_tod * projection

naive = mapper_naive(tod, model)
naive[np.isnan(naive)] = 0

prior = BlockColumnOperator([DiscreteDifference(axis, shapein=(103,97)) for axis in (0,1)])

dli = DoubleLoopAlgorithm(model, tod, prior)
map_dli = dli()
