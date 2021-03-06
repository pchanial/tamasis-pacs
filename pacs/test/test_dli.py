import numpy as np
import os
import pyoperators
import tamasis

from pyoperators import BlockColumnOperator, MaskOperator
from pyoperators.iterative.algorithms import StopCondition
from pyoperators.iterative.dli import DoubleLoopAlgorithm
from tamasis import PacsObservation, DiscreteDifferenceOperator, mapper_naive

pyoperators.memory.verbose=False
tamasis.var.verbose = True
data_dir = os.path.dirname(__file__) + '/data/'
obs = PacsObservation(filename=data_dir+'frames_blue.fits')
obs.pointing.chop[:] = 0
tod = obs.get_tod()

projection = obs.get_projection_operator(resolution=3.2,
                                         downsampling=True,
                                         npixels_per_sample=6)
masking_tod = MaskOperator(tod.mask)
model = masking_tod * projection

naive = mapper_naive(tod, model)
naive[np.isnan(naive)] = 0

prior = BlockColumnOperator([DiscreteDifferenceOperator(axis, shapein=(103,97)) for axis in (0,1)], new_axisout=0)

stop_condition = StopCondition(maxiter=2)
dli = DoubleLoopAlgorithm(model, tod, prior, stop_condition=stop_condition, lanczos={'maxiter':5}, fmin_args={'maxiter':2})

map_dli = dli.run()

def test():
    pass
