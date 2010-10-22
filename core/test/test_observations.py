import numpy
from tamasis.observations import *

flags = ['bad', 'u1', 'u2']
good_policy = ['kEep', 'removE', 'MASK']
bad_policy = ['remove', 'KKeep']

mask_policy = MaskPolicy(flags, good_policy)
if numpy.any(numpy.array(mask_policy) != (0,2,1)): raise TestFailure()
if mask_policy.bad != 'keep' or mask_policy.u1 != 'remove' or mask_policy.u2 != 'mask': raise TestFailure()
try:
    junk = MaskPolicy(flags, bad_policy)
except ValueError:
    pass
try:
    junk = MaskPolicy(flags[0:2], bad_policy)
except KeyError:
    pass

print 'OK.'
