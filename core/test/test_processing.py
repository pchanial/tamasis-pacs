import numpy as np
from numpy.testing import assert_array_equal
from tamasis.processing import remove_nonfinite
from tamasis.datatypes import Tod

def test_remove_non_finite1():
    t = Tod([1,np.nan,5,np.inf,-np.inf,8])
    remove_nonfinite(t)
    assert_array_equal(t, [1,0,5,0,0,8])
    assert_array_equal(t.mask, [False,True,False,True,True,False])

def test_remove_non_finite2():
    t = Tod([1,np.nan,5,np.inf,-np.inf,8],
            mask=[False,False,False,True,False,True])
    remove_nonfinite(t)
    assert_array_equal(t, [1,0,5,0,0,8])
    assert_array_equal(t.mask, [False,True,False,True,True,True])
