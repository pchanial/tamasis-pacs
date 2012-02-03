import numpy as np
from numpy.testing import assert_array_equal
from tamasis.processing import filter_polynomial, interpolate_linear, remove_nonfinite
from tamasis.datatypes import Tod


def test_filter_polynomial():
    y = [1.,0,3,4]
    m = [False, True, False, True]

    def func(out):
        out_ = out.reshape((-1,4))
        for i in range(out_.shape[0]):
            assert np.allclose(out_[i][~np.asarray(m)], 0, atol=1e-14)

    for y_, m_ in zip((y, [y,y], [[y],[y]]), (m, [m,m], [[m],[m]])):
        o_ = filter_polynomial(y_, 1, m_)
        y_ = np.asarray(y_)
        assert o_.shape == y_.shape
        yield func, o_
        
    y2 = np.concatenate([y,y])
    m2 = np.concatenate([m,m])
    o2 = filter_polynomial(y2, 1, m2, partition=(4,4))
    assert o2.shape == y2.shape
    yield func, o2

    #interpolate_linear(y2, m2, partition=(4,4), out=y2)
    #yield func, y2

def test_interpolate_linear():
    y = [1.,0,3,4]
    m = [False, True, False, True]
    o = interpolate_linear(y, m)
    assert np.all(o == [1,2,3,3])
    def func(out):
        out_ = out.reshape((-1,4))
        for i in range(out_.shape[0]):
            assert np.all(out_[i] == o)
    for y_, m_ in zip((y, [y,y], [[y],[y]]), (m, [m,m], [[m],[m]])):
        o_ = interpolate_linear(y_,m_)
        y_ = np.asarray(y_)
        assert o_.shape == y_.shape
        yield func, o_

        y2 = np.ma.masked_array(y_, mask=m_)
        o_ = interpolate_linear(y2)
        yield func, o_

        interpolate_linear(y_,m_,out=y_)
        yield func, y_

    y2 = np.concatenate([y,y])
    m2 = np.concatenate([m,m])
    o2 = interpolate_linear(y2, m2, partition=(4,4))
    assert o2.shape == y2.shape
    yield func, o2

    interpolate_linear(y2, m2, partition=(4,4), out=y2)
    yield func, y2
    
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
