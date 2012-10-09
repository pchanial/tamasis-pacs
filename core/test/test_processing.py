import numpy as np

from pysimulators.datatypes import Tod
from tamasis.processing import filter_polynomial, interpolate_linear, filter_nonfinite
from tamasis.utils import assert_all_eq


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
        
    for y_, m_ in zip((y, [y,y], [[y],[y]]), (m, [m,m], [[m],[m]])):
        y_ = np.asarray(y_)
        filter_polynomial(y_, 1, m_, out=y_)
        yield func, y_
        
    y2 = np.concatenate([y,y])
    m2 = np.concatenate([m,m])
    o2 = filter_polynomial(y2, 1, m2, partition=(4,4))
    assert o2.shape == y2.shape
    yield func, o2

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
    
def test_filter_nonfinite():

    class ndarraywrap(np.ndarray):
        pass

    def func(x):
        y = filter_nonfinite(x)
        y_expected = np.asanyarray(x).copy()
        y_expected[~np.isfinite(x)] = 0
        if hasattr(x, 'mask') and x.mask is not None:
            y_expected.mask = x.mask.copy()
            y_expected.mask[~np.isfinite(x)] = True
        assert_all_eq(y, y_expected)
            
        x = np.asanyarray(x)
        filter_nonfinite(x, out=x)
        assert_all_eq(x, y_expected)

    mask = [False,False,False,True,False,True]
    data1 = [1,np.nan,5,np.inf,-np.inf,8]
    data2 = np.array(data1).view(ndarraywrap); data2.mask = np.array(mask)
    data3 = Tod(data1)
    data4 = Tod(data1, mask=mask)

    for x in (data1, data2, data3, data4):
        yield func, x

    mask = [[False,False,False],[True,False,True]]
    data1 = [[1,np.nan,5], [np.inf,-np.inf,8]]
    data2 = np.array(data1).view(ndarraywrap); data2.mask = np.array(mask)
    data3 = Tod(data1)
    data4 = Tod(data1, mask=mask)

    for x in (data1, data2, data3, data4):
        yield func, x

