import numpy as np
from numpy.testing import assert_array_equal, assert_almost_equal
from tamasis.utils import distance, gaussian, airy_disk, profile, _distance_slow, diff, shift

def test_distance1():
    origin = (1.,)
    d0 = np.arange(5.)*0.5
    d1 = distance(5, origin=origin, resolution=0.5)
    d2 = _distance_slow((5,), origin, [0.5], None)
    assert_array_equal(d0, d1)
    assert_array_equal(d0, d2)

def test_distance2():
    origin = (1.,2.)
    d0 = np.array([[2., np.sqrt(5), np.sqrt(8)], [0,1,2]])
    d1 = distance((2,3), origin=origin, resolution=(1.,2.))
    d2 = _distance_slow((2,3), origin, [1.,2.], None)
    assert_array_equal(d0, d1)
    assert_array_equal(d0, d2)

def test_fwhm():
    def func(f, fwhm, resolution, n):
        m = f((1000,1000), fwhm=fwhm, resolution=resolution)
        assert np.sum(m[500,:] > np.max(m)/2) == n
    for f in [gaussian, airy_disk]:
        for fwhm, resolution, n in zip([10,10,100], [0.1,1,10], [100,10,10]):
            yield func, f, fwhm, resolution, n

def test_diff():
    def func(naxis, axis):
        shape = [i+8 for i in range(naxis)]
        a=np.random.random_integers(1, 10, size=shape).astype(float)
        ref = -np.diff(a, axis=axis)
        diff(a, a, axis=axis)
        s = tuple([slice(0,s) for s in ref.shape])
        print s, a.shape, ref.shape
        a[s] -= ref
        assert_array_equal(a, 0)
    for naxis in range(5):
        for axis in range(naxis):
            yield func, naxis, axis

def test_shift1():
    def func(a, s):
        b = np.empty_like(a)
        shift(b, b, s, axis=-1)
        assert_array_equal(b, 0)
    for a in (np.ones(10),np.ones((12,10))):
        for s in (10,11,100,-10,-11,-100):
            yield func, a, s

def test_shif2():
    a = np.array([[1.,1.,1.,1.],[2.,2.,2.,2.]])
    shift(a, a, [1,-1], axis=1)
    assert_array_equal(a, [[0,1,1,1],[2,2,2,0]])

def test_shif3():
    a = np.array([[0.,0,0],[0,1,0],[0,0,0]])
    b = np.empty_like(a)
    shift(a, b, 1, axis=0)
    assert_array_equal(b, np.roll(a,1,axis=0))
    shift(b, b, -2, axis=0)
    assert_array_equal(b, np.roll(a,-1,axis=0))

def test_shif4():
    a = np.array([[0.,0,0],[0,1,0],[0,0,0]])
    b = np.empty_like(a)
    shift(a, b, 1, axis=1)
    assert_array_equal(b, np.roll(a,1,axis=1))
    shift(b, b, -2, axis=1)
    assert_array_equal(b, np.roll(a,-1,axis=1))

def test_profile():
    def profile_slow(input, origin=None, bin=1.):
        d = distance(input.shape, origin=origin)
        d /= bin
        d = d.astype(int)
        m = np.max(d)
        p = np.ndarray(int(m+1))
        n = np.zeros(m+1, int)
        for i in range(m+1):
            p[i] = np.mean(input[d == i])
            n[i] = np.sum(d==i)
        return p, n
    d = distance((10,20))
    x, y, n = profile(d, origin=(4,5), bin=2., histogram=True)
    y2, n2 = profile_slow(d, origin=(4,5), bin=2.)
    assert_almost_equal(y, y2[0:y.size])
    assert_almost_equal(n, n2[0:n.size])
