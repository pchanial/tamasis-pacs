import numpy as np

from numpy.testing import assert_array_equal, assert_almost_equal
from tamasis.acquisitionmodels import CompressionAverage, DownSampling, Padding, ResponseTruncatedExponential

def test_compression_average1():
    data = np.array([1., 2., 2., 3.])
    compression = CompressionAverage(2)
    compressed = compression(data)
    assert_array_equal(compressed, [1.5, 2.5])
    assert_array_equal(compression.T(compressed), [0.75, 0.75, 1.25, 1.25])

def test_compression_average2():
    nsamples = (10,5)
    tod = np.empty((2,15), float)
    tod[0,:] = [1,1,1,1,1,3,3,3,3,3,4,4,4,4,4]
    tod[1,:] = [1,2,1.,0.5,0.5,5,0,0,0,0,1,2,1,2,1.5]
    compression = CompressionAverage(5, nsamples=nsamples)
    tod2 = compression(tod)
    assert tod2.shape == (2,3)
    assert_array_equal(tod2, [[1.,3.,4.],[1.,1.,1.5]])

    tod3 = compression.T(tod2)
    assert tod3.shape == (2,15)
    assert_almost_equal(tod3[0,:], (0.2,0.2,0.2,0.2,0.2,0.6,0.6,0.6,0.6,0.6,0.8,0.8,0.8,0.8,0.8))
    
    tod = np.array([1,2,2,3,3,3,4,4,4,4])
    compression = CompressionAverage([1,2,3,4], nsamples=[1,2,3,4])
    tod2 = compression(tod)
    assert_almost_equal(tod2, [1,2,3,4])
    tod3 = compression.T(tod2)
    assert_almost_equal(tod3, 10*[1])

    a = CompressionAverage(3)
    assert_almost_equal(a.todense().T, a.T.todense())


def test_downsampling():
    nsamples = (1,2,3,4)
    tod = np.array([1,2,1,3,1,1,4,1,1,1])
    compression=DownSampling([1,2,3,4], nsamples=nsamples)
    tod2 = compression(tod)
    assert_array_equal(tod2, [1,2,3,4])
    tod3 = compression.T(tod2)
    assert_array_equal(tod3, [1,2,0,3,0,0,4,0,0,0])


def test_padding():
    padding = Padding(left=1,right=20)
    a = np.arange(10*15).reshape((10,15))
    b = padding(a)
    assert b.shape == (10,36)
    assert_array_equal(b[:,0:1], 0)
    assert_array_equal(b[:,1:16], a)
    assert_array_equal(b[:,16:], 0)

    c = padding.T(b)
    assert c.shape == a.shape
    assert_array_equal(c, a)

    padding = Padding(left=1,right=(4,20), nsamples=(12,3))
    b = padding(a)
    assert b.shape == (10,41)
    assert_array_equal(b[:,0:1], 0)
    assert_array_equal(b[:,1:13], a[:,0:12])
    assert_array_equal(b[:,13:17], 0)
    assert_array_equal(b[:,17:18], 0)
    assert_array_equal(b[:,18:21], a[:,12:])
    assert_array_equal(b[:,21:], 0)

    c = padding.T(b)
    assert c.shape == a.shape
    assert_array_equal(c, a)


def test_response_truncated_exponential():
    r = ResponseTruncatedExponential(1., shapein=(1,10))
    a = np.ones((1,10))
    b = r(a)
    assert np.allclose(a, b)

    a[0,1:]=0
    b = r(a)
    assert_almost_equal(b[0,:], [np.exp(-t/1.) for t in range(0,10)])
    assert_almost_equal(r.T.todense(), r.todense().T)

test_compression_average1()
test_compression_average2()
test_downsampling()
test_padding()
test_response_truncated_exponential()
