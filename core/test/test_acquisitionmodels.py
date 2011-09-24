import numpy as np

from operators import Operator, PartitionOperator
from operators.utils import isscalar
from numpy.testing import assert_array_equal, assert_almost_equal, assert_raises
from tamasis.acquisitionmodels import CompressionAverage, DownSampling, FftHalfComplex, Padding, ResponseTruncatedExponential, partitioned


def test_partitioning():

    @partitioned('value', 'mykey')
    class MyOp(Operator):
        """ Bla1. """
        def __init__(self, arg1, value, arg3, mykey=None, **keywords):
            """ Bla2. """
            Operator.__init__(self, **keywords)
            self.arg1 = arg1
            self.value = value
            self.arg3 = arg3
            self.mykey = mykey
        def direct(self, input, output):
            output[:] = self.value * input
        __str__ = Operator.__repr__

    arg1 = [1,2,3,4,5]
    arg3 = ['a', 'b', 'c', 'd']

    def func(op, n, v, k):
        input = np.ones(np.sum(n))
        if not isinstance(n, tuple) or len(n) == 1:
            if n is not None:
                if not isscalar(v): v = v[0]
                if not isscalar(k): k = k[0]
            assert op.__class__ is MyOp
            assert op.arg1 is arg1
            assert op.value is v
            assert op.arg3 is arg3
            assert op.mykey is k
            if not isinstance(v, tuple):
                output = op(input)
                assert_array_equal(output, v)
        else:
            assert op.__class__ is PartitionOperator
            v = len(n) * [v] if isscalar(v) else v[0:len(n)]
            k = len(n) * [k] if isscalar(k) else k[0:len(n)]
            for i, p in enumerate(n):
                func(op.operands[i], None, v[i], k[i])
            expected = np.hstack(n_*[v_] for n_, v_ in zip(n,v))
            output = op(input)
            assert_array_equal(output, expected)

    for n in (None, 2, (2,), (4,2), (5,4,2)):
        for v in (2., (2.,), (2., 3)):
            for k in (0., (0.,), (0., 1.)):
                if isinstance(n, tuple) and len(n) > 1 and \
                  (isinstance(v, tuple) and len(v) < len(n) or \
                   isinstance(k, tuple) and len(k) < len(n)):
                    yield assert_raises, IndexError, \
                          lambda: MyOp(arg1, v, arg3, mykey=k, partition=n)
                else:
                    myop = MyOp(arg1, v, arg3, mykey=k, partition=n)
                    yield func, myop, n, v, k


def test_compression_average1():
    data = np.array([1., 2., 2., 3.])
    compression = CompressionAverage(2)
    compressed = compression(data)
    assert_array_equal(compressed, [1.5, 2.5])
    assert_array_equal(compression.T(compressed), [0.75, 0.75, 1.25, 1.25])

def test_compression_average2():
    partition = (10,5)
    tod = np.empty((2,15), float)
    tod[0,:] = [1,1,1,1,1,3,3,3,3,3,4,4,4,4,4]
    tod[1,:] = [1,2,1.,0.5,0.5,5,0,0,0,0,1,2,1,2,1.5]
    compression = CompressionAverage(5, partition=partition)
    tod2 = compression(tod)
    assert tod2.shape == (2,3)
    assert_array_equal(tod2, [[1.,3.,4.],[1.,1.,1.5]])

    tod3 = compression.T(tod2)
    assert tod3.shape == (2,15)
    assert_almost_equal(tod3[0,:], (0.2,0.2,0.2,0.2,0.2,0.6,0.6,0.6,0.6,0.6,0.8,0.8,0.8,0.8,0.8))
    
    tod = np.array([1,2,2,3,3,3,4,4,4,4])
    compression = CompressionAverage([1,2,3,4], partition=[1,2,3,4])
    tod2 = compression(tod)
    assert_almost_equal(tod2, [1,2,3,4])
    tod3 = compression.T(tod2)
    assert_almost_equal(tod3, 10*[1])

def test_compression_average3():
    a = CompressionAverage(3)
    assert_almost_equal(a.todense(9).T, a.T.todense(3))

def test_downsampling1():
    partition = (1,2,3,4)
    tod = np.array([1,2,1,3,1,1,4,1,1,1])
    compression=DownSampling([1,2,3,4], partition=partition)
    tod2 = compression(tod)
    assert_array_equal(tod2, [1,2,3,4])
    tod3 = compression.T(tod2)
    assert_array_equal(tod3, [1,2,0,3,0,0,4,0,0,0])

def test_downsampling2():
    a = CompressionAverage(3)
    assert_almost_equal(a.todense(9).T, a.T.todense(3))

def test_padding1():
    padding = Padding(left=1,right=20)
    a = np.arange(10*15).reshape((10,15))
    b = padding(a)
    assert b.shape == (10,36)
    assert_array_equal(b[:,0:1], 0)
    assert_array_equal(b[:,1:16], a)
    assert_array_equal(b[:,16:], 0)
    shapein = (10,15)
    shapeout = (10,15+1+20)
    assert_array_equal(padding.T.todense(shapeout), padding.todense(shapein).T)

def test_padding2():
    padding = Padding(left=1,right=(4,20), partition=(12,3))
    a = np.arange(10*15).reshape((10,15))
    b = padding(a)
    assert b.shape == (10,41)
    assert_array_equal(b[:,0:1], 0)
    assert_array_equal(b[:,1:13], a[:,0:12])
    assert_array_equal(b[:,13:17], 0)
    assert_array_equal(b[:,17:18], 0)
    assert_array_equal(b[:,18:21], a[:,12:])
    assert_array_equal(b[:,21:], 0)
    shapein = (10,15)
    shapeout = (10,(12+1+4)+(3+1+20))
    assert_array_equal(padding.T.todense(shapeout), padding.todense(shapein).T)

def test_response_truncated_exponential():
    r = ResponseTruncatedExponential(1., shapein=(1,10))
    a = np.ones((1,10))
    b = r(a)
    assert np.allclose(a, b)

    a[0,1:]=0
    b = r(a)
    assert_almost_equal(b[0,:], [np.exp(-t/1.) for t in range(0,10)])
    assert_almost_equal(r.T.todense(), r.todense().T)

def test_ffthalfcomplex1():
    n = 100
    fft = FftHalfComplex(n)
    a = np.random.random(n)+1
    b = fft.T(fft(a))
    assert np.allclose(a, b)
    a = np.random.random((3,n))+1
    b = fft.T(fft(a))
    assert np.allclose(a, b)

def test_ffthalfcomplex2():
    nsamples = 1000
    fft = FftHalfComplex(nsamples)
    a = np.random.random((10,nsamples))+1
    b = fft.T(fft(a))
    assert np.allclose(a, b)

def test_ffthalfcomplex3():
    partition = (100,300,5,1000-100-300-5)
    ffts = [FftHalfComplex(p) for p in partition]
    fft = PartitionOperator(ffts, partitionin=partition, axisin=-1)
    a = np.random.random((10,np.sum(partition)))+1
    b = fft(a)
    b_ = np.hstack([ffts[0](a[:,:100]), ffts[1](a[:,100:400]), ffts[2](a[:,400:405]), ffts[3](a[:,405:])])
    assert np.allclose(b, b_)
    b = fft.T(fft(a))
    assert np.allclose(a, b)
