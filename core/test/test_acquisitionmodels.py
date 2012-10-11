import numpy as np

from numpy.testing import assert_equal, assert_almost_equal
from pyoperators import AdditionOperator, CompositionOperator, BlockDiagonalOperator, IdentityOperator, MaskOperator
from pyoperators.utils.testing import assert_is_instance
from tamasis.acquisitionmodels import CompressionAverageOperator, ConvolutionTruncatedExponentialOperator, DdTddOperator, DiscreteDifferenceOperator, DownSamplingOperator, FftHalfComplexOperator, InvNttUncorrelatedOperator,  InvNttUncorrelatedPythonOperator, PackOperator, PadOperator, ShiftOperator, UnpackOperator
from tamasis.utils import all_eq
    
def test_compression_average1():
    data = np.array([1., 2., 2., 3.])
    compression = CompressionAverageOperator(2)
    compressed = compression(data)
    assert_equal(compressed, [1.5, 2.5])
    assert_equal(compression.T(compressed), [0.75, 0.75, 1.25, 1.25])

def test_compression_average2():
    partition = (10,5)
    tod = np.empty((2,15), float)
    tod[0,:] = [1,1,1,1,1,3,3,3,3,3,4,4,4,4,4]
    tod[1,:] = [1,2,1.,0.5,0.5,5,0,0,0,0,1,2,1,2,1.5]
    compression = CompressionAverageOperator(5, partitionin=partition)
    tod2 = compression(tod)
    assert tod2.shape == (2,3)
    assert_equal(tod2, [[1.,3.,4.],[1.,1.,1.5]])

    tod3 = compression.T(tod2)
    assert tod3.shape == (2,15)
    assert_almost_equal(tod3[0,:], (0.2,0.2,0.2,0.2,0.2,0.6,0.6,0.6,0.6,0.6,0.8,0.8,0.8,0.8,0.8))
    
    tod = np.array([1,2,2,3,3,3,4,4,4,4])
    compression = CompressionAverageOperator([1,2,3,4], partitionin=[1,2,3,4])
    tod2 = compression(tod)
    assert_almost_equal(tod2, [1,2,3,4])
    tod3 = compression.T(tod2)
    assert_almost_equal(tod3, 10*[1])

def test_compression_average3():
    a = CompressionAverageOperator(3)
    assert_almost_equal(a.todense(9).T, a.T.todense(3))

def test_downsampling1():
    partition = (1,2,3,4)
    tod = np.array([1,2,1,3,1,1,4,1,1,1])
    compression=DownSamplingOperator([1,2,3,4], partitionin=partition)
    tod2 = compression(tod)
    assert_equal(tod2, [1,2,3,4])
    tod3 = compression.T(tod2)
    assert_equal(tod3, [1,2,0,3,0,0,4,0,0,0])

def test_downsampling2():
    a = CompressionAverageOperator(3)
    assert_almost_equal(a.todense(9).T, a.T.todense(3))

def test_padding1():
    padding = PadOperator(left=1,right=20)
    a = np.arange(10*15).reshape((10,15))
    b = padding(a)
    assert b.shape == (10,36)
    assert_equal(b[:,0:1], 0)
    assert_equal(b[:,1:16], a)
    assert_equal(b[:,16:], 0)
    shapein = (10,15)
    shapeout = (10,15+1+20)
    assert_equal(padding.T.todense(shapeout), padding.todense(shapein).T)

def test_padding2():
    padding = PadOperator(left=1,right=(4,20), partitionin=(12,3))
    a = np.arange(10*15).reshape((10,15))
    b = padding(a)
    assert b.shape == (10,41)
    assert_equal(b[:,0:1], 0)
    assert_equal(b[:,1:13], a[:,0:12])
    assert_equal(b[:,13:17], 0)
    assert_equal(b[:,17:18], 0)
    assert_equal(b[:,18:21], a[:,12:])
    assert_equal(b[:,21:], 0)
    shapein = (10,15)
    shapeout = (10,(12+1+4)+(3+1+20))
    assert_equal(padding.T.todense(shapeout), padding.todense(shapein).T)

def test_convolution_truncated_exponential():
    r = ConvolutionTruncatedExponentialOperator(1., shapein=(1,10))
    a = np.ones((1,10))
    b = r(a)
    assert np.allclose(a, b)

    a[0,1:]=0
    b = r(a)
    assert_almost_equal(b[0,:], [np.exp(-t/1.) for t in range(0,10)])
    assert_almost_equal(r.T.todense(), r.todense().T)

def test_ffthalfcomplex1():
    n = 100
    fft = FftHalfComplexOperator(n)
    a = np.random.random(n)+1
    b = fft.T(fft(a))
    assert np.allclose(a, b)
    a = np.random.random((3,n))+1
    b = fft.T(fft(a))
    assert np.allclose(a, b)

def test_ffthalfcomplex2():
    nsamples = 1000
    fft = FftHalfComplexOperator(nsamples)
    a = np.random.random((10,nsamples))+1
    b = fft.T(fft(a))
    assert np.allclose(a, b)

def test_ffthalfcomplex3():
    partition = (100,300,5,1000-100-300-5)
    ffts = [FftHalfComplexOperator(p) for p in partition]
    fft = BlockDiagonalOperator(ffts, partitionin=partition, axisin=-1)
    a = np.random.random((10,np.sum(partition)))+1
    b = fft(a)
    b_ = np.hstack([ffts[0](a[:,:100]), ffts[1](a[:,100:400]), ffts[2](a[:,400:405]), ffts[3](a[:,405:])])
    assert np.allclose(b, b_)
    b = fft.T(fft(a))
    assert np.allclose(a, b)

def test_invntt_uncorrelated():
    filter = np.array([0., 0.2, 0.3, 0.4, 0.4, 0.3, 0.2, 0]).reshape((1,-1))
    ncorrelations = 2
    invntt = InvNttUncorrelatedOperator(filter, ncorrelations, 3)
    invntt_todense = invntt.todense()
    assert all_eq(invntt_todense, invntt.todense(inplace=True))
    invntt2 = InvNttUncorrelatedPythonOperator(filter, ncorrelations, 3)
    invntt2_todense = invntt2.todense()
    assert all_eq(invntt2_todense, invntt2.todense(inplace=True))
    assert all_eq(invntt_todense, invntt2_todense)

def test_addition():
    def func(nops):
        ops = [ DiscreteDifferenceOperator(axis=axis,shapein=(2,3,4,5)) \
                for axis in range(nops) ]
        model = AdditionOperator(ops)
        v = np.arange(2*3*4*5.).reshape(2,3,4,5)
        a = model(v)
        b = sum([op(v) for op in ops])
        c = model(v,v)
        assert np.all(a == b)
        assert np.all(b == c)
    for nops in range(1, 5):
        yield func, nops

def test_additionT():
    def func(nops):
        ops = [ DiscreteDifferenceOperator(axis=axis,shapein=(2,3,4,5)) \
                for axis in range(nops) ]
        model = AdditionOperator(ops)
        assert model.T.T is model
        v = np.arange(2*3*4*5.).reshape(2,3,4,5)
        a = model.T(v)
        b = sum([op.T(v) for op in ops])
        c = model.T(v,v)
        assert np.all(a == b)
        assert np.all(b == c)
    for nops in range(1, 5):
        yield func, nops

def test_composition():
    def func(nops):
        ops = [ DiscreteDifferenceOperator(axis=axis,shapein=(2,3,4,5)) \
                for axis in range(nops) ]
        model = CompositionOperator(ops)
        v = np.arange(2*3*4*5.).reshape(2,3,4,5)
        a = model(v)
        b = v.copy()
        for m in reversed(ops):
            b = m(b)
        c = model(v,v)
        assert np.all(a == b) and np.all(b == c)
    for nops in range(1, 5):
        yield func, nops

def test_compositionT():
    def func(nops):
        ops = [ DiscreteDifferenceOperator(axis=axis,shapein=(2,3,4,5)) \
                for axis in range(nops) ]
        model = CompositionOperator(ops)
        assert model.T.T is model
        v = np.arange(2*3*4*5.).reshape(2,3,4,5)
        a = model.T(v)
        b = v.copy()
        for m in ops:
            b = m.T(b)
        c = model.T(v,v)
        assert np.all(a == b) and np.all(b == c)
    for nops in range(1, 5):
        yield func, nops


def test_packing():

    p = PackOperator([False, True, True, False])
    assert all_eq(p([1,2,3,4]), [1,4])
    assert all_eq(p.T([1,4]), [1,0,0,4])

    u = UnpackOperator([False, True, True, False])
    assert all_eq(u([1,4]), [1,0,0,4])
    assert all_eq(u.T([1,2,3,4]), [1,4])

    pdense = p.todense()
    udense = u.todense()
    assert all_eq(pdense, p.todense(inplace=True))
    assert all_eq(udense, u.todense(inplace=True))
    assert all_eq(pdense, udense.T)

    assert_is_instance(p*u, IdentityOperator)
    assert_is_instance(u*p, MaskOperator)
    m = u * p
    assert all_eq(np.dot(udense, pdense), m.todense())

def test_diff():
    def func(shape, axis):
        dX = DiscreteDifferenceOperator(axis=axis, shapein=shape)
        dTd = DdTddOperator(axis=axis, shapein=shape)
        dX_dense = dX.todense()
        assert_equal(dX_dense, dX.todense(inplace=True))

        dXT_dense = dX.T.todense()
        assert_equal(dXT_dense, dX.T.todense(inplace=True))
        assert_equal(dX_dense.T, dXT_dense)

        dtd_dense = dTd.todense()
        assert_equal(dtd_dense, dTd.todense(inplace=True))
        assert_equal(np.matrix(dXT_dense) * np.matrix(dX_dense), dtd_dense)

    for shape in ((3,), (3,4), (3,4,5), (3,4,5,6)):
        for axis in range(len(shape)):
            yield func, shape, axis

def test_shift1():
    for axis in range(4):
        shift = ShiftOperator(1, axis=axis, shapein=(3,4,5,6))
        yield assert_equal, shift.todense().T, shift.T.todense()

def test_shift2():
    for axis in range(1,4):
        shift = ShiftOperator(((1,2,3),), axis=axis, shapein=(3,4,5,6))
        yield assert_equal, shift.todense().T, shift.T.todense()

def test_shift3():
    for offset in ( (3,), (3,4), (3,4,5) ):
        for axis in range(len(offset),4):
            s = np.random.random_integers(-2,2,offset)
            shift = ShiftOperator([s], axis=axis, shapein=(3,4,5,6))
            yield assert_equal, shift.todense().T, shift.T.todense()
