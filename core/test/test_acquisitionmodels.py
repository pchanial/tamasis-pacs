import numpy
import scipy.sparse.linalg

from tamasis import *

class TestFailure(Exception): pass


#--------------------
# CompressionAverage
#--------------------

tod = Tod(numpy.ndarray((2,15)), nsamples=(10,5))
tod[0,:] = [1,1,1,1,1,3,3,3,3,3,4,4,4,4,4]
tod[1,:] = [1,2,1.,0.5,0.5,5,0,0,0,0,1,2,1,2,1.5]
compression = CompressionAverage(5)
tod2 = compression(tod)
if tod2.shape != (2,3): raise TestFailure('compAvg1')
if tuple(tod2.ravel()) != (1.,3.,4.,1.,1.,1.5): raise TestFailure('compAvg2')
if tod2.nsamples != (2,1): raise TestFailure('compAvg3')
tod3 = compression.T(tod2)
if tod3.shape != (2,15): raise TestFailure('compAvg1t')
if any_neq(tod3[0,:], (0.2,0.2,0.2,0.2,0.2,0.6,0.6,0.6,0.6,0.6,0.8,0.8,0.8,0.8,0.8)): raise TestFailure('compAvg2t')
if tod3.nsamples != (10,5): raise TestFailure('compAvg3t')

tod = Tod(numpy.ndarray((2,15)), nsamples=(11,4))
try:
    tod2 = compression(tod)
except ValidationError:
    pass

#---------
# Masking
#---------

try:
    mask = Masking('lkj')
except TypeError:
    pass
else: 
    raise TestFailure('mask str')
try:
    mask = Masking(3)
except TypeError:
    pass
else: 
    raise TestFailure('mask 3')
try:
    mask = Masking(object)
except TypeError:
    pass
else: 
    raise TestFailure('mask object')
try:
    mask = Masking(numpy.array(3))
except TypeError:
    pass
else: 
    raise TestFailure('mask ndarray rank=0')

mask = Masking(None)
if mask.shapein is not None or mask.shapeout is not None or mask.mask is not None:
    raise TestFailure('Masking(None)')
if mask(3) != 3:
    raise TestFailure('Masking(None)(data)')

try:
    junk = mask('str')
except TypeError:
    pass
else: 
    raise TestFailure('mask None wrong type str')

junk = mask(32)
try:
    junk = mask(object)
except TypeError:
    pass
else: 
    raise TestFailure('mask None wrong type object')

b = numpy.array([3., 4., 1., 0., 3., 2.])
c = numpy.array([3., 4., 0., 0., 3., 0.])

mask = Masking(numpy.array([0, 0., 1., 1., 0., 1], dtype='int8'))
if numpy.any(mask(b) != c): raise TestFailure('mask 1d direct1')
mask = Masking(numpy.array([1, 1., 0., 0., 1., 0]))
if numpy.any(mask(b) != c): raise TestFailure('mask 1d direct2')
mask = Masking(numpy.array([False, False, True, True, False, True]))
if numpy.any(mask(b) != c): raise TestFailure('mask 1d direct3')

b = numpy.array([[3., 4.], [1., 0.], [3., 2.]])
c = numpy.array([[3., 4.], [0., 0.], [3., 0.]])

mask = Masking(numpy.array([[0, 0.], [1., 1.], [0., 1.]], dtype='int8'))
if numpy.any(mask(b) != c): raise TestFailure('mask 2d direct1')
mask = Masking(numpy.array([[1, 1.], [0., 0.], [1., 0.]]))
if numpy.any(mask(b) != c): raise TestFailure('mask 2d direct2')
mask = Masking(numpy.array([[False, False], [True, True], [False, True]]))
if numpy.any(mask(b) != c): raise TestFailure('mask 2d direct3')

b = numpy.array([[[3, 4.], [1., 0.]], [[3., 2], [-1, 9]]])
c = numpy.array([[[3, 4.], [0., 0.]], [[3., 0], [0, 0]]])

mask = Masking(numpy.array([[[0, 0.], [1., 1.]], [[0., 1], [1, 1]]], dtype='int8'))
if numpy.any(mask(b) != c): raise TestFailure('mask 3d direct')
mask = Masking(numpy.array([[[1, 1], [0., 0]], [[1, 0], [0, 0]]]))
if numpy.any(mask(b) != c): raise TestFailure('mask 3d direct')
mask = Masking(numpy.array([[[False, False], [True, True]], [[False, True], [True, True]]]))
if numpy.any(mask(b) != c): raise TestFailure('mask 3d direct')

c = mask(b, reusein=True, reuseout=True)
if id(b) != id(c):
    raise TestFailure('mask no copy')


#---------
# Padding
#---------

padding = Padding(left=1,right=20)
tod = Tod(numpy.ndarray((10,15)))
tod2 = padding(tod)
if tod2.shape != (10,36):          raise TestFailure('padding shapea')
if any_neq(tod2[:,0:1], 0):    raise TestFailure('padding left')
if any_neq(tod2[:,1:16], tod): raise TestFailure('padding middle')
if any_neq(tod2[:,16:], 0):    raise TestFailure('padding right')
tod3 = padding.transpose(tod2)
if tod3.shape != tod.shape: raise TestFailure('padding shapeb')
if any_neq(tod3,tod): raise TestFailure('padding end')

padding = Padding(left=1,right=(4,20))
tod = Tod.ones((10,(12,3)))
tod2 = padding(tod)
if tod2.shape != (10,41):          raise TestFailure('padding shapea')
if tod2.nsamples != (17,24):       raise TestFailure('padding nsamples')
if any_neq(tod2[:,0:1], 0):    raise TestFailure('padding left2')
if any_neq(tod2[:,1:13], tod[:,0:12]): raise TestFailure('padding middle2')
if any_neq(tod2[:,13:17], 0):  raise TestFailure('padding right2')
if any_neq(tod2[:,17:18], 0):  raise TestFailure('padding left2')
if any_neq(tod2[:,18:21], tod[:,12:]): raise TestFailure('padding middle2')
if any_neq(tod2[:,21:], 0):    raise TestFailure('padding right2')
tod3 = padding.transpose(tod2)
if tod3.shape != tod.shape: raise TestFailure('padding shapeb')
if any_neq(tod3,tod): raise TestFailure('padding end')


#-----
# Fft
#-----

n = 100
fft = FftHalfComplex(n)
a = numpy.random.random(n)+1
b = fft.transpose(fft.direct(a))
if any_neq(a, b): raise TestFailure('fft1')

a = numpy.random.random((3,n))+1
b = fft.transpose(fft.direct(a))
if any_neq(a, b): raise TestFailure('fft2')

tod = Tod(numpy.random.random((10,1000))+1)
fft = FftHalfComplex(tod.nsamples)
tod2 = fft.transpose(fft.direct(tod))
if any_neq(tod, tod2): raise TestFailure('fft3')

tod = Tod(numpy.random.random((10,1000))+1, nsamples=(100,300,5,1000-100-300-5))
fft = FftHalfComplex(tod.nsamples)
tod2 = fft.transpose(fft.direct(tod))
if any_neq(tod, tod2): raise TestFailure('fft4')


#--------------------------------
# Operations on AcquisitionModel
#--------------------------------

model = 1. + Identity()
if model((3,)) != [6] or model(3) != 6: raise TestFailure('Addition')

model = - Identity()
if model((3,)) != [-3] or model(3) != -3: raise TestFailure('-I')

model = (-2) * Identity()
if model((3,)) != [-6] or model(3) != -6: raise TestFailure('-2 I')

model = -(2 * Identity())
if model((3,)) != [-6] or model(3) != -6: raise TestFailure('-(2I)')

model = 1. - Identity()
if model((3,)) != [0] or model(3) != 0: raise TestFailure('Subtraction1')

model = 1. - 2 * Identity()
if model((3,)) != [-3] or model(3) != -3: raise TestFailure('Subtraction2')

model += 2*Identity()
if model((3,)) != [3] or model(3) != 3: raise TestFailure('+=')

model -= 2*Identity()
if model((3,)) != [-3] or model(3) != -3: raise TestFailure('-=')

model *= 2
if model((3,)) != [-6] or model(3) != -6: raise TestFailure('*=')


#------------------------------------
# Linearoperator -> AcquisitionModel
#------------------------------------

diagonal = numpy.arange(10.)
M = scipy.sparse.dia_matrix((diagonal, 0), shape=2*diagonal.shape)
model = asacquisitionmodel(M)
vec = numpy.ones(10)
if numpy.any(model(vec)   != diagonal): raise TestFailure()
if numpy.any(model.T(vec) != diagonal): raise TestFailure()


#------------------------------------
# AcquisitionModel -> Linearoperator
#------------------------------------

M2 = model.aslinearoperator()
if M2.__class__ != scipy.sparse.linalg.interface.LinearOperator: raise TestFailure()
if numpy.any(M2.matvec(vec)  != diagonal): raise TestFailure()
if numpy.any(M2.rmatvec(vec) != diagonal): raise TestFailure()


#---------------------------
# AcquisitionModelTranspose
#---------------------------

a = CompressionAverage(3)
b = a.T.T
if id(a) != id(b): raise TestFailure()

a = Identity()
if id(a) != id(a.T): raise TestFailure()

a = Diagonal([1,2,3])
if id(a) != id(a.T): raise TestFailure()

a = Masking([0,1,1,1])
if id(a) != id(a.T): raise TestFailure()

# dtype
a = Fft((3,4))
if a.dtype is not numpy.dtype('complex128'): raise TestFailure()
a = Scalar(complex(0,1))
if a.dtype is not numpy.dtype('complex128'): raise TestFailure()
a = Scalar(complex(1,0))
if a.dtype is not numpy.dtype('float64'): raise TestFailure()
a = Scalar(3) + Scalar(1)
if a.dtype is not numpy.dtype('float64'): raise TestFailure()
a = Diagonal(numpy.array([2,3,4]))
if a.dtype is not numpy.dtype('float64'): raise TestFailure()
a = Diagonal(numpy.array([2,complex(3,1),4]))
if a.dtype is not numpy.dtype('complex128'): raise TestFailure()
a = Masking(numpy.array([1, complex(2,2)]))
if a.dtype is not numpy.dtype('complex128'): raise TestFailure()
a = Masking([True, False])
if a.dtype is not numpy.dtype('float64'): raise TestFailure()
a = Masking(numpy.array([0,1,0], dtype='int8'))
if a.dtype is not numpy.dtype('float64'): raise TestFailure()

print 'OK.'

