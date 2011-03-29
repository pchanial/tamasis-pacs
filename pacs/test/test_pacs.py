import numpy as np
import pyfits
import os
import tamasis

from tamasis import *
from uuid import uuid1

class TestFailure(Exception): pass

tamasis.var.verbose = True

data_dir = os.path.dirname(__file__) + '/data/'

# slice[10:19]
obs = PacsObservation(data_dir+'frames_blue.fits[11:20]')
tod = obs.get_tod(flatfielding=False)

filename = 'obs-'+str(uuid1())+'.fits'
try:
    obs.save(filename, tod)
    obs2 = PacsObservation(filename)
    tod2 = obs2.get_tod(raw=True)
    status2 = obs2.status
finally:
    try:
        os.remove(filename)
    except:
        pass
for field in obs.status.dtype.names:
    status = obs.status[10:20]
    if isinstance(status[field][0], str):
        if np.any(status[field] != status2[field]): raise TestFailure('Status problem with: '+field)
    elif not np.allclose(status[field], status2[field]): raise TestFailure('Status problem with: '+field)
if not np.allclose(tod, tod2): raise TestFailure()
if not np.all(tod.mask == tod2.mask): raise TestFailure()

# all observation
obs = PacsObservation(data_dir+'frames_blue.fits')

# get_mask
proj = Projection(obs, npixels_per_sample=6, oversampling=False)
o = Tod.ones(proj.shapeout)
nocoverage = mapper_naive(o, proj).coverage == 0
if any_neq(nocoverage, proj.get_mask().magnitude): raise TestFailure()
tod = obs.get_tod()

# packed projection
proj2 = Projection(obs, npixels_per_sample=6, packed=True, oversampling=False)
proj3 = proj2 * Unpacking(proj2.get_mask()).T

if any_neq(proj.T(tod), proj3.T(tod)): raise TestFailure()


filename = 'obs-'+str(uuid1())+'.fits'
try:
    obs.save(filename, tod)
    obs2 = PacsObservation(filename)
    tod2 = obs2.get_tod(raw=True)
    status2 = obs2.status
finally:
    try:
        os.remove(filename)
    except:
        pass
for field in obs.status.dtype.names:
    if isinstance(obs.status[field][0], str):
        if np.any(obs.status[field] != status2[field]): raise TestFailure('Status problem with: '+field)
    elif not np.allclose(obs.status[field], status2[field]): raise TestFailure('Status problem with: '+field)
if not np.allclose(tod, tod2): raise TestFailure()
if not np.all(tod.mask == tod2.mask): raise TestFailure()

telescope    = Identity(description='Telescope PSF')
projection   = Projection(obs, resolution=3.2, oversampling=False, npixels_per_sample=6)
multiplexing = CompressionAverage(obs.instrument.fine_sampling_factor, description='Multiplexing')
crosstalk    = Identity(description='Crosstalk')
compression  = CompressionAverage(obs.slice.compression_factor)
masking      = Masking(tod.mask)

model = masking * crosstalk * multiplexing * projection * telescope
print model

# naive map
tod.inunit('Jy/arcsec^2')
backmap = model.T(tod.magnitude)
unity = Tod.ones(tod.shape, nsamples=tod.nsamples)
weights = model.T(unity)
map_naive = Map(backmap / weights, unit='Jy/arcsec^2')

header = projection.header
header2 = header.copy()
header2['NAXIS1'] += 500
header2['CRPIX1'] += 250
projection2 = Projection(obs, header=header2, oversampling=False)
map_naive2 = mapper_naive(tod, projection2)
map_naive2.inunit('Jy/arcsec^2')
map_naive3 = map_naive2[:,250:header['NAXIS1']+250]
if any_neq(map_naive, map_naive3, 2.e-7): raise TestFailure('mapper_naive, with custom header')

# test compatibility with photproject
tod = obs.get_tod('Jy/arcsec^2', flatfielding=False, subtraction_mean=False)
map_naive4 = mapper_naive(tod, projection)
hdu_ref = pyfits.open(data_dir + 'frames_blue_map_hcss_photproject.fits')[1]
map_ref = Map(hdu_ref.data, hdu_ref.header, unit=hdu_ref.header['qtty____']+'/pixel')
std_naive = np.std(map_naive4[40:60,40:60])
std_ref = np.std(map_ref[40:60,40:60])
relerror = abs(std_naive-std_ref) / std_ref
if relerror > 0.025: raise TestFailure('Tncompatibility with HCSS photproject: ' + str(relerror*100)+'%.')

map_naive_ref = Map(data_dir + 'frames_blue_map_naive.fits')
obs = PacsObservation(data_dir + 'frames_blue.fits')
obs.pointing.chop[:] = 0
projection = Projection(obs, header=map_naive_ref.header, oversampling=False, npixels_per_sample=6)
tod = obs.get_tod(flatfielding=False)
masking = Masking(tod.mask)
model = masking * projection
map_naive = mapper_naive(tod, model)
if any_neq(map_naive, map_naive_ref, 1.e-11): raise TestFailure()

obs_rem = PacsObservation(data_dir + 'frames_blue.fits', policy_bad_detector='remove')
obs_rem.pointing.chop[:] = 0
projection_rem = Projection(obs_rem, header=map_naive.header, oversampling=False, npixels_per_sample=7)
tod_rem = obs_rem.get_tod(flatfielding=False)
masking_rem = Masking(tod_rem.mask)
model_rem = masking_rem * projection_rem
map_naive_rem = mapper_naive(tod_rem, model_rem)
if any_neq(map_naive, map_naive_rem, 1.e-11): raise TestFailure()

# pack/unpack
for channel, nrows, ncolumns in ('red',16,32), ('blue',32,64):
    obs = PacsSimulation(Pointing(0., 0., 0., 0.), channel)
    for dt in (np.uint8, np.uint16, np.uint32, np.uint64):
        a = np.arange(nrows*ncolumns*3, dtype=dt).reshape((nrows,ncolumns,-1))
        p = obs.pack(a)
        if np.any(a[1,0:16,:] != p[16:32,:]):
            raise TestFailure()

        u = obs.unpack(p)
        if np.any(a != u):
            raise TestFailure()

# pmatrix with nsamples_per_pixel == 0
proj = Projection(obs, npixels_per_sample=2)
header = proj.header.copy()
header['crval1']=245.998427916727+1
proj2 = Projection(obs, header=header)
m = np.ones((header['NAXIS2'],header['NAXIS1']))
t = proj2(m)
if any_neq(minmax(t), [0,0]): raise TestFailure()
t[:] = 1
if any_neq(minmax(proj2.T(t)), [0,0]): raise TestFailure()


