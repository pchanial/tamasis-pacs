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
tod = obs.get_tod()
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
tod.mask[:] = 0
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
obs = PacsObservation(data_dir + 'frames_blue.fits', policy_detector='mask')
obs.pointing.chop[:] = 0
projection = Projection(obs, header=map_naive_ref.header, oversampling=False, npixels_per_sample=6)
tod = obs.get_tod(flatfielding=False)
masking = Masking(tod.mask)
model = masking * projection
map_naive = mapper_naive(tod, model)
if any_neq(map_naive, map_naive_ref, 1.e-11): raise TestFailure()

obs_rem = PacsObservation(data_dir + 'frames_blue.fits', policy_detector='remove')
obs_rem.pointing.chop[:] = 0
projection_rem = Projection(obs_rem, header=map_naive.header, oversampling=False, npixels_per_sample=7)
tod_rem = obs_rem.get_tod(flatfielding=False)
masking_rem = Masking(tod_rem.mask)
model_rem = masking_rem * projection_rem
map_naive_rem = mapper_naive(tod_rem, model_rem)
if any_neq(map_naive, map_naive_rem, 1.e-11): raise TestFailure()

