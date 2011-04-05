#-------------------------------------------------------------
# Creation of a map using the regularised least square method
#-------------------------------------------------------------
import os
import numpy as np
from tamasis import *

# Specify the Frames observations as FITS files
path = 'mydir/'
frames_files = [path+'1342202088_red_level1Frames.fits[1501:7000]',
                path+'1342202089_red_level1Frames.fits[1501:6500]']

# Setup the instrument and pointings for these files
# The policy for each read-out can be 'mask', 'remove' or 'keep'
obs = PacsObservation(frames_files,
                      policy_inscan='keep',
                      policy_turnaround='keep',
                      policy_other='remove',
                      policy_invalid='mask')

# Read the Time Ordered Data: the signal and mask
tod = obs.get_tod(flatfielding=False,
                  subtraction_mean=True,
		  masks='badpixels, blindpixels, nonscience, sa_badpix, saturation, saturation_high, saturation_low, uncleanchop')

# Get the projection matrix used for the map-level deglitching
# 'oversampling=False' means that the acquisition model will not
# try to sample at a frequency higher than that of the observation
# (10Hz for prime mode, 5Hz for parallel mode)
projection = Projection(obs,
                        method='sharp',
                        oversampling=False,
                        npixels_per_sample=6)

# Remove low frequency drifts. The specified window length is the
# number of samples used to compute the median (unlike HCSS, where
# half the length should be specified).
# If these are all masked, an interpolated value will be computed 
# once the whole filtered timeline is computed.
hpf_length = 10000
tod = filter_median(tod, hpf_length)

# Map-level deglitching using the MAD (median absolute deviation to
# the mean). We highpass filter the Tod using a short filtering 
# window, to remove the low-frequency offsets
tod_glitch = filter_median(tod, 100)
tod.mask = deglitch_l2mad(tod_glitch, projection, nsigma=25.)

# save the Tod as a FITS file
tod.save('tod_preprocessed.fits')

# We solve the equation y = H x, where y is the Tod, x the unknown map
# and H the acquisition model.
# To take into account bad samples such as glitches, we solve
# M y == M H x, M is the mask operator which sets bad samples values to 0
projection = Projection(obs,
                        method='sharp',
                        npixels_per_sample=5)
response = ResponseTruncatedExponential(obs.pack(
    obs.instrument.detector_time_constant) / obs.SAMPLING_PERIOD)
compression = CompressionAverage(obs.slice.compression_factor)
masking = Masking(tod.mask)
model = masking * compression * response * projection

# The naive map is given by
map_naive = mapper_naive(tod, model)

length = np.asarray(2**np.ceil(np.log2(np.array(tod.nsamples) + 200)), dtype='int')
invNtt = InvNtt(length, obs.get_filter_uncorrelated())
fft = FftHalfComplex(length)
padding = Padding(left=invNtt.ncorrelations, right=length-tod.nsamples-invNtt.ncorrelations)
weight = padding.T * fft.T * invNtt * fft * padding

# The regularised least square map is obtained by minimising the criterion
# J(x) = ||y-Hx||^2 + hyper ||Dx||^2, the first ||.||^2 being the N^-1 norm
# it is equivalent to solving the equation (H^T H + hyper D^T D ) x = H^T y
hyper = 1
map_tamasis = mapper_rls(tod, model,
                         weight=weight,
                         unpacking=Unpacking(map_naive.coverage==0),
                         tol=1.e-5,
                         hyper=hyper)

# save the map a as fits file
map_tamasis.save('map_tamasis_h'+str(hyper)+'_hpf'+str(hpf_length)+'.fits')
