!#----------------------------------------------
# Creation of a map using the MADmap algorithm
#----------------------------------------------
import numpy as np
import os
from tamasis import *

# Specify the Frames observations as FITS files
path = os.getenv('PACS_DATA')+'transpScan/'
frames_files = [path+'1342185454_red_PreparedFrames.fits[10001:]',
                path+'1342185455_red_PreparedFrames.fits[10001:]']

# Setup the instrument and pointings for these files
# The policy for each read-out can be 'mask', 'remove' or 'keep'
obs = PacsObservation(frames_files,
                      policy_inscan='keep',
                      policy_turnaround='keep',
                      policy_other='remove',
                      policy_invalid='mask')

# Read the Time Ordered Data: the signal and mask
tod = obs.get_tod(flatfielding=True,
                  subtraction_mean=True,
                  unit='Jy/detector')

# Get the projection matrix used for the map-level deglitching
# 'downsampling=True' means that the acquisition model will not
# sample at the instrument frequency of 40Hz, but at the compressed frequency
# (10Hz for prime mode, 5Hz for parallel mode)
projection = obs.get_projection_operator(
                 downsampling=True,
                 method='sharp',
                 npixels_per_sample=5)

# Remove low frequency drifts. The specified window length is the
# number of samples used to compute the median (unlike HCSS, where
# half the length should be specified).
# If these are all masked, an interpolated value will be computed 
# once the whole filtered timeline is computed.
tod = filter_median(tod, 1000)

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
projection = obs.get_projection_operator(
                 downsampling=True,
                 method='nearest',
                 npixels_per_sample=5)
masking = MaskOperator(tod.mask)
model = masking * projection

# Get the filter operator N^-1
invntt = InvNttOperator(obs)

# The naive map is given by
map_naive = mapper_naive(tod, model)

# The MADmap map is obtained by minimising the criterion
# J(x) = ||y-Hx||^2, ||.||^2 being the N^-1 norm
# it is equivalent to solving the equation H^T N^-1 H x = H^T N^-1 y
map_madmap = mapper_ls(tod, model,
                       invntt=invntt,
                       unpacking=UnpackOperator(projection.get_mask()),
                       M=DiagonalOperator(1/map_naive.coverage),
                       tol=1.e-5)
map_madmap.save('map_madmap.fits')
