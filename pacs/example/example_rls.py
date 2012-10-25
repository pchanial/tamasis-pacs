#-------------------------------------------------------------
# Creation of a map using the regularised least square method
#-------------------------------------------------------------
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
                  subtraction_mean=True)

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
tod.mask = deglitch_l2mad(tod_glitch, projection, nsigma=5.)

# save the Tod as a FITS file
tod.save('tod_preprocessed.fits')

# We solve the equation y = H x, where y is the Tod, x the unknown map
# and H the acquisition model.
# To take into account bad samples such as glitches, we solve
# M y == M H x, M is the mask operator which sets bad samples values to 0
compression = CompressionAverageOperator(obs.slice.compression_factor)
masking = MaskOperator(tod.mask)
model = masking * compression * projection

# The regularised least square map is obtained by minimising the criterion
# J(x) = ||y-Hx||^2 + hyper ||Dx||^2, the first ||.||^2 being euclidian norm
# it is equivalent to solving the equation (H^T H + hyper D^T D ) x = H^T y
hyper = 0.1
map_rls = mapper_rls(tod, model,
                     unpacking=UnpackOperator(projection.get_mask()),
                     tol=1.e-6,
                     hyper=hyper)
map_rls.save('map_rls.fits')


