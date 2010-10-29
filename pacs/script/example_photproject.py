#---------------------------------------------------
# Creation of a map using the photproject algorithm
#---------------------------------------------------
from tamasis import *

# Specify the Frames observations as FITS files
frames_files = ['myFrame_inscan.fits', 'myFrame_crossscan.fits']

# Setup the instrument and pointings for these files
# The policy for each read-out can be 'mask', 'remove' or 'keep'
obs = PacsObservation(frames_files,
                      policy_inscan='keep',
                      policy_turnaround='remove',
                      policy_other='remove',
                      policy_invalid='mask')

# Get the projection matrix
# 'oversampling=False' means that the acquisition model will not
# try to sample at a frequency higher than that of the observation
# (10Hz for prime mode, 5Hz for parallel mode)
projection = Projection(obs, 
                        oversampling=False,
                        npixels_per_sample=6)

# Read the Time Ordered Data: the signal and mask
tod = obs.get_tod(flatfielding=True,
                  subtraction_mean=True,
                  unit='Jy/detector')

# Remove very low frequency drift by removing a fitted polynomial
# of degree 6
tod_filtered = filter_polynomial(tod, 6)

# Remove low frequency drifts. The specified window length is the
# number of samples used to compute the median (unlike HCSS, where
# half the length should be specified).
# If these are all masked, an interpolated value will be computed 
# once the whole filtered timeline is computed.
tod = filter_median(tod, 10000)

# Map-level deglitching using the MAD (median absolute deviation to
# the mean). We highpass filter the Tod using a short filtering 
# window, to remove the low-frequency offsets
tod_glitch = filter_median(tod, 100)
tod.mask = deglitch_l2mad(tod_glitch, projection, nsigma=5.)

# We solve the equation y = H x, where y is the Tod, x the unknown map
# and H the acquisition model.
# To take into account bad samples such as glitches, we solve
# M y == M H x, M is the mask operator which sets bad samples values to 0
masking = Masking(tod.mask)
model = masking * projection
tod = masking(tod)

# Get the photproject map, which is equal to:
#    H^T y / H^T 1 (^T is tranpose operator)
map_naive = mapper_naive(tod, projection)

# and display it in a matplotlib's figure
map_naive.imshow()

# or using ds9
map_naive.ds9()

# inspect the glitch spatial distribution on the map:
map_mask = projection.T(numpy.asarray(tod.mask, float))
map_mask.imshow()

