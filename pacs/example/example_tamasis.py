#-------------------------------------------------------------
# Creation of a map using the regularised least square method
#-------------------------------------------------------------
import os
import scipy
from tamasis import *

# Specify the Frames observations as FITS files
# The optional brackets can be used to remove the beginning or the end of an
# an observation
path = os.getenv('PACS_DATA')+'transpScan/'
frames_files = [path+'1342184598_blue_PreparedFrames.fits[10001:]',
                path+'1342184599_blue_PreparedFrames.fits[10001:]']

# Setup the instrument and pointings for these files
# The policy for each read-out can be 'mask', 'remove' or 'keep'
obs = PacsObservation(frames_files,
                      policy_inscan='keep',
                      policy_turnaround='keep',
                      policy_other='mask',
                      policy_invalid='mask')

# Read the Time Ordered Data: the signal and mask
# Set flatfielding to False if the observation has already been flat-fielded.
# Set subtraction_mean to False to keep the offsets. The mask keyword can be
# a list of the names of the masks to be combined. By default, the activated
# masks are selected.
tod = obs.get_tod(flatfielding=True,
                  subtraction_mean=True,
		  masks='saturation')

# Get the projection matrix used for the map-level deglitching
# 'downsampling=True' means that the acquisition model will not
# sample at the instrument frequency of 40Hz, but at the compressed frequency
# (10Hz for prime mode, 5Hz for parallel mode)
projection = ProjectionOperator(obs.get_pointing_matrix(
                 header=obs.get_map_header(),
                 npixels_per_sample=5,
                 method='sharp',
                 downsampling=True))

# Map-level deglitching using the MAD (median absolute deviation to
# the mean). We highpass filter the Tod using a short filtering 
# window, to remove the low-frequency offsets
tod_glitch = filter_median(tod, 100)
tod.mask = deglitch_l2mad(tod_glitch, projection, nsigma=20.)

# Inspect that the glitches have correctly been removed in the timeline of
# detector #40
plot_tod(tod[40])

# Backproject the Tod's mask on a map, to ensure that glitches are uniformly
# distributed, i.e. that the bright sources have not been impacted
projection.T(tod.mask-tod_glitch.mask).imshow()

# we don't need the coarse pointing matrix anymore
del projection

# remove the very low frequency drift by subtraction a polynomial of degree 3
tod = filter_polynomial(tod, 3)

# Remove low frequency drifts. The specified window length is the
# number of samples used to compute the median (unlike HCSS, where
# half the length should be specified).
# If these are all masked, an interpolated value will be computed 
# once the whole filtered timeline is computed.
# Remove very low frequency drifts

hpf_length = 1000
tod = filter_median(tod, hpf_length)

# Inspect the tod by displaying the Tod as an image: X for time and Y for
# the detector number. In the following a example a stride of 10 is used
# to display one sample out of 10, because the display is memory intensive
tod[:,::10].imshow()

# Save the Tod as a FITS file
tod.save('tod_preprocessed.fits')

# Save the Observation and the Tod as a FITS file, that can later be read
# as a PacsObservation
obs.save('obs_preprocessed.fits', tod)

# We solve the equation y = H x, where y is the Tod, x the unknown map
# and H the acquisition model.
# To take into account bad samples such as glitches, we solve
# M y == M H x, M is the mask operator which sets bad samples values to 0
projection = ProjectionOperator(obs.get_pointing_matrix(
                 header=obs.get_map_header(resolution=3.2),
                 npixels_per_sample=5,
                 method='sharp'))
response = ConvolutionTruncatedExponentialOperator(obs.pack(
                 obs.instrument.detector.time_constant) / obs.SAMPLING_PERIOD)
compression = CompressionAverageOperator(obs.slice.compression_factor)
masking = MaskOperator(tod.mask)
model = masking * compression * response * projection

# The naive map is given by
map_naive = mapper_naive(tod, model)

# Inspect the naive map
ds9(map_naive)

# The regularised least square map is obtained by minimising the criterion
# J(x) = ||y-Hx||^2 + hyper ||Dx||^2, the first ||.||^2 being the N^-1 norm
# it is equivalent to solving the equation (H^T H + hyper D^T D ) x = H^T y
hyper = 0.1
map_tamasis = mapper_rls(tod, model,
                         invntt=InvNttOperator(obs),
                         unpacking=UnpackOperator(projection.get_mask()),
                         tol=1.e-5,
                         solver=scipy.sparse.linalg.bicgstab,
                         hyper=hyper)

# save the map as FITS file
map_tamasis.save('map_tamasis_h'+str(hyper)+'_hpf'+str(hpf_length)+'.fits')
