#-------------------------------------------------------------
# Creation of a map using the double loop inference algorithm
#-------------------------------------------------------------
import os
from tamasis import *
import pyoperators

pyoperators.memory.verbose = False

# Specify the Frames observations as FITS files
#path = os.getenv('PACS_DATA')+'transpScan/'
#frames_files = [path+'1342185454_red_PreparedFrames.fits[10001:]',
#                path+'1342185455_red_PreparedFrames.fits[10001:]']

path = os.path.dirname(__file__) + '/../test/data/'
frames_files = [os.path.join(path, 'frames_blue.fits'),]

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
tod = filter_median(tod, 1000)

# We solve the equation y = H x, where y is the Tod, x the unknown map
# and H the acquisition model.
# To take into account bad samples such as glitches, we solve
# M y == M H x, M is the mask operator which sets bad samples values to 0
model = projection

# the prior tells the algorithm to look for a sparse solution in
# the corresponding space. For wavelets, the solution would have
# few non-zero coefficients, for FFT, it would have few non zero
# Fourier coefficient, etc ...
prior = Wavelet2("haar", shapein=model.shapein)

# tau drives the "level of sparseness" of the solution
# it can be a vector of size the number of prior coefficients
tau = 1e3
# generate an algorithm instance, it does not perform estimation
# but set all the information required by the algorithm to be run.
# A lot of information is accessible as attributes of algo.
algo = DoubleLoopAlgorithm(model, tod, prior, tau=tau)
# Calling DoubleLoopAlgorithm to estimate a map
map_dli = algo()
# save the resulting map
map_dli.save('map_dli.fits')
