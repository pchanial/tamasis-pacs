import numpy
import pyfits

from tamasis import *

class TestFailure(Exception): pass

datadir = tamasis_dir + 'tests/'
obs = PacsObservation(datadir+'frames_blue.fits',
                      fine_sampling_factor=1, 
                      detector_policy='remove')

tod = obs.get_tod()

telescope    = Identity('Telescope PSF')
projection   = Projection(obs, resolution=3.2, oversampling=False, npixels_per_sample=6)
multiplexing = CompressionAverage(obs.fine_sampling_factor, 'Multiplexing')
crosstalk    = Identity('Crosstalk')
compression  = CompressionAverage(obs.compression_factor)
masking      = Masking(tod.mask)

model = masking * crosstalk * multiplexing * projection * telescope
print model
model = projection
print model

# naive map
backmap = model.transpose(tod)
unity = Tod.ones(tod.shape, nsamples=tod.nsamples)
weights = model.transpose(unity)
map_naive = backmap / weights

header = projection.header
header2 = header.copy()
header2['NAXIS1'] += 500
header2['CRPIX1'] += 250
projection2 = Projection(obs, header=header2, oversampling=False)
map_naive2 = mapper_naive(tod, projection2)
map_naive3 = map_naive2[:,250:header['NAXIS1']+250]
if any_neq(map_naive, map_naive3, 7): raise TestFailure('mapper_naive, with custom header')

# test compatibility with photproject
tod = obs.get_tod('Jy/arcsec^2')
map_naive4 = mapper_naive(tod, projection, unit='Jy/pixel')
hdu_ref = pyfits.open(datadir + 'frames_blue_map_hcss_photproject.fits')[1]
map_ref = Map(hdu_ref.data, hdu_ref.header, unit=hdu_ref.header['qtty____']+'/pixel')
std_naive = numpy.std(map_naive4[40:60,40:60])
std_ref = numpy.std(map_ref[40:60,40:60])
relerror = abs(std_naive-std_ref) / std_ref
if relerror > 0.025: raise TestFailure('Tncompatibility with HCSS photproject: ' + str(relerror*100)+'%.')

print 'OK.'

#ra0  = 20.
#dec0 = 0.1
#time = numpy.arange(0.,100., 1./40)
#simulation = PacsSimulation(inputmap           = 
#                            time               = time \
#                            ra                 = numpy.linspace(ra0, ra0+0.1, nsamples)   \
#                            dec                = numpy.linspace(dec0, dec0+0.1, nsamples) \
#                            pa                 = numpy.zeros(nsamples) \
#                            chop               = numpy.zeros(nsamples) \
#                            array              = 'blue'        \
#                            npixelsPerSample   = 9             \
#                            observingMode      = 'transparent' \
#                            fineSamplingFactor = 1             \
#                            compressionFactor  = 1             \
#                            keepBadDetectors   = True)

