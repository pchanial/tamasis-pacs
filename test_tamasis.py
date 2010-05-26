from tamasis import *
import numpy

class TestFailure(Exception): pass

pacs = PacsObservation(filename=tamasis_dir+'tests/frames_blue.fits',
                       fine_sampling_factor=1, 
                       keep_bad_detectors=False)

tod = pacs.get_tod()

telescope    = Identity('Telescope PSF')
projection   = Projection(pacs, resolution=3.2, finer_sampling=False, npixels_per_sample=6)
multiplexing = CompressionAverage(pacs.fine_sampling_factor, 'Multiplexing')
crosstalk    = Identity('Crosstalk')
compression  = CompressionAverage(pacs.compression_factor)
masking      = Masking(tod.mask)

model = masking * crosstalk * multiplexing * projection * telescope
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
projection2 = Projection(pacs, header=header2, finer_sampling=False)
map_naive2 = mapper_naive(tod, projection2)
map_naive3 = map_naive2[:,250:header['NAXIS1']+250]
if any_neq(map_naive, map_naive3, 7): raise TestFailure('mapper_naive, with custom header')

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

