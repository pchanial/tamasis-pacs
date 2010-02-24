from tamasis import Identity, CompressionAverage, PacsProjectionSharpEdges, PacsMultiplexing, Masking, PacsObservation, PacsSimulation, Map, Tod, hcss_photproject
from copy import copy

datadir = '/home/pchanial/work/pacs/data/'
pacs = PacsObservation(filename=datadir+'transparent/NGC6946/1342184520_blue', \
                       first=20000,                 \
                       last=86000,                  \
                       resolution=3.,               \
                       fine_sampling_factor=1,      \
                       keep_bad_detectors=False,    \
                       npixels_per_sample=6)
tod = pacs.get_tod()

telescope    = Identity('Telescope PSF')
projection   = PacsProjectionSharpEdges(pacs)
#multiplexing = PacsMultiplexing(pacs)
multiplexing = CompressionAverage(pacs.fine_sampling_factor, 'Multiplexing')
crosstalk    = Identity('Crosstalk')
compression  = CompressionAverage(pacs.compression_factor)
masking      = Masking(tod.mask)

model = masking * compression * crosstalk * multiplexing * projection * telescope
print model

mymap = copy(model.transpose(tod))
tod[:] = 1.
weights = model.transpose(tod)
mymap /= weights
mymap.writefits('/home/pchanial/work/tamasis/ngc6946.fits')

#map = hcss_photproject(pacs)

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
#                            keepBadPixels      = True)

