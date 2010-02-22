from tamasis import Identity, CompressionAverage, PacsProjectionSharpEdges, PacsMultiplexing, PacsObservation, applymask

def hcssPhotProject(pacs):
    """
    Returns a map, as calculated by HCSS's PhotProject
    Inplace substitutions can affect inputs.
    """
    if pacs.fineSamplingFactor != 1:
        raise ValueError('Fine sampling factor should be 1 for hcssPhotProject.') # or add decimation
    model = PacsProjectionSharpEdges(pacs) 
    signal, mask = pacs.get_timeline()
    backmap = model.transpose(applymask(signal, mask)).copy()
    signal[:] = 1.
    weights = model.transpose(applymask(signal, mask))
    return backmap/weights

datadir = '/home/pchanial/work/pacs/data/'
pacs = PacsObservation(filename=datadir+'transparent/NGC6946/1342184520_blue', \
                       first=20000,                 \
                       last=86000,                  \
                       fineSamplingFactor=1,        \
                       npixelsPerSample=9)

telescope    = Identity()
projection   = PacsProjectionSharpEdges(pacs)
#multiplexing = PacsMultiplexing(pacs)
multiplexing = CompressionAverage(pacs.fineSamplingFactor)
crosstalk    = Identity()
compression  = CompressionAverage(pacs.compressionFactor)

model = compression * crosstalk * multiplexing * projection * telescope
print model

map = hcssPhotProject(pacs)

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

