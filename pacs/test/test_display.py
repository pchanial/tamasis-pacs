#!/usr/bin/env python
import numpy as np
import os
import tamasis

from   matplotlib.pyplot import clim, figure, plot, show, ioff
from   tamasis import *

tamasis.var.verbose = False

datadir  = os.getenv('PACS_DATA', '')+'/transpScan/'
datafile = [datadir+'1342184598_blue_PreparedFrames.fits',
            datadir+'1342184599_blue_PreparedFrames.fits']

if not all(map(os.path.exists, datafile)):
    print('The data files are not found: ' + ', '.join(datafile))
    exit(0)

pacs = PacsObservation([datafile[0]+'[6065:20000]', datafile[1]+'[6066:20001]'],
                       fine_sampling_factor=1, calblock_extension_time=0.)

telescope    = Identity(description='Telescope PSF')
projection   = Projection(pacs, resolution=3.2, npixels_per_sample=5)
multiplexing = CompressionAverage(1, description='Multiplexing')
crosstalk    = Identity(description='Crosstalk')
compression  = CompressionAverage(4)

model = compression * crosstalk * multiplexing * projection * telescope

# read the Tod off the disk
tod40Hz = pacs.get_tod()

# remove drift
tod40Hz_filtered = filter_polynomial(tod40Hz, 6)
drift = tod40Hz - tod40Hz_filtered

tod40Hz = filter_median(tod40Hz_filtered, 10000)

# second level deglitching
tod40Hz.mask = deglitch_l2mad(tod40Hz, projection)

idetector=5
plot_tod((tod40Hz+drift)[idetector])
plot(drift[idetector], 'r')

masking   = Masking(tod40Hz.mask)
model40Hz = masking * projection
map_naive40Hz = mapper_naive(tod40Hz, model40Hz)

map_naive40Hz.imshow()
# clim doesn't work anymore with AnnotatedImage
#clim(-0.00002,0.00002)
