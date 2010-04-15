from   matplotlib.pyplot import clim, figure, plot, show
from   tamasis import *
import os

do_plot = True

datadir = os.getenv('PACS_DATA')+'/transpScan/'
pacs = PacsObservation(filename=[datadir+'1342184598_blue_PreparedFrames.fits',
                                 datadir+'1342184599_blue_PreparedFrames.fits'],
                       resolution=3.2,
                       fine_sampling_factor=1,
                       npixels_per_sample=5, 
                       keep_bad_detectors=False)

telescope    = Identity('Telescope PSF')
projection   = Projection(pacs)
multiplexing = CompressionAverage(1, 'Multiplexing')
crosstalk    = Identity('Crosstalk')
compression  = CompressionAverage(4)

model = compression * crosstalk * multiplexing * projection * telescope

# read the Tod off the disk
tod40Hz = pacs.get_tod()

# remove drift
x = numpy.arange(tod40Hz.shape[0])
slope = scipy.polyfit(x, tod40Hz, deg=6)
drift = tod40Hz.copy()
for i in xrange(tod40Hz.shape[1]):
    drift[:, i] = scipy.polyval(slope[:, i], x)

idetector=5
if do_plot:
    figure()
    plot(tod40Hz[:,idetector])
    plot(drift[:,idetector],'r')
    show()

tod40Hz -= drift

mask_before = tod40Hz.mask.copy('a')

# second level deglitching
deglitch_l2mad(tod40Hz, projection)

if do_plot:
    figure()
    plot(tod40Hz[:,idetector])
    index=numpy.where(tod40Hz.mask[:,idetector])
    plot(index,tod40Hz.data[index,idetector],'ro')
    show()

masking   = Masking(tod40Hz.mask)
model40Hz = masking * projection
map_naive40Hz = naive_mapper(tod40Hz, model40Hz)

if do_plot:
    map_naive40Hz.imshow()
    clim(-0.00002,0.00002)

# compressed TOD
tod = compression.direct(tod40Hz).copy()

# naive map
map_naive = naive_mapper(tod, model)

