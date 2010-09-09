import os
from tamasishcss import *

fa = FitsArchive(reader=FitsArchive.HCSS_READER)
path = '/mnt/herschel1/mapmaking/data/pacs/M81prime/'
#frames = (fa.load(path+'1342186085_red_PreparedFrames.fits'), \
#          fa.load(path+'1342186086_red_PreparedFrames.fits'))
frames = fa.load(tamasis_dir+'pacs/test/data/frames_blue.fits')
map = tamasisPhotProject(frames, \
                         updateFrames=True, \
                         deglitching='l2mad', \
                         medianFiltering=10000, \
                         flatfielding=True, \
                         framePolicyTurnaround='keep', \
                         npixelsPerSample=6, \
                         ds9=False)
Display(map)
print 'Done.'
