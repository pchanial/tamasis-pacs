import os
from tamasishcss import *

fa = FitsArchive(reader=FitsArchive.HCSS_READER)
frames = fa.load(tamasis_dir+'pacs/test/data/frames_blue.fits')
map = tamasisPhotProject(frames,
                         updateFrames=True,
                         deglitching='l2mad',
                         medianFiltering=10000,
                         flatfielding=True,
                         policyTurnaround='keep',
                         npixelsPerSample=6,
                         ds9=False)
Display(map)
print 'Done.'
