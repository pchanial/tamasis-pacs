import os
#from tamasis import *

fa =FitsArchive(reader=FitsArchive.HCSS_READER)
#frames = fa.load(os.getenv('PACS_DATA')+'/M81prime/1342186085_red_PreparedFrames.fits')
frames = fa.load(tamasis_dir+'tests/frames_blue.fits')
map = tamasisPhotProject(frames, updateFrames=True, deglitching='l2mad', medianFiltering=10000, flatfielding=True, ds9=True)
print 'Done.'
