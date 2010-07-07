import os
import tamasishcss as tm

tamasis_dir = '/home/pchanial/work/tamasis/tamasis-1.0.4/'

fa =FitsArchive(reader=FitsArchive.HCSS_READER)
#frames = fa.load(os.getenv('PACS_DATA')+'/transpScan/1342184598_blue_PreparedFrames.fits')
frames = fa.load(tamasis_dir+'tests/frames_blue.fits')
map = tm.tamasis_mapper(frames)
