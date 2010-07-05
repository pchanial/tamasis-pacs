import os

tamasis_dir = '/home/pchanial/work/tamasis/tamasis-1.0.4/'

def get_filename():
    filename = '/tmp/'+os.getenv('USER')+long(time.time()*1000)
    
def tamasis_preprocessor(frames, flatfielding=True, subtraction_mean=True, median_filtering=200, deglitching=None, nsigma=5.):
    if frames.class is herschel.pacs.signal.Frames:
        frames = (frames,)

    cmd = tamasis_dir + 'pacs/scripts/hcss_interface.py'
    
    # write frames to disk
    filename = get_filename()+'_tod_'

    for ifile in range(len(frames)):
        file = file + str(ifile) + '.fits'
        fa.save(file, frames[ifile])
        cmd = cmd + ' ' + file

    print cmd
        
def tamasis_mapper(frames, header=None, resolution=None, npixels_per_sample=None):
    pass
