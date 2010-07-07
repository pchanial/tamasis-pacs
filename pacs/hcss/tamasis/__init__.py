import os
import herschel
import time

__all__ = ['tamasis_dir', 'tamasis_preprocessor', 'tamasis_mapper']
tamasis_dir = '/home/pchanial/work/tamasis/tamasis-1.0.4/'
False = 0

fa = herschel.ia.io.fits.FitsArchive()
simpleFitsReader = herschel.ia.toolbox.util.SimpleFitsReaderTask()

def _get_filename():
    return '/tmp/' + os.getenv('USER') + '-' + str(long(time.time()*1000))

def _get_cmd_preprocessor(frames, flatfielding=False, subtraction_mean=False, median_filtering=None, deglitching=None, nsigma=5., npixels_per_sample=None):

    if frames.class is herschel.pacs.signal.Frames:
        frames = (frames,)


    cmd = tamasis_dir + 'pacs/scripts/hcss_interface.py'
    
    # write frames to disk
    filename = _get_filename()+'_tod_'
    files = []

    for ifile in range(len(frames)):
        file = filename + str(ifile) + '.fits'
        fa.save(file, frames[ifile])
        files.append(file)
        cmd += ' ' + file

    if flatfielding:
        cmd += ' --flatfielding'

    if subtraction_mean:
        cmd += ' --subtraction-mean'

    if median_filtering is not None:
        cmd += ' --median-filtering ' + str(median_filtering)

    if deglitching is not None:
        cmd += ' --deglitching ' + deglitching + ' --nsigma ' + str(nsigma)
    
    if npixels_per_sample is not None:
        cmd += ' --npixels-per-sample ' + str(npixels_per_sample)
    
    return cmd, files


def tamasis_preprocessor(frames, flatfielding=False, subtraction_mean=False, median_filtering=None, deglitching=None, nsigma=5., npixels_per_sample=None):
    cmd, files = _get_cmd_preprocessor(frames, flatfielding=flatfielding, subtraction_mean=subtraction_mean, median_filtering=median_filtering, deglitching=deglitching, nsigma=nsigma, npixels_per_sample=npixels_per_sample)
    cmd += ' --output-tod'
    os.system(cmd)
    
    for file, frame in zip(files, frames):
        tod = simpleFitsReader(file+'_tod.fits')
        print file, ': shape=', tod['PrimaryImage'].data.dimensions
        frame.signal = tod['PrimaryImage'].data


def tamasis_mapper(frames, flatfielding=False, subtraction_mean=False, median_filtering=None, deglitching=None, nsigma=5., header=None, resolution=None, npixels_per_sample=None, ds9=False, update_frames=False):
    cmd, files = _get_cmd_preprocessor(frames, flatfielding=flatfielding, subtraction_mean=subtraction_mean, median_filtering=median_filtering, deglitching=deglitching, nsigma=nsigma)

    print 'check tod units'
    if header is not None:
        cmd += ' --header '+header

    if resolution is not None:
        cmd += ' --resolution ' + str(resolution)

    if ds9:
        cmd += ' --ds9'

    if update_frames:
        cmd += ' --output-tod'

    cmd += ' --output-map'

    print cmd
    os.system(cmd)

    if update_frames:
        for file, frame in zip(files, frames):
            tod = simpleFitsReader(file+'_tod.fits')
            print file, ': shape=', tod['PrimaryImage'].data.dimensions
            frame.signal = tod['PrimaryImage'].data

    m = simpleFitsReader(files[0]+'_map.fits')
    image = herschel.ia.dataset.image.SimpleImage()
    image.setImage(m['PrimaryImage'].data, herschel.share.unit.Unit.parse(m.meta['BUNIT'].string))
    image.coverage = m['COVERAGE'].data

    return image
