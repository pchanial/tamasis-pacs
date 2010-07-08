import os
import herschel
import time

from org.python.core import PyFloat, PyInteger, PyList, PyString, PyTuple

from herschel.ia.dataset.image import SimpleImage
from herschel.ia.dataset.image.wcs import Wcs
from herschel.ia.task.all import *
from herschel.ia.gui.kernel import ParameterValidator, ParameterValidationException
from herschel.pacs.signal import Frames

__all__ = ['tamasis_dir', 'tamasisPreprocessor', 'tamasisPhotProject']
tamasis_dir = '/home/pchanial/work/tamasis/tamasis-1.0.4/'
False = 0

fa = herschel.ia.io.fits.FitsArchive()
simpleFitsReader = herschel.ia.toolbox.util.SimpleFitsReaderTask()


def _get_filename():
    return '/tmp/' + os.getenv('USER') + '-' + str(long(time.time()*1000))


def _get_cmd_preprocessor(frames, flatfielding=False, subtraction_mean=False, median_filtering=None, deglitching=None, nsigma=5., npixels_per_sample=None):

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


def _import_tod(file, frame):
    tod = simpleFitsReader(file+'_tod.fits')
    signal = frame.get('Signal')
    d1 = tod['PrimaryImage'].data.dimensions
    d2 = signal.getData().dimensions
    if not d1 == d2:
        raise ValueError('The dimension of the cube '+str(tuple(d1))+' imported from Tamasis is incompatible with the input one '+str(tuple(d2))+" in '"+file+"'.")
    signal.setData(tod['PrimaryImage'].data)
    signal.setUnit(herschel.share.unit.Unit.parse(tod.meta['BUNIT'].string))
    frame.addMaskType('Tamasis', 'Mask imported from TAMASIS')
    mask = herschel.ia.numeric.Bool3d(tod['MASK'].data)
    frame.setMask('Tamasis', mask)


def _import_map(file):
    m = simpleFitsReader(file)
    image = SimpleImage()
    image.setImage(m['PrimaryImage'].data, herschel.share.unit.Unit.parse(m.meta['BUNIT'].string))
    image.coverage = m['COVERAGE'].data
    image.wcs = Wcs(m.meta)
    return image


class _InputFramesValidator(ParameterValidator):
    def validate(self, value):
        if isinstance(value, Frames):
            return
        if not isinstance(value, PyList) and not isinstance(value, PyTuple) or len(value) == 0:
            raise ParameterValidationException()
        for f in value:
            if not isinstance(f, Frames):
                raise ParameterValidationException()
    

class TamasisPreprocessor(JTask):

    def __init__(self, name='Tamasis Preprocessor'):

        print 'tamasispreprocessor __init__'

        p = TaskParameter("frames", valueType=java.lang.Object, mandatory=1, parameterValidator=_InputFramesValidator())
        self.addTaskParameter(p)

        p = TaskParameter('flatfielding', valueType=PyInteger, defaultValue=False)
        self.addTaskParameter(p)

        p = TaskParameter('subtractionMean', valueType=PyInteger, defaultValue=False)
        self.addTaskParameter(p)

        p = TaskParameter('medianFiltering', valueType=PyInteger)
        self.addTaskParameter(p)

        p = TaskParameter('deglitching', valueType=PyString)
        self.addTaskParameter(p)

        p = TaskParameter('nsigma', valueType=PyFloat, default=5.)
        self.addTaskParameter(p)

        p = TaskParameter('npixelsPerSample', valueType=PyInteger, default=6)
        self.addTaskParameter(p)

    def execute(self):

        frames = self.frames
        if isinstance(frames, Frames):
            frames = (frames,)

        cmd, files = _get_cmd_preprocessor(frames, flatfielding=self.flatfielding, subtraction_mean=self.subtractionMean, median_filtering=self.medianFiltering, deglitching=self.deglitching, nsigma=self.nsigma, npixels_per_sample=self.npixelsPerSample)

        cmd += ' --output-tod'
                
        print
        print 'Running command:'
        print cmd
        print
        os.system(cmd)

        map(_import_tod, files, frames)
    


class TamasisPhotProject(TamasisPreprocessor):

    def __init__(self, name='Tamasis PhotProject'):

        TamasisPreprocessor__init__(self, name=name)

        p = TaskParameter('result', valueType=SimpleImage, type=OUT)
        self.addTaskParameter(p)

        p = TaskParameter('updateFrames', valueType=PyInteger, default=False)
        self.addTaskParameter(p)

        p = TaskParameter('header', valueType=PyString)
        self.addTaskParameter(p)

        p = TaskParameter('resolution', valueType=PyFloat)
        self.addTaskParameter(p)

        p = TaskParameter('ds9', valueType=PyInteger, default=0)
        self.addTaskParameter(p)


    def execute(self):

        frames = self.frames
        if isinstance(frames, Frames):
            frames = (frames,)

        cmd, files = _get_cmd_preprocessor(frames, flatfielding=self.flatfielding, subtraction_mean=self.subtractionMean, median_filtering=self.medianFiltering, deglitching=self.deglitching, nsigma=self.nsigma, npixels_per_sample=self.npixelsPerSample)

        if self.header is not None:
            cmd += ' --header ' + self.header

        if self.resolution is not None:
            cmd += ' --resolution ' + str(self.resolution)

        if self.ds9:
            cmd += ' --ds9'

        if self.updateFrames:
            cmd += ' --output-tod'

        cmd += ' --output-map'

        print
        print 'Running command:'
        print cmd
        print
        os.system(cmd)

        if self.update_frames:
            map(_import_tod, files, frames)

        self.result = _import_map(files[0] + '_map.fits')

tamasisPreprocessor = TamasisPreprocessor()
tamasisPhotProject = TamasisPhotProject()
