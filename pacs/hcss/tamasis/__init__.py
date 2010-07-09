import herschel
import java
import os
import time

from java.lang import Double, Integer, String
from org.python.core import PyList, PyTuple
from herschel.ia.dataset.image import SimpleImage
from herschel.ia.dataset.image.wcs import Wcs
from herschel.ia.task.all import *
from herschel.ia.task.views import TaskToolRegistry
from herschel.ia.gui.kernel import ParameterValidator, ParameterValidationException
from herschel.ia.gui.kernel.Tool import Category
from herschel.pacs.signal import Frames

__all__ = ['tamasis_dir', 'tamasisPreprocessor', 'tamasisPhotProject']
tamasis_dir = '/home/pchanial/work/tamasis/tamasis-latest/'

ld_library_path = '/mnt/local/opt/core-3.1-amd64/python/2.5/lib:/opt/core-3.1-amd64/pgplot/5.2/:/opt/core-3.1-amd64/gcc/4.5/lib64:/opt/core-3.1-amd64/libelf/0.8.13-gcc45/lib:/opt/core-3.1-amd64/cloog-ppl/0.15.7-gcc45/lib:/opt/core-3.1-amd64/polylib/5.22.4-gcc45/lib:/opt/core-3.1-amd64/ppl/0.10.2-gcc45/lib:/opt/core-3.1-amd64/mpc/0.8.1-gcc45/lib:/opt/core-3.1-amd64/mpfr/2.4.2-gcc45/lib:/opt/core-3.1-amd64/gmp/4.3.2-gcc45/lib:/opt/core-3.1-amd64/fftw3/3.2.2/lib:/opt/core-3.1-amd64/python/2.6/lib:/opt/core-3.1-amd64/ifc/11.1/lib/intel64:/opt/core-3.1-amd64/ifc/11.1/mkl/lib/em64t:/opt/core-3.1-amd64/cfitsio/3.23/lib'
python_path = '/home/pchanial/software/lib/python2.6/site-packages/:/home/pchanial/work/tamasis/tamasis-latest/:/opt/core-3.1-amd64/python/2.6'
False = 0

fa = herschel.ia.io.fits.FitsArchive()
simpleFitsReader = herschel.ia.toolbox.util.SimpleFitsReaderTask()


def _get_filename():
    return '/tmp/' + os.getenv('USER') + '-' + str(long(time.time()*1000))


def _get_cmd_preprocessor(frames, options):

    cmd = tamasis_dir + 'pacs/scripts/hcss_interface.py'
    
    # write frames to disk
    filename = _get_filename()+'_tod_'
    files = []

    for ifile in range(len(frames)):
        file = filename + str(ifile) + '.fits'
        fa.save(file, frames[ifile])
        files.append(file)
        cmd += ' ' + file

    if options.flatfielding:
        cmd += ' --flatfielding'

    if options.subtractionMean:
        cmd += ' --subtraction-mean'

    if options.medianFiltering is not None:
        cmd += ' --median-filtering ' + str(options.medianFiltering)

    if options.deglitching is not None:
        cmd += ' --deglitching ' + options.deglitching + ' --nsigma ' + str(options.nsigma)
    
    if options.npixelsPerSample is not None:
        cmd += ' --npixels-per-sample ' + str(options.npixelsPerSample)

    for policy in (options.framePolicyInscan, options.framePolicyTurnaround, options.framePolicyOther, options.framePolicyInvalid):
        if policy.lower() not in ('keep', 'mask'):
            raise ValueError("In HCSS, a policy must be 'keep' or 'mask'.")
        
    cmd += ' --frame-policy-inscan ' + options.framePolicyInscan
    cmd += ' --frame-policy-turnaround ' + options.framePolicyTurnaround
    cmd += ' --frame-policy-other ' + options.framePolicyOther
    cmd += ' --frame-policy-invalid ' + options.framePolicyInvalid

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


def _run_command(cmd):
    try:
        ld_tmp = os.getenv('LD_LIBRARY_PATH')
    except KeyError:
        ld_tmp = None
    os.putenv('LD_LIBRARY_PATH', ld_library_path)
    try:
        py_tmp = os.getenv('PYTHONPATH')
    except:
        py_tmp = None
    os.putenv('PYTHONPATH', python_path)
    print
    print 'Running command:'
    print cmd
    print
    os.system(cmd)
    if ld_tmp is not None:
        os.putenv('LD_LIBRARY_PATH', ld_tmp)
    if py_tmp is not None:
        os.putenv('PYTHONPATH', py_tmp)

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
    """
    The following frames preprocessing step can be done with TAMASIS.
       - V-to-Jy conversion
       - subtraction of mean timeline value
       - flat-fielding
       - second level deglitching using or standard deviation or MAD
       - median filtering
    """
    def __init__(self, name='Tamasis Preprocessor'):
        p = TaskParameter("frames", valueType=java.lang.Object, mandatory=1, parameterValidator=_InputFramesValidator(), description='Input Frames, or array of Frames')
        self.addTaskParameter(p)

        p = TaskParameter('flatfielding', valueType=Integer, defaultValue=False, description='Apply flat-field correction')
        self.addTaskParameter(p)

        p = TaskParameter('subtractionMean', valueType=Integer, defaultValue=False, description="Subtract detector timeline's mean value for each detector")
        self.addTaskParameter(p)

        p = TaskParameter('medianFiltering', valueType=Integer, description='Length of the median filtering window')
        self.addTaskParameter(p)

        p = TaskParameter('deglitching', valueType=String, description="Deglitching method 'l2std' or 'l2mad'")
        self.addTaskParameter(p)

        p = TaskParameter('nsigma', valueType=Double, defaultValue=5., description='N-sigma value for deglitching')
        self.addTaskParameter(p)

        p = TaskParameter('npixelsPerSample', valueType=Integer, defaultValue=6, description='Maximum number of sky pixels seen by a detector.')
        self.addTaskParameter(p)

        p = TaskParameter('framePolicyInscan', valueType=String, defaultValue='keep', description="Policy for in-scan frames: 'keep' or 'mask'")
        self.addTaskParameter(p)

        p = TaskParameter('framePolicyTurnaround', valueType=String, defaultValue='keep', description="Policy for turn-around frames: 'keep' or 'mask'")
        self.addTaskParameter(p)

        p = TaskParameter('framePolicyOther', valueType=String, defaultValue='mask', description="Policy for other frames: 'keep' or 'mask'")
        self.addTaskParameter(p)

        p = TaskParameter('framePolicyInvalid', valueType=String, defaultValue='mask', description="Policy for invalid frames: 'keep' or 'mask'")
        self.addTaskParameter(p)

    def execute(self):

        frames = self.frames
        if isinstance(frames, Frames):
            frames = (frames,)

        cmd, files = _get_cmd_preprocessor(frames, self)
        
        cmd += ' --output-tod'

        _run_command(cmd)

        map(_import_tod, files, frames)
    


class TamasisPhotProject(TamasisPreprocessor):

    def __init__(self, name='Tamasis PhotProject'):

        TamasisPreprocessor.__init__(self, name=name)

        p = TaskParameter('result', valueType=SimpleImage, type=OUT, description='Output map')
        self.addTaskParameter(p)

        p = TaskParameter('updateFrames', valueType=Integer, defaultValue=False, description='Update input Frames with Tamasis preprocessed data')
        self.addTaskParameter(p)

        p = TaskParameter('header', valueType=String, description='FITS filename containing the header to be used for the map creation')
        self.addTaskParameter(p)

        p = TaskParameter('resolution', valueType=Double, description='Map pixel size in arcsec')
        self.addTaskParameter(p)

        p = TaskParameter('ds9', valueType=Integer, defaultValue=False, description='Load map in ds9')
        self.addTaskParameter(p)


    def execute(self):

        frames = self.frames
        if isinstance(frames, Frames):
            frames = (frames,)

        cmd, files = _get_cmd_preprocessor(frames, self)

        if self.header is not None:
            cmd += ' --header ' + self.header

        if self.resolution is not None:
            cmd += ' --resolution ' + str(self.resolution)

        if self.ds9:
            cmd += ' --ds9'

        if self.updateFrames:
            cmd += ' --output-tod'

        cmd += ' --output-map'

        _run_command(cmd)
        
        if self.updateFrames:
            map(_import_tod, files, frames)

        self.result = _import_map(files[0] + '_map.fits')

TaskToolRegistry.getInstance().register(TamasisPreprocessor(), [Category.PACS])
TaskToolRegistry.getInstance().register(TamasisPhotProject(), [Category.PACS])
