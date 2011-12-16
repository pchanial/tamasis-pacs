import herschel
import java
import os
import time

from java.lang import Double, Integer, String, System
from java.util import Date
from org.python.core import PyList, PyTuple
from herschel.ia.dataset.image import SimpleImage
from herschel.ia.dataset.image.wcs import Wcs
from herschel.ia.gui.kernel import ParameterValidator, ParameterValidationException
from herschel.ia.gui.kernel.Tool import Category
from herschel.ia.numeric import Bool2d
from herschel.ia.task.all import *
from herschel.ia.task.views import TaskToolRegistry
from herschel.pacs.cal import *
from herschel.pacs.spg import *
from herschel.pacs.spg.common import *
from herschel.pacs.spg.phot import *


__all__ = [ 'tamasisPrepareFrames', 'tamasisPreprocessor', 'tamasisPhotProject' ]

python_path = '/home/pchanial/software/lib/python2.6/site-packages'
tamasis_path = python_path + '/tamasis/'
ld_library_path = '/opt/core-3.1-amd64/ifc/2011.0.084/composerxe-2011.0.084/composerxe-2011.1.107/compiler/lib/intel64:/opt/core-3.1-amd64/ifc/2011.0.084/composerxe-2011.0.084/composerxe-2011.1.107/mpirt/lib/intel64:/opt/core-3.1-amd64/ifc/2011.0.084/composerxe-2011.0.084/composerxe-2011.1.107/compiler/lib/intel64:/opt/core-3.1-amd64/ifc/2011.0.084/composerxe-2011.0.084/composerxe-2011.1.107/mkl/lib/intel64:/opt/core-3.1-amd64/openmpi/1.4.1-gcc45/lib:/opt/core-3.1-amd64/cfitsio/3.25/lib:/opt/core-3.1-amd64/pgplot/5.2/:/opt/core-3.1-amd64/gcc/4.5/lib64:/opt/core-3.1-amd64/libelf/0.8.13-gcc45/lib:/opt/core-3.1-amd64/cloog-ppl/0.15.7-gcc45/lib:/opt/core-3.1-amd64/polylib/5.22.4-gcc45/lib:/opt/core-3.1-amd64/ppl/0.10.2-gcc45/lib:/opt/core-3.1-amd64/mpc/0.8.1-gcc45/lib:/opt/core-3.1-amd64/mpfr/2.4.2-gcc45/lib:/opt/core-3.1-amd64/gmp/4.3.2-gcc45/lib:/opt/core-3.1-amd64/fftw3/3.2.2/lib:/opt/core-3.1-amd64/python/2.6/lib:/opt/core-3.1-amd64/ifc/12.0/lib/intel64:/opt/core-3.1-amd64/ifc/12.0/mkl/lib/em64t:/opt/core-3.1-amd64/cfitsio/3.23/lib'

False = 0
True  = 1

fa               = herschel.ia.io.fits.FitsArchive()
simpleFitsReader = herschel.ia.toolbox.util.SimpleFitsReaderTask()

class TamasisPrepareFrames(JTask):
    """
    This task converts a level-0 Frames into a level-0.5 Frames which
    can be used by the Tamasis package.
    The following processing steps are operated:
        - findBlocks
        - photFlagBadPixels
        - set BADPIXELS mask
        - photFlagSaturation
        - photConvDigit2Volts
        - convertChopper2Angle
        - photAddInstantPointing
        - cleanPlateauFrames
    """

    def __init__(self, name='Frames Preparation for Tamasis'):

        p = TaskParameter("filename", valueType=java.lang.String, mandatory=1, description='Input prefix for the level-0 set of files (*_ObsCont.fits, *_HPPHK.fits, *_Pointing.fits, *_OrbitEphem.fits, *_TimeCorr.fits)')
        self.addTaskParameter(p)

        p = TaskParameter("path", valueType=java.lang.String, mandatory=0, defaultValue='', description='Path of the output level-0.5 frames')
        self.addTaskParameter(p)

    def execute(self):
        import herschel, os

        obsCont    = fa.load(self.filename + '_ObsCont.fits')
        hk         = fa.load(self.filename + '_HPPHK.fits')['hk']
        pointing   = fa.load(self.filename + '_Pointing.fits')
        orbitEphem = fa.load(self.filename + '_OrbitEphem.fits')
        timeCorr   = fa.load(self.filename + '_TimeCorr.fits')
        
        calTree  = getCalTree("FM")
    
        for band in ('blue', 'green', 'red'):
            
            try:
                frames = fa.load(self.filename + '_' + band + '_level0Frames.fits')
            except java.io.IOException:
                continue

            print "running findBlocks", Date()
            frames = findBlocks(frames)

            print "photFlagBadPixels", Date()
            frames = photFlagBadPixels(frames, calTree=calTree)

            myMask = _mkMyMask(band)
            badPix = frames.getMask('BADPIXELS')
            dims = badPix.dimensions
            for ipix in range(dims[0]):
                for jpix in range(dims[1]):
                    if (myMask[ipix,jpix] == True):
                        badPix[ipix,jpix,:] = True
            frames.setMask('BADPIXELS',badPix)

            # This checks for ADC & CL saturation
            print "photFlagSaturation", Date()
            frames = photFlagSaturation(frames, calTree=calTree, hkdata=hk, check='full')

            print "photConvDigit2Volts", Date()
            frames = photConvDigit2Volts(frames, calTree=calTree)

            print "Convert chopper position to angle", Date()
            frames = convertChopper2Angle(frames,calTree=calTree)
            System.gc()

            print "Add pointing", Date()
            frames = photAddInstantPointing(frames, pointing, calTree=calTree, orbitEphem = orbitEphem)
            System.gc()
            # we shall never use the pointing and ephemris data again

            print "cleanPlateauFrames ; " , Date()
            frames = cleanPlateauFrames(frames, calTree=calTree)
            System.gc()

            # now save the file
            print 'Write prepared frames to disk', Date()
            prepared_filename = os.path.join(self.path, os.path.basename(self.filename))
            fa.save(prepared_filename + '_' + band + '_PreparedFrames.fits', frames)

            print


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

        p = TaskParameter('policyInscan', valueType=String, defaultValue='keep', description="Policy for in-scan frames: 'keep' or 'mask'")
        self.addTaskParameter(p)

        p = TaskParameter('policyTurnaround', valueType=String, defaultValue='keep', description="Policy for turn-around frames: 'keep' or 'mask'")
        self.addTaskParameter(p)

        p = TaskParameter('policyOther', valueType=String, defaultValue='mask', description="Policy for other frames: 'keep' or 'mask'")
        self.addTaskParameter(p)

        p = TaskParameter('policyInvalid', valueType=String, defaultValue='mask', description="Policy for invalid frames: 'keep' or 'mask'")
        self.addTaskParameter(p)

    def execute(self):

        frames = self.frames
        if isinstance(frames, Frames):
            frames = (frames,)

        cmd, files = _get_cmd_preprocessor(frames, self)
        
        cmd += ' --output-tod'

        try:
            _run_command(cmd)
            for file, frame in zip(files, frames):
                _import_tod(file, frame)
        except:
            pass

        _cleanup(files)
    

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

        try:
            _run_command(cmd)
            if self.updateFrames:
                for file, frame in zip(files, frames):
                    _import_tod(file, frame)
            self.result = _import_map(files[0] + '_map.fits')
        except:
            pass

        _cleanup(files)


def _mkMyMask(band):
    """     
    This function creates a dead pixel map to be multiplied to the data
    Directly imported from Nicolas
    """
    # HISTORY
    # 05-Feb-09 MS I remove the pixels that are already in the calibration
    #              mask
    # 01-Jun-10 MS Again I remove 2 pixels that made it to the official mask
    
    if band != 'red':
        bMask = Bool2d(32,64)
        bMask[:,:] = False
        bMask[3,47] = True 
        # slow pixels
        bMask[9,0] = True
        bMask[4,20] = True
        bMask[4,24] = True
        bMask[13,51] = True
        return bMask
        
    rMask = Bool2d(16,32)
    rMask[:,:] = False
    rMask[13,15] = True
    rMask[14,15] = True
    # slow pixels
    rMask[14,3] = True
    rMask[3,6] = True
    rMask[4,5] = True
    return rmask

def _get_filename():
    return '/tmp/' + os.getenv('USER') + '-' + str(long(time.time()*1000))

def _get_cmd_preprocessor(frames, options):

    cmd = 'python ' + tamasis_path + 'pacs_photproject_deprecated.py'
    
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

    for policy in (options.policyInscan, options.policyTurnaround, options.policyOther, options.policyInvalid):
        if policy.lower() not in ('keep', 'mask'):
            raise ValueError("In HCSS, a policy must be 'keep' or 'mask'.")
        
    cmd += ' --policy-inscan ' + options.policyInscan
    cmd += ' --policy-turnaround ' + options.policyTurnaround
    cmd += ' --policy-other ' + options.policyOther
    cmd += ' --policy-invalid ' + options.policyInvalid

    return cmd, files

def _import_tod(file, frame):
    tod = simpleFitsReader(file+'_tod.fits')
    signal = frame.get('Signal')
    d1 = tod['PrimaryImage'].data.dimensions
    d2 = signal.getData().dimensions
    if d1 != d2:
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

def _cleanup(files):
    for file in files:
        try:
            os.remove(file)
            os.remove(file + '_tod.fits')
        except:
            pass
    try:
        os.remove(files[0] + '_map.fits')
    except:
        pass
        
class _InputFramesValidator(ParameterValidator):
    def validate(self, value):
        if isinstance(value, Frames):
            return
        if not isinstance(value, PyList) and not isinstance(value, PyTuple) or len(value) == 0:
            raise ParameterValidationException()
        for f in value:
            if not isinstance(f, Frames):
                raise ParameterValidationException()


TaskToolRegistry.getInstance().register(TamasisPrepareFrames(), [Category.PACS])
TaskToolRegistry.getInstance().register(TamasisPreprocessor(), [Category.PACS])
TaskToolRegistry.getInstance().register(TamasisPhotProject(), [Category.PACS])
