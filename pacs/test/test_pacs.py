import astropy.io.fits as pyfits
import numpy as np
import os
import tamasis

from glob import glob
from pyoperators import MaskOperator
from pyoperators.utils.testing import assert_eq
from pysimulators import Map, Pointing
from tamasis import PacsInstrument, PacsObservation, PacsSimulation, CompressionAverageOperator, mapper_naive
from tamasis.utils import all_eq, minmax
from uuid import uuid1

tamasis.var.verbose = True
data_dir = os.path.dirname(__file__) + '/data/'
uuid = str(uuid1())
tol = 1.e-10
PacsInstrument.info.CALFILE_BADP = tamasis.var.path + '/pacs/PCalPhotometer_Ba'\
                                   'dPixelMask_FM_v5.fits'
PacsInstrument.info.CALFILE_RESP = tamasis.var.path + '/pacs/PCalPhotometer_Re'\
                                   'sponsivity_FM_v5.fits'

def test_save():
    obs = PacsObservation(data_dir+'frames_blue.fits', reject_bad_line=False)
    obs.pointing.ra += 0.1
    obs.pointing.dec -= 0.1
    obs.pointing.pa += 20
    obs.pointing.chop = 0
    tod = obs.get_tod()

    filename = 'obs-' + uuid + '.fits'
    obs.save(filename, tod)

    obs2 = PacsObservation(filename, reject_bad_line=False)
    tod2 = obs2.get_tod()

    assert all_eq(obs.pointing, obs2.pointing)

    obs.status.RaArray = obs.pointing.ra
    obs.status.DecArray = obs.pointing.dec
    obs.status.PaArray = obs.pointing.pa
    obs.status.CHOPFPUANGLE = obs.pointing.chop
    assert all_eq(obs.status, obs2.status)
    assert all_eq(tod, tod2)

def test_header():
    obs = PacsObservation(data_dir+'frames_blue.fits', reject_bad_line=False)
    obs.pointing.chop = 0
    header = obs.get_map_header()
    projection = obs.get_projection_operator(header=header, downsampling=True,
                                             npixels_per_sample=6)
    tod = obs.get_tod()

    map_naive = mapper_naive(tod, projection)

    header2 = header.copy()
    header2['NAXIS1'] += 500
    header2['CRPIX1'] += 250
    projection2 = obs.get_projection_operator(header=header2, downsampling=True)
    map_naive2 = mapper_naive(tod, MaskOperator(tod.mask) * projection2)
    map_naive2.inunit('Jy/arcsec^2')
    map_naive3 = map_naive2[:,250:header['NAXIS1']+250]
    assert all_eq(map_naive.tounit('Jy/arcsec^2'), map_naive3, 2.e-7)

    # test compatibility with photproject
    tod = obs.get_tod()
    map_naive4 = mapper_naive(tod, projection)
    hdu_ref = pyfits.open(data_dir + 'frames_blue_map_hcss_photproject.fits')[1]
    map_ref = Map(hdu_ref.data, hdu_ref.header, unit=hdu_ref.header['qtty____']+'/pixel')
    std_naive = np.std(map_naive4[40:60,40:60])
    std_ref = np.std(map_ref[40:60,40:60])
    assert abs(std_naive-std_ref) / std_ref < 0.025

def test_detector_policy():
    map_naive_ref = Map(data_dir + '../../../core/test/data/frames_blue_map_naive.fits')
    obs = PacsObservation(data_dir + 'frames_blue.fits', reject_bad_line=False)
    obs.pointing.chop[:] = 0
    projection = obs.get_projection_operator(header=map_naive_ref.header,
                                             downsampling=True,
                                             npixels_per_sample=6)
    tod = obs.get_tod()
    masking = MaskOperator(tod.mask)
    model = masking * projection
    map_naive = mapper_naive(tod, model)
    assert all_eq(map_naive, map_naive_ref, tol)

    obs_rem = PacsObservation(data_dir + 'frames_blue.fits',
                              policy_bad_detector='remove',
                              reject_bad_line=False)
    obs_rem.pointing.chop[:] = 0
    projection_rem = obs_rem.get_projection_operator(header=map_naive.header,
                                                     downsampling=True,
                                                     npixels_per_sample=7)
    tod_rem = obs_rem.get_tod()
    masking_rem = MaskOperator(tod_rem.mask)
    model_rem = masking_rem * projection_rem
    map_naive_rem = mapper_naive(tod_rem, model_rem)
    assert all_eq(map_naive, map_naive_rem, tol)

def test_pack():
    for channel, nrows, ncolumns in ('red',16,32), ('blue',32,64):
        obs = PacsSimulation(Pointing((0., 0., 0.), 0.), channel)
        for dt in (np.uint8, np.uint16, np.uint32, np.uint64):
            a = np.arange(nrows*ncolumns*3, dtype=dt) \
                  .reshape((nrows,ncolumns,-1))
            p = obs.pack(a)
            assert all_eq(a[1,0:16,:], p[16:32,:])
            u = obs.unpack(p)
            assert all_eq(a, u)

def test_npixels_per_sample_is_zero():
    obs = PacsObservation(data_dir + 'frames_blue.fits')
    header = obs.get_map_header()
    header['crval1'] += 1
    proj2 = obs.get_projection_operator(header=header)
    assert proj2.matrix.shape[-1] == 0
    t = proj2(np.ones((header['NAXIS2'],header['NAXIS1'])))
    assert all_eq(minmax(t), [0,0])
    t[:] = 1
    assert all_eq(minmax(proj2.T(t)), [0,0])

def test_slice1():
    obs = PacsObservation(data_dir+'frames_blue.fits[11:20]')
    tod = obs.get_tod()

    filename = 'obs-' + uuid + '.fits'
    obs.save(filename, tod)
    obs2 = PacsObservation(filename)
    tod2 = obs2.get_tod()
    assert all_eq(obs.status[10:20], obs2.status)
    assert all_eq(tod, tod2)

def test_slice2():
    obs1 = PacsObservation(data_dir + 'frames_blue.fits')
    obs2 = PacsObservation([data_dir + 'frames_blue.fits[1:41]',
                            data_dir + 'frames_blue.fits[42:43]',
                            data_dir + 'frames_blue.fits[44:360]'])
    assert all_eq(obs1.pointing, obs2.pointing[~obs2.pointing.removed])
    obs1.pointing.chop = 0
    obs2.pointing.chop = 0

    tod1 = obs1.get_tod()
    tod2 = obs2.get_tod()
    assert all_eq(tod1, tod2)

    header = obs1.get_map_header()

    proj1 = obs1.get_projection_operator(header=header)
    proj2 = obs2.get_projection_operator(header=header)
    proj3 = obs2.get_projection_operator(header=header, storage='on fly')
    assert all_eq(proj1.get_mask(), proj2.get_mask())
    assert all_eq(proj1.get_mask(), proj3.get_mask())
    assert all_eq(proj1.matrix, np.concatenate([p.matrix for p in \
                  proj2.operands], axis=1))
    assert all_eq(proj1.matrix, np.concatenate([p.matrix for p in \
                  proj3.operands], axis=1))

    model1 = CompressionAverageOperator(obs1.slice.compression_factor) * proj1
    model2 = CompressionAverageOperator(obs2.slice.compression_factor) * proj2
    model3 = CompressionAverageOperator(obs2.slice.compression_factor) * proj3
    
    m1 = model1.T(tod1)
    m2 = model2.T(tod2)
    m3 = model3.T(tod2)
    assert all_eq(m1, m2, tol)
    assert all_eq(m1, m3, tol)
    assert all_eq(model1(m1), model2(m1))
    assert all_eq(model1(m1), model3(m1))

def test_pTx_pT1():
    obs1 = PacsObservation(data_dir + 'frames_blue.fits')
    obs2 = PacsObservation([data_dir + 'frames_blue.fits[1:41]',
                            data_dir + 'frames_blue.fits[42:43]',
                            data_dir + 'frames_blue.fits[44:360]'])
    obs1.pointing.chop = 0
    obs2.pointing.chop = 0
    header = obs1.get_map_header()

    tod = obs1.get_tod()

    model1 = obs1.get_projection_operator(downsampling=True, header=header)
    ref = mapper_naive(tod, model1, unit='Jy/arcsec^2')

    model1.apply_mask(tod.mask)
    tod.inunit('Jy/arcsec^2')
    b1, w1 = model1.get_pTx_pT1(tod)
    m1 = b1 / w1
    assert all_eq(ref, m1, tol)

    model2 = obs2.get_projection_operator(downsampling=True, header=header)
    model2.apply_mask(tod.mask)
    
    b2, w2 = model2.get_pTx_pT1(tod)
    m2 = b2 / w2
    assert all_eq(ref, m2, tol)

    model3 = obs2.get_projection_operator(downsampling=True, header=header,
                                          storage='on fly')
    MaskOperator(tod.mask)(tod, tod)
    b3, w3 = model3.get_pTx_pT1(tod)
    m3 = b3 / w3
    assert all_eq(ref, m3, tol)

    
def teardown():
    files = glob('*' + uuid + '.fits')
    for f in files:
        try:
            os.remove(f)
        except IOError:
            pass
