import numpy as np
import pyfits
import os
import tamasis

from glob import glob
from tamasis import PacsObservation, PacsSimulation, Pointing, CompressionAverageOperator, ProjectionOperator, Map, Tod, IdentityOperator, MaskOperator, mapper_naive
from tamasis.numpyutils import all_eq, minmax
from uuid import uuid1

tamasis.var.verbose = True
data_dir = os.path.dirname(__file__) + '/data/'
uuid = str(uuid1())

def test():
    # all observation
    obs = PacsObservation(data_dir+'frames_blue.fits', reject_bad_line=False)
    obs.pointing.chop[:] = 0
    header = obs.get_map_header()

    # get mask
    proj = ProjectionOperator(obs, npixels_per_sample=6, downsampling=True,
                              header=header)
    o = Tod.ones(proj.shapeout)
    nocoverage = mapper_naive(o, proj).coverage == 0
    assert all_eq(nocoverage, proj.get_mask())
    tod = obs.get_tod()

    # packed projection
    proj2 = ProjectionOperator(obs, npixels_per_sample=6, packed=True,
                               downsampling=True, header=header)
    m1 = proj.matrix
    m2 = proj2.operands[0].matrix
    assert all_eq(m1.value, m2.value)
    assert all_eq(proj.T(tod), proj2.T(tod), 1.e-11)

    filename = 'obs-' + uuid + '.fits'
    obs.save(filename, tod)
    obs2 = PacsObservation(filename, reject_bad_line=False)
    obs2.pointing.chop[:] = 0
    tod2 = obs2.get_tod(raw=True)
    assert all_eq(obs.status, obs2.status)
    assert all_eq(tod, tod2)

    telescope  = IdentityOperator()
    projection = ProjectionOperator(obs, header=header, downsampling=True,
                                    npixels_per_sample=6)
    crosstalk  = IdentityOperator()
    masking    = MaskOperator(tod.mask)

    model = masking * crosstalk * projection * telescope
    print model

    # naive map
    tod_sb = tod.tounit('Jy/arcsec^2')
    backmap = model.T(tod_sb.magnitude)
    unity = Tod.ones(tod_sb.shape)
    weights = model.T(unity)
    map_naive = Map(backmap / weights, unit='Jy/arcsec^2')

    header2 = header.copy()
    header2['NAXIS1'] += 500
    header2['CRPIX1'] += 250
    projection2 = ProjectionOperator(obs, header=header2, downsampling=True)
    map_naive2 = mapper_naive(tod, MaskOperator(tod.mask) * projection2)
    map_naive2.inunit('Jy/arcsec^2')
    map_naive3 = map_naive2[:,250:header['NAXIS1']+250]
    assert all_eq(map_naive.tounit('Jy/arcsec^2'), map_naive3, 2.e-7)

    # test compatibility with photproject
    tod = obs.get_tod(flatfielding=False, subtraction_mean=False)
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
    projection = ProjectionOperator(obs, header=map_naive_ref.header,
                                    downsampling=True, npixels_per_sample=6)
    tod = obs.get_tod(flatfielding=False, subtraction_mean=False)
    masking = MaskOperator(tod.mask)
    model = masking * projection
    map_naive = mapper_naive(tod, model)
    assert all_eq(map_naive, map_naive_ref, 1.e-11)

    obs_rem = PacsObservation(data_dir + 'frames_blue.fits',
                              policy_bad_detector='remove',
                              reject_bad_line=False)
    obs_rem.pointing.chop[:] = 0
    projection_rem = ProjectionOperator(obs_rem, header=map_naive.header,
                                        downsampling=True, npixels_per_sample=7)
    tod_rem = obs_rem.get_tod(flatfielding=False, subtraction_mean=False)
    masking_rem = MaskOperator(tod_rem.mask)
    model_rem = masking_rem * projection_rem
    map_naive_rem = mapper_naive(tod_rem, model_rem)
    assert all_eq(map_naive, map_naive_rem, 1.e-11)

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
    proj2 = ProjectionOperator(obs, header=header)
    assert proj2.matrix.shape[-1] == 0
    t = proj2(np.ones((header['NAXIS2'],header['NAXIS1'])))
    assert all_eq(minmax(t), [0,0])
    t[:] = 1
    assert all_eq(minmax(proj2.T(t)), [0,0])

def test_slice1():
    obs = PacsObservation(data_dir+'frames_blue.fits[11:20]')
    tod = obs.get_tod(flatfielding=False)

    filename = 'obs-' + uuid + '.fits'
    obs.save(filename, tod)
    obs2 = PacsObservation(filename)
    tod2 = obs2.get_tod(raw=True)
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

    tod1 = obs1.get_tod(subtraction_mean=False)
    tod2 = obs2.get_tod(subtraction_mean=False)
    assert all_eq(tod1, tod2)

    header = obs1.get_map_header()

    proj1 = ProjectionOperator(obs1, header=header)
    proj2 = ProjectionOperator(obs2, header=header)
    assert all_eq(proj1.get_mask(), proj2.get_mask())
    assert all_eq(proj1.matrix, np.concatenate([p.matrix for p in \
                  proj2.operands], axis=1))

    model1 = CompressionAverageOperator(obs1.slice.compression_factor) * proj1
    model2 = CompressionAverageOperator(obs2.slice.compression_factor) * proj2
    
    m1 = model1.T(tod1)
    m2 = model2.T(tod2)
    assert all_eq(m1, m2, 1e-10)
    assert all_eq(model1(m1), model2(m1))

def test_pTx_pT1():
    obs1 = PacsObservation(data_dir + 'frames_blue.fits')
    obs2 = PacsObservation([data_dir + 'frames_blue.fits[1:41]',
                            data_dir + 'frames_blue.fits[42:43]',
                            data_dir + 'frames_blue.fits[44:360]'])
    obs1.pointing.chop = 0
    obs2.pointing.chop = 0

    tod = obs1.get_tod(subtraction_mean=False)

    model1 = ProjectionOperator(obs1, downsampling=True)
    m1 = mapper_naive(tod, model1, unit='Jy/arcsec^2')

    model1.apply_mask(tod.mask)
    tod.inunit('Jy/arcsec^2')
    b,w = model1.get_pTx_pT1(tod)
    m2 = (b / w)
    assert all_eq(m1, m2)

    model3 = ProjectionOperator(obs2, downsampling=True)
    model3.apply_mask(tod.mask)
    
    b,w = model3.get_pTx_pT1(tod)
    m3 = (b / w)
    assert all_eq(m1, m3, 1e-10)

    
def teardown():
    files = glob('*' + uuid + '.fits')
    for f in files:
        try:
            os.remove(f)
        except IOError:
            pass
