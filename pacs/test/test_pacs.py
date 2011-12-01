import numpy as np
import pyfits
import os
import tamasis

from glob import glob
from tamasis import PacsObservation, PacsSimulation, Pointing, Projection, Map, Tod, IdentityOperator, MaskOperator, mapper_naive
from tamasis.numpyutils import all_eq, minmax
from uuid import uuid1

tamasis.var.verbose = True
data_dir = os.path.dirname(__file__) + '/data/'
uuid = str(uuid1())

def test():
    # all observation
    obs = PacsObservation(data_dir+'frames_blue.fits', reject_bad_line=False)
    obs.pointing.chop[:] = 0

    # get mask
    proj = Projection(obs, npixels_per_sample=6, oversampling=False)
    o = Tod.ones(proj.shapeout)
    nocoverage = mapper_naive(o, proj).coverage == 0
    assert all_eq(nocoverage, proj.get_mask())
    tod = obs.get_tod()

    # packed projection
    print 'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX'
    print 'XXX FIX ME                                              XXX'
    print 'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX'
    print

    #proj2 = Projection(obs, npixels_per_sample=6, packed=True, oversampling=False)
    #if any_neq(proj.T(tod), proj2.T(tod), 1.e-12): raise TestFailure()

    filename = 'obs-' + uuid + '.fits'
    obs.save(filename, tod)
    obs2 = PacsObservation(filename, reject_bad_line=False)
    tod2 = obs2.get_tod(raw=True)
    assert all_eq(obs.status, obs2.status)
    assert all_eq(tod, tod2)

    telescope  = IdentityOperator()
    projection = Projection(obs, resolution=3.2, oversampling=False,
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

    header = projection.header
    header2 = header.copy()
    header2['NAXIS1'] += 500
    header2['CRPIX1'] += 250
    projection2 = Projection(obs, header=header2, oversampling=False)
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
    projection = Projection(obs, header=map_naive_ref.header,
                            oversampling=False, npixels_per_sample=6)
    tod = obs.get_tod(flatfielding=False)
    masking = MaskOperator(tod.mask)
    model = masking * projection
    map_naive = mapper_naive(tod, model)
    assert all_eq(map_naive, map_naive_ref, 1.e-11)

    obs_rem = PacsObservation(data_dir + 'frames_blue.fits',
                              policy_bad_detector='remove',
                              reject_bad_line=False)
    obs_rem.pointing.chop[:] = 0
    projection_rem = Projection(obs_rem, header=map_naive.header,
                                oversampling=False, npixels_per_sample=7)
    tod_rem = obs_rem.get_tod(flatfielding=False)
    masking_rem = MaskOperator(tod_rem.mask)
    model_rem = masking_rem * projection_rem
    map_naive_rem = mapper_naive(tod_rem, model_rem)
    assert all_eq(map_naive, map_naive_rem, 1.e-11)

def test_pack():
    for channel, nrows, ncolumns in ('red',16,32), ('blue',32,64):
        obs = PacsSimulation(Pointing(0., 0., 0., 0.), channel)
        for dt in (np.uint8, np.uint16, np.uint32, np.uint64):
            a = np.arange(nrows*ncolumns*3, dtype=dt) \
                  .reshape((nrows,ncolumns,-1))
            p = obs.pack(a)
            assert all_eq(a[1,0:16,:], p[16:32,:])
            u = obs.unpack(p)
            assert all_eq(a, u)

def test_npixels_per_sample_is_zero():
    obs = PacsObservation(data_dir + 'frames_blue.fits')
    proj = Projection(obs, npixels_per_sample=2)
    header = proj.header.copy()
    header['crval1'] += 1
    proj2 = Projection(obs, header=header)
    assert proj2.npixels_per_sample == 0
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

    proj1 = Projection(obs1, header=header, oversampling=False)
    proj2 = Projection(obs2, header=header, oversampling=False)
    assert all_eq(proj1.get_mask(), proj2.get_mask())
    assert all_eq(proj1.pmatrix,
                  np.concatenate([p.pmatrix for p in proj2.operands], axis=1))

    m1 = proj1.T(tod1)
    m2 = proj2.T(tod2)
    assert all_eq(m1, m2, 1e-12)
    assert all_eq(proj1(m1), proj2(m1))
    
def tearDown():
    files = glob('*' + uuid + '.fits')
    for f in files:
        try:
            os.remove(f)
        except IOError:
            pass
