import numpy as np
import os
import tamasis

from tamasis import *
from tamasis.numpyutils import all_eq
from uuid import uuid1

tamasis.var.verbose = False

def test1():
    # creation of the sky map
    msize = 50
    mymap = gaussian(2*(msize*2+1,), 10, unit='Jy/pixel')
    cd = np.array([[-1., 0.],[0., 1.]]) / 3600.
    header = create_fitsheader(fromdata=mymap, crval=[53.,27.], cd=cd)
    mymap.header = header

    # creation of the simulation
    scan = pacs_create_scan(header['CRVAL1'], header['CRVAL2'], cam_angle=0., scan_length=60, scan_nlegs=1, scan_angle=20.)
    simul = PacsSimulation(scan, 'red', policy_bad_detector='keep')

    # build the acquisition model
    model = CompressionAverage(simul.slice.compression_factor) * \
            Projection(simul, header=header, oversampling=True, npixels_per_sample=49)

    # get the noiseless tod
    tod = model(mymap)

    filename = 'simul-'+str(uuid1())+'.fits'

    try:
        simul.save(filename, tod)
        simul2 = PacsObservation(filename, policy_bad_detector='keep')
        status2 = simul2.status
        tod2 = simul2.get_tod(raw=True)
    finally:
        try:
            os.remove(filename)
        except:
            pass

    for field in simul.status.dtype.names:
        if field == 'BAND': continue
        assert all_eq(simul.status[field], status2[field])

    assert all_eq(tod, tod2)
    fields = [x for x in simul.slice.dtype.names if x not in ('filename','unit')]
    for field in fields:
        if getattr(simul.slice[0], field) != getattr(simul2.slice[0], field):
            msg = "Field '" + field + "'"
            if field == 'scan_step':
                print(msg + ' not implemented.')
            else:
                assert False

def test2():
    # test
    data_dir = os.path.dirname(__file__) + '/data/'
    obs = PacsObservation(filename=data_dir+'frames_blue.fits') 
    pointing = obs.pointing.copy()
    pointing.header = obs.pointing.header.copy()

    # check mode
    time = pointing.time / 4

    result = {
        (1,'blue'): 'calibration', (2,'blue'):'calibration', (4,'blue'):'prime', (8,'blue'):'parallel',
        (1,'green'): 'calibration', (2,'green'):'calibration', (4,'green'):'prime', (8,'green'):'parallel',
        (1,'red'): 'calibration', (2,'red'):'calibration', (4,'red'):'prime', (8,'red'):'calibration',
    }

    def func1(c, b):
        pointing.time = time * c
        simul = PacsSimulation(pointing, b)
        assert simul.slice[0].mode == result[(c,b)]
        assert simul.slice[0].compression_factor == c

    for c in (1,2,4,8):
        for b in ('blue', 'green', 'red'):
            yield func1, c, b

    ok = [ (1,'blue','transparent'), (1,'green','transparent'),
           (1,'red','transparent'), (4,'blue','prime'), (4,'green','prime'),
           (4,'red','prime'), (4,'red','parallel'), (8,'blue','parallel'),
           (8,'green','parallel'),
         ]

    def func2(c, b, m):
        pointing.time = time * c
        if (c,b,m) in ok:
            simul = PacsSimulation(pointing, b, mode=m)
            assert simul.slice[0].compression_factor == c
            assert simul.slice[0].mode == m
        else:
            try:
                simul = PacsSimulation(pointing, b, mode=m)
            except ValueError:
                pass
            else:
                print c,b,m
                assert False

    for c in (1,2,4,8):
        for b in ('blue', 'green', 'red'):
            for m in ('prime', 'parallel', 'transparent'):
                yield func2, c, b, m

    def func3(c, b):
        pointing.time = time * c
        simul = PacsSimulation(pointing, b, mode='calibration')
        assert simul.slice[0].compression_factor == c

    for c in (1,2,4,8):
        for b in ('blue', 'green', 'red'):
            yield func3, c, b

def test_multiple_pointings():
    ra = (23, 24)
    dec = (50, 51)
    cam_angle=(10, 11)
    scan_angle = (0, -90)
    scan_length = (10, 20)
    scan_nlegs = (2, 3)
    scan_step = (147, 149)
    scan_speed = (20, 60)
    compression_factor = (4, 4)

    pointings = [pacs_create_scan(r,d,ca,sa,sl,sn,sst,ssp,cf) for r,d,ca,sa,sl,
                 sn,sst,ssp,cf in zip(ra,dec,cam_angle,scan_angle, scan_length,
                 scan_nlegs,scan_step,scan_speed,compression_factor)]
    simul = PacsSimulation(pointings, 'blue')
    s = simul.slice
    assert len(s) == 2
    assert all_eq(s.ra, ra)
    assert all_eq(s.dec, dec)
    assert all_eq(s.cam_angle, cam_angle)
    assert all_eq(s.scan_angle, scan_angle)
    assert all_eq(s.scan_length, scan_length)
    assert all_eq(s.scan_nlegs, scan_nlegs)
    assert all_eq(s.scan_step, scan_step)
    assert all_eq(s.scan_speed, scan_speed)
    assert all_eq(s.compression_factor, compression_factor)
