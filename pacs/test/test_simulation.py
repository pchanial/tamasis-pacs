import numpy as np
import os
import tamasis

from tamasis import *
from uuid import uuid1

class TestFailure(Exception):
    pass

tamasis.var.verbose = False

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
    if not np.allclose(simul.status[field], status2[field]): raise TestFailure('Status problem with: '+field)

if not np.allclose(tod, tod2): raise TestFailure()
fields = [x for x in simul.slice.dtype.names if x not in ('filename','unit')]
for field in fields:
    if getattr(simul.slice[0], field) != getattr(simul2.slice[0], field):
        str = "Field '" + field + "'"
        if field == 'scan_step':
            print(str + ' not implemented.')
        else:
            raise TestFailure(str)

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
for c in (1,2,4,8):
    pointing.time = time * c
    for b in ('blue', 'green', 'red'):
        simul = PacsSimulation(pointing, b)
        if simul.slice[0].mode != result[(c,b)]: raise TestFailure()
        if simul.slice[0].compression_factor != c: raise TestFailure()

ok = [ (1,'blue','transparent'), (1,'green','transparent'), (1,'red','transparent'),
       (4,'blue','prime'), (4,'green','prime'), (4,'red','prime'), (4,'red','parallel'),
       (8,'blue','parallel'), (8,'green','parallel'),
     ]
for c in (1,2,4,8):
    pointing.time = time * c
    for b in ('blue', 'green', 'red'):
        for m in ('prime', 'parallel', 'transparent'):
            if (c,b,m) in ok:
                simul = PacsSimulation(pointing, b, mode=m)
                if simul.slice[0].compression_factor != c: raise TestFailure()
                if simul.slice[0].mode != m: raise TestFailure()
            else:
                try:
                    simul = PacsSimulation(pointing, b, mode=m)
                except ValueError:
                    pass
                else:
                    print c,b,m
                    raise TestFailure()

for c in (1,2,4,8):
    pointing.time = time * c
    for b in ('blue', 'green', 'red'):
        simul = PacsSimulation(pointing, b, mode='calibration')
        if simul.slice[0].compression_factor != c: raise TestFailure()

# simulation #>1
ra = (23, 24)
dec = (50, 51)
cam_angle=(10, 11)
scan_angle = (0, -90)
scan_length = (10, 20)
scan_nlegs = (2, 3)
scan_step = (147, 149)
scan_speed = (20, 60)
compression_factor = (4, 4)

pointings = [pacs_create_scan(r,d,ca,sa,sl,sn,sst,ssp,cf) for \
             r,d,ca,sa,sl,sn,sst,ssp,cf in zip(ra,dec,cam_angle,scan_angle,
             scan_length,scan_nlegs,scan_step,scan_speed,compression_factor)]
simul = PacsSimulation(pointings, 'blue')
s = simul.slice

if any_neq(s.ra, ra) or \
   any_neq(s.dec, dec) or \
   any_neq(s.cam_angle, cam_angle) or \
   any_neq(s.scan_angle, scan_angle) or \
   any_neq(s.scan_length, scan_length) or \
   any_neq(s.scan_nlegs, scan_nlegs) or \
   any_neq(s.scan_step, scan_step) or \
   any_neq(s.scan_speed, scan_speed) or \
   any_neq(s.compression_factor, compression_factor): raise TestFailure()
