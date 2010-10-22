import numpy
import os
from tamasis import *
from uuid import uuid1

class TestFailure(Exception):
    pass

# creation of the sky map
msize = 50
mymap = gaussian(2*(msize*2+1,), 10)
cd = numpy.array([[-1., 0.],[0., 1.]]) / 3600.
header = create_fitsheader(mymap, crval=[53.,27.], cd=cd)
mymap.header = header

# creation of the simulation
simul = PacsSimulation('red', (header['CRVAL1'],header['CRVAL2']), scan_length=60, scan_nlegs=1, scan_angle=20., detector_mask=None)

# build the acquisition model
model = CompressionAverage(simul.slice.compression_factor) * \
        Projection(simul, header=header, oversampling=True, npixels_per_sample=49)

# get the noiseless tod
tod = model(mymap)

filename = 'simul-'+str(uuid1())+'.fits'

try:
    simul.save(filename, tod)
    simul2 = PacsObservation(filename, detector_mask=None)
    status2 = simul2.status
    tod2 = simul2.get_tod(raw_data=True)
finally:
    try:
        os.remove(filename)
    except:
        pass

for field in simul.status.dtype.names:
    if field == 'BAND': continue
    if not numpy.allclose(simul.status[field], status2[field]): raise TestFailure('Status problem with: '+field)

if not numpy.allclose(tod, tod2): raise TestFailure()
fields='nsamples,nfinesamples,mode,compression_factor,unit,ra,dec,cam_angle,scan_angle,scan_length,scan_nlegs,scan_step,scan_speed,ninscans,nturnarounds,nothers,ninvalids'.split(',')
ok = True
for field in fields:
    if getattr(simul.slice[0], field) != getattr(simul2.slice[0], field):
        ok = False
        print("Field '" + field + "' is wrong.")
if not ok: raise TestFailure()

