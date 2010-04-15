from   tamasis import *
import pyfits

observation = MadMap1Observation('tests/madmap1/todSpirePsw_be', 'tests/madmap1/invnttSpirePsw_be', 
                                 'tests/madmap1/madmapSpirePsw.fits[coverage]', 'big_endian', 135)

projection = Projection(observation)
pack = Pack(observation.mapmask)
model = projection * pack

map_naive = naive_mapper(observation.get_tod(), model)
map_ref = pyfits.fitsopen('tests/madmap1/naivemapSpirePsw.fits')['image'].data.T
if any_neq(map_naive,map_ref,15): print 'FAILED: naive_mapper madcap'

print 'OK.'
