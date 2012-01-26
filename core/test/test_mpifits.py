from __future__ import division

import os
import numpy as np
import tamasis
from glob import glob
from numpy.testing import assert_equal
from tamasis import Map, MPI
from tamasis.numpyutils import assert_all_eq
from uuid import uuid1

rank = MPI.COMM_WORLD.rank
size = MPI.COMM_WORLD.size
path_data = os.path.dirname(__file__) + '/data/'
id = MPI.COMM_WORLD.bcast(str(uuid1()))

filename = path_data + 'frames_blue_map_naive.fits'
ref = Map(filename, comm=MPI.COMM_SELF)
lmap = Map(filename, comm=MPI.COMM_WORLD)

def test_read():
    shape_global = ref.shape
    shape_local = tamasis.mpiutils.distribute_shape(shape_global)
    assert_equal(shape_local, lmap.shape)

    lmaps = MPI.COMM_WORLD.allgather(lmap)
    gmap = np.concatenate(lmaps, axis=0).magnitude

    assert_all_eq(ref, gmap)

def test_write():
    filename2 = 'test-' + id + '.fits'
    lmap.save(filename2)
    ref2 = Map(filename2, comm=MPI.COMM_SELF)
    lmap2 = Map(filename2, comm=MPI.COMM_WORLD)
    lmap2.save('test'+str(rank)+'fits', comm=MPI.COMM_SELF)
    assert_all_eq(lmap, lmap2)
    assert_all_eq(ref, ref2)

#    local2 = ref.tolocal()

#    assert local.shape[0] == int(np.ceil(ref.shape[0] / size))
#    assert_all_eq(local, local2)
#
#    ref2 = local.toglobal()
#    assert_all_eq(ref, ref2)

def teardown():
    if rank > 0:
        return
    files = glob('*' + id + '.fits')
    for f in files:
        os.remove(f)
