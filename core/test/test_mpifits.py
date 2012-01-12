from __future__ import division

import os
import numpy as np
from glob import glob
from tamasis import Map, MPI
from tamasis.numpyutils import assert_all_eq
from uuid import uuid1

rank = MPI.COMM_WORLD.Get_rank()
size = MPI.COMM_WORLD.Get_size()
path_data = os.path.dirname(__file__) + '/data/'
id = MPI.COMM_WORLD.bcast(str(uuid1()))

def test():
    filename = path_data + 'frames_blue_map_naive.fits'
    ref = Map(filename, comm=MPI.COMM_SELF)
    local = Map(filename, comm=MPI.COMM_WORLD)
    local2 = ref.tolocal()

    assert local.shape[0] == int(np.ceil(ref.shape[0] / size))
    assert_all_eq(local, local2)

    ref2 = local.toglobal()
    assert_all_eq(ref, ref2)

    filename2 = 'test-' + id + '.fits'
    local.save(filename2)
    ref3 = Map(filename2, comm=MPI.COMM_SELF)
    local3 = Map(filename2, comm=MPI.COMM_WORLD)

    assert_all_eq(local, local3)
    assert_all_eq(ref, ref3)

def teardown():
    if rank > 0:
        return
    files = glob('*' + id + '.fits')
    for f in files:
        os.remove(f)
