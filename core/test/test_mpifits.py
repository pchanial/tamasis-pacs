from __future__ import division

import os
import pyfits
import tamasis
from mpi4py import MPI
from tamasis import *
from uuid import uuid1

class TestFailure(Exception): pass

tamasis.var.verbose=True

rank = MPI.COMM_WORLD.Get_rank()
size = MPI.COMM_WORLD.Get_size()

filename = os.path.dirname(__file__) + '/data/' + 'frames_blue_map_naive.fits'
map_ref = pyfits.open(filename)[0].data

ref = Map(filename, comm=MPI.COMM_SELF)
local = Map(filename, comm=MPI.COMM_WORLD)
local2 = ref.tolocal()

if local.shape[0] != int(np.ceil(ref.shape[0] / size)): raise TestFailure()
if any_neq(local, local2): raise TestFailure()

ref2 = local.toglobal()
if any_neq(ref, ref2): raise TestFailure()

if rank == 0:
    filename2 = 'test-'+str(uuid1())+'.fits'
else:
    filename2 = None
filename2 = MPI.COMM_WORLD.bcast(filename2)

try:
    local.save(filename2)
    ref2 = Map(filename2, comm=MPI.COMM_SELF)
    local2 = Map(filename2, comm=MPI.COMM_WORLD)
finally:
    try:
        os.remove(filename2)
    except:
        pass

if any_neq(ref, ref2):
    raise TestFailure()

if any_neq(local, local2):
    raise TestFailure()

