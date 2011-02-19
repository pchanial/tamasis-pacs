import numpy as np
from kapteyn import wcs
from mpi4py import MPI
from tamasis import any_neq, mpiutils as mu, var, wcsutils as wu

rank = MPI.COMM_WORLD.Get_rank()
size = MPI.COMM_WORLD.Get_size()

class TestFailure(Exception):
    def __init__(self, *args):
        print('Test failure: rank ' + str(rank) + '\n' + ',\n'.join([repr(a) for a in args]))
        raise Exception()

if var.mpi_comm.Get_rank() != MPI.COMM_WORLD.Get_rank(): raise TestFailure()
