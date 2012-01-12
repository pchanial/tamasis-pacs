"""
MPI-wrapper module for non-MPI enabled platforms.
"""

try:

    from mpi4py.MPI import *

except:

    _g = globals()
    constants = 'SUM', 'MIN', 'MAX', 'BYTE', 'FLOAT', 'DOUBLE', 'IN_PLACE'
    for i, c in enumerate(constants):
        _g[c] = i

    class FakeComm(object):
        def __init__(self, rank, size):
            self.rank = rank
            self.size = size
        def Get_rank(self):
            return self.rank
        def Get_size(self):
            return self.size
        def allgather(self,sendobj=None, recvobj=None):
            return [sendobj]
        def allreduce(self, sendobj=None, recvobj=None, op=SUM):
            return sendobj
        def AllGatherv(self, i, o, op=None):
            if i == IN_PLACE:
                return
            o[0][...] = i[0]
        def AllReduce(self, i, o, op=None):
            if i == IN_PLACE:
                return
            o[0][...] = i[0]
        def Barrier(self):
            pass
        @staticmethod
        def f2py(fcomm):
             return COMM_SELF
        def py2f(self):
            return 0

    def Get_processor_name():
        return ''

    COMM_NULL = FakeComm(0, 0)
    COMM_SELF = FakeComm(0, 1)
    COMM_WORLD = COMM_SELF

# mpiutils
# datatypes.py
