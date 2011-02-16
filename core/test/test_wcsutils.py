import numpy as np
from tamasis import *

class TestFailure(Exception):
    pass

# mean_degrees
if any_neq(mean_degrees([1,2]), 1.5): raise TestFailure()
if any_neq(mean_degrees([1,359.1]), 0.05, rtol=1.e-12): raise TestFailure()
if any_neq(mean_degrees([0.1,359.1]), 359.6, rtol=1.e-12): raise TestFailure()


# angle_lonlat
if any_neq(angle_lonlat(30, 0, 40, 0), 10.): raise TestFailure()

input = ( ((30,0), (40,0), 10),
          ((39,0), (92, 90), 90),
          ((39,0), (37,-90), 90),
          ((37,-90), (39,0), 90),
          ((100,90),(134,-90), 180),
          ((24,30),(24,32), 2),
        )
for c1, c2, angle in input:
    if any_neq(angle_lonlat(c1, c2), angle, 1.e-10): raise TestFailure()

# barycenter_lonlat
if any_neq(barycenter_lonlat([30,40], [0, 0]), [35,0]): raise TestFailure()
if any_neq(barycenter_lonlat([20,20,20], [-90,0,90]), [20,0]): raise TestFailure()
if any_neq(barycenter_lonlat([20,20,20], [0,45,90]), [20,45]): raise TestFailure()
