import numpy as np
import tamasis.wcsutils as wu
from tamasis import *
from kapteyn import wcs

class TestFailure(Exception):
    pass

# mean_degrees
if any_neq(wu.mean_degrees([1,2]), 1.5): raise TestFailure()
if any_neq(wu.mean_degrees([1,359.1]), 0.05, rtol=1.e-12): raise TestFailure()
if any_neq(wu.mean_degrees([0.1,359.1]), 359.6, rtol=1.e-12): raise TestFailure()


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
if any_neq(wu.barycenter_lonlat([30,40], [0, 0]), [35,0]): raise TestFailure()
if any_neq(wu.barycenter_lonlat([20,20,20], [-90,0,90]), [20,0]): raise TestFailure()
if any_neq(wu.barycenter_lonlat([20,20,20], [0,45,90]), [20,45]): raise TestFailure()

# get_cdelt_crota2
header = create_fitsheader(np.ones((10,10)), cdelt=(-1.2,3))
cdelt, rot = wu.get_cdelt_crota2(header)
if any_neq(cdelt, (-1.2,3)): raise TestFailure()
if any_neq(rot, 0): raise TestFailure()

cdelt  = (-1.5, 3)
crota2 = 25.
header = create_fitsheader(np.ones((10,10)), cdelt=cdelt, crota2=crota2)
cdelt_, crota2_ = wu.get_cdelt_crota2(header)
if any_neq(cdelt, cdelt_): raise TestFailure()
if any_neq(crota2, crota2_): raise TestFailure()

# combine_fitsheader
headers = [
    wu.create_fitsheader(np.ones((1,1)), cdelt=3., crval=(0,0)),
    wu.create_fitsheader(np.ones((3,3)), cdelt=1., crval=(1,1)),
    wu.create_fitsheader(np.ones((5,5)), cdelt=1., crval=(3,3)),
    wu.create_fitsheader(np.ones((2,2)), cdelt=1., crval=(5,5)),
]
header0 = wu.combine_fitsheader(headers)
proj0 = wcs.Projection(header0)

for iheader, header in enumerate(headers):
    nx = header['NAXIS1']
    ny = header['NAXIS2']
    x = header['CRPIX1']
    y = header['CRPIX2']
    edges = (np.array(3*(0.5,x,nx+0.5)), np.array((0.5,0.5,0.5,y,y,y,ny+0.5,ny+0.5,ny+0.5)))

    a,d = wcs.Projection(header).toworld(edges)
    x0,y0 = proj0.topixel((a,d))

    epsilon = 1.e-10
    if np.any(x0 < 0.5-epsilon) or np.any(y0 < 0.5-epsilon): raise TestFailure()
    if np.any(x0 > header0['NAXIS1'] + 0.5 + epsilon): raise TestFailure()
    if np.any(y0 > header0['NAXIS2'] + 0.5 + epsilon): raise TestFailure()
    if iheader == 0:
        if y0[0] > 1.5 or x0[-1] < header0['NAXIS1']-0.5:
            raise TestFailure()
    if iheader == len(headers)-1:
        if x0[0] > 1.5 or y0[-1] < header0['NAXIS2']-0.5:
            raise TestFailure()
