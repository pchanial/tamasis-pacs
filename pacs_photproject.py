import tamasisfortran as tmf
from tamasis import Map, PacsObservation, PacsProjectionSharpEdges
import numpy
from matplotlib.pyplot import figure, plot, show

def pacs_photproject(filename='/home/pchanial/work/pacs/data/transparent/NGC6946/1342184520_blue', first=20000, last=86000, header=None, resolution=3., npixels_per_sample=6):

    pacs = PacsObservation(filename=filename,
                           first=first,
                           last=last,
                           header=header,
                           resolution=resolution,
                           fine_sampling_factor=1,
                           keep_bad_detectors=False,
                           npixels_per_sample=npixels_per_sample)

    tod = pacs.get_tod()
    nx = pacs.header['naxis1']
    ny = pacs.header['naxis2']

    # compute pointing matrix
    projection = PacsProjectionSharpEdges(pacs)

    # deglitching
    tod.mask = tmf.deglitch_l2b_std(projection.pmatrix, nx, ny, tod, tod.mask, 5., pacs.npixels_per_sample)

    # back projection
    mymap = Map.zeros([nx, ny], order='f')
    tmf.backprojection_weighted(projection.pmatrix, tod, tod.mask, mymap, pacs.npixels_per_sample)

    return mymap

filename='/home/pchanial/work/pacs/data/transparent/NGC6946/1342184520_blue'
first=20000
last=86000
header=None
resolution=3.
npixels_per_sample=6
pacs = PacsObservation(filename=filename,
                       first=first,
                       last=last,
                       header=header,
                       resolution=resolution,
                       fine_sampling_factor=1,
                       keep_bad_detectors=False,
                       npixels_per_sample=npixels_per_sample)

projection = PacsProjectionSharpEdges(pacs)

nx = pacs.header['naxis1']
ny = pacs.header['naxis2']

# deglitching
tod = pacs.get_tod()
nbads = numpy.sum(tod.mask)
tod.mask = tmf.deglitch_l2b_std(projection.pmatrix, nx, ny, tod, tod.mask, 10., pacs.npixels_per_sample)
print 'Number of glitches detected:', numpy.sum(tod.mask) - nbads

mymap = Map.zeros([nx, ny], order='f')
tmf.backprojection_weighted(projection.pmatrix, tod, tod.mask, mymap, pacs.npixels_per_sample)

idetector = 0
figure()
plot(tod[:,idetector])
index=numpy.where(tod.mask[:,idetector])
plot(index, tod.data[index,idetector], 'ro')
show()

#! command parsing
#allocate(parser)
#call parser%init('pacs_photproject [options] fitsfile', 0, 1)
#call parser%add_option('', 'o', 'Filename of the output map (FITS format)',&
#                       has_value=.true., default='photproject.fits')
#call parser%add_option('header', 'h', 'Input FITS header of the map',   &
#                       has_value=.true.)
#call parser%add_option('resolution', '', 'Input pixel size of the map', &
#                       has_value=.true., default='3.')
#call parser%add_option('npixels-per-sample', 'n', 'Maximum number of sky &
#                       &pixels intercepted by a PACS detector',          &
#                       has_value=.true., default='6')
#call parser%add_option('no-flatfield','','Do not divide by calibration flat field')
#call parser%add_option('filtering', 'f', 'Timeline filtering (mean|none)', &
#                       has_value=.true., default='mean')
#call parser%add_option('deglitching', 'd', 'Timeline deglitching &
#                       &(l2bl2bc|none)', has_value=.true., default='none')
#call parser%add_option('nsigma-deglitching', '', 'N-sigma for deglitching',&
#                       has_value=.true., default='5.')
#call parser%add_option('first', '', 'First sample in timeline', &
#                       has_value=.true., default='12001')
#call parser%add_option('last', '', 'Last sample in timeline',   &
#                       has_value=.true., default='86000')
#
#call parser%parse(status)
#if (status == -1) stop
#if (status /=  0) stop 'Aborting.'
#
#call parser%print_options()
#
#if (parser%get_argument_count() == 1) then
#   infile = parser%get_argument(1,status)
#else
#   infile = filename
#end if
#outfile = parser%get_option('o', status) !XXX status gfortran bug
#headerfile = parser%get_option('header', status)
#resolution = parser%get_option_as_real('resolution', status)
#npixels_per_sample = parser%get_option_as_integer('npixels-per-sample', status)
#first = parser%get_option_as_integer('first', status)
#last  = parser%get_option_as_integer('last', status)
#do_flatfield = .not. parser%get_option_as_logical('no-flatfield', status)
#do_meansubtraction = parser%get_option('filtering', status) == 'mean'
#do_deglitching_l2b = parser%get_option('deglitching', status) == 'l2b'
#do_deglitching_l2c = parser%get_option('deglitching', status) == 'l2c'
#nsigma_deglitching = parser%get_option_as_real('nsigma-deglitching', status)
