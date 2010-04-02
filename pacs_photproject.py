#!/usr/bin/python
from matplotlib.pyplot import figure, plot, show
import numpy
from optparse import OptionParser
import sys
import tamasisfortran as tmf
from tamasis import Map, PacsObservation, PacsProjectionSharpEdges

parser = OptionParser('Usage: %prog [options] fitsfile...')
parser.add_option('-o', help='write output map to FILE in FITS format [default: %default]', metavar='FILE', dest='outputfile', default='photproject.fits')
parser.add_option('--header', help='use FITS header in FILE to specify the map projection [default: automatic]', metavar='FILE')
parser.add_option('--resolution', help='input pixel size of the map in arcseconds [default: %default]', default=3.2)
parser.add_option('-n', '--npixels-per-sample', help='Maximum number of sky pixels intercepted by a PACS detector [default: %default]', default=6)
#parser.add_option('--no-flatfield', help='do not divide by calibration flat-field [default: False]', dest='do_flatfield', action='store_false', default=True)
#parser.add_option('-f', '--filtering', help='method for timeline filtering: mean or none [default: %default]', metavar='METHOD', default='mean')
parser.add_option('-d', '--deglitching', help='method for timeline deglitching: l2std, l2mad or none [default: %default]', metavar='METHOD', default='none')
parser.add_option('--nsigma-deglitching', help='N-sigma deglitching value [default: %default]', dest='nsigma', default=5.)
parser.add_option('--plot', help='plot the map', action='store_true', dest='do_plot', default=False)

(options, filename) = parser.parse_args(sys.argv[1:])

if len(filename) == 0:
    raise SystemExit(parser.print_help() or 1)

# check options
options.deglitching = options.deglitching.lower()
if options.deglitching not in ('none', 'l2std', 'l2mad'):
    raise ValueError("Invalid deglitching method '"+options.deglitching+"'. Valid methods are 'l2std', 'l2mad' or 'none'.")

#options.filtering = options.filtering.lower()
#if options.filtering not in ('none', 'mean'):
#    raise ValueError("Invalid filtering method '"+options.filtering+"'. Valid methods are 'mean', or 'none'.")

pacs = PacsObservation(filename=filename,
                       header=options.header,
                       resolution=options.resolution,
                       fine_sampling_factor=1,
                       keep_bad_detectors=False,
                       npixels_per_sample=options.npixels_per_sample)

projection = PacsProjectionSharpEdges(pacs, finer_sampling=False)

nx = pacs.header['naxis1']
ny = pacs.header['naxis2']

tod = pacs.get_tod()

# deglitching
if options.deglitching != 'none':
    nbads = numpy.sum(tod.mask)
    if options.deglitching == 'l2std':
        tod.mask = tmf.deglitch_l2b_std(projection.pmatrix, nx, ny, tod, tod.mask.astype('int8'), options.nsigma, pacs.npixels_per_sample)
    else:
        tod.mask = tmf.deglitch_l2b_mad(projection.pmatrix, nx, ny, tod, tod.mask.astype('int8'), options.nsigma, pacs.npixels_per_sample)
    print 'Number of glitches detected:', numpy.sum(tod.mask) - nbads

mymap = Map.zeros((nx, ny), order='f', header=pacs.header)
tmf.backprojection_weighted(projection.pmatrix, tod, tod.mask.astype('int8'), mymap, pacs.npixels_per_sample)


mymap.writefits(options.outputfile)

if options.do_plot:
     mymap.imshow()
     show()
#    idetector = 0
#    figure()
#    plot(tod[:,idetector])
#    index=numpy.where(tod.mask[:,idetector])
#    plot(index, tod.data[index,idetector], 'ro')
#    show()



