#!/usr/bin/env python2.6
#
# NAME: pacs_photproject
# 
# DESCRIPTION: create a map from a set of PACS observations, by backprojecting 
# the timelines onto the sky map, and by dividing it by the weights, which are
# the backprojection of 1.
# The routine is meant to replicate HCSS' photproject
#
# Author: P. Chanial

from   matplotlib.pyplot import figure, plot, show
import numpy
from   optparse import OptionParser
import sys
import tamasisfortran as tmf
from   tamasis import *

# set up the option parser
parser = OptionParser('Usage: %prog [options] fitsfile...')
parser.add_option('-o', help='write output map to FILE in FITS format [default:'
                  ' %default]', metavar='FILE', dest='outputfile',     default=
                  'photproject.fits')
parser.add_option('--header', help='use FITS header in FILE to specify the map '
                  'projection [default: automatic]', metavar='FILE')
parser.add_option('--resolution', help='input pixel size of the map in arcsecon'
                  'ds [default: %default]', default=3.2)
parser.add_option('-n', '--npixels-per-sample', help='maximum number of sky pix'
                  'els intercepted by a PACS detector [default: 5 for the blue '
                  'channel and 11 for the red one]')
parser.add_option('--no-flatfield', help='do not divide by calibration flat-fie'
                  'ld [default: False]', dest='do_flatfielding', action='store_'
                  'false', default=True)
parser.add_option('-f', '--filtering', help='method for timeline filtering: mea'
                  'n or none [default: %default]', metavar='METHOD', default='n'
                  'one')
parser.add_option('-l', '--length', help='filtering length [default: %default]',
                  default=200)
parser.add_option('-d', '--deglitching', help='method for timeline deglitching:'
                  ' l2std, l2mad or none [default: %default]', metavar='METHOD',
                  default='none')
parser.add_option('--nsigma', help='N-sigma deglitching value [default: %defaul'
                  't]', default=5.)
parser.add_option('--plot', help='plot the map', action='store_true', dest='do_'
                  'plot', default=False)

(options, filename) = parser.parse_args(sys.argv[1:])

if len(filename) == 0:
    raise SystemExit(parser.print_help() or 1)

# Check options
options.deglitching = options.deglitching.lower()
if options.deglitching not in ('none', 'l2std', 'l2mad'):
    raise ValueError("Invalid deglitching method '"+options.deglitching+"'. Val"
                     "id methods are 'l2std', 'l2mad' or 'none'.")

options.filtering = options.filtering.lower()
if options.filtering not in ('none', 'median'):
    raise ValueError("Invalid filtering method '"+options.filtering+"'. Valid m"
                     "ethods are 'median', or 'none'.")

if options.npixels_per_sample is not None:
    options.npixels_per_sample = int(options.npixels_per_sample)
options.length = int(options.length)

print

bad_detector_mask = None
# uncomment the following lines to make a map with fewer detectors
# 1 means bad detector
#
#bad_detector_mask = numpy.ones([32,64], dtype='int8')
#bad_detector_mask[0,0] = 0

# Set up the PACS observation(s)
obs = PacsObservation(filename=filename,
                      header=options.header,
                      resolution=options.resolution,
                      fine_sampling_factor=1,
                      bad_detector_mask = bad_detector_mask,
                      keep_bad_detectors=False)

# Get map dimensions
nx = obs.header['naxis1']
ny = obs.header['naxis2']

# Set up the acquisition model. oversampling is set to False because
# photproject does not attempt to sample better than what is transmitted
projection = Projection(obs, oversampling=False, npixels_per_sample=options.npixels_per_sample)

# Read the timeline
tod = obs.get_tod(do_flatfielding=options.do_flatfielding, do_subtraction_mean=options.filtering == 'none')

if options.filtering == 'median':
    tod = filter_median(tod, options.length)

# Deglitch
if options.deglitching != 'none':
    nbads = numpy.sum(tod.mask % 2)
    if options.deglitching == 'l2std':
        tod.mask = deglitch_l2std(tod, projection, nsigma=options.nsigma)
    else:
        tod.mask = deglitch_l2mad(tod, projection, nsigma=options.nsigma)

# Backproject the timeline and divide it by the weight
print 'Computing the map...'
mymap = Map.zeros((ny,nx), header=obs.header)
tmf.backprojection_weighted(projection.pmatrix, tod.T, tod.mask.T, 
                            mymap.T, projection.npixels_per_sample)

# Write resulting map as a FITS file
print 'Writing the map...'
mymap.writefits(options.outputfile)

# Plot the map
if options.do_plot:
     mymap.imshow()
     show()
#    idetector = 0
#    figure()
#    plot(tod[:,idetector])
#    index=numpy.where(tod.mask[:,idetector])
#    plot(index, tod.data[index,idetector], 'ro')
#    show()



