# Copyrights 2010-2011 Pierre Chanial
# All rights reserved
#
#!/usr/bin/env python2.6
#
# NAME: tamasis_photproject
# 
# DESCRIPTION: create a map from a set of PACS observations, by backprojecting 
# the timelines onto the sky map, and by dividing it by the weights, which are
# the backprojection of 1.
# The routine is meant to replicate HCSS' photproject using TAMASIS' tools
#
# Author: P. Chanial

import numpy as np
from   optparse import OptionParser
import sys
from   tamasis import *

# preprocessor options
parser = OptionParser('Usage: %prog [options] fitsfile...')
parser.add_option('--output-tod', help='write tod to disk in FITS format [default:'
                  ' %default]', action='store_true', dest='do_outputtod', default=False)
parser.add_option('--flatfielding', help='divide by calibration flat-fie'
                  'ld [default: %default]', dest='do_flatfielding', action='store_'
                  'true', default=False)
parser.add_option('--subtraction-mean', help='subtract mean value [default: %default]', dest='do_subtraction_mean', action='store_'
                  'true', default=False)
parser.add_option('--median-filtering', help='window length for timeline median filtering.', metavar='LENGTH')
parser.add_option('-d', '--deglitching', help='method for timeline deglitching:'
                  ' l2std, l2mad or none [default: %default]', metavar='METHOD',
                  default='none')
parser.add_option('--nsigma', help='N-sigma deglitching value [default: %defaul'
                  't]', default=5.)
parser.add_option('-n', '--npixels-per-sample', help='maximum number of sky pix'
                  'els intercepted by a PACS detector [default: 6]', default=6)

parser.add_option('--policy-inscan', help='Policy for in-scan frames [default: %default]', default='keep')

parser.add_option('--policy-turnaround', help='Policy for turn-around frames [default: %default]', default='keep')

parser.add_option('--policy-other', help='Policy for non in-scan and non turn-around frames [default: %default]', default='remove')

parser.add_option('--policy-invalid', help='Policy for invalid frames [default: %default]', default='mask')

# mapper options
parser.add_option('--output-map', help='write output map to disk in FITS format [default:'
                  ' %default]', action='store_true', dest='do_outputmap', default=True)
parser.add_option('--header', help='use FITS header in FILE to specify the map '
                  'projection [default: automatic]', metavar='FILE')
parser.add_option('--resolution', help='input pixel size of the map in arcsecon'
                  'ds [default: 3.2 for the blue channel, 6.4 for the red one]')
parser.add_option('--ds9', help='display the map using ds9', action='store_true', dest='do_ds9', default=False)

(options, filename) = parser.parse_args(sys.argv[1:])

if len(filename) == 0:
    raise SystemExit(parser.print_help() or 1)

# Check options
options.deglitching = options.deglitching.lower()
if options.deglitching not in ('none', 'l2std', 'l2mad'):
    raise ValueError("Invalid deglitching method '"+options.deglitching+"'. Val"
                     "id methods are 'l2std', 'l2mad' or 'none'.")

if options.median_filtering is not None:
    try:
        length = int(options.median_filtering)
    except:
        raise ValueError("Invalid filtering length '"+options.median_filtering+"'.")

if options.npixels_per_sample is not None:
    options.npixels_per_sample = int(options.npixels_per_sample)

# Set up the PACS observation(s)
obs = PacsObservation(filename,
                      policy_inscan=options.policy_inscan,
                      policy_turnaround=options.policy_turnaround,
                      policy_other=options.policy_other,
                      policy_invalid=options.policy_invalid)

# Read the timeline
tod = obs.get_tod(flatfielding=options.do_flatfielding,
                  subtraction_mean=options.do_subtraction_mean,
                  unit='Jy/arcsec^2')

if options.median_filtering is not None:
    tod = filter_median(tod, length)

# Set up the acquisition model. oversampling is set to False because
# photproject does not attempt to sample better than what is transmitted
if options.deglitching != 'none' or options.do_outputmap:
    projection = Projection(obs,
                            header=options.header,
                            resolution=options.resolution,
                            oversampling=False,
                            npixels_per_sample=options.npixels_per_sample)

# Deglitch
if options.deglitching != 'none':
    nbads = np.sum(tod.mask % 2)
    if options.deglitching == 'l2std':
        tod.mask = deglitch_l2std(tod, projection, nsigma=options.nsigma)
    else:
        tod.mask = deglitch_l2mad(tod, projection, nsigma=options.nsigma)

if options.do_outputtod:
    if len(tod.nsamples) != len(filename):
        raise ValueError('The number of tod slices is not the number of input filenames.')
    dest = 0
    tod_ = obs.unpack(tod)
    for nsamples, f in zip(tod.nsamples, filename):
        tod_[:,:,dest:dest+nsamples].save(f+'_tod.fits')
        dest += nsamples

if not options.do_outputmap:
    exit()

# Get map dimensions
nx = projection.header['naxis1']
ny = projection.header['naxis2']

# Backproject the timeline and divide it by the weight
print('Computing the map...')
mymap = mapper_naive(tod, projection, unit='Jy/pixel')

# Write resulting map as a FITS file
print('Writing the map...')
mymap.save(filename[0] + '_map.fits')

# Display map
if options.do_ds9:
     mymap.ds9()
