
############################
# inputs are files
############################

from tamasis import *
import numpy as np
exit()
print
print 'Input: file'
print '-----------'
filename = '/home/pchanial/work/pacs/data/transparent/NGC6946/1342184520_blue'
ndetectors = tmf.pacs_info_ndetectors(filename, True)
ij = tmf.pacs_info_ij(filename, True, ndetectors)
ntotsamples = tmf.pacs_info_nsamples(filename)
print 'ntotsamples: %i' % ntotsamples
print 'ndetectors: %i' % ndetectors
first = 20000
last  = 86000
npixels_per_sample = 9
nsamples = last - first + 1
fine_sampling_factor = 1
compression_factor = 1

header = tmf.pacs_map_header_file(filename, True, first, last, fine_sampling_factor, compression_factor)
signal, mask = tmf.pacs_timeline(filename, True, first, last, ndetectors)
sizeofPmatrix = npixels_per_sample*nsamples*ndetectors
print 'Allocating '+str(sizeofPmatrix/2.**17)+' MiB for the pointing matrix.'
pmatrix = np.zeros(sizeofPmatrix, dtype=np.int64)
tmf.pacs_pointing_matrix_file(filename, True, npixels_per_sample, first, last, ndetectors, fine_sampling_factor, compression_factor, header, pmatrix)

# nx and ny should be read from the header. They are hardcoded for now
nx = 379
ny = 422
print 'Allocating '+str(nx*ny/2.**17)+' MiB for the map.'
map2d = np.zeros((nx,ny), dtype=np.float64, order='fortran')
tmf.pacs_projection_sharp_edges_transpose(pmatrix, signal, map2d, npixels_per_sample)
tmf.pacs_projection_sharp_edges_direct(pmatrix, map2d, signal, npixels_per_sample)
tmf.apply_mask(signal, mask)


############################
# inputs are arrays
############################

print
print 'Input: arrays'
print '-------------'
time = np.arange(0.,100., 1./40)
nsamples = np.size(time)
ra0  = 20.
dec0 = 0.1
ra = np.linspace(ra0, ra0+0.1, nsamples)
dec = np.linspace(dec0, dec0+0.1, nsamples)
pa = np.zeros(nsamples)
chop = np.zeros(nsamples)
npixels_per_sample = 9
ndetectors = 250
sizeofPmatrix = npixels_per_sample * nsamples * ndetectors
print 'Allocating '+str(sizeofPmatrix/2.**17)+' MiB for the pointing matrix.'
pmatrix = np.zeros(sizeofPmatrix, dtype=np.int64)

header = tmf.pacs_map_header_array(time, ra, dec, pa, chop, 'blue', True)
tmf.pacs_pointing_matrix_array(time, ra, dec, pa, chop, npixels_per_sample, ndetectors, 'blue', True, header, pmatrix)

print 'Allocating '+str(np.product((nsamples,ndetectors))/2.**17)+' MiB for the 1 timeline.'
signal=np.ones((nsamples,ndetectors), order='fortran')

# nx and ny should be read from the header. They are hardcoded for now
nx = 139
ny = 139
print 'Allocating '+str(np.product((nx,ny))/2.**17)+' MiB for the map.'
map2d = np.zeros((nx,ny), dtype=np.float64, order='fortran')
tmf.pacs_projection_sharp_edges_transpose(pmatrix, signal, map2d, npixels_per_sample)

