module pointing

    ! This module contains routines about pointing manipulation.

    use module_fitstools,      only : ft_read_keyword
    use module_math,           only : DEG2RAD, mInf, pInf
    use module_pointingmatrix, only : PointingElement, xy2pmatrix, xy2roi, roi2pmatrix
    use module_projection,     only : convex_hull
    use module_tamasis,        only : p
    use module_wcs,            only : ad2xy_gnomonic, ad2xys_gnomonic, init_astrometry

    implicit none
    private

    public :: NEAREST_NEIGHBOUR, SHARP_EDGES

    public :: instrument2ad
    public :: instrument2xy_minmax
    public :: instrument2pmatrix_nearest_neighbour
    public :: instrument2pmatrix_sharp_edges

    integer, parameter :: NEAREST_NEIGHBOUR = 0
    integer, parameter :: SHARP_EDGES = 1

contains

    subroutine instrument2ad(input, output, ncoords, ra, dec, pa)
        ! Convert coordinates in the instrument frame (focal plane) into celestial coordinates,
        ! assuming a pointing direction and a position angle.
        ! The routine is not accurate at the poles.
        ! Input units are in arc seconds, output units in degrees.

        integer, intent(in)    :: ncoords           ! number of coordinates
        real(p), intent(in)    :: input(2,ncoords)  ! input coordinates in instrument frame
        real(p), intent(inout) :: output(2,ncoords) ! output in celestial coordinates
        real(p), intent(in)    :: ra, dec, pa       ! pointing direction is (0,0) in the local frame

        real(p) :: cospa, sinpa, c1, c2
        integer :: i

        cospa = cos(pa * DEG2RAD)
        sinpa = sin(pa * DEG2RAD)

        do i = 1, ncoords
            c1 = input(1,i) / 3600._p
            c2 = input(2,i) / 3600._p
            output(2,i) = dec + (-c1 * sinpa + c2 * cospa)
            output(1,i) = ra  + ( c1 * cospa + c2 * sinpa) / cos(output(2,i) * DEG2RAD)
        end do

    end subroutine instrument2ad


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine instrument2xy_minmax(coords, ncoords, ra, dec, pa, npointings, header, xmin, ymin, xmax, ymax, status)
        ! Return the minimum and maximum sky pixel coordinate values of coordinates in the instrument frame.

        integer*8, intent(in)                      :: ncoords, npointings     ! #coordinates, #pointings
        real(p), intent(in)                        :: coords(2,ncoords)       ! instrument frame coordinates
        real(p), intent(in), dimension(npointings) :: ra, dec, pa             ! input pointings in celestial coordinates
        character(len=2880), intent(in)            :: header                  ! input FITS header
        real(p), intent(out)                       :: xmin, ymin, xmax, ymax  ! min and max values of the map coordinates
        integer, intent(out)                       :: status                  ! status flag

        real(p), allocatable :: hull_instrument(:,:), hull(:,:)
        integer, allocatable :: ihull(:)
        integer*8            :: ipointing

        call init_astrometry(header, status=status)
        if (status /= 0) return

        call convex_hull(coords, ihull)
        allocate (hull_instrument(2,size(ihull)), hull(2,size(ihull)))
        hull_instrument = coords(:,ihull)

        xmin = pInf
        xmax = mInf
        ymin = pInf
        ymax = mInf

#ifndef IFORT
        !$omp parallel do reduction(min:xmin,ymin) reduction(max:xmax,ymax) private(hull)
#endif
        do ipointing = 1, npointings

            call instrument2ad(hull_instrument, hull, size(ihull), ra(ipointing), dec(ipointing), pa(ipointing))
            hull = ad2xy_gnomonic(hull)
            xmin = min(xmin, minval(hull(1,:)))
            ymin = min(ymin, minval(hull(2,:)))
            xmax = max(xmax, maxval(hull(1,:)))
            ymax = max(ymax, maxval(hull(2,:)))

        end do
#ifndef IFORT
        !$omp end parallel do
#endif

    end subroutine instrument2xy_minmax


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine instrument2pmatrix_nearest_neighbour(coords, ncoords, area, ra, dec, pa, masked, npointings, header, pmatrix, out,  &
                                                    status)
        !f2py integer*8, depend(npointings,ncoords) :: pmatrix(npointings*ncoords)
        integer*8, intent(in)                        :: ncoords, npointings     ! #coordinates, #pointings
        real(p), intent(in)                          :: coords(2,ncoords)       ! instrument frame coordinates
        real(p), intent(in)                          :: area(ncoords)           ! detector area / reference_area
        real(p), intent(in), dimension(npointings)   :: ra, dec, pa             ! input pointings in celestial coordinates
        logical*1, intent(in), dimension(npointings) :: masked                  ! pointing flags: true if masked, removed
        character(len=*), intent(in)                 :: header
        type(PointingElement), intent(inout)         :: pmatrix(1,npointings,ncoords)
        logical, intent(out)                         :: out
        integer, intent(out)                         :: status

        real(p)   :: coords2(2,ncoords), x(ncoords), y(ncoords), s(ncoords)
        integer*8 :: isample
        integer   :: nx, ny

        out = .false.
        call init_astrometry(header, status=status)
        if (status /= 0) return
        ! get the size of the map
        call ft_read_keyword(header, 'naxis1', nx, status=status)
        if (status /= 0) return
        call ft_read_keyword(header, 'naxis2', ny, status=status)
        if (status /= 0) return

        !$omp parallel do private(isample, coords2, x, y, s) reduction(.or. : out)
        ! loop over the samples which have not been removed
        do isample = 1, npointings

            if (masked(isample)) then
                pmatrix(1,isample,:)%pixel = -1
                pmatrix(1,isample,:)%weight = 0
                cycle
            end if

            call instrument2ad(coords, coords2, int(ncoords), ra(isample), dec(isample), pa(isample))
            
            call ad2xys_gnomonic(coords2, x, y, s)

            ! the input map has flux densities, not surface brightness
            ! f_idetector = f_imap * weight
            ! with weight = detector_area / pixel_area
            ! and pixel_area = reference_area / s
            call xy2pmatrix(x, y, nx, ny, out, pmatrix(1,isample,:))
            pmatrix(1,isample,:)%weight = real(s * area, kind=kind(pmatrix%weight))

        end do
        !$omp end parallel do

    end subroutine instrument2pmatrix_nearest_neighbour


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine instrument2pmatrix_sharp_edges(coords, ncoords, ra, dec, pa, masked, npointings, header, pmatrix,                   &
                                              npixels_per_sample, new_npixels_per_sample, out, status)
        !f2py integer*8, depend(npixels_per_sample,npointings,ncoords) :: pmatrix(npixels_per_sample*npointings*ncoords/4)
        integer*8, intent(in)                        :: ncoords, npointings ! #coordinates, #pointings
        real(p), intent(in)                          :: coords(2,ncoords)   ! instrument frame coordinates
        real(p), intent(in), dimension(npointings)   :: ra, dec, pa         ! input pointings in celestial coordinates
        logical*1, intent(in), dimension(npointings) :: masked              ! pointing flags: true if masked, removed
        character(len=*), intent(in)                 :: header              ! sky map FITS header
        type(PointingElement), intent(inout)         :: pmatrix(npixels_per_sample,npointings,ncoords/4) ! the pointing matrix
        integer, intent(in)  :: npixels_per_sample     ! input maximum number of sky pixels intersected by a detector
        integer, intent(out) :: new_npixels_per_sample ! actual maximum number of sky pixels intersected by a detector
        logical, intent(out) :: out                    ! true if some coordinates fall outside of the map
        integer, intent(out) :: status

        real(p)   :: coords2(2,ncoords)
        integer   :: roi(2,2,ncoords/4)
        integer*8 :: isample
        integer   :: nx, ny

        new_npixels_per_sample = 0
        out = .false.
        call init_astrometry(header, status=status)
        if (status /= 0) return
        ! get the size of the map
        call ft_read_keyword(header, 'naxis1', nx, status=status)
        if (status /= 0) return
        call ft_read_keyword(header, 'naxis2', ny, status=status)
        if (status /= 0) return

        !$omp parallel do private(isample, coords2, roi) &
        !$omp reduction(max : npixels_per_sample) reduction(.or. : out)
        do isample = 1, npointings

            if (masked(isample)) then
                pmatrix(:,isample,:)%pixel = -1
                pmatrix(:,isample,:)%weight = 0
                cycle
            end if

            call instrument2ad(coords, coords2, int(ncoords), ra(isample), dec(isample), pa(isample))
            
            coords2 = ad2xy_gnomonic(coords2)
            roi = xy2roi(coords2, 4)
            call roi2pmatrix(roi, 4, coords2, nx, ny, new_npixels_per_sample, out, pmatrix(:,isample,:))

        end do
        !$omp end parallel do

    end subroutine instrument2pmatrix_sharp_edges


end module pointing
