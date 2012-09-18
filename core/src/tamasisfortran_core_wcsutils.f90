module wcsutils

  use module_math,           only : DEG2RAD, NaN
  use module_pointingmatrix, only : PointingElement, xy2roi_ => xy2roi, roi2pmatrix_ => roi2pmatrix
  use module_tamasis,        only : p
  implicit none

contains

    subroutine create_grid_square(nx, ny, size, filling_factor, xreflection, yreflection, rotation, xcenter, ycenter, coords)
        !f2py threadsafe
        integer, intent(in)    :: nx                 ! number of detectors along the x axis
        integer, intent(in)    :: ny                 ! number of detectors along the y axis
        real(p), intent(in)    :: size               ! size of the detector placeholder
        real(p), intent(in)    :: filling_factor     ! fraction of transmitting detector area
        logical, intent(in)    :: xreflection        ! reflection along the x-axis (before rotation)
        logical, intent(in)    :: yreflection        ! reflection along the y-axis (before rotation)
        real(p), intent(in)    :: rotation           ! counter-clockwise rotation in degrees (before translation)
        real(p), intent(in)    :: xcenter, ycenter   ! coordinates of the grid center
        real(p), intent(inout) :: coords(2,4,nx,ny)  ! output coordinates of the detector corners (first dimension is x and y)

        integer :: i, j, k
        real(p) :: x, y, x0, y0, size_eff, r11, r12, r21, r22, tmp

        size_eff = size * sqrt(filling_factor)
        r11 = cos(DEG2RAD * rotation)
        r21 = sin(DEG2RAD * rotation)
        r12 = -r21
        r22 = r11
        if (xreflection) then
            r11 = -r11
            r21 = -r21
        end if
        if (yreflection) then
            r12 = -r12
            r22 = -r22
        end if
        x0 = -0.5_p * ((nx + 1) * size + size_eff)
        y0 = -0.5_p * ((ny + 1) * size + size_eff)
        
        do j = 1, ny
            y = y0 + size * j
            do i = 1, nx 
                x = x0 + size * i
                coords(1,1,i,j) = x
                coords(2,1,i,j) = y
                coords(1,2,i,j) = x + size_eff
                coords(2,2,i,j) = y
                coords(1,3,i,j) = x + size_eff
                coords(2,3,i,j) = y + size_eff
                coords(1,4,i,j) = x
                coords(2,4,i,j) = y + size_eff
                do k = 1, 4
                    tmp = coords(1,k,i,j)
                    coords(1,k,i,j) = xcenter + r11 * coords(1,k,i,j) + r12 * coords(2,k,i,j)
                    coords(2,k,i,j) = ycenter + r21 * tmp             + r22 * coords(2,k,i,j)
                end do
            end do
        end do
        
    end subroutine create_grid_square

    subroutine coordinate_boundary_inplace(input, boundary, ndims, ncoords)
        !f2py threadsafe
        real(p), intent(inout) :: input(ndims, ncoords)
        real(p), intent(in)    :: boundary(ndims, 2)
        integer                :: ndims
        integer*8              :: ncoords

        integer*8 :: i, j

        do j = 1, ncoords
            do i = 1, ndims
                if (input(i,j) < boundary(i,1) .or. input(i,j) > boundary(i,2)) then
                    input(i,j) = NaN
                end if
            end do
        end do

    end subroutine coordinate_boundary_inplace

    subroutine coordinate_boundary_outplace(input, output, boundary, ndims, ncoords)
        !f2py threadsafe
        real(p), intent(in)  :: input(ndims, ncoords)
        real(p), intent(out) :: output(ndims, ncoords)
        real(p), intent(in)  :: boundary(ndims, 2)
        integer              :: ndims
        integer*8            :: ncoords

        integer*8 :: i, j

        do j = 1, ncoords
            do i = 1, ndims
                if (input(i,j) < boundary(i,1) .or. input(i,j) > boundary(i,2)) then
                    output(i,j) = NaN
                else
                    output(i,j) = input(i,j)
                end if
            end do
        end do

    end subroutine coordinate_boundary_outplace

    subroutine xy2roi(xy, ncoords, nvertices, roi)
        !f2py threadsafe
        integer*8, intent(in) :: ncoords
        integer, intent(in)   :: nvertices
        real(p), intent(in)   :: xy(2, ncoords*nvertices)
        integer, intent(out)  :: roi(2, 2, ncoords)

        roi = xy2roi_(xy, nvertices)

    end subroutine xy2roi


    subroutine roi2pmatrix(roi, coords, ncoords, nvertices, npixels_per_sample, nx, ny, pmatrix, new_npixels_per_sample)
        !f2py threadsafe
        !f2py integer*8, depend(ncoords,npixels_per_sample) :: pmatrix(ncoords*npixels_per_sample)
        integer*8, intent(in)                :: ncoords
        integer, intent(in)                  :: nvertices
        integer, intent(in)                  :: npixels_per_sample
        integer, intent(in)                  :: roi(2, 2, ncoords)
        real(p), intent(in)                  :: coords(2, nvertices*ncoords)
        integer, intent(in)                  :: nx, ny
        type(PointingElement), intent(inout) :: pmatrix(npixels_per_sample, ncoords)
        integer, intent(out)                 :: new_npixels_per_sample
        
        logical :: out
        call roi2pmatrix_(roi, nvertices, coords, nx, ny, new_npixels_per_sample, out, pmatrix)

    end subroutine roi2pmatrix


end module wcsutils
