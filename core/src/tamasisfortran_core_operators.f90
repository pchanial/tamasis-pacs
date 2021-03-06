module operators

    ! This module contains optimised operator routines

    use module_tamasis, only : p
    implicit none

    public :: diff
    public :: diffT
    public :: diffTdiff
    public :: invntt_uncorrelated_inplace
    public :: invntt_uncorrelated_outplace
    public :: pack_inplace
    public :: pack_outplace
    public :: unpack_inplace
    public :: unpack_outplace

contains

    subroutine diff(input, output, asize, dim, ashape, arank, inplace)

        integer*8, intent(in)  :: asize
        real(p), intent(inout) :: input(asize)
        real(p), intent(inout) :: output(asize)
        integer, intent(in)    :: dim
        integer, intent(in)    :: arank
        integer*8, intent(in)  :: ashape(arank)
        logical, intent(in)    :: inplace

        integer*8 :: nfast, ndiff, nslow

        nfast = product(ashape(1:dim-1))
        ndiff = ashape(dim)
        nslow = product(ashape(dim+1:arank))

        if (dim == 1) then
            if (inplace) then
                call xxx_diff_fast_inplace(input(1), ndiff, nslow)
            else
                call xxx_diff_fast_outplace(input(1), output(1), ndiff, nslow)
            end if
        else
            if (inplace) then
                call xxx_diff_medium_inplace(input(1), nfast, ndiff, nslow)
            else
                call xxx_diff_medium_outplace(input(1), output(1), nfast, ndiff, nslow)
            end if
        end if

    end subroutine diff


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine diffT(input, output, asize, dim, ashape, arank, inplace)

        integer*8, intent(in)  :: asize
        real(p), intent(inout) :: input(asize)
        real(p), intent(inout) :: output(asize)
        integer, intent(in)    :: dim
        integer, intent(in)    :: arank
        integer*8, intent(in)  :: ashape(arank)
        logical, intent(in)    :: inplace

        integer*8 :: nfast, ndiff, nslow

        nfast = product(ashape(1:dim-1))
        ndiff = ashape(dim)
        nslow = product(ashape(dim+1:arank))

        if (dim == 1) then
            if (inplace) then
                call xxx_diffT_fast_inplace(input(1), ndiff, nslow)
            else
                call xxx_diffT_fast_outplace(input(1), output(1), ndiff, nslow)
            end if
        else
            if (inplace) then
                call xxx_diffT_medium_inplace(input(1), nfast, ndiff, nslow)
            else
                call xxx_diffT_medium_outplace(input(1), output(1), nfast, ndiff, nslow)
            end if
        end if

    end subroutine diffT


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine diffTdiff(input, output, asize, dim, ashape, arank, scalar, inplace)

        integer*8, intent(in)  :: asize
        real(p), intent(inout) :: input(asize)
        real(p), intent(inout) :: output(asize)
        integer, intent(in)    :: dim
        integer, intent(in)    :: arank
        integer*8, intent(in)  :: ashape(arank)
        real(p), intent(in)    :: scalar
        logical, intent(in)    :: inplace

        integer*8 :: nfast, ndiff, nslow

        nfast = product(ashape(1:dim-1))
        ndiff = ashape(dim)
        nslow = product(ashape(dim+1:arank))

        if (dim == 1) then
            if (inplace) then
                call xxx_diffTdiff_fast_inplace(input(1), ndiff, nslow, scalar)
            else
                call xxx_diffTdiff_fast_outplace(input(1), output(1), ndiff, nslow, scalar)
            end if
        else
            if (inplace) then
                call xxx_diffTdiff_medium_inplace(input(1), nfast, ndiff, nslow, scalar)
            else
                call xxx_diffTdiff_medium_outplace(input(1), output(1), nfast, ndiff, nslow, scalar)
            end if
        end if

    end subroutine diffTdiff


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine xxx_diff_fast_inplace(input, m, n)

        integer*8, intent(in)  :: m, n
        real(p), intent(inout) :: input(m,n)

        integer*8 :: i, j

        if (m <= 1) then
            input = 0
            return
        end if

        !$omp parallel do private(i,j)
        do j = 1, n
            do i = 1, m-1
                input(i,j) = input(i,j) - input(i+1,j)
            end do
            input(m,j) = 0
        end do
        !$omp end parallel do

    end subroutine xxx_diff_fast_inplace


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine xxx_diff_fast_outplace(input, output, m, n)

        integer*8, intent(in)  :: m, n
        real(p), intent(in)    :: input(m,n)
        real(p), intent(inout) :: output(m,n)

        integer*8 :: i, j

        if (m <= 1) then
            output = 0
            return
        end if

        !$omp parallel do private(i,j)
        do j = 1, n
            do i = 1, m-1
                output(i,j) = input(i,j) - input(i+1,j)
            end do
            output(m,j) = 0
        end do
        !$omp end parallel do

    end subroutine xxx_diff_fast_outplace


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine xxx_diff_medium_inplace(input, m, n, o)

        integer*8, intent(in)  :: m, n, o
        real(p), intent(inout) :: input(m,n,o)

        integer, parameter     :: block = 4096
        integer*8              :: i, j, k, a, z

        if (n <= 1) then
            input = 0
            return
        end if

        !$omp parallel do private(i,j,k,a,z)
        do k = 1, o
            ! block computation to avoid cache exhaustion
            do i = 1, (m-1) / block + 1
                a = (i-1) * block + 1
                z = min(m, i * block)
                do j = 1, n - 1
                    input(a:z,j,k) = input(a:z,j,k) - input(a:z,j+1,k)
                end do
            end do
            input(:,n,k) = 0
        end do
        !$omp end parallel do

    end subroutine xxx_diff_medium_inplace


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine xxx_diff_medium_outplace(input, output, m, n, o)

        integer*8, intent(in)  :: m, n, o
        real(p), intent(in)    :: input(m,n,o)
        real(p), intent(inout) :: output(m,n,o)

        integer, parameter     :: block = 4096
        integer*8              :: i, j, k, a, z

        if (n <= 1) then
            output = 0
            return
        end if

        !$omp parallel do private(i,j,k,a,z)
        do k = 1, o
            ! block computation to avoid cache exhaustion
            do i = 1, (m-1) / block + 1
                a = (i-1) * block + 1
                z = min(m, i * block)
                do j = 1, n - 1
                    output(a:z,j,k) = input(a:z,j,k) - input(a:z,j+1,k)
                end do
            end do
            output(:,n,k) = 0
        end do
        !$omp end parallel do

    end subroutine xxx_diff_medium_outplace


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine xxx_diff_slow_inplace(input, m, n, boundary, islastrank)

        integer*8, intent(in)         :: m, n
        real(p), intent(inout)        :: input(m,n)
        real(p), intent(in), optional :: boundary(m)
        logical, intent(in), optional :: islastrank

        integer, parameter :: block = 4096
        integer*8          :: i, j, a, z
        logical            :: islastrank_

        if (n == 0) return

        islastrank_ = .true.
        if (present(islastrank)) islastrank_ = islastrank

        !$omp parallel do private(i,j,a,z)
        do i = 1, (m-1) / block + 1
            a = (i-1) * block + 1
            z = min(m, i * block)
            do j = 1, n - 1
                input(a:z,j) = input(a:z,j) - input(a:z,j+1)
            end do
            if (present(boundary) .and. .not. islastrank_) then
                input(a:z,n) = input(a:z,n) - boundary(a:z)
            else
                input(a:z,n) = 0
            end if
        end do
        !$omp end parallel do

    end subroutine xxx_diff_slow_inplace


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine xxx_diff_slow_outplace(input, output, m, n, boundary, islastrank)

        integer*8, intent(in)         :: m, n
        real(p), intent(in)           :: input(m,n)
        real(p), intent(inout)        :: output(m,n)
        real(p), intent(in), optional :: boundary(m)
        logical, intent(in), optional :: islastrank

        integer, parameter :: block = 4096
        integer*8          :: i, j, a, z
        logical            :: islastrank_

        if (n == 0) return

        islastrank_ = .true.
        if (present(islastrank)) islastrank_ = islastrank

        !$omp parallel do private(i,j,a,z)
        do i = 1, (m-1) / block + 1
            a = (i-1) * block + 1
            z = min(m, i * block)
            do j = 1, n - 1
                output(a:z,j) = input(a:z,j) - input(a:z,j+1)
            end do
            if (present(boundary) .and. .not. islastrank_) then
                output(a:z,n) = input(a:z,n) - boundary(a:z)
            else
                output(a:z,n) = 0
            end if
        end do
        !$omp end parallel do

    end subroutine xxx_diff_slow_outplace


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine xxx_diffT_fast_inplace(input, m, n)

        integer*8, intent(in)  :: m, n
        real(p), intent(inout) :: input(m,n)

        integer*8 :: i, j

        if (m <= 1) then
            input = 0
            return
        end if

        !$omp parallel do private(i,j)
        do j = 1, n
            input(m,j) = -input(m-1,j)
            do i = m-1, 2, -1
                input(i,j) = input(i,j) - input(i-1,j)
            end do
            input(1,j) = input(1,j)
        end do
        !$omp end parallel do

    end subroutine xxx_diffT_fast_inplace


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine xxx_diffT_fast_outplace(input, output, m, n)

        integer*8, intent(in)  :: m, n
        real(p), intent(in)    :: input(m,n)
        real(p), intent(inout) :: output(m,n)

        integer*8 :: i, j

        if (m <= 1) then
            output = 0
            return
        end if

        !$omp parallel do private(i,j)
        do j = 1, n
            output(m,j) = -input(m-1,j)
            do i = m-1, 2, -1
                output(i,j) = input(i,j) - input(i-1,j)
            end do
            output(1,j) = input(1,j)
        end do
        !$omp end parallel do

    end subroutine xxx_diffT_fast_outplace


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine xxx_diffT_medium_inplace(input, m, n, o)

        integer*8, intent(in)  :: m, n, o
        real(p), intent(inout) :: input(m,n,o)

        integer, parameter     :: block = 4096
        integer*8              :: i, j, k, a, z

        if (n <= 1) then
            input = 0
            return
        end if

        !$omp parallel do private(i,j,k,a,z)
        do k = 1, o
            do i = 1, (m-1) / block + 1
                a = (i-1) * block + 1
                z = min(m, i * block)
                input(a:z,n,k) = -input(a:z,n-1,k)
                do j = n-1, 2, -1
                    input(a:z,j,k) = input(a:z,j,k) - input(a:z,j-1,k)
                end do
                input(a:z,1,k) = input(a:z,1,k)
            end do
        end do
        !$omp end parallel do

    end subroutine xxx_diffT_medium_inplace


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine xxx_diffT_medium_outplace(input, output, m, n, o)

        integer*8, intent(in)  :: m, n, o
        real(p), intent(in)    :: input(m,n,o)
        real(p), intent(inout) :: output(m,n,o)

        integer, parameter     :: block = 4096
        integer*8              :: i, j, k, a, z

        if (n <= 1) then
            output = 0
            return
        end if

        !$omp parallel do private(i,j,k,a,z)
        do k = 1, o
            do i = 1, (m-1) / block + 1
                a = (i-1) * block + 1
                z = min(m, i * block)
                output(a:z,n,k) = -input(a:z,n-1,k)
                do j = n-1, 2, -1
                    output(a:z,j,k) = input(a:z,j,k) - input(a:z,j-1,k)
                end do
                output(a:z,1,k) = input(a:z,1,k)
            end do
        end do
        !$omp end parallel do

    end subroutine xxx_diffT_medium_outplace


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine xxx_diffT_slow_inplace(input, m, n, boundary, isfirstrank, islastrank)

        integer*8, intent(in)         :: m, n
        real(p), intent(inout)        :: input(m,n)
        real(p), intent(in), optional :: boundary(m)
        logical, intent(in), optional :: isfirstrank
        logical, intent(in), optional :: islastrank

        integer, parameter :: block = 4096
        integer*8          :: i, j, a, z
        logical            :: isfirstrank_, islastrank_

        if (n == 0) return

        isfirstrank_ = .true.
        islastrank_ = .true.
        if (present(isfirstrank)) isfirstrank_ = isfirstrank
        if (present(islastrank)) islastrank_ = islastrank

        if (n == 1 .and. isfirstrank_ .and. islastrank_) then
            input = 0
            return
        end if

        !$omp parallel do private(i,j,a,z)
        do i = 1, (m-1) / block + 1
            a = (i-1) * block + 1
            z = min(m, i * block)
            if (n == 1) then
                if (isfirstrank_) then
                    continue
                else if (islastrank_) then
                    input(a:z,1) = -boundary(a:z)
                else
                    input(a:z,1) = input(a:z,1) - boundary(a:z)
                end if
                cycle
            end if                    
            if (islastrank_) then
                input(a:z,n) = -input(a:z,n-1)
            else
                input(a:z,n) = input(a:z,n) - input(a:z,n-1)
            end if
            do j = n-1, 2, -1
                input(a:z,j) = input(a:z,j) - input(a:z,j-1)
            end do
            if (.not. isfirstrank_) then
                input(a:z,1) = input(a:z,1) - boundary(a:z)
            end if
        end do
        !$omp end parallel do

    end subroutine xxx_diffT_slow_inplace


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine xxx_diffT_slow_outplace(input, output, m, n, boundary, isfirstrank, islastrank)

        integer*8, intent(in)         :: m, n
        real(p), intent(in)           :: input(m,n)
        real(p), intent(inout)        :: output(m,n)
        real(p), intent(in), optional :: boundary(m)
        logical, intent(in), optional :: isfirstrank
        logical, intent(in), optional :: islastrank

        integer, parameter :: block = 4096
        integer*8          :: i, j, a, z
        logical            :: isfirstrank_, islastrank_

        if (n == 0) return

        isfirstrank_ = .true.
        islastrank_ = .true.
        if (present(isfirstrank)) isfirstrank_ = isfirstrank
        if (present(islastrank)) islastrank_ = islastrank

        if (n == 1 .and. isfirstrank_ .and. islastrank_) then
            output = 0
            return
        end if

        !$omp parallel do private(i,j,a,z)
        do i = 1, (m-1) / block + 1
            a = (i-1) * block + 1
            z = min(m, i * block)
            if (n == 1) then
                if (isfirstrank_) then
                    output(a:z,1) = input(a:z,1)
                else if (islastrank_) then
                    output(a:z,1) = -boundary(a:z)
                else
                    output(a:z,1) = input(a:z,1) - boundary(a:z)
                end if
                cycle
            end if                    
            if (islastrank_) then
                output(a:z,n) = -input(a:z,n-1)
            else
                output(a:z,n) = input(a:z,n) - input(a:z,n-1)
            end if
            do j = n-1, 2, -1
                output(a:z,j) = input(a:z,j) - input(a:z,j-1)
            end do
            if (.not. isfirstrank_) then
                output(a:z,1) = input(a:z,1) - boundary(a:z)
            else
                output(a:z,1) = input(a:z,1)
            end if
        end do
        !$omp end parallel do

    end subroutine xxx_diffT_slow_outplace


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine xxx_diffTdiff_fast_inplace(input, m, n, scalar)

        integer*8, intent(in)         :: m, n
        real(p), intent(inout)        :: input(m,n)
        real(p), intent(in), optional :: scalar

        integer*8 :: i, j
        real(p)   :: v, w, s

        if (m <= 1) then
            input = 0
            return
        end if

        if (present(scalar)) then
            s = scalar
        else
            s = 1._p
        end if

        !$omp parallel do private(i,j,v,w)
        do j = 1, n
            v = input(1,j)
            input(1,j) = s * (input(1,j) - input(2,j))
            do i = 2, m-1
                w = input(i,j)
                input(i,j) = s * (2 * w - v - input(i+1,j))
                v = w
            end do
            input(m,j) = s * (input(m,j) - v)
        end do
        !$omp end parallel do

    end subroutine xxx_diffTdiff_fast_inplace


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine xxx_diffTdiff_fast_outplace(input, output, m, n, scalar)

        integer*8, intent(in)         :: m, n
        real(p), intent(in)           :: input(m,n)
        real(p), intent(inout)        :: output(m,n)
        real(p), intent(in), optional :: scalar

        integer*8 :: i, j
        real(p)   :: v, w, s

        if (m <= 1) then
            output = 0
            return
        end if

        if (present(scalar)) then
            s = scalar
        else
            s = 1._p
        end if

        !$omp parallel do private(i,j,v,w)
        do j = 1, n
            v = input(1,j)
            output(1,j) = s * (input(1,j) - input(2,j))
            do i = 2, m-1
                w = input(i,j)
                output(i,j) = s * (2 * w - v - input(i+1,j))
                v = w
            end do
            output(m,j) = s * (input(m,j) - v)
        end do
        !$omp end parallel do

    end subroutine xxx_diffTdiff_fast_outplace


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine xxx_diffTdiff_medium_inplace(input, m, n, o, scalar)

        integer*8, intent(in)         :: m, n, o
        real(p), intent(inout)        :: input(m,n,o)
        real(p), intent(in), optional :: scalar

        integer, parameter :: block = 4096
        integer*8 :: h, i, j, k, a, z
        real(p)   :: w, v(block), s
        
        if (n <= 1) then
            input = 0
            return
        end if

        if (present(scalar)) then
            s = scalar
        else
            s = 1._p
        end if

        !$omp parallel do private(i,j,k,a,z,v,w)
        do k = 1, o
            do i = 1, (m-1) / block + 1
                a = (i-1) * block + 1
                z = min(i * block, m)
                v(1:z-a+1) = input(a:z,1,k)
                input(a:z,1,k) = s * (v(1:z-a+1) - input(a:z,2,k))
                do j = 2, n-1
                    do h = a, z
                        w = input(h,j,k)
                        input(h,j,k) = s * (2 * w - v(h-a+1) - input(h,j+1,k))
                        v(h-a+1) = w
                    end  do
                end do
                input(a:z,n,k) = s * (input(a:z,n,k) - v(1:z-a+1))
            end do
        end do
        !$omp end parallel do

    end subroutine xxx_diffTdiff_medium_inplace


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine xxx_diffTdiff_medium_outplace(input, output, m, n, o, scalar)

        integer*8, intent(in)         :: m, n, o
        real(p), intent(in)           :: input(m,n,o)
        real(p), intent(inout)        :: output(m,n,o)
        real(p), intent(in), optional :: scalar

        integer, parameter :: block = 4096
        integer*8 :: h, i, j, k, a, z
        real(p)   :: w, v(block), s
        
        if (n <= 1) then
            output = 0
            return
        end if

        if (present(scalar)) then
            s = scalar
        else
            s = 1._p
        end if

        !$omp parallel do private(i,j,k,a,z,v,w)
        do k = 1, o
            do i = 1, (m-1) / block + 1
                a = (i-1) * block + 1
                z = min(i * block, m)
                v(1:z-a+1) = input(a:z,1,k)
                output(a:z,1,k) = s * (v(1:z-a+1) - input(a:z,2,k))
                do j = 2, n-1
                    do h = a, z
                        w = input(h,j,k)
                        output(h,j,k) = s * (2 * w - v(h-a+1) - input(h,j+1,k))
                        v(h-a+1) = w
                    end  do
                end do
                output(a:z,n,k) = s * (input(a:z,n,k) - v(1:z-a+1))
            end do
        end do
        !$omp end parallel do

    end subroutine xxx_diffTdiff_medium_outplace


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine xxx_diffTdiff_slow_inplace(input, m, n, scalar, boundary1, boundary2, isfirstrank, islastrank)

        integer*8, intent(in)         :: m, n
        real(p), intent(inout)        :: input(m,n)
        real(p), intent(in), optional :: scalar
        real(p), intent(in), optional :: boundary1(m), boundary2(m)
        logical, intent(in), optional :: isfirstrank
        logical, intent(in), optional :: islastrank

        integer, parameter :: block = 4096
        integer*8          :: h, i, j, a, z
        logical            :: isfirstrank_, islastrank_
        real(p)            :: w, v(block), s

        if (n == 0) return

        if (present(scalar)) then
            s = scalar
        else
            s = 1._p
        end if

        isfirstrank_ = .true.
        islastrank_ = .true.
        if (present(isfirstrank)) isfirstrank_ = isfirstrank
        if (present(islastrank)) islastrank_ = islastrank

        if (n == 1 .and. isfirstrank_ .and. islastrank_) then
            input = 0
            return
        end if

        !$omp parallel do private(i,j,a,z,v,w)
        do i = 1, (m-1) / block + 1
            a = (i-1) * block + 1
            z = min(i * block, m)
            if (n == 1) then
                if (isfirstrank_) then
                    input(a:z,1) = s * (input(a:z,1) - boundary2(a:z))
                else if (islastrank_) then
                    input(a:z,1) = s * (input(a:z,1) - boundary1(a:z))
                else
                    input(a:z,1) = s * (2 * input(a:z,1) - boundary1(a:z) - boundary2(a:z))
                end if
                cycle
            end if
            v(1:z-a+1) = input(a:z,1)
            if (isfirstrank_) then
                input(a:z,1) = s * (v(1:z-a+1) - input(a:z,2))
            else
                input(a:z,1) = s * (2 * v(1:z-a+1) - boundary1(a:z) - input(a:z,2))
            end if
            do j = 2, n-1
                do h = a, z
                    w = input(h,j)
                    input(h,j) = s * (2 * w - v(h-a+1) - input(h,j+1))
                    v(h-a+1) = w
                end  do
            end do
            if (islastrank_) then
                input(a:z,n) = s * (input(a:z,n) - v(1:z-a+1))
            else
                input(a:z,n) = s * (2 * input(a:z,n) - v(1:z-a+1) - boundary2(a:z))
            end if
        end do
        !$omp end parallel do

    end subroutine xxx_diffTdiff_slow_inplace


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine xxx_diffTdiff_slow_outplace(input, output, m, n, scalar, boundary1, boundary2, isfirstrank, islastrank)

        integer*8, intent(in)         :: m, n
        real(p), intent(in)           :: input(m,n)
        real(p), intent(inout)        :: output(m,n)
        real(p), intent(in), optional :: scalar
        real(p), intent(in), optional :: boundary1(m), boundary2(m)
        logical, intent(in), optional :: isfirstrank
        logical, intent(in), optional :: islastrank

        integer, parameter :: block = 4096
        integer*8          :: h, i, j, a, z
        logical            :: isfirstrank_, islastrank_
        real(p)            :: w, v(block), s

        if (n == 0) return

        if (present(scalar)) then
            s = scalar
        else
            s = 1._p
        end if

        isfirstrank_ = .true.
        islastrank_ = .true.
        if (present(isfirstrank)) isfirstrank_ = isfirstrank
        if (present(islastrank)) islastrank_ = islastrank

        if (n == 1 .and. isfirstrank_ .and. islastrank_) then
            output = 0
            return
        end if

        !$omp parallel do private(i,j,a,z,v,w)
        do i = 1, (m-1) / block + 1
            a = (i-1) * block + 1
            z = min(i * block, m)
            if (n == 1) then
                if (isfirstrank_) then
                    output(a:z,1) = s * (input(a:z,1) - boundary2(a:z))
                else if (islastrank_) then
                    output(a:z,1) = s * (input(a:z,1) - boundary1(a:z))
                else
                    output(a:z,1) = s * (2 * input(a:z,1) - boundary1(a:z) - boundary2(a:z))
                end if
                cycle
            end if
            v(1:z-a+1) = input(a:z,1)
            if (isfirstrank_) then
                output(a:z,1) = s * (v(1:z-a+1) - input(a:z,2))
            else
                output(a:z,1) = s * (2 * v(1:z-a+1) - boundary1(a:z) - input(a:z,2))
            end if
            do j = 2, n-1
                do h = a, z
                    w = input(h,j)
                    output(h,j) = s * (2 * w - v(h-a+1) - input(h,j+1))
                    v(h-a+1) = w
                end  do
            end do
            if (islastrank_) then
                output(a:z,n) = s * (input(a:z,n) - v(1:z-a+1))
            else
                output(a:z,n) = s * (2 * input(a:z,n) - v(1:z-a+1) - boundary2(a:z))
            end if
        end do
        !$omp end parallel do

    end subroutine xxx_diffTdiff_slow_outplace


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine invntt_uncorrelated_inplace(input, ninputs, nsamples, istride,  fft_filter, filter_length, ndetectors, fplan, bplan,&
                                           left, right)

        real(p), intent(inout) :: input(ninputs)
        integer, intent(in)    :: ninputs
        integer, intent(in)    :: nsamples, istride
        integer, intent(in)    :: filter_length, ndetectors
        real(p), intent(in)    :: fft_filter(filter_length, ndetectors)
        integer*8, intent(in)  :: fplan, bplan
        integer, intent(in)    :: left, right

        real(p) :: buffer(filter_length)
        integer :: i

        !$omp parallel do private(buffer)
        do i = 1, ndetectors
            buffer(:left) = 0
            buffer(left+1:filter_length-right) = input((i-1)*istride+1:(i-1)*istride+nsamples)
            buffer(filter_length-right+1:) = 0
            call dfftw_execute_r2r(fplan, buffer, buffer)
            buffer = buffer * fft_filter(:,i)
            call dfftw_execute_r2r(bplan, buffer, buffer)
            input((i-1)*istride+1:(i-1)*istride+nsamples) = buffer(left+1:filter_length-right) / filter_length
        end do
        !$omp end parallel do

    end subroutine invntt_uncorrelated_inplace


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine invntt_uncorrelated_outplace(input, ninputs, nsamples, istride, output, noutputs, ostride, fft_filter,              &
                                            filter_length, ndetectors, fplan, bplan, left, right)

        real(p), intent(in)    :: input(ninputs)
        integer, intent(in)    :: ninputs
        integer, intent(in)    :: nsamples, istride
        real(p), intent(inout) :: output(noutputs)
        integer, intent(in)    :: noutputs
        integer, intent(in)    :: ostride
        integer, intent(in)    :: filter_length, ndetectors
        real(p), intent(in)    :: fft_filter(filter_length, ndetectors)
        integer*8, intent(in)  :: fplan, bplan
        integer, intent(in)    :: left, right

        real(p) :: buffer(filter_length)
        integer :: i

        !$omp parallel do private(buffer)
        do i = 1, ndetectors
            buffer(:left) = 0
            buffer(left+1:filter_length-right) = input((i-1)*istride+1:(i-1)*istride+nsamples)
            buffer(filter_length-right+1:) = 0
            call dfftw_execute_r2r(fplan, buffer, buffer)
            buffer = buffer * fft_filter(:,i)
            call dfftw_execute_r2r(bplan, buffer, buffer)
            output((i-1)*ostride+1:(i-1)*ostride+nsamples) = buffer(left+1:filter_length-right) / filter_length
        end do
        !$omp end parallel do

    end subroutine invntt_uncorrelated_outplace


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine pack_inplace(array, mask, ninputs)

        real(p), intent(inout) :: array(ninputs)
        logical*1, intent(in)  :: mask(ninputs)
        integer, intent(in)    :: ninputs

        integer :: ii, io

        io = 1
        do ii = 1, ninputs
            if (mask(ii)) cycle
            array(io) = array(ii)
            io = io + 1
        end do

    end subroutine pack_inplace


    !-------------------------------------------------------------------------------------------------------------------------------
    
    
    subroutine pack_outplace(input, mask, ninputs, output, noutputs)

        real(p), intent(in)    :: input(ninputs)
        logical*1, intent(in)  :: mask(ninputs)
        integer, intent(in)    :: ninputs
        real(p), intent(inout) :: output(noutputs)
        integer, intent(in)    :: noutputs

        output = pack(input, .not. mask)

    end subroutine pack_outplace


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine unpack_inplace(array, mask, noutputs, ninputs)

        real(p), intent(inout) :: array(noutputs)
        logical*1, intent(in)  :: mask(noutputs)
        integer, intent(in)    :: noutputs, ninputs

        integer :: ii, io

        ! unlike unpack, the following works in-place
        ii = ninputs
        do io = noutputs, 1, -1
            if (mask(io)) then
                array(io) = 0
            else
                array(io) = array(ii)
                ii = ii - 1
            end if
        end do

    end subroutine unpack_inplace


    !-------------------------------------------------------------------------------------------------------------------------------
    
    
    subroutine unpack_outplace(input, ninputs, mask, output, noutputs)

        real(p), intent(in)    :: input(ninputs)
        integer, intent(in)    :: ninputs
        logical*1, intent(in)  :: mask(noutputs)
        real(p), intent(inout) :: output(noutputs)
        integer, intent(in)    :: noutputs

        output = unpack(input, .not. mask, 0._p)

    end subroutine unpack_outplace


end module operators
