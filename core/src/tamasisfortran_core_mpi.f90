! Copyright 2012 Pierre Chanial
! All rights reserved
!
! MPI routines
!
! Author: P. Chanial

module operators_mpi

    use module_tamasis, only : p
    use operators,      only : xxx_diff_slow_inplace, xxx_diff_slow_outplace, xxx_diffT_slow_inplace, xxx_diffT_slow_outplace,     &
                               xxx_diffTdiff_slow_inplace, xxx_diffTdiff_slow_outplace
#ifdef HAVE_MPI_MODULE
        use mpi
#endif
    implicit none

#ifdef HAVE_MPI_HEADER
        include 'mpif.h'
#endif

    public :: diff
    public :: diffT
    public :: diffTdiff
    public :: allreducelocal
    public :: allscatterlocal

contains

    subroutine diff(input, output, asize, dim, ashape, arank, inplace, comm, status)

        integer*8, intent(in)  :: asize
        real(p), intent(inout) :: input(asize)
        real(p), intent(inout) :: output(asize)
        integer, intent(in)    :: dim
        integer, intent(in)    :: arank
        integer*8, intent(in)  :: ashape(arank)
        logical, intent(in)    :: inplace
        integer, intent(in)    :: comm
        integer, intent(out)   :: status

        integer                :: rank, size, nfast, ndiff, nslow, mpistatus(6)
        real(p), allocatable   :: boundary(:)

        nfast = product(ashape(1:dim-1))
        ndiff = ashape(dim)
        nslow = product(ashape(dim+1:arank))

        call MPI_Comm_rank(comm, rank, status)
        if (status /= 0) return
        call MPI_Comm_size(comm, size, status)
        if (status /= 0) return
        allocate (boundary(nfast))

        call MPI_Sendrecv(input(1:nfast), nfast, MPI_DOUBLE_PRECISION, modulo(rank-1, size), 0,                                    &
             boundary, nfast, MPI_DOUBLE_PRECISION, modulo(rank+1, size), 0, comm, mpistatus, status)
        if (status /= 0) return

        if (inplace) then
            call xxx_diff_slow_inplace(input(1), nfast, ndiff, boundary, rank == size - 1)
        else
            call xxx_diff_slow_outplace(input(1), output(1), nfast, ndiff, boundary, rank == size - 1)
        end if

    end subroutine diff


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine diffT(input, output, asize, dim, ashape, inplace, arank, comm, status)

        integer*8, intent(in)  :: asize
        real(p), intent(inout) :: input(asize)
        real(p), intent(inout) :: output(asize)
        integer, intent(in)    :: dim
        integer, intent(in)    :: arank
        integer*8, intent(in)  :: ashape(arank)
        logical, intent(in)    :: inplace
        integer, intent(in)    :: comm
        integer, intent(out)   :: status

        integer                :: rank, size, nfast, ndiff, nslow, mpistatus(6)
        real(p), allocatable   :: boundary(:)

        nfast = product(ashape(1:dim-1))
        ndiff = ashape(dim)
        nslow = product(ashape(dim+1:arank))

        call MPI_Comm_rank(comm, rank, status)
        if (status /= 0) return
        call MPI_Comm_size(comm, size, status)
        if (status /= 0) return
        allocate (boundary(nfast))

        call MPI_Sendrecv(input(asize-nfast+1:asize), nfast, MPI_DOUBLE_PRECISION, modulo(rank+1, size), 0,                        &
             boundary, nfast, MPI_DOUBLE_PRECISION, modulo(rank-1, size), 0, comm, mpistatus, status)
        if (status /= 0) return

        if (inplace) then
            call xxx_diffT_slow_inplace(input(1), nfast, ndiff, boundary, rank == 0, rank == size - 1)
        else
            call xxx_diffT_slow_outplace(input(1), output(1), nfast, ndiff, boundary, rank == 0, rank == size - 1)
        end if

    end subroutine difft


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine diffTdiff(input, output, asize, dim, ashape, arank, scalar, inplace, comm, status)

        integer*8, intent(in)  :: asize
        real(p), intent(inout) :: input(asize)
        real(p), intent(inout) :: output(asize)
        integer, intent(in)    :: dim
        integer, intent(in)    :: arank
        integer*8, intent(in)  :: ashape(arank)
        real(p), intent(in)    :: scalar
        logical, intent(in)    :: inplace
        integer, intent(in)    :: comm
        integer, intent(out)   :: status

        integer                :: rank, size, nfast, ndiff, nslow, mpistatus(6)
        real(p), allocatable   :: boundary1(:), boundary2(:)

        nfast = product(ashape(1:dim-1))
        ndiff = ashape(dim)
        nslow = product(ashape(dim+1:arank))

        call MPI_Comm_rank(comm, rank, status)
        if (status /= 0) return
        call MPI_Comm_size(comm, size, status)
        if (status /= 0) return
        allocate (boundary1(nfast), boundary2(nfast))
        call MPI_Sendrecv(input(asize-nfast+1:asize), nfast, MPI_DOUBLE_PRECISION, modulo(rank+1, size), 0,                        &
             boundary1, nfast, MPI_DOUBLE_PRECISION, modulo(rank-1, size), 0, comm, mpistatus, status)
        if (status /= 0) return
        call MPI_Sendrecv(input(1:nfast), nfast, MPI_DOUBLE_PRECISION, modulo(rank-1, size), 0,                                    &
             boundary2, nfast, MPI_DOUBLE_PRECISION, modulo(rank+1, size), 0, comm, mpistatus, status)
        if (status /= 0) return

        if (inplace) then
            call xxx_diffTdiff_slow_inplace(input(1), nfast, ndiff, scalar, boundary1, boundary2, rank == 0, rank == size - 1)
        else
            call xxx_diffTdiff_slow_outplace(input(1), output(1), nfast, ndiff, scalar, boundary1, boundary2, rank==0, rank==size-1)
        end if

    end subroutine difftdiff


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine allreducelocal(input, ninput, mask, nmask, output, noutput, chunk_size, op, comm, status)

        use module_math, only : mInf, pInf

        real(p), intent(in)    :: input(ninput)   ! packed local input
        logical*1, intent(in)  :: mask(nmask)     ! global mask controlling how the input is packed
        real(p), intent(inout) :: output(noutput) ! local output (section of the global output) that will host the reduction
        integer, intent(in)    :: ninput, nmask, noutput, chunk_size, op, comm
        integer, intent(out)   :: status ! >0:MPI, -1: invalid mask, 

        integer                :: size, root, iinput, imask, nglobal, nlocal, a, z
        real(p)                :: field
        real(p), allocatable   :: buffer(:)

        call MPI_Comm_size(comm, size, status)
        if (status /= 0) return

        ! get the global number of work items
        if (modulo(nmask, chunk_size) /= 0) then
            status = -1
            return
        end if
        nglobal = nmask / chunk_size

        select case (op)
        case (MPI_PROD)
            field = 1
        case (MPI_MIN)
            field = pInf
        case (MPI_MAX)
            field = mInf
        case default
            field = 0
        end select

        allocate (buffer((nglobal / size + 1) * chunk_size))

        iinput = 1
        ! loop over output sections
        do root = 0, size - 1

            call distribute(nglobal, root, size, nlocal)
            call distribute_slice(nglobal, root, size, a, z)
            a = (a - 1) * chunk_size + 1
            z = z * chunk_size

            ! unpack input as a local section corresponding to the running output section, under the control of the mask
            do imask = a, z
                if (mask(imask)) then
                    buffer(imask-a+1) = field
                else
                    buffer(imask-a+1) = input(iinput)
                    iinput = iinput + 1
                end if
            end do

            ! reduce on node root
            call MPI_Reduce(buffer, output, nlocal * chunk_size, MPI_DOUBLE_PRECISION, op, root, comm, status)
            if (status /= 0) return

        end do

    end subroutine allreducelocal


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine allscatterlocal(input, ninput, mask, nmask, output, noutput, chunk_size, comm, status)

        real(p), intent(in)    :: input(ninput)   ! local input (section of the global input)
        logical*1, intent(in)  :: mask(nmask)     ! global mask controlling how to pack the output
        real(p), intent(inout) :: output(noutput) ! the packed output
        integer, intent(in)    :: ninput, nmask, noutput, chunk_size, comm
        integer, intent(out)   :: status          ! >0:MPI, -1: invalid mask, -2: invalid input, -3: transmission check

        integer :: a_unpacked, z_unpacked, tag, nrecv, nsend, mpistatus(6), ioutput, imask, nglobal, nlocal
        integer :: source, dest, dp, rank, size
        integer, allocatable :: a_packed(:), z_packed(:)

        logical*1 :: maskbuffer(ninput)
        real(p)   :: databuffer(ninput)

        call MPI_Comm_size(comm, size, status)
        if (status /= 0) return
        call MPI_Comm_rank(comm, rank, status)
        if (status /= 0) return

        ! get the global number of work items
        if (modulo(nmask, chunk_size) /= 0) then
            status = -1
            return
        end if
        nglobal = nmask / chunk_size

        call distribute(nglobal, rank, size, nlocal)
        if (ninput /= nlocal * chunk_size) then
            status = -2
            return
        end if

        ! before exchanging data, compute the local map bounds in the packed output
        allocate (a_packed(0:size-1), z_packed(0:size-1))
        a_packed(0) = 1
        z_packed(0) = 0
        source = 0
        do
            call distribute_slice(nglobal, source, size, a_unpacked, z_unpacked)
            a_unpacked = (a_unpacked - 1) * chunk_size + 1
            z_unpacked = z_unpacked * chunk_size

            ! unpack input as a local section corresponding to the running output section, under the control of the mask
            do imask = a_unpacked, z_unpacked
                if (mask(imask)) cycle
                z_packed(source) = z_packed(source) + 1
            end do
            if (source == size - 1) exit
            source = source + 1
            a_packed(source) = min(z_packed(source-1) + 1, noutput)
            z_packed(source) = a_packed(source) - 1
        end do

        ! loop over the distances between mpi nodes, logically organised in a ring
        ioutput = 1
        tag = 99
        do dp = 0, size - 1

            dest = modulo(rank + dp, size)
            source = modulo(rank - dp, size)

            ! send mask to 'dest' and receive it from 'source'
            call distribute_slice(nglobal, dest, size, a_unpacked, z_unpacked)
            a_unpacked = (a_unpacked - 1) * chunk_size + 1
            z_unpacked = z_unpacked * chunk_size
            call MPI_Sendrecv(mask(a_unpacked:z_unpacked), z_unpacked-a_unpacked+1, MPI_BYTE, dest, tag,                           &
                 maskbuffer, ninput, MPI_BYTE, source, tag, comm, mpistatus, status)
            if (status /= 0) return
            call MPI_Get_count(mpistatus, MPI_BYTE, nrecv, status)
            if (status /= 0) return

            if (nrecv /= ninput) then
                status = -3
                return
            end if

            ! pack pixels observed by 'source' from local map of the current processor
            nsend = 0
            do imask = 1, ninput
                if (maskbuffer(imask)) cycle
                nsend = nsend + 1
                databuffer(nsend) = input(imask)
            end do

            ! send data to 'source' and receive them from 'dest'
            call MPI_Sendrecv(databuffer, nsend, MPI_DOUBLE_PRECISION, source, tag + 1, output(a_packed(dest):z_packed(dest)),     &
                 z_packed(dest) - a_packed(dest) + 1, MPI_DOUBLE_PRECISION, dest, tag + 1, comm, mpistatus, status)
             if (status /= 0) return
            call MPI_Get_count(mpistatus, MPI_DOUBLE_PRECISION, nrecv, status)
            if (status /= 0) return

        end do

    end subroutine allscatterlocal


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine distribute(nglobal, rank, size, nlocal)

        integer, intent(in)  :: nglobal, rank, size
        integer, intent(out) :: nlocal

        nlocal = nglobal / size
        if (modulo(nglobal, size) > rank) then
            nlocal = nlocal + 1
        end if

    end subroutine distribute


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine distribute_slice(nglobal, rank, size, a, z)

        integer, intent(in)  :: nglobal, rank, size
        integer, intent(out) :: a, z
        integer              :: nlocal

        call distribute(nglobal, rank, size, nlocal)
        a = nglobal / size * rank + min(rank, modulo(nglobal, size)) + 1
        z = a + nlocal - 1

    end subroutine distribute_slice


end module operators_mpi
