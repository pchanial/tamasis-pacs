! Copyright 2012 Pierre Chanial
! All rights reserved
!
! MPI routines
!
! Author: P. Chanial

subroutine diff_mpi(input, output, asize, dim, ashape, arank, comm, status)

    use module_math,    only : diff_fast, diff_medium, diff_slow
    use module_tamasis, only : p
#ifdef HAVE_MPI_MODULE
    use mpi
#endif
    implicit none

    integer*8, intent(in)  :: asize
    real(p), intent(in)    :: input(asize)
    real(p), intent(inout) :: output(asize)
    integer, intent(in)    :: dim
    integer, intent(in)    :: arank
    integer*8, intent(in)  :: ashape(arank)
    integer, intent(in)    :: comm
    integer, intent(out)   :: status

    integer                :: rank, size, nfast, ndiff, nslow, mpistatus(6)
    real(p), allocatable   :: boundary(:)

#ifdef HAVE_MPI_HEADER
    include 'mpif.h'
#endif

    nfast = product(ashape(1:dim-1))
    ndiff = ashape(dim)
    nslow = product(ashape(dim+1:arank))

    if (dim == arank) then
        call MPI_Comm_rank(comm, rank, status)
        if (status /= 0) return
        call MPI_Comm_size(comm, size, status)
        if (status /= 0) return
        allocate (boundary(nfast))
        call MPI_Sendrecv(input(1:nfast), nfast, MPI_DOUBLE_PRECISION, modulo(rank-1, size), 0,                                    &
                          boundary, nfast, MPI_DOUBLE_PRECISION, modulo(rank+1, size), 0, comm, mpistatus, status)
        if (status /= 0) return
        call diff_slow(input(1), output(1), nfast, ndiff, boundary, rank == size - 1)
    else if (dim > 1) then
        call diff_medium(input(1), output(1), nfast, ndiff, nslow)
    else
        call diff_fast(input(1), output(1), ndiff, nslow)
    end if

end subroutine diff_mpi


!-----------------------------------------------------------------------------------------------------------------------------------


subroutine diffT_mpi(input, output, asize, dim, ashape, arank, comm, status)

    use module_math,    only : diffT_fast, diffT_medium, diffT_slow
    use module_tamasis, only : p
#ifdef HAVE_MPI_MODULE
    use mpi
#endif
    implicit none

    integer*8, intent(in)  :: asize
    real(p), intent(in)    :: input(asize)
    real(p), intent(inout) :: output(asize)
    integer, intent(in)    :: dim
    integer, intent(in)    :: arank
    integer*8, intent(in)  :: ashape(arank)
    integer, intent(in)    :: comm
    integer, intent(out)   :: status

    integer                :: rank, size, nfast, ndiff, nslow, mpistatus(6)
    real(p), allocatable   :: boundary(:)

#ifdef HAVE_MPI_HEADER
    include 'mpif.h'
#endif

    nfast = product(ashape(1:dim-1))
    ndiff = ashape(dim)
    nslow = product(ashape(dim+1:arank))

    if (dim == arank) then
        call MPI_Comm_rank(comm, rank, status)
        if (status /= 0) return
        call MPI_Comm_size(comm, size, status)
        if (status /= 0) return
        allocate (boundary(nfast))
        call MPI_Sendrecv(input(asize-nfast+1:asize), nfast, MPI_DOUBLE_PRECISION, modulo(rank+1, size), 0,                        &
                          boundary, nfast, MPI_DOUBLE_PRECISION, modulo(rank-1, size), 0, comm, mpistatus, status)
        if (status /= 0) return
        call diffT_slow(input(1), output(1), nfast, ndiff, boundary, rank == 0, rank == size - 1)
    else if (dim > 1) then
        call diffT_medium(input(1), output(1), nfast, ndiff, nslow)
    else
        call diffT_fast(input(1), output(1), ndiff, nslow)
    end if

end subroutine difft_mpi


!-----------------------------------------------------------------------------------------------------------------------------------


subroutine diffTdiff_mpi(input, output, asize, dim, ashape, arank, scalar, comm, status)

    use module_math,    only : diffTdiff_fast, diffTdiff_medium, diffTdiff_slow
    use module_tamasis, only : p
#ifdef HAVE_MPI_MODULE
    use mpi
#endif
    implicit none

    integer*8, intent(in)  :: asize
    real(p), intent(in)    :: input(asize)
    real(p), intent(inout) :: output(asize)
    integer, intent(in)    :: dim
    integer, intent(in)    :: arank
    integer*8, intent(in)  :: ashape(arank)
    real(p), intent(in)    :: scalar
    integer, intent(in)    :: comm
    integer, intent(out)   :: status

    integer                :: rank, size, nfast, ndiff, nslow, mpistatus(6)
    real(p), allocatable   :: boundary1(:), boundary2(:)

#ifdef HAVE_MPI_HEADER
    include 'mpif.h'
#endif

    nfast = product(ashape(1:dim-1))
    ndiff = ashape(dim)
    nslow = product(ashape(dim+1:arank))

    if (dim == arank) then
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
        call diffTdiff_slow(input(1), output(1), nfast, ndiff, scalar, boundary1, boundary2, rank == 0, rank == size - 1)
    else if (dim > 1) then
        call diffTdiff_medium(input(1), output(1), nfast, ndiff, nslow, scalar)
    else
        call diffTdiff_fast(input(1), output(1), ndiff, nslow, scalar)
    end if

end subroutine difftdiff_mpi


!-----------------------------------------------------------------------------------------------------------------------------------


subroutine mpi_allreducelocal(input, ninputs, mask, nmasks, output, noutputs, op, comm, status)
    
    use module_math,    only : mInf, pInf
    use module_tamasis, only : p
#ifdef HAVE_MPI_MODULE 
    use mpi
#endif
    implicit none
 
    real(p), intent(in)    :: input(ninputs)
    logical*1, intent(in)  :: mask(nmasks)
    real(p), intent(inout) :: output(noutputs)
    integer, intent(in)    :: ninputs, nmasks, noutputs, op, comm
    integer, intent(out)   :: status

    integer :: size, root, a, z, iinput, imask
    real(p) :: input_(noutputs), field

#ifdef HAVE_MPI_HEADER
    include 'mpif.h'
#endif

    call MPI_Comm_size(comm, size, status)
    if (status /= 0) return

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
 
    output = 0
    iinput = 1
    do root = 0, size - 1

        ! unpack the input under the control of mask
        a = root * noutputs + 1
        z = (root + 1) * noutputs
        do imask = a, z
            if (imask > nmasks) then
                input_(imask-a+1:) = field
                exit
            end if
            if (mask(imask)) then
                input_(imask-a+1) = field
            else
                input_(imask-a+1) = input(iinput)
                iinput = iinput + 1
            end if
        end do

        ! reduce on node root
        call MPI_Reduce(input_, output, noutputs, MPI_DOUBLE_PRECISION, op, root, comm, status)
        if (status /= 0) return

    end do

end subroutine mpi_allreducelocal


!-----------------------------------------------------------------------------------------------------------------------------------


subroutine mpi_allscatterlocal(input, ninputs, mask, nmasks, output, noutputs, comm, status)
    
    use module_tamasis, only : p
#ifdef HAVE_MPI_MODULE
    use mpi
#endif
    implicit none

    real(p), intent(in)    :: input(ninputs) ! local image of the distributed map
    logical*1, intent(in)  :: mask(nmasks)   ! global mask for the locally observed pixels
    real(p), intent(inout) :: output(noutputs)
    integer, intent(in)    :: ninputs, nmasks, noutputs, comm
    integer, intent(out)   :: status

    integer :: a_unpacked, z_unpacked, tag, nrecvs, nsends, mpistatus(6), ioutput, imask
    integer :: source, dest, dp, rank, size
    integer, allocatable :: a_packed(:), z_packed(:)

    logical*1 :: maskbuffer(ninputs)
    real(p)   :: databuffer(min(ninputs,noutputs))

#ifdef HAVE_MPI_HEADER
    include 'mpif.h'
#endif

    call MPI_Comm_size(comm, size, status)
    if (status /= 0) return
    call MPI_Comm_rank(comm, rank, status)
    if (status /= 0) return

    ! before exchanging data, compute the local map bounds in the packed output
    allocate (a_packed(0:size-1), z_packed(0:size-1))
    a_packed(0) = 1
    z_packed(0) = 0
    source = 0
    do
        do imask = source * ninputs + 1, min((source+1) * ninputs, nmasks)
            if (mask(imask)) cycle
            z_packed(source) = z_packed(source) + 1
        end do
        if (source == size - 1) exit
        source = source + 1
        a_packed(source) = min(z_packed(source-1) + 1, noutputs)
        z_packed(source) = a_packed(source) - 1
    end do

    ! loop over the distances between mpi nodes, logically organised in a ring
    ioutput = 1
    tag = 99
    do dp = 0, size - 1

        dest = modulo(rank + dp, size)
        source = modulo(rank - dp, size)

        ! send mask to 'dest' and receive it from source
        a_unpacked = min(dest * ninputs + 1, nmasks+1)
        z_unpacked = min((dest+1) * ninputs, nmasks)
        call MPI_Sendrecv(mask(a_unpacked:z_unpacked), z_unpacked-a_unpacked+1, MPI_BYTE, dest, tag,                               &
                          maskbuffer, ninputs, MPI_BYTE, source, tag, comm, mpistatus, status)
        if (status /= 0) return
        call MPI_Get_count(mpistatus, MPI_BYTE, nrecvs, status)
        if (status /= 0) return

        ! pack pixels observed by 'source' from local map of the current processor
        nsends = 0
        do imask = 1, nrecvs
            if (maskbuffer(imask)) cycle
            nsends = nsends + 1
            databuffer(nsends) = input(imask)
        end do

        ! send data to 'source' and receive them from 'dest'
        call MPI_Sendrecv(databuffer, nsends, MPI_DOUBLE_PRECISION, source, tag + 1, output(a_packed(dest)),                       &
             z_packed(dest) - a_packed(dest) + 1, MPI_DOUBLE_PRECISION, dest, tag + 1, comm, mpistatus, status)
        if (status /= 0) return
        call MPI_Get_count(mpistatus, MPI_DOUBLE_PRECISION, nrecvs, status)
        if (status /= 0) return

    end do

end subroutine mpi_allscatterlocal
