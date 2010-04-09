program test_deglitching

    use iso_fortran_env,       only : ERROR_UNIT, OUTPUT_UNIT
    use module_deglitching
    use module_math,           only : linspace
    use module_pointingmatrix, only : pointingelement, pmatrix_direct, xy2roi, roi2pmatrix
    use module_precision,      only : p
    implicit none

    type(pointingelement),allocatable :: pmatrix(:,:,:)
    integer                           :: i, j, itime, idetector, ipixel
    integer                           :: ntimes, ndetectors, npixels_per_sample
    integer                           :: nrepeats
    integer                           :: nx, ny
    real(kind=p), allocatable         :: xc(:), yc(:)
    real(kind=p), allocatable         :: map(:), signal(:,:)
    logical(kind=1), allocatable      :: mask(:,:)
    
    nrepeats = 6
    ndetectors = 3
    npixels_per_sample = 2
    nx = 6
    ny = 1
    ntimes = (nx-3)*nrepeats

    allocate (map(0:nx*ny-1))
    allocate (xc(ntimes))
    allocate (yc(ntimes))
    allocate (signal(ntimes,ndetectors))
    allocate (mask(ntimes,ndetectors))

    map = linspace(1._p, 2.25_p, nx*ny)
    do i=1,nrepeats
        xc((i-1)*(nx-3)+1:i*(nx-3)) = linspace(1._p, real(nx-3, kind=p), nx-3)
    end do
    yc = 1._p


    !----------------------------------------------------------------------
    ! test 1: detectors do not overlap map pixels (npixels_per_sample = 1)
    !----------------------------------------------------------------------
    !
    !    .-----------------------.
    !    | 0 | 1 | 2 | 3 | 4 | 5 |
    !    .-----------------------.
    ! D1: xxx ooo ooo
    !     xxx ooo ooo
    !     ooo ooo ooo
    !     ooo ooo ooo
    !     ooo ooo ooo
    !     ooo xxx ooo
    !
    ! D2:     ooo ooo xxx
    !         ooo ooo ooo
    !         ooo xxx xxx
    !         ooo ooo ooo
    !         ooo ooo ooo
    !         ooo ooo ooo
    !
    ! xxx : glitch
    ! ooo : no glitch

    allocate (pmatrix(1, ntimes, 2))
    mask = .false.
    call get_pmatrix(pmatrix(:,:,1:1), nx, ny, xc,      yc, 1._p)
    call get_pmatrix(pmatrix(:,:,2:2), nx, ny, xc+1._p, yc, 1._p)

    ! read timeline from map
    call pmatrix_direct(pmatrix(:,:,1:2), map, signal(:,1:2))
    ! add noise
    signal(:,1) = signal(:,1) + [([(0.001_p, j=1,nx-3)]*(i-1), i=1, nrepeats)]
    signal(:,2) = signal(:,2) - [([(0.001_p, j=1,nx-3)]*(i-1), i=1, nrepeats)]
    ! add glitches
    signal([1,4,17],1) = 100._p
    signal([3,8,9 ],2) = 10._p

    call deglitch_l2b(pmatrix(:,:,1:2), nx, ny, signal(:,1:2), mask(:,1:2), 5._p, .true.)
    if (count(mask) /= 6 .or. any(.not. mask([1,4,17],1)) .or. any(.not. mask([3,8,9],2))) stop 'FAILED: deglitch_l2b 1'

    deallocate (pmatrix)


    !-------------------------------------
    ! test 2 detectors overlap map pixels
    !-------------------------------------
    !    .-----------------------.
    !    | 0 | 1 | 2 | 3 | 4 | 5 |
    !    .-----------------------.
    ! D1:   xxx ooo oox
    !       xxx ooo ooo
    !       ooo oox oox
    !       ooo ooo ooo
    !       ooo ooo ooo
    !       oox xxx ooo
    !
    ! D2:       oox ooo xxx
    !           oox ooo ooo
    !           ooo xxx xxx
    !           ooo ooo ooo
    !           ooo ooo ooo
    !           xxx oox ooo
    !
    ! xxx : glitch
    ! oox : masked by contamination
    ! ooo : no glitch

    allocate (pmatrix(npixels_per_sample, ntimes, ndetectors))
    mask = .false.
    call get_pmatrix(pmatrix(:,:,1:1), nx, ny, xc+0.5_p, yc, 1._p)
    call get_pmatrix(pmatrix(:,:,2:2), nx, ny, xc+1.5_p, yc, 1._p)

    ! read map
    call pmatrix_direct(pmatrix(:,:,1:2), map, signal(:,1:2))
    ! add noise
    signal(:,1) = signal(:,1) + [([(0.001_p, j=1,nx-3)]*(i-1), i=1, nrepeats)]
    signal(:,2) = signal(:,2) - [([(0.001_p, j=1,nx-3)]*(i-1), i=1, nrepeats)]
    ! add glitches
    signal([1,4,17  ],1) = 100._p
    signal([3,8,9,16],2) = 100._p

    call deglitch_l2b(pmatrix(:,:,1:2), nx, ny, signal(:,1:2), mask(:,1:2), 5._p, .true.)
    if (count(mask) /= 14 .or. any(.not. mask([1,3,4,8,9,16,17],1)) .or. &
        any(.not. mask([1,3,4,8,9,16,17],2))) stop 'FAILED: deglitch_l2b 2'


    !-------------------------------------
    ! test 3: 3 detectors overlap map pixels
    !-------------------------------------
    !    .-----------------------.
    !    | 0 | 1 | 2 | 3 | 4 | 5 |
    !    .-----------------------.
    ! D1:   xxx ooo oox
    !       xxx ooo ooo
    !       ooo oox oox
    !       ooo ooo ooo
    !       ooo ooo ooo
    !       oox xxx ooo
    !
    ! D2:       oox ooo xxx
    !           oox ooo ooo
    !           ooo xxx xxx
    !           oox ooo ooo
    !           ooo ooo oox
    !           xxx oox ooo
    !
    ! D3:           ooo ooo oox
    !               ooo ooo ooo
    !               ooo oox oox
    !               xxx ooo ooo
    !               ooo ooo xxx
    !               oox ooo ooo
    !
    ! xxx : glitch
    ! oox : masked by contamination
    ! ooo : no glitch


    mask = .false.
    call get_pmatrix(pmatrix(:,:,1:1), nx, ny, xc+0.5_p, yc, 1._p)
    call get_pmatrix(pmatrix(:,:,2:2), nx, ny, xc+1.5_p, yc, 1._p)
    call get_pmatrix(pmatrix(:,:,3:3), nx, ny, xc+2.5_p, yc, 1._p)

    ! read map
    call pmatrix_direct(pmatrix, map, signal)
    ! add noise
    signal(:,1) = signal(:,1) + [([(0.001_p, j=1,nx-3)]*(i-1), i=1, nrepeats)]
    signal(:,2) = signal(:,2) - [([(0.001_p, j=1,nx-3)]*(i-1), i=1, nrepeats)]
    signal(:,3) = signal(:,3) + [([(0.004_p, j=1,nx-3)]*(i-1), i=1, nrepeats)]
    ! add glitches
    signal([1,4,17   ],1) = 100._p
    signal([3,8,9,16 ],2) = 100._p
    signal([10,15],    3) = 100._p

    call deglitch_l2b(pmatrix, nx, ny, signal, mask, 20._p, .true.)
    if (count(mask) /= 22 .or. any(.not. mask([1,3,4,8,9,16,17],1)) .or. &
        any(.not. mask([1,3,4,8,9,10,15,16,17],2)) .or. &
        any(.not. mask([3,8,9,10,15,16],3))) stop 'FAILED: deglitch_l2b 3'

    stop 'OK.'

    do idetector = 1, ndetectors
        write (*,'(a,i0,a)') 'Detector ', idetector, ': '
        do itime=1, ntimes
            write (*,'(a,i0,a,$)') '   t', itime, ':'
            do ipixel = 1, npixels_per_sample
               write (*,'(x,i0,":",f8.5,$)') pmatrix(ipixel,itime,idetector)%pixel, pmatrix(ipixel,itime,idetector)%weight
            end do
            if (mask(itime,idetector)) then
                write (*,'(a,$)') ' masked'
            end if
            write (*,*)
        end do
    end do


contains


    subroutine get_pmatrix(pmatrix, nx, ny, xc, yc, detectorsize)
        type(pointingelement), intent(inout) :: pmatrix(:,:,:)
        integer, intent(in)                  :: nx, ny
        real(kind=p), intent(in)             :: xc(:), yc(:)
        real(kind=p), intent(in)             :: detectorsize
        real(kind=p) :: h, xy(2,4)
        integer      :: roi(2,2,1)
        integer      :: ntimes, npixels_per_sample, itime
        logical      :: out

        npixels_per_sample = 0
        ntimes = size(xc)
        if (size(yc) /= ntimes) stop 'Error: x and y do not have the same size.'

        out = .false.
        h = detectorsize / 2._p
        do itime=1, ntimes
            xy(1,:) = [xc(itime)-h, xc(itime)+h, xc(itime)+h, xc(itime)-h]
            xy(2,:) = [yc(itime)-h, yc(itime)-h, yc(itime)+h, yc(itime)+h]
!IFORT BUG: print *, shape([yc(itime)-h, yc(itime)-h, yc(itime)+h, yc(itime)+h])
            xy(2,1) = yc(itime)-h
            xy(2,2) = yc(itime)-h
            xy(2,3) = yc(itime)+h
            xy(2,4) = yc(itime)+h

            roi = xy2roi(xy,4)
            call roi2pmatrix(roi, 4, xy, nx, ny, itime, npixels_per_sample, out, pmatrix)
        end do

        if (npixels_per_sample > size(pmatrix,1)) then
            write(ERROR_UNIT,'(a,i0,a)') 'Error: Please update npixels_per_sample to ', npixels_per_sample, '.'
        else if (npixels_per_sample < size(pmatrix,1)) then
            write(OUTPUT_UNIT,'(a,i0,a)') 'Warning: You may update npixels_per_sample to ', npixels_per_sample, '.'
        end if

        if (out) then
            write (OUTPUT_UNIT,'(a)') 'Warning: Some detectors fall outside the map.'
        end if

    end subroutine get_pmatrix


end program test_deglitching

