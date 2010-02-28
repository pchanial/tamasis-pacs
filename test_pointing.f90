program test_pointing
 use precision, only : test_real_eq
 use module_pacspointing
 implicit none

 class(pacspointing), allocatable :: ptg
 character(len=*), parameter :: filename = '/home/pchanial/work/pacs/data/transparent/NGC6946/1342184520_blue'
 integer                     :: i, index
 real*8                      :: ra, dec, pa, chop, nan, zero
 real*8, allocatable         :: times(:), ras(:), decs(:), pas(:)
 real*8, parameter           :: timetest(5) = [1632946094.374053d0, 1632946094.399055d0, 1632946094.424049d0, 1632946094.449043d0, &
                                               1632946094.474037d0]
 integer                     :: status


 zero = 0.d0
 nan = zero / zero

 ! check values in files: time, ra, dec, pa
 !      1632946094374053       258.50445       66.964140       252.34927
 !      1632946094399055       258.51125       66.964432       252.35550
 !      1632946094424049       258.51805       66.964723       252.36173
 !      1632946094449043       258.52485       66.965014       252.36796
 !      1632946094474037       258.53164       66.965305       252.37419

 allocate(ptg)
 call ptg%load(filename, status)
 if (status /= 0) stop 'ptg%load: FAILED.'

 index = 2
 do i = 1, 5
     call ptg%get_position(timetest(i), ra, dec, pa, chop, index)
     write (*,*) 'sample: ', i, timetest(i), ra, dec, pa, chop
 end do
 call ptg%destructor()
 deallocate(ptg)

 allocate(times(10))
 allocate(ras  (10))
 allocate(decs (10))
 allocate(pas  (10))

 times = [(i * 0.5d0, i = -2, 7)]

 ! test fast interpolation (evenly spaced sampling)
 allocate(ptg)
 call ptg%load([1.d0, 2.d0], [0.d0, 1.d0], [3.d0, 3.d0], [2.d0, 1.d0], [0.d0, 0.d0], status)
 if (status /= 0) stop 'ptg%load: FAILED.'

 index = 2
 do i = 1, 10
     call ptg%get_position(times(i), ra, dec, pa, chop, index)
     ras (i) = ra
     decs(i) = dec
     pas (i) = pa
 end do
 if (any(.not. test_real_eq(ras,  [nan, nan,-1.0d0,-0.5d0, 0.0d0, 0.5d0, 1.0d0, 1.5d0, 2.0d0, nan], 10) .or.                       &
         .not. test_real_eq(decs, [nan, nan, 3.0d0, 3.0d0, 3.0d0, 3.0d0, 3.0d0, 3.0d0, 3.0d0, nan], 10) .or.                       &
         .not. test_real_eq(pas,  [nan, nan, 3.0d0, 2.5d0, 2.0d0, 1.5d0, 1.0d0, 0.5d0, 0.0d0, nan], 10))) then
     write (*,*) 'Test get_position_ev: FAILED'
 else
     write (*,*) 'Test get_position_ev: OK'
 end if
 call ptg%destructor()
 deallocate(ptg)

 ! test slow interpolation (unevenly spaced sampling)
 allocate(ptg)
 call ptg%load([1.d0, 2.d0, 2.5d0], [0.d0, 1.d0, 1.5d0], [3.d0, 3.d0, 3.d0], [2.d0, 1.d0, 0.5d0], [0.d0, 0.d0, 0.d0], status)
 if (status /= 0) stop 'ptg%load: FAILED.'

 index = 2
 do i = 1, 10
     call ptg%get_position(times(i), ra, dec, pa, chop, index)
     ras (i) = ra
     decs(i) = dec
     pas (i) = pa
 end do
 if (any(.not. test_real_eq(ras,  [nan, nan,-1.0d0,-0.5d0, 0.0d0, 0.5d0, 1.0d0, 1.5d0, 2.0d0, 2.5d0], 10) .or.                     &
         .not. test_real_eq(decs, [nan, nan, 3.0d0, 3.0d0, 3.0d0, 3.0d0, 3.0d0, 3.0d0, 3.0d0, 3.0d0], 10) .or.                     &
         .not. test_real_eq(pas,  [nan, nan, 3.0d0, 2.5d0, 2.0d0, 1.5d0, 1.0d0, 0.5d0, 0.0d0,-0.5d0], 10))) then
     write (*,*) 'Test get_position_gen failed'
 else
     write (*,*) 'Test get_position_gen: OK'
 end if
 call ptg%destructor()
 deallocate(ptg)

 stop "OK."

end program test_pointing
