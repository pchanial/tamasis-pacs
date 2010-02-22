module module_sort

    implicit none

    integer, pointer, private :: array_int(:)
    real*8,  pointer, private :: array_double(:)

    interface qsorti
        module procedure qsorti_int, qsorti_double
    end interface qsorti

    public :: qsortid

contains

! From HDK@psuvm.psu.edu Thu Dec  8 15:27:16 MST 1994
!
! The following was converted from Algol recursive to Fortran iterative
! by a colleague at Penn State (a long time ago - Fortran 66, please
! excuse the GoTo's). The following code also corrects a bug in the
! Quicksort algorithm published in the ACM (see Algorithm 402, CACM,
! Sept. 1970, pp 563-567; also you younger folks who weren't born at
! that time might find interesting the history of the Quicksort
! algorithm beginning with the original published in CACM, July 1961,
! pp 321-322, Algorithm 64). Note that the following algorithm sorts
! integer data; actual data is not moved but sort is affected by sorting
! a companion index array (see leading comments). The data type being
! sorted can be changed by changing one line; see comments after
! declarations and subsequent one regarding comparisons(Fortran
! 77 takes care of character comparisons of course, so that comment
! is merely historical from the days when we had to write character
! compare subprograms, usually in assembler language for a specific
! mainframe platform at that time). But the following algorithm is
! good, still one of the best available.


      SUBROUTINE QSORTID (A, ORD)
      implicit none
      integer :: i
      integer :: ip
      integer :: iq
      integer :: ix
      integer :: iz
      integer :: l
      integer :: l1
      integer :: ndeep
      integer :: p
      integer :: q
      integer :: u
      integer :: u1
      integer :: yp
!
      real*8, intent(in) :: A(:)
      integer, intent(out) :: ORD(size(A))
      integer :: POPLST(2,1000)
      integer :: N
      real*8 X,XX,Z,ZZ,Y
!
      n = size(A)
!     TO SORT DIFFERENT INPUT TYPES, CHANGE THE FOLLOWING
!     SPECIFICATION STATEMENTS; FOR EXAMPLE, FOR FORTRAN CHARACTER
!     USE THE FOLLOWING:  CHARACTER *(*) A(N)
!
      NDEEP=0
      U1=N
      L1=1
      DO 1  I=1,N
    1 ORD(I)=I
    2 IF (U1.LE.L1) RETURN
!
    3 L=L1
      U=U1
!
! PART
!
    4 P=L
      Q=U
!     FOR CHARACTER SORTS, THE FOLLOWING 3 STATEMENTS WOULD BECOME
!     X = ORD(P)
!     Z = ORD(Q)
!     IF (A(X) .LE. A(Z)) GO TO 2
!
!     WHERE "CLE" IS A LOGICAL FUNCTION WHICH RETURNS "TRUE" IF THE
!     FIRST ARGUMENT IS LESS THAN OR EQUAL TO THE SECOND, BASED ON "LEN"
!     CHARACTERS.
!
      X=A(ORD(P))
      Z=A(ORD(Q))
      IF (X.LE.Z) GO TO 5
      Y=X
      X=Z
      Z=Y
      YP=ORD(P)
      ORD(P)=ORD(Q)
      ORD(Q)=YP
    5 IF (U-L.LE.1) GO TO 15
      XX=X
      IX=P
      ZZ=Z
      IZ=Q
!
! LEFT
!
    6 P=P+1
      IF (P.GE.Q) GO TO 7
      X=A(ORD(P))
      IF (X.GE.XX) GO TO 8
      GO TO 6
    7 P=Q-1
      GO TO 13
!
! RIGHT
!
    8 Q=Q-1
      IF (Q.LE.P) GO TO 9
      Z=A(ORD(Q))
      IF (Z.LE.ZZ) GO TO 10
      GO TO 8
    9 Q=P
      P=P-1
      Z=X
      X=A(ORD(P))
!
! DIST
!
   10 IF (X.LE.Z) GO TO 11
      Y=X
      X=Z
      Z=Y
      IP=ORD(P)
      ORD(P)=ORD(Q)
      ORD(Q)=IP
   11 IF (X.LE.XX) GO TO 12
      XX=X
      IX=P
   12 IF (Z.GE.ZZ) GO TO 6
      ZZ=Z
      IZ=Q
      GO TO 6
!
! OUT
!
   13 CONTINUE
      IF (.NOT.(P.NE.IX.AND.X.NE.XX)) GO TO 14
      IP=ORD(P)
      ORD(P)=ORD(IX)
      ORD(IX)=IP
   14 CONTINUE
      IF (.NOT.(Q.NE.IZ.AND.Z.NE.ZZ)) GO TO 15
      IQ=ORD(Q)
      ORD(Q)=ORD(IZ)
      ORD(IZ)=IQ
   15 CONTINUE
      IF (U-Q.LE.P-L) GO TO 16
      L1=L
      U1=P-1
      L=Q+1
      GO TO 17
   16 U1=U
      L1=Q+1
      U=P-1
   17 CONTINUE
      IF (U1.LE.L1) GO TO 18
!
! START RECURSIVE CALL
!
      NDEEP=NDEEP+1
      POPLST(1,NDEEP)=U
      POPLST(2,NDEEP)=L
      GO TO 3
   18 IF (U.GT.L) GO TO 4
!
! POP BACK UP IN THE RECURSION LIST
!
      IF (NDEEP.EQ.0) GO TO 2
      U=POPLST(1,NDEEP)
      L=POPLST(2,NDEEP)
      NDEEP=NDEEP-1
      GO TO 18
!
! END SORT
! END QSORT
!
      end subroutine QSORTID

      SUBROUTINE QSORTGI (nitems, compare, ORD)
      implicit none
      integer, intent(in) :: nitems
      integer, intent(out) :: ORD(nitems)
      interface
          integer function compare(first, second)
              integer, intent(in) :: first, second
          end function compare
      end interface

      integer :: i
      integer :: ip
      integer :: iq
      integer :: ix
      integer :: iz
      integer :: l
      integer :: l1
      integer :: ndeep
      integer :: p
      integer :: q
      integer :: u
      integer :: u1
      integer :: yp
!
      integer :: POPLST(2,1000)
      integer X,XX,Z,ZZ,Y
!
!     TO SORT DIFFERENT INPUT TYPES, CHANGE THE FOLLOWING
!     SPECIFICATION STATEMENTS; FOR EXAMPLE, FOR FORTRAN CHARACTER
!     USE THE FOLLOWING:  CHARACTER *(*) A(N)
!
      NDEEP=0
      U1=NITEMS
      L1=1
      DO I=1,NITEMS
          ORD(I)=I
      END DO

    2 IF (U1.LE.L1) RETURN
!
    3 L=L1
      U=U1
!
! PART
!
    4 P=L
      Q=U
!     FOR CHARACTER SORTS, THE FOLLOWING 3 STATEMENTS WOULD BECOME
!     X = ORD(P)
!     Z = ORD(Q)
!     IF (A(X) .LE. A(Z)) GO TO 2
!
!     WHERE "CLE" IS A LOGICAL FUNCTION WHICH RETURNS "TRUE" IF THE
!     FIRST ARGUMENT IS LESS THAN OR EQUAL TO THE SECOND, BASED ON "LEN"
!     CHARACTERS.
!
      X=ORD(P)
      Z=ORD(Q)
      IF (compare(X, Z) >= 0) GO TO 5
      Y = X
      X = Z
      Z = Y
      YP=ORD(P)
      ORD(P)=ORD(Q)
      ORD(Q)=YP
    5 IF (U - L <= 1) GO TO 15
      XX=X
      IX=P
      ZZ=Z
      IZ=Q
!
! LEFT
!
    6 P=P+1
      IF (P.GE.Q) GO TO 7
      X=ORD(P)
      IF (compare(X, XX) <= 0) GO TO 8
      GO TO 6
    7 P=Q-1
      GO TO 13
!
! RIGHT
!
    8 Q=Q-1
      IF (Q.LE.P) GO TO 9
      Z=ORD(Q)
      IF (compare(Z, ZZ) >= 0) GO TO 10
      GO TO 8
    9 Q=P
      P=P-1
      Z=X
      X=ORD(P)
!
! DIST
!
   10 IF (compare(X, Z) >= 0) GO TO 11
      Y = X
      X = Z
      Z = Y
      IP=ORD(P)
      ORD(P)=ORD(Q)
      ORD(Q)=IP
   11 IF (compare(X, XX) >= 0) GO TO 12
      XX=X
      IX=P
   12 IF (compare(Z, ZZ) <=0) GO TO 6
      ZZ=Z
      IZ=Q
      GO TO 6
!
! OUT
!
   13 CONTINUE
      IF (.NOT.(P /= IX .AND. compare(X,XX) /= 0)) GO TO 14
      IP=ORD(P)
      ORD(P)=ORD(IX)
      ORD(IX)=IP
   14 CONTINUE
      IF (.NOT.(Q /= IZ .AND. compare(Z,ZZ) /= 0)) GO TO 15
      IQ=ORD(Q)
      ORD(Q)=ORD(IZ)
      ORD(IZ)=IQ
   15 CONTINUE
      IF (U-Q <= P-L) GO TO 16
      L1=L
      U1=P-1
      L=Q+1
      GO TO 17
   16 U1=U
      L1=Q+1
      U=P-1
   17 CONTINUE
      IF (U1 <= L1) GO TO 18
!
! START RECURSIVE CALL
!
      NDEEP=NDEEP+1
      POPLST(1,NDEEP)=U
      POPLST(2,NDEEP)=L
      GO TO 3
   18 IF (U.GT.L) GO TO 4
!
! POP BACK UP IN THE RECURSION LIST
!
      IF (NDEEP.EQ.0) GO TO 2
      U=POPLST(1,NDEEP)
      L=POPLST(2,NDEEP)
      NDEEP=NDEEP-1
      GO TO 18
!
! END SORT
! END QSORT
!
    end subroutine qsortgi


    !---------------------------------------------------------------------------
    ! INTEGER
    !---------------------------------------------------------------------------

    subroutine qsorti_int(array, index)

        integer, intent(in), target :: array(:)
        integer, intent(out) :: index(size(array))

        array_int => array
        call qsortgi(size(array_int), compare_int, index)
        array_int => null()

    end subroutine qsorti_int


    !---------------------------------------------------------------------------


    function compare_int(first, second)

        integer, intent(in) :: first, second
        integer             :: compare_int

        if (array_int(first) < array_int(second)) then
           compare_int = 1
        else if (array_int(first) /= array_int(second)) then
           compare_int = -1
        else
           compare_int = 0
        end if

    end function compare_int


    !---------------------------------------------------------------------------
    ! DOUBLE
    !---------------------------------------------------------------------------

    subroutine qsorti_double(array, index)

        real*8, intent(in), target :: array(:)
        integer, intent(out) :: index(size(array))

        array_double => array
        call qsortgi(size(array_double), compare_double, index)
        array_double => null()

    end subroutine qsorti_double


    !---------------------------------------------------------------------------


    function compare_double(first, second)

        integer, intent(in) :: first, second
        integer             :: compare_double

        if (array_double(first) < array_double(second)) then
           compare_double = 1
        else if (array_double(first) /= array_double(second)) then
           compare_double = -1
        else
           compare_double = 0
        end if

    end function compare_double


end module module_sort
