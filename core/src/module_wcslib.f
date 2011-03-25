! Copyright 2010-2011 Pierre Chanial
! All rights reserved
!
! interface to WCSLIB
module module_wcslib

    use, intrinsic :: ISO_C_BINDING
    implicit none

    include 'cel.inc'
    include 'prj.inc'
    include 'wcs.inc'
    include 'wcsfix.inc'
    include 'wcshdr.inc'

#if WCSLIB_MISSING_EXTERNAL
    external wcss2p, wcspih, wcsvcopy, wcsput, wcscopy, wcsfree, wcsvfree, wcsfix, wcsset
#endif

    ! routines from cel.inc
!!$    interface
!!$
!!$        function celini(cel) bind(C, name='celini_')
!!$            import C_INT, CELLEN
!!$            implicit none
!!$            integer(kind=C_INT) :: cel(CELLEN), celini
!!$        end function celini
!!$
!!$        function celput(cel, what, value, i) bind(C, name='celput_')
!!$            import C_INT, C_PTR, CELLEN
!!$            implicit none
!!$            integer(kind=C_INT) :: cel(CELLEN), what, i, celput
!!$            type(C_PTR), value :: value
!!$        end function celput
!!$
!!$        !int celget_(const int *cel, const int *what, void *value)
!!$        function celget(cel, what, value) bind(C, name='celget_')
!!$            import C_INT, C_PTR, CELLEN
!!$            implicit none
!!$            integer(kind=C_INT) :: cel(CELLEN), what, celget
!!$            type(C_PTR), value :: value
!!$        end function celget
!!$
!!$        !int celprt_(int *cel)
!!$        function celprt(cel) bind(C, name='celprt_')
!!$            import C_INT, CELLEN
!!$            implicit none
!!$            integer(kind=C_INT) :: cel(CELLEN), celprt
!!$        end function celprt
!!$
!!$        !int celset_(int *cel)
!!$        function celset(cel) bind(C, name='celset_')
!!$            import C_INT, CELLEN
!!$            implicit none
!!$            integer(kind=C_INT) :: cel(CELLEN), celset
!!$        end function celset
!!$
!!$        !int celx2s_(int *cel,const int *nx,const int *ny,const int *sxy,const int *sll,const double x[],const double y[],double phi[],double theta[],double lng[],double lat[],int stat[])
!!$        function celx2s(cel, nx, ny, sxy, sll, x, y, phi, theta, lng, lat, stat) bind(C, name='celx2s_')
!!$            import C_DOUBLE, C_INT, CELLEN
!!$            implicit none
!!$            integer(kind=C_INT) :: cel(CELLEN), nx, ny, sxy, sll, stat(*), celx2s
!!$            real(kind=C_DOUBLE) :: x(*), y(*), phi(*), theta(*), lng(*), lat(*)
!!$        end function celx2s
!!$
!!$        !int cels2x_(int *cel,const int *nlng,const int *nlat,const int *sll,const int *sxy,const double x[],const double y[],double phi[],double theta[],double lng[],double lat[],int stat[])
!!$        function cels2x(cel, nlng, nlat, sll, sxy, lng, lat, phi, theta, x, y, stat) bind(C, name='cels2x_')
!!$            import C_DOUBLE, C_INT, CELLEN
!!$            implicit none
!!$            integer(kind=C_INT) :: cel(CELLEN), nlng, nlat, sll, sxy, stat(*), cels2x
!!$            real(kind=C_DOUBLE) :: lng(*), lat(*), phi(*), theta(*), x(*), y(*)
!!$        end function cels2x
!!$
!!$    end interface
!!$
!!$    ! routines from wcs.inc
!!$    interface
!!$        function wcscopy(wcssrc, wcsdst) bind(C, name='wcscopy_')
!!$            import C_INT, WCSLEN
!!$            implicit none
!!$            integer(kind=C_INT), intent(in)    :: wcssrc(WCSLEN)
!!$            integer(kind=C_INT), intent(inout) :: wcsdst(WCSLEN)
!!$            integer(kind=C_INT)                :: wcscopy
!!$        end function wcscopy
!!$
!!$        function wcsfree(wcs) bind(C, name='wcsfree_')
!!$            import C_INT, WCSLEN
!!$            implicit none
!!$            integer(kind=C_INT) :: wcs(WCSLEN), wcsfree
!!$        end function wcsfree
!!$
!!$        function wcsprt(wcs) bind(C, name='wcsprt_')
!!$            import C_INT, WCSLEN
!!$            implicit none
!!$            integer(kind=C_INT) :: wcs(WCSLEN), wcsprt
!!$        end function wcsprt
!!$
!!$        function wcsset(wcs) bind(C, name='wcsset_')
!!$            import C_INT, WCSLEN
!!$            implicit none
!!$            integer(kind=C_INT) :: wcs(WCSLEN), wcsset
!!$        end function wcsset
!!$
!!$        !int wcsget_(const int *wcs, const int *what, void *value)
!!$        !function wcsget(wcs, what, value) bind(C, name='wcsget_')
!!$        !    import C_INT, C_PTR, WCSLEN
!!$        !    integer(kind=C_INT) :: wcs(WCSLEN), what, wcsget
!!$        !    type(C_PTR), value  :: value
!!$        !end function wcsget
!!$
!!$        !int wcsput_(int *wcs,const int *what,const void *value,const int *i,const int *j)
!!$        !function wcsput(wcs, what, value, i, j) bind(C, name='wcsput_')
!!$        !    import C_INT, C_PTR, WCSLEN
!!$        !    integer(kind=C_INT) :: wcs(WCSLEN), what, i, j, wcsput
!!$        !    type(C_PTR), value  :: value
!!$        !end function wcsput
!!$
!!$    end interface
!!$
!!$    interface
!!$
!!$        function wcspih(header, nkeys, relax, ctrl, nreject, nwcs, wcsp) bind(C, name='wcspih_')
!!$            import C_CHAR, C_INT
!!$            implicit none
!!$            character(kind=C_CHAR), intent(in) :: header(*)
!!$            integer(kind=C_INT), intent(in)    :: nkeys, relax, ctrl
!!$            integer(kind=C_INT), intent(out)   :: nreject, nwcs, wcsp
!!$            integer(kind=C_INT)                :: wcspih
!!$        end function wcspih
!!$
!!$         function wcsidx(nwcs, wcsp, alts) bind(C, name='wcsidx_')
!!$            import C_INT
!!$            implicit none
!!$            integer(kind=C_INT) :: nwcs, wcsp, alts(27), wcsidx
!!$        end function wcsidx
!!$
!!$        function wcsvfree(nwcs, wcsp) bind(C, name='wcsvfree_')
!!$             import C_INT
!!$             implicit none
!!$             integer(kind=C_INT) :: nwcs, wcsp, wcsvfree
!!$        end function wcsvfree
!!$
!!$        function wcsvcopy(wcspp, i, wcs) bind(C, name='wcsvcopy_')
!!$             import C_INT, WCSLEN
!!$             implicit none
!!$             integer(kind=C_INT) :: wcspp, i, wcs(WCSLEN), wcsvcopy
!!$        end function wcsvcopy
!!$
!!$    end interface
!!$
!!$    !routines from wcsfix.inc
!!$    interface
!!$        function wcsfix(ctrl, naxis, wcs, stat) bind(C, name='wcsfix_')
!!$            import C_INT, C_PTR, WCSLEN, WCSFIX_NWCS
!!$            implicit none
!!$            integer(kind=C_INT) :: ctrl, wcs(WCSLEN), stat(WCSFIX_NWCS), wcsfix
!!$            type(C_PTR) :: naxis
!!$        end function wcsfix
!!$    end interface
!!$
!!$    ! routines from ./wcshdr.c
!!$    interface
!!$        !int wcshdo(int relax, struct wcsprm *wcs, int *nkeyrec, char **header)
!!$        function wcshdo(relax, wcs, nkeyrec, header) bind(C)
!!$            import C_INT, C_PTR, WCSLEN
!!$            implicit none
!!$            integer(kind=C_INT), intent(in), value :: relax
!!$            integer(kind=C_INT), intent(in)        :: wcs(WCSLEN), nkeyrec
!!$            type(C_PTR), intent(out)               :: header
!!$            integer(kind=C_INT)                    :: wcshdo
!!$        end function wcshdo
!!$    end interface

end module module_wcslib
