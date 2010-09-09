! interface to WCSLIB 4.4
module module_wcslibc

    use, intrinsic :: ISO_C_BINDING
    implicit none

!-------------------------------------------------------------------------------
! cel.h
!-------------------------------------------------------------------------------

    !struct celprm {
    !  /* Initialization flag (see the prologue above).                          */
    !  /*------------------------------------------------------------------------*/
    !  int    flag;                  /* Set to zero to force initialization.     */
    !
    !  /* Parameters to be provided (see the prologue above).                    */
    !  /*------------------------------------------------------------------------*/
    !  int    offset;                /* Force (x,y) = (0,0) at (phi_0,theta_0).  */
    !  double phi0, theta0;          /* Native coordinates of fiducial point.    */
    !  double ref[4];                /* Celestial coordinates of fiducial        */
    !                                /* point and native coordinates of          */
    !                                /* celestial pole.                          */
    !
    !  struct prjprm prj;            /* Projection parameters (see prj.h).       */
    !
    !  /* Information derived from the parameters supplied.                      */
    !  /*------------------------------------------------------------------------*/
    !  double euler[5];              /* Euler angles and functions thereof.      */
    !  int    latpreq;               /* LATPOLEa requirement.                    */
    !  int    isolat;                /* True if |latitude| is preserved.         */
    !};

    !int celini(struct celprm *cel);
    !int celprt(const struct celprm *cel);
    !int celset(struct celprm *cel);
    !int celx2s(struct celprm *cel, int nx, int ny, int sxy, int sll,
    !           const double x[], const double y[],
    !           double phi[], double theta[], double lng[], double lat[],
    !           int stat[]);
    !int cels2x(struct celprm *cel, int nlng, int nlat, int sll, int sxy,
    !           const double lng[], const double lat[],
    !           double phi[], double theta[], double x[], double y[],
    !           int stat[]);

!-------------------------------------------------------------------------------
! prj.h
!-------------------------------------------------------------------------------

    !struct prjprm {
    !  /* Initialization flag (see the prologue above).                          */
    !  /*------------------------------------------------------------------------*/
    !  int   flag;                   /* Set to zero to force initialization.     */
    !
    !  /* Parameters to be provided (see the prologue above).                    */
    !  /*------------------------------------------------------------------------*/
    !  char   code[4];               /* Three-letter projection code.            */
    !  double r0;                    /* Radius of the generating sphere.         */
    !  double pv[PVN];               /* Projection parameters.                   */
    !  double phi0, theta0;          /* Fiducial native coordinates.             */
    !  int   bounds;                 /* Enable strict bounds checking.           */
    !
    !  /* Information derived from the parameters supplied.                      */
    !  /*------------------------------------------------------------------------*/
    !  char   name[40];              /* Projection name.                         */
    !  int   category;               /* Projection category.                     */
    !  int   pvrange;                /* Range of projection parameter indices.   */
    !  int   simplezen;              /* Is it a simple zenithal projection?      */
    !  int   equiareal;              /* Is it an equal area projection?          */
    !  int   conformal;              /* Is it a conformal projection?            */
    !  int   global;                 /* Can it map the whole sphere?             */
    !  int   divergent;              /* Does the projection diverge in latitude? */
    !  double x0, y0;                /* Fiducial offsets.                        */
    !
    !  double w[10];                 /* Intermediate values.                     */
    !  int   n;                      /* Intermediate value.                      */
    !  int   padding;                /* (Dummy inserted for alignment purposes.) */
    !
    !  int (*prjx2s)(PRJX2S_ARGS);   /* Pointers to the spherical projection and */
    !  int (*prjs2x)(PRJS2X_ARGS);   /* deprojection functions.                  */
    !};

    !int prjini(struct prjprm *prj);
    !int prjprt(const struct prjprm *prj);
    !int prjset(struct prjprm *prj);
    !int prjx2s(PRJX2S_ARGS);
    !int prjs2x(PRJS2X_ARGS);
    !... all projection routines

!-------------------------------------------------------------------------------
! sph.h
!-------------------------------------------------------------------------------

    !int sphx2s(const double eul[5], int nphi, int ntheta, int spt, int sxy,
    !           const double phi[], const double theta[],
    !           double lng[], double lat[]);
    !int sphs2x(const double eul[5], int nlng, int nlat, int sll , int spt,
    !           const double lng[], const double lat[],
    !           double phi[], double theta[]);
    !int sphdpa(int nfield, double lng0, double lat0,
    !           const double lng[], const double lat[],
    !           double dist[], double pa[]);

!-------------------------------------------------------------------------------
! wcs.h
!-------------------------------------------------------------------------------

    integer(kind=C_INT), parameter :: WCSSUB_LONGITUDE = Z'1001'
    integer(kind=C_INT), parameter :: WCSSUB_LATITUDE  = Z'1002'
    integer(kind=C_INT), parameter :: WCSSUB_CUBEFACE  = Z'1004'
    integer(kind=C_INT), parameter :: WCSSUB_CELESTIAL = Z'1007'
    integer(kind=C_INT), parameter :: WCSSUB_SPECTRAL  = Z'1008'
    integer(kind=C_INT), parameter :: WCSSUB_STOKES    = Z'1010'

    ! pscard Struct Reference
    !• int i
    !• int m
    !• char value [72]

    type, BIND(C) :: pscard
        integer(kind=C_INT)           :: i, m
        character(len=72,kind=C_CHAR) :: value
    end type pscard

    ! pvcard Struct Reference
    !• int i
    !• int m
    !• double value

    type, BIND(C) :: pvcard
        integer(kind=C_INT) :: i, m
        real(kind=C_DOUBLE) :: value
    end type pvcard

    !struct wtbarr {
    !  int  i;                       /* Image axis number.                       */
    !  int  m;                       /* Array axis number for index vectors.     */
    !  int  kind;                    /* wcstab array type.                       */
    !  char extnam[72];              /* EXTNAME of binary table extension.       */
    !  int  extver;                  /* EXTVER  of binary table extension.       */
    !  int  extlev;                  /* EXTLEV  of binary table extension.       */
    !  char ttype[72];               /* TTYPEn of column containing the array.   */
    !  long row;                     /* Table row number.                        */
    !  int  ndim;                    /* Expected wcstab array dimensionality.    */
    !  int  *dimlen;                 /* Where to write the array axis lengths.   */
    !  double **arrayp;              /* Where to write the address of the array  */
    !                                /* allocated to store the wcstab array.     */
    !};

    type, bind(C) :: wtbarr
        integer(kind=C_INT)           :: i, m, kind
        character(len=72,kind=C_CHAR) :: extnam
        integer(kind=C_INT)           :: extver, extlev
        character(len=72,kind=C_CHAR) :: ttype
        integer(kind=C_LONG)          :: row
        integer(kind=C_INT)           :: ndim
        type(C_PTR)                   :: dimlen, arrayp
    end type wtbarr

    !struct wcsprm {
    !  /* Initialization flag (see the prologue above).                          */
    !  /*------------------------------------------------------------------------*/
    !  int    flag;                  /* Set to zero to force initialization.     */
    !
    !  /* FITS header keyvalues to be provided (see the prologue above).         */
    !  /*------------------------------------------------------------------------*/
    !  int    naxis;                 /* Number of axes (pixel and coordinate).   */
    !  double *crpix;                /* CRPIXja keyvalues for each pixel axis.   */
    !  double *pc;                   /* PCi_ja  linear transformation matrix.    */
    !  double *cdelt;                /* CDELTia keyvalues for each coord axis.   */
    !  double *crval;                /* CRVALia keyvalues for each coord axis.   */
    !
    !  char   (*cunit)[72];          /* CUNITia keyvalues for each coord axis.   */
    !  char   (*ctype)[72];          /* CTYPEia keyvalues for each coord axis.   */
    !
    !  double lonpole;               /* LONPOLEa keyvalue.                       */
    !  double latpole;               /* LATPOLEa keyvalue.                       */
    !
    !  double restfrq;               /* RESTFRQa keyvalue.                       */
    !  double restwav;               /* RESTWAVa keyvalue.                       */
    !
    !  int    npv;                   /* Number of PVi_ma keywords, and the       */
    !  int    npvmax;                /* number for which space was allocated.    */
    !  struct pvcard *pv;            /* PVi_ma keywords for each i and m.        */
    !
    !  int    nps;                   /* Number of PSi_ma keywords, and the       */
    !  int    npsmax;                /* number for which space was allocated.    */
    !  struct pscard *ps;            /* PSi_ma keywords for each i and m.        */
    !
    !  /* Alternative header keyvalues (see the prologue above).                 */
    !  /*------------------------------------------------------------------------*/
    !  double *cd;                   /* CDi_ja linear transformation matrix.     */
    !  double *crota;                /* CROTAia keyvalues for each coord axis.   */
    !  int    altlin;                /* Alternative representations              */
    !                                /*   Bit 0: PCi_ja  is present,             */
    !                                /*   Bit 1: CDi_ja  is present,             */
    !                                /*   Bit 2: CROTAia is present.             */
    !  int    padding;               /* (Dummy inserted for alignment purposes.) */
    !
    !  /* Auxiliary coordinate system information, not used by WCSLIB.           */
    !  char   alt[4];
    !  int    colnum;
    !  int    *colax;
    !
    !  char   (*cname)[72];
    !  double *crder;
    !  double *csyer;
    !  char   dateavg[72];
    !  char   dateobs[72];
    !  double equinox;
    !  double mjdavg;
    !  double mjdobs;
    !  double obsgeo[3];
    !  char   radesys[72];
    !  char   specsys[72];
    !  char   ssysobs[72];
    !  double velosys;
    !  double zsource;
    !  char   ssyssrc[72];
    !  double velangl;
    !  char   wcsname[72];
    !
    !  /* Coordinate lookup tables (see the prologue above).                     */
    !  /*------------------------------------------------------------------------*/
    !  int    ntab;                  /* Number of separate tables.               */
    !  int    nwtb;                  /* Number of wtbarr structs.                */
    !  struct tabprm *tab;           /* Tabular transformation parameters.       */
    !  struct wtbarr *wtb;           /* Array of wtbarr structs.                 */
    !
    !  /* Information derived from the FITS header keyvalues by wcsset().        */
    !  /*------------------------------------------------------------------------*/
    !  int    *types;                /* Coordinate type codes for each axis.     */
    !  char   lngtyp[8], lattyp[8];  /* Celestial axis types, e.g. RA, DEC.      */
    !  int    lng, lat, spec;        /* Longitude, latitude and spectral axis    */
    !                                /* indices (0-relative).                    */
    !  int    cubeface;              /* True if there is a CUBEFACE axis.        */
    !
    !  struct linprm lin;            /* Linear    transformation parameters.     */
    !  struct celprm cel;            /* Celestial transformation parameters.     */
    !  struct spcprm spc;            /* Spectral  transformation parameters.     */
    !
    !  int    m_flag, m_naxis;       /* The remainder are for memory management. */
    !  double *m_crpix, *m_pc, *m_cdelt, *m_crval;
    !  char  (*m_cunit)[72], (*m_ctype)[72];
    !  struct pvcard *m_pv;
    !  struct pscard *m_ps;
    !  double *m_cd, *m_crota;
    !  int    *m_colax;
    !  char  (*m_cname)[72];
    !  double *m_crder, *m_csyer;
    !  struct tabprm *m_tab;
    !  struct wtbarr *m_wtb;
    !};

    type, BIND(C) :: wcsprm
        integer(kind=C_INT)           :: flag, naxis
        real(kind=C_FLOAT)            :: r
        type(C_PTR)                   :: crpix, pc, cdelt, crval
        type(C_PTR)                   :: cunit, ctype
        real(kind=C_DOUBLE)           :: lonpole, latpole, restfrq, restwav
        integer(C_INT)                :: npv, npvmax
        type(C_PTR)                   :: pv
        integer(C_INT)                :: nps, npsmax
        type(C_PTR)                   :: ps, cd, crota
        integer(C_INT)                :: altlin, padding
        character(len=4,kind=C_CHAR)  :: alt
        integer(C_INT)                :: colnum
        type(C_PTR)                   :: colax, cname, crder, csyer
        character(len=72,kind=C_CHAR) :: dateavg, dateobs
        real(kind=C_DOUBLE)           :: equinox, mjdavg, mjdobs, obsgeo(3)
        character(len=72,kind=C_CHAR) :: radesys, specsys, ssysobs
        real(kind=C_DOUBLE)           :: velosys, zsource
        character(len=72,kind=C_CHAR) :: ssyssrc
        real(kind=C_DOUBLE)           :: velangl
        character(len=72,kind=C_CHAR) :: wcsname
        integer(C_INT)                :: ntab, nwtb
        type(C_PTR)                   :: tab, wtb, types
        character(len=8,kind=C_CHAR)  :: lngtyp, lattyp
        integer(C_INT)                :: lng, lat, spec, cubeface
        type(C_PTR)                   :: lin, cel, spc
        integer(C_INT)                :: m_flag, m_naxis
        type(C_PTR)                   :: m_crpix, m_pc, m_cdelt, m_crval, m_cunit, m_ctype, m_pv, m_ps, m_cd, m_crota, m_colax
        type(C_PTR)                   :: m_cname, m_crder, m_csyer, m_tab, m_wtb
    end type wcsprm

    interface

       !int wcsnpv(int n);
       function wcsnpv(n) bind(C)
           import C_INT
           integer(kind=C_INT), value :: n
           integer(kind=C_INT)        :: wcsnpv
       end function wcsnpv

       !int wcsnps(int n);
       function wcsnps(n) bind(C)
           import C_INT
           integer(kind=C_INT), value :: n
           integer(kind=C_INT)        :: wcsnps
       end function wcsnps

       !int wcsini(int alloc, int naxis, struct wcsprm *wcs);
       function wcsini(alloc, naxis, wcs) bind(C)
           import C_INT, C_PTR
           integer(kind=C_INT), value :: alloc, naxis
           type(C_PTR), value         :: wcs
           integer(kind=C_INT)        :: wcsini
       end function wcsini

       !int wcssub(int alloc, const struct wcsprm *wcssrc, int *nsub, int axes[], struct wcsprm *wcsdst);
       function wcssub(alloc, wcssrc, nsub, axes, wcsdst) bind(C)
           import C_INT, C_PTR
           integer(kind=C_INT), value :: alloc
           type(C_PTR), value         :: wcssrc
           integer(kind=C_INT)        :: nsub
           integer(kind=C_INT)        :: axes(*)
           type(C_PTR), value         :: wcsdst
           integer(kind=C_INT)        :: wcssub
       end function wcssub

       !int wcsfree(struct wcsprm *wcs);
       function wcsfree(wcs) bind(C)
           import C_INT, C_PTR
           type(C_PTR), value, intent(in) :: wcs
           integer(kind=C_INT)            :: wcsfree
       end function wcsfree

       !int wcsprt(const struct wcsprm *wcs);
       function wcsprt(wcs) bind(C)
           import C_INT, C_PTR
           type(C_PTR), value, intent(in) :: wcs
           integer(kind=C_INT)            :: wcsprt
       end function wcsprt

       !int wcsset(struct wcsprm *wcs);
       function wcsset(wcs) bind(C)
           import C_INT, C_PTR
           type(C_PTR), value, intent(in) :: wcs
           integer(kind=C_INT)            :: wcsset
       end function wcsset

   end interface
   !int wcsp2s(struct wcsprm *wcs, int ncoord, int nelem, const double pixcrd[],
   !           double imgcrd[], double phi[], double theta[], double world[],
   !           int stat[]);
   !int wcss2p(struct wcsprm *wcs, int ncoord, int nelem, const double world[],
   !           double phi[], double theta[], double imgcrd[], double pixcrd[],
   !           int stat[]);
   !int wcsmix(struct wcsprm *wcs, int mixpix, int mixcel, const double vspan[],
   !           double vstep, int viter, double world[], double phi[],
   !           double theta[], double imgcrd[], double pixcrd[]);
   !int wcssptr(struct wcsprm *wcs, int *i, char ctype[9]);


!-------------------------------------------------------------------------------
! wcsfix.h
!-------------------------------------------------------------------------------

   integer(kind=C_INT), parameter :: CDFIX   = 0
   integer(kind=C_INT), parameter :: DATFIX  = 1
   integer(kind=C_INT), parameter :: UNITFIX = 2
   integer(kind=C_INT), parameter :: CELFIX  = 3
   integer(kind=C_INT), parameter :: SPCFIX  = 4
   integer(kind=C_INT), parameter :: CYLFIX  = 5
   integer(kind=C_INT), parameter :: NWCSFIX = 6

   !extern const char *wcsfix_errmsg[];
   !#define cylfix_errmsg wcsfix_errmsg

   interface

       !int wcsfix(int ctrl, const int naxis[], struct wcsprm *wcs, int stat[]);
       function wcsfix(ctrl, naxis, wcs, stat) bind(C)
           import C_INT, C_PTR
           integer(kind=C_INT), value :: ctrl
           type(C_PTR)                :: naxis
           type(C_PTR), value         :: wcs
           integer(kind=C_INT)        :: stat(*)
           integer(kind=C_INT)        :: wcsfix
       end function wcsfix

   end interface

   !int cdfix(struct wcsprm *wcs);
   !int datfix(struct wcsprm *wcs);
   !int unitfix(int ctrl, struct wcsprm *wcs);
   !int celfix(struct wcsprm *wcs);
   !int spcfix(struct wcsprm *wcs);
   !int cylfix(const int naxis[], struct wcsprm *wcs);

!-------------------------------------------------------------------------------
! wcshdr.h
!-------------------------------------------------------------------------------

    integer(kind=C_INT), parameter :: WCSHDR_none      = Z'00000000'
    integer(kind=C_INT), parameter :: WCSHDR_all       = Z'000FFFFF'
    integer(kind=C_INT), parameter :: WCSHDR_reject    = Z'10000000'

    integer(kind=C_INT), parameter :: WCSHDR_CROTAia   = Z'00000001'
    integer(kind=C_INT), parameter :: WCSHDR_EPOCHa    = Z'00000002'
    integer(kind=C_INT), parameter :: WCSHDR_VELREFa   = Z'00000004'
    integer(kind=C_INT), parameter :: WCSHDR_CD00i00j  = Z'00000008'
    integer(kind=C_INT), parameter :: WCSHDR_PC00i00j  = Z'00000010'
    integer(kind=C_INT), parameter :: WCSHDR_PROJPn    = Z'00000020'
    integer(kind=C_INT), parameter :: WCSHDR_RADECSYS  = Z'00000040'
    integer(kind=C_INT), parameter :: WCSHDR_VSOURCE   = Z'00000080'
    integer(kind=C_INT), parameter :: WCSHDR_DOBSn     = Z'00000100'
    integer(kind=C_INT), parameter :: WCSHDR_LONGKEY   = Z'00000200'
    integer(kind=C_INT), parameter :: WCSHDR_CNAMn     = Z'00000400'
    integer(kind=C_INT), parameter :: WCSHDR_AUXIMG    = Z'00000800'
    integer(kind=C_INT), parameter :: WCSHDR_ALLIMG    = Z'00001000'

    integer(kind=C_INT), parameter :: WCSHDR_IMGHEAD   = Z'00010000'
    integer(kind=C_INT), parameter :: WCSHDR_BIMGARR   = Z'00020000'
    integer(kind=C_INT), parameter :: WCSHDR_PIXLIST   = Z'00040000'

    integer(kind=C_INT), parameter :: WCSHDO_none      = Z'00'
    integer(kind=C_INT), parameter :: WCSHDO_all       = Z'FF'
    integer(kind=C_INT), parameter :: WCSHDO_safe      = Z'0F'
    integer(kind=C_INT), parameter :: WCSHDO_DOBSn     = Z'01'
    integer(kind=C_INT), parameter :: WCSHDO_TPCn_ka   = Z'02'
    integer(kind=C_INT), parameter :: WCSHDO_PVn_ma    = Z'04'
    integer(kind=C_INT), parameter :: WCSHDO_CRPXna    = Z'08'
    integer(kind=C_INT), parameter :: WCSHDO_CNAMna    = Z'10'
    integer(kind=C_INT), parameter :: WCSHDO_WCSNna    = Z'20'

    interface

        !int wcspih (char *header, int nkeyrec, int relax, int ctrl, int *nreject, int *nwcs, struct wcsprm **wcs)
        function wcspih(header, nkeyrec, relax, ctrl, nreject, nwcs, wcs) bind(C)
            import C_CHAR, C_INT, C_PTR
            implicit none
            character(kind=C_CHAR)     :: header(*)
            integer(kind=C_INT), value :: nkeyrec, relax, ctrl
            integer(kind=C_INT)        :: nreject, nwcs
            type(C_PTR), intent(out)   :: wcs
            integer(kind=C_INT)        :: wcspih
        end function wcspih

        !int wcsidx(int nwcs, struct wcsprm **wcs, int alts[27]);
        function wcsidx(nwcs, wcs, alts) bind(C)
            import C_INT, C_PTR
            implicit none
            integer(kind=C_INT), value :: nwcs
            type(C_PTR), intent(out)   :: wcs
            integer(kind=C_INT)        :: alts(27)
            integer(kind=C_INT)        :: wcsidx
        end function wcsidx

       !int wcsvfree(int *nwcs, struct wcsprm **wcs);
       function wcsvfree(nwcs, wcs) bind(C)
            import C_INT, C_PTR
            implicit none
            integer(kind=C_INT) :: nwcs
            type(C_PTR)         :: wcs
            integer(kind=C_INT) :: wcsvfree
       end function wcsvfree

    end interface

    !int wcsbth(char *header, int nkeyrec, int relax, int ctrl, int keysel,
    !           int *colsel, int *nreject, int *nwcs, struct wcsprm **wcs);
    !int wcstab(struct wcsprm *wcs);
    !int wcsbdx(int nwcs, struct wcsprm **wcs, int type, short alts[1000][28]);
    !int wcshdo(int relax, struct wcsprm *wcs, int *nkeyrec, char **header);

end module module_wcslibc
