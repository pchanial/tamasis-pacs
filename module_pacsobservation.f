!-------------------------------------------------------------------------------
!
! This module defines the pacsobservation derived type.
! It contains all the information of a PACS observation relevant to its data
! processing. In particular:
!    - file name of the observation
!    - channel ('b', 'g' or 'r')
!    - mode ('prime', 'parallel', 'transparent')
!    - section to be processed
!
! Author: Pierre Chanial
!
!-------------------------------------------------------------------------------

module module_pacsobservation

    use iso_fortran_env,    only : ERROR_UNIT, OUTPUT_UNIT
    use module_fitstools,   only : ft_open, ft_open_bintable, ft_read_column, ft_read_extension, ft_close
    use module_observation, only : maskarray, observation, pointing, pacsobservationslice, pacsobservation
    use module_precision,   only : dp, p
    use module_string,      only : strinteger, strreal, strsection, strternary
    implicit none
    private

    public :: maskarray
    public :: pacsobservationslice
    public :: pacsobservation

 end module module_pacsobservation
