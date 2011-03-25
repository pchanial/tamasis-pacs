! Copyright 2010-2011 Pierre Chanial
! All rights reserved
!
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
    use module_fitstools,   only : ft_open, ft_open_bintable, ft_read_column, ft_read_image, ft_close
    use module_observation, only : MaskPolicy, Observation, Pointing, PacsObservationSlice, PacsObservation
    use module_string,      only : strinteger, strreal, strsection, strternary
    implicit none
    private

    public :: MaskPolicy
    public :: PacsObservationSlice
    public :: PacsObservation

 end module module_pacsobservation
