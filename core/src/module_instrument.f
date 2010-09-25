module module_instrument
 
    use module_string,  only : strcompress, strinteger, strlowcase, strsplit
    use module_tamasis, only : p
    implicit none
    private

    public :: read_instrument_configfile
    public :: write_configuration

    ! valid parameters
    character(len=*), parameter :: c_receiver_name        = 'Receiver name'
    character(len=*), parameter :: c_dimension_number     = 'Dimension number'
    character(len=*), parameter :: c_detector_geometry    = 'Detector geometry'
    character(len=*), parameter :: c_detector_number      = 'Detector number'
    character(len=*), parameter :: c_offset_unit          = 'Detector offset unit'
    character(len=*), parameter :: c_detector_description = 'Detector description'
    character(len=*), parameter :: c_level_number         = 'Level number'
    character(len=*), parameter :: c_level_names          = 'Level names'

    ! enumerations
    integer, parameter          :: unknown = 0
    integer, parameter          :: uniform_square = 1
    integer, parameter          :: arcsec = 100
    character(len=*), parameter :: c_unknown = 'unknown'
    character(len=*), parameter :: c_uniform_square = 'uniform square'
    character(len=*), parameter :: c_arcsec = 'arcsec'

    ! misc
    integer, parameter          :: nlevels_max = 10

    ! description of the instrument
    type, public :: receiver
        character(len=80)                            :: name
        integer                                      :: ndimensions
        integer                                      :: geometry
        integer                                      :: ndetectors
        integer                                      :: nvertices
        integer                                      :: offset_unit
        integer, dimension(:,:), allocatable         :: id
        real(p), dimension(:,:,:), allocatable       :: vertex
        character(len=80), dimension(:), allocatable :: flag
        integer                                      :: nlevels
        character(len=80), dimension(nlevels_max)    :: level_name
        character(len=1024)                          :: configfile
    end type receiver


contains


    function read_instrument_configfile(configfile) result(output)

        type(receiver)               :: output
        character(len=*), intent(in) :: configfile

        character(len=1024)          :: line, param, value
        integer                      :: unit, iostat, iline
        integer                      :: idetector
   
        ! initialisation
        output%name = ' '
        output%ndimensions = 2
        output%geometry = unknown
        output%ndetectors = 0
        output%nvertices = 0
        output%offset_unit = arcsec
        output%nlevels = 1
        output%level_name = ' '
        output%configfile = configfile
        
        ! set a unit to read the config file
        unit = 1
        
        open(unit=unit, file=trim(configfile), status='old', action='read', iostat=iostat)
        if (iostat /= 0) then
           write(*,'(a)') "The instrument configuration file '" // trim(configfile) // "' can not be read."
           stop
        end if
        
        iline = 0
        do

            iline = iline + 1
            read(unit,'(a)', iostat=iostat) line
            if (iostat /= 0) exit

            call read_configfile_one_line(line, param, value)
 
            if (param(1:1) == '#' .or. len_trim(line) == 0) cycle

            ! receiver name
            if (strcomp_config(param, c_receiver_name)) then
                output%name = value(1:80)

            ! number of levels
            else if (strcomp_config(param, c_dimension_number)) then
                read(value,*) output%ndimensions

            ! detector geometry
            else if (strcomp_config(param, c_detector_geometry)) then
                if (strcomp_config(value, c_uniform_square)) then
                    output%geometry = uniform_square
                    output%nvertices = 4
                else
                    write(*,'(a)') "Value '" // trim(value) // "' for parameter '" //  c_detector_geometry // "' is unknown. " //  &
                         "It can only be '" // c_uniform_square // "'."
                end if

            ! number of detectors
            else if (strcomp_config(param, c_detector_number)) then
                read (value,*) output%ndetectors

            ! number of levels
            else if (strcomp_config(param, c_level_number)) then
                read (value,'(i2)') output%nlevels

            ! name of levels
            else if (strcomp_config(param, c_level_names)) then
                if (output%nlevels == 0) then
                    write (*,'(a)') "In the configuration file, the parameter '" // c_level_names // "' should be after '" //      &
                          c_level_number // "'."
                    stop
                end if
         
                call strsplit(value, output%level_name, delimiter=',')
                output%level_name = adjustl(output%level_name)

            ! offset unit
            else if (strcomp_config(param, c_offset_unit)) then
                if (strcomp_config(value, c_arcsec)) then
                    output%offset_unit = arcsec
                else
                    write (*,'(a)') "Value '" // trim(value) // "' for detector offset unit is unknown. It can only be '" //       &
                          c_arcsec // "'"
                end if

            ! description of detectors
            else if (strcomp_config(param, c_detector_description)) then
         
                if (output%nvertices == 0) then
                    write (*,'(a)') "In the configuration file, the parameter '" // c_detector_description // "' should be after '"&
                          // c_detector_geometry // "'."
                    stop
                end if
                if (output%ndetectors == 0) then
                    write (*,'(a)') "In the configuration file, the parameter '" // c_detector_description // "' should be after '"&
                          // c_detector_number // "'."
                    stop
                end if

                allocate(output%id(0:output%nlevels,output%ndetectors))
                allocate(output%flag(output%ndetectors))
                allocate(output%vertex(output%ndimensions, output%nvertices, output%ndetectors))

                ! loop over the detectors: each detector has a description
                do idetector = 1, output%ndetectors
                    do
                        iline = iline + 1
                        read (unit,'(a)',iostat=iostat) line
                        if (iostat /= 0) then
                            write (*,'(a)') 'Read error, make sure the number of detector description records matches the number of&
                                & detectors.'
                            stop
                        end if
                        if ((line(1:1) /= '#' .and. len_trim(line) /= 0)) exit
                    end do
            
                    ! read one description line
                    read(line,*, iostat=iostat) output%id(:,idetector), output%vertex(:,:,idetector), output%flag(idetector)
                    if (iostat /= 0) then
                        ! the description has no flag
                        read(line,*, iostat=iostat) output%id(:,idetector), output%vertex(:,:,idetector)
                        if (iostat == 0) then
                            output%flag(idetector) = ' '
                        else
                            ! the description can not be understood
                            write(*,'(a)') 'Detector description improperly formatted:'
                            write(*,'(a)') trim(line)
                            stop
                        end if
                    end if
                    !write(*,*) output%id(0:output%nlevels, idetector), output%vertex(:,:,idetector), output%flag(idetector)
                end do

                ! line is not understood
            else
                write(*,'(a)') 'Line ' // strinteger(iline) // ' is not understood: ' // trim(line)            
            end if

        end do

        close(unit=unit, iostat=iostat)
        if (iostat /= 0) then
             write(*,'(a)') 'The unit of the configuration file cannot be closed'
        end if

    end function read_instrument_configfile


    !-------------------------------------------------------------------------------------------------------------------------------


    subroutine read_configfile_one_line(line, param, value)

        character(len=*), intent(in)  :: line
        character(len=*), intent(out) :: param, value

        integer                       :: idelim

        idelim = scan(line, ':')
        if (idelim == 0) then
            param = line
            value = ' '
        else
            param = line(:idelim-1)
            value = line(idelim+1:)
        end if
        
        param = adjustl(param)
        value = adjustl(value)

      end subroutine read_configfile_one_line

  
      !-----------------------------------------------------------------------------------------------------------------------------


      function strcomp_config(string1, string2)

          logical                      :: strcomp_config
          character(len=*), intent(in) :: string1, string2

          strcomp_config = strlowcase(strcompress(string1)) == strlowcase(strcompress(string2))

      end function strcomp_config


      !-----------------------------------------------------------------------------------------------------------------------------


      subroutine write_configuration(input)

          type(receiver), intent(in) :: input

          character(len=80)          :: string
          integer                    :: i

          select case(input%geometry)
          case (unknown)
             string = c_unknown
          case (uniform_square)
             string = c_uniform_square
          case default
             string = 'Error'
          end select

          write(*,'(a)') 'Configuration file: ' // trim(input%configfile)
          write(*,'(a)') 'Receiver name     : ' // trim(input%name)
          write(*,'(a)') 'Dimension number  : ' // strinteger(input%ndimensions)
          write(*,'(a)') 'Detector geometry : ' // trim(string)
          write(*,'(a)') 'Detector number   : ' // strinteger(input%ndetectors)
          write(*,'(a)') 'Level number      : ' // strinteger(input%nlevels)
          string = ' '
          write(string,  '(a,i2,a)') '(a,', input%nlevels-1, '(a,", "),a)'
          write(*,string)'Level names       : ', (trim(input%level_name(i)), i=1,input%nlevels)

      end subroutine write_configuration


end module module_instrument
