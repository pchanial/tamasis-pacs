! Copyright 2010-2011 Pierre Chanial
! All rights reserved
!
! interface to C 'macros' stdin, stdout, stderr
module module_stdio

    use iso_c_binding
    implicit none
    private

    public :: init_stdio
    public :: stderr
    public :: stdin
    public :: stdout

    ! macros from unistd.h
    integer(C_INT), parameter :: STDIN_FILENO  = 0 ! Standard input.
    integer(C_INT), parameter :: STDOUT_FILENO = 1 ! Standard output.
    integer(C_INT), parameter :: STDERR_FILENO = 2 ! Standard error output.

    type(C_PTR) :: stderr
    type(C_PTR) :: stdin
    type(C_PTR) :: stdout

    interface
        function fdopen(filedes, mode) bind(C)
            import C_CHAR, C_INT, C_PTR
            integer(kind=C_INT), value, intent(in) :: filedes
            character(kind=C_CHAR), intent(in)     :: mode(*)
            type(C_PTR)                            :: fdopen
        end function fdopen
    end interface


contains


    subroutine init_stdio()
        stdin  = fdopen(STDIN_FILENO,  "r" // C_NULL_CHAR)
        stdout = fdopen(STDOUT_FILENO, "w" // C_NULL_CHAR)
        stderr = fdopen(STDERR_FILENO, "w" // C_NULL_CHAR)
    end subroutine init_stdio


end module module_stdio
