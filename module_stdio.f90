! interface to C 'macros' stdin, stdout, stderr
module module_stdio
    use, intrinsic :: ISO_C_BINDING
    implicit none

    ! macros from unistd.h
    integer, parameter :: STDIN_FILENO  = 0 ! Standard input.
    integer, parameter :: STDOUT_FILENO = 1 ! Standard output.
    integer, parameter :: STDERR_FILENO = 2 ! Standard error output.

    type(C_PTR) :: stdin
    type(C_PTR) :: stdout
    type(C_PTR) :: stderr

    public :: init_stdio

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
