program test_stack

    use iso_fortran_env, only : ERROR_UNIT
    use module_stack
    implicit none

    class(stack_int), allocatable :: stack
    integer, allocatable          :: array(:)

    allocate(stack)

    call stack%push(10)
    if (stack%pop() /= 10) call failure('pop.')
    if (stack%nelements /= 0) call failure('nelements 1')
    call stack%push(10)
    call stack%push(20)
    call stack%push(30)
    call stack%to_array(array)
    if (array(1) /= 30 .or. array(2) /= 20 .or. array(3) /= 10) call failure('push')
    if (stack%nelements /= 3) call failure('nelements 2')
    if (stack%pop() /= 30) call failure('pop 2')
    if (stack%pop() /= 20) call failure('pop 3')
    if (stack%pop() /= 10) call failure('pop 4')
    if (stack%nelements /= 0) call failure('nelements 3')

contains

    subroutine failure(errmsg)
        character(len=*), intent(in) :: errmsg
        write (ERROR_UNIT,'(a)'), 'FAILED: ' // errmsg
        stop 1
    end subroutine failure

end program test_stack
