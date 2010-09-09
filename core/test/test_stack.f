program test_stack

    use module_stack
    implicit none

    class(stack_int), allocatable :: stack
    integer, allocatable          :: array(:)

    allocate(stack)

    call stack%push(10)
    if (stack%pop() /= 10) stop 'FAILED: pop.'
    if (stack%nelements /= 0) stop 'FAILED: nelements 1'
    call stack%push(10)
    call stack%push(20)
    call stack%push(30)
    call stack%to_array(array)
    if (array(1) /= 30 .or. array(2) /= 20 .or. array(3) /= 10) stop 'FAILED: push'
    if (stack%nelements /= 3) stop 'FAILED: nelements 2'
    if (stack%pop() /= 30) stop 'FAILED: pop 2'
    if (stack%pop() /= 20) stop 'FAILED: pop 3'
    if (stack%pop() /= 10) stop 'FAILED: pop 4'
    if (stack%nelements /= 0) stop 'FAILED: nelements 3'

    stop 'OK.'

end program test_stack
