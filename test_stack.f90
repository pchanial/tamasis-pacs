program test_stack

    use module_stack
    implicit none

    class(stack_int), allocatable :: stack
    integer, allocatable          :: array(:)

    allocate(stack)

    call stack%print()
    call stack%push(10)
    call stack%print()
    write (*,*) stack%pop() == 10
    call stack%push(10)
    call stack%push(20)
    call stack%push(30)
    call stack%to_array(array)
    write (*,*) array(size(array):1:-1)
    call stack%print()
    write (*,*) stack%pop() == 30
    call stack%print()
    write (*,*) stack%pop() == 20
    call stack%print()
    write (*,*) stack%pop() == 10
    call stack%print()

end program test_stack
