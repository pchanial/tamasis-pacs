program test_string

    use iso_fortran_env, only : ERROR_UNIT
    use module_string
    implicit none

    character(len=5), allocatable :: output(:)

    if (strsection(0,0)     /= ':'      ) call failure('strsection 1')
    if (strsection(1,0)     /= '1:'     ) call failure('strsection 2')
    if (strsection(-234,0)  /= '-234:'  ) call failure('strsection 3')
    if (strsection(0,3)     /= ':3'     ) call failure('strsection 4')
    if (strsection(0,-32)   /= ':-32'   ) call failure('strsection 5')
    if (strsection(30,43)   /= '30:43'  ) call failure('strsection 6')
    if (strsection(-30,321) /= '-30:321') call failure('strsection 7')
    if (strsection(32,2324) /= '32:2324') call failure('strsection 8')

    call strsplit('', ',', output)
    if (size(output) /= 1) call failure('strsplit 1')
    if (output(1) /= '') call failure('strsplit 1b')
    deallocate (output)
    call strsplit(',', ',', output)
    if (size(output) /= 2) call failure('strsplit 2')
    if (output(1) /= '' .or. output(2) /= '') call failure('strsplit 2b')
    deallocate (output)
    call strsplit('a,b,c', ',', output)
    if (size(output) /= 3) call failure('strsplit 3')
    if (output(1) /= 'a' .or. output(2) /= 'b' .or. output(3) /= 'c') call failure('strsplit 3b')
    deallocate (output)
    if (strternary(.true., 'a', 'bc') /= 'a')     call failure('strternary 1')
    if (strternary(.false., 'a', 'bc') /= 'bc')   call failure('strternary 2')
    if (len(strternary(.true., 'a', 'bc')) /= 1)  call failure('strternary 3')
    if (len(strternary(.false., 'a', 'bc')) /= 2) call failure('strternary 4')

contains

    subroutine failure(errmsg)
        character(len=*), intent(in) :: errmsg
        write (ERROR_UNIT,'(a)'), 'FAILED: ' // errmsg
        stop 1
    end subroutine failure

end program test_string
