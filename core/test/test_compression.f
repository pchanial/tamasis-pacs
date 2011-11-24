program test_compression

    use iso_fortran_env, only : ERROR_UNIT
    use module_compression
    use module_math,    only : neq_real
    use module_tamasis, only : p
    implicit none

    integer, parameter :: factor = 5
    real(p) :: tod(10)
    real(p) :: compressed(size(tod,1)/factor)
    integer :: i

    tod = [(real( i,p), i=1, size(tod))]

    call downsampling_direct(tod, compressed, factor)
    if (any(neq_real(compressed(:),  [1._p,6._p]))) call failure('downsampling_direct')
    
    call downsampling_transpose(compressed, tod, factor)
    if (any(neq_real(tod(:),  [1._p,0._p,0._p,0._p,0._p,6._p,0._p,0._p,0._p,0._p]))) call failure('downsampling_transpose')

contains

    subroutine failure(errmsg)
        character(len=*), intent(in) :: errmsg
        write (ERROR_UNIT,'(a)'), 'FAILED: ' // errmsg
        stop 1
    end subroutine failure

end program test_compression
