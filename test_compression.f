program test_compression

    use module_compression
    use module_math,      only : neq_real
    use module_precision, only : p
    implicit none

    integer, parameter :: factor = 5
    real(p) :: tod(10,2)
    real(p) :: compressed(size(tod,1)/factor,size(tod,2))
    integer :: i

    tod(:,1) = [(real( i,p), i=1, size(tod,1))]
    tod(:,2) = [(real(-i,p), i=1, size(tod,1))]

    call downsampling_direct(tod, compressed, factor)
    if (any(neq_real(compressed(:,1),  [1._p,6._p], 15)) .or.                                                                      &
        any(neq_real(compressed(:,2), -[1._p,6._p], 15))) stop 'FAILED: downsampling_direct'
    
    call downsampling_transpose(compressed, tod, factor)
    if (any(neq_real(tod(:,1),  [1._p,0._p,0._p,0._p,0._p,6._p,0._p,0._p,0._p,0._p], 15)) .or.                                     &
        any(neq_real(tod(:,2), -[1._p,0._p,0._p,0._p,0._p,6._p,0._p,0._p,0._p,0._p], 15))) stop 'FAILED: downsampling_transpose'
    
    stop 'OK.'

end program test_compression
