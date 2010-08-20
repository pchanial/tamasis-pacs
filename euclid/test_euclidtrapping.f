program test_euclidtrapping

    use module_euclidtrapping
    use module_precision, only : p => dp

    integer, parameter :: ntirages = 100000
    integer, parameter :: nlines = 1024, ncolumns = 1
    real(p) :: proba(nlines, ncolumns)
    integer :: output(nlines, ncolumns, ntirages)
    
    !call init
    proba = 0.5_p
    call tirage_proba(proba, ntirages, output)

end program test_euclidtrapping
