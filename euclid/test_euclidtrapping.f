program test_euclidtrapping

    use module_euclidtrapping
    use module_fitstools, only : ft_write_image
    use module_precision, only : p => dp

    integer, parameter :: ntirages = 100000
    integer :: status

    call run(ntirages)

contains

    subroutine run(ntirages)

        integer, intent(in) :: ntirages
        
        integer :: output(ntirages, nlines, ncolumns)
        real(p) :: proba(nlines, ncolumns)

        print *, 'starting with ', shape(output)
        proba = 0.5_p
        call tirage_proba(proba, ntirages, output)
        call ft_write_image('/mnt/herschel1/mapmaking/pchanial/work/euclid/simul_trapping/tirage_proba_ifort.fits', output,        &
             status=status)

    end subroutine run

end program test_euclidtrapping
