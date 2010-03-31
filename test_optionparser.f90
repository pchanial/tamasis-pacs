program test_optionparser

    use module_optionparser
    implicit none

    class(optionparser), allocatable :: parser
    integer                          :: status

    allocate(parser)
    call parser%init('test_optionparser [options] arg1 [arg2]', 1, 2)

    call parser%add_option('recursive','r','Recursive loop',status=status)
    if (status /= 0) stop 'FAILED: parser%add_option -r'

    call parser%add_option('recursivity','R','Recursivity index',.true., status=status)
    if (status /= 0) stop 'FAILED: parser%add_option -R'

    call parser%add_option('xoffset','x','X offset (arcsec)',.true., default="1.5", status=status)
    if (status /= 0) stop 'FAILED: parser%add_option -x'
    
    call parser%add_option('input','i','input file',.true., status=status)
    if (status /= 0) stop 'FAILED: parser%add_option -i'
    
    call parser%add_option('output','o','output file',.true., status=status)
    if (status /= 0) stop 'FAILED: parser%add_option -o'
    
    call parser%add_option('one','1','first guess',action='store true', status=status)
    if (status /= 0) stop 'FAILED: parser%add_option -1'
    
    call parser%add_option('two','2','second guess',action='store false', status=status)
    if (status /= 0) stop 'FAILED: parser%add_option -2'
    
    call parser%parse(status)
    if (status == 0) stop 'FAILED: parser%parse 1'

    call parser%parse(status, '')
    if (status == 0) stop 'FAILED: parser%parse 2'

    call parser%parse(status, 'my\ arg1 "my arg2" --recursive')
    if (status /= 0) stop 'FAILED: parser%parse 3'
    if (parser%get_option('1', status) /= 'False') stop 'FAILED: parser%parse 3b'
    if (parser%get_option('two', status) /= 'True') stop 'FAILED: parser%parse 3c'
    if (parser%get_argument_count() /= 2) stop 'FAILED: parser%parse 3d'
    call parser%reset()

    call parser%parse(status, 'arg -R 100')
    if (status /= 0) stop 'FAILED: parser%parse 4'
    if (parser%get_option_as_integer('R', status) /= 100) stop 'FAILED: parser%parse 4b'
    call parser%reset()

    call parser%parse(status, 'arg -i my\ input -o "my o\ utput"')
    if (status /= 0) stop 'FAILED: parser%parse 5'
    if (parser%get_option('i', status) /= 'my input') stop 'FAILED: parser%parse 5b'
    if (parser%get_option('o', status) /= 'my o\ utput') stop 'FAILED: parser%parse 5c'
    call parser%reset()

    call parser%parse(status, 'arg -Rx')
    if (status == 0) stop 'FAILED: parser%parse 6'

    call parser%parse(status, 'arg -1 --two')
    if (status /= 0) stop 'FAILED: parser%parse 7'
    if (parser%get_option('1', status) /= 'True') stop 'FAILED: parser%parse 7b'
    if (parser%get_option_as_logical('two', status)) stop 'FAILED: parser%parse 7c'

    call parser%reset()


    stop 'OK.'

end program test_optionparser
