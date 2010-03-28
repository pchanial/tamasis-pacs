program test_string

    use string
    implicit none

    if (strsection(0,0)     /= ':'      ) stop 'FAILED strsection 1'
    if (strsection(1,0)     /= '1:'     ) stop 'FAILED strsection 2'
    if (strsection(-234,0)  /= '-234:'  ) stop 'FAILED strsection 3'
    if (strsection(0,3)     /= ':3'     ) stop 'FAILED strsection 4'
    if (strsection(0,-32)   /= ':-32'   ) stop 'FAILED strsection 5'
    if (strsection(30,43)   /= '30:43'  ) stop 'FAILED strsection 6'
    if (strsection(-30,321) /= '-30:321') stop 'FAILED strsection 7'
    if (strsection(32,2324) /= '32:2324') stop 'FAILED strsection 8'

    if (strternary(.true., 'a', 'bc') /= 'a') stop 'FAILED strternary 1'
    if (strternary(.false., 'a', 'bc') /= 'bc') stop 'FAILED strternary 2'
    if (len(strternary(.true., 'a', 'bc')) /= 1) stop 'FAILED strternary 3'
    if (len(strternary(.false., 'a', 'bc')) /= 2) stop 'FAILED strternary 4'
    stop 'OK.'

end program test_string
