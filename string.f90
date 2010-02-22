module string

 implicit none

 private
 public :: strlowcase
 public :: strupcase
 public :: strcompress
 public :: strsplit
 public :: strinteger
 public :: strjoin

 interface strjoin
     module procedure strjoin_trim, strjoin_opt
 end interface strjoin

contains

 ! convert a word to lower case
 elemental function strlowcase(input)

  character(len=*), intent(in) :: input
  character(len=len(input))    :: strlowcase
  integer                      :: i,ic

  strlowcase = input
  do i=1, len(input)
     ic = iachar(input(i:i))
     if (ic >= 65 .and. ic < 90) strlowcase(i:i) = achar(ic+32)
  end do

 end function strlowcase

 ! convert a word to upper case
 elemental function strupcase(input)

  character(len=*), intent(in) :: input
  character(len=len(input))    :: strupcase
  integer                      :: i,ic

  strupcase = input
  do i=1, len(input)
     ic = iachar(input(i:i))
     if (ic >= 97 .and. ic < 122) strupcase(i:i) = achar(ic-32)
  end do

 end function strupcase

 ! remove blanks
 ! strcompress should have an allocatable length, but not implemented in gfortran 4.4
 elemental function strcompress(input)

  character(len=*), intent(in) :: input
  character(len=len(input))    :: strcompress
  integer                      :: i, j
  
  strcompress = ' '
  j = 1
  do i=1, len_trim(input)
     if (iachar(input(i:i)) == 32) cycle
     strcompress(j:j) = input(i:i)
     j = j + 1
  end do
      
 end function strcompress


 subroutine strsplit(input, output, delimiter)

  character(len=*), intent(in)                  :: input
  character(len=*), dimension(:), intent(inout) :: output
  character, intent(in), optional               :: delimiter
  character                                     :: delim
  integer                                       :: i, ic, j, len_in, len_out

  if (present(delimiter)) then
     delim = delimiter
  else
     delim = ','
  end if

  len_in  = len_trim(input)
  len_out = len(output)

  do i=1, size(output, 1)
     output(i) = ' '
  enddo

  ic = 1
  i  = 1
  do j=1, len_in
     if (input(j:j) == delim) then
        i = i + 1
        if (i > size(output,1)) exit
        ic = 1
        cycle
     endif
     if (ic <= len_out) output(i)(ic:ic) = input(j:j)
     ic = ic + 1
  end do

 end subroutine strsplit

 pure function strinteger(input)

  integer, intent(in)                  :: input
  character(len=strinteger_len(input)) :: strinteger
  character(len=80)                    :: string

  string = ' '
  write(string, '(i80)') input
  strinteger = adjustl(string)
      
 end function strinteger

 pure function strinteger_len(input)
 
  integer, intent(in) :: input
  integer             :: strinteger_len

  if (input == 0)  then
     strinteger_len = 1
     return
  end if

  strinteger_len = floor(log10(dble(abs(input))))+1

  if (input < 0) strinteger_len = strinteger_len + 1

 end function strinteger_len

 pure function strjoin_trim(input)
     character(len=*), intent(in)           :: input(:)
     character(len=strjoin_trim_len(input)) :: strjoin_trim
     integer                                :: i, k

     k = 1
     do i = 1, size(input)
         strjoin_trim(k:k+len_trim(input(i))-1) = trim(input(i))
         k = k + len_trim(input(i))
     end do

 end function strjoin_trim

 pure function strjoin_trim_len(input)
     character(len=*), intent(in)  :: input(:)
     integer                       :: strjoin_trim_len
     integer                       :: i
     strjoin_trim_len = sum([(len_trim(input(i)), i=1, size(input))])
 end function strjoin_trim_len

 pure function strjoin_opt(input, dotrim)
     character(len=*), intent(in)                 :: input(:)
     logical, intent(in)                          :: dotrim
     character(len=strjoin_opt_len(input,dotrim)) :: strjoin_opt
     integer                                      :: i, k

     if (dotrim) then
         k = 1
         do i = 1, size(input)
            strjoin_opt(k:k+len_trim(input(i))-1) = trim(input(i))
            k = k + len_trim(input(i))
         end do
     else
         k = len(input)
         do i = 1, size(input)
            strjoin_opt((i-1)*k+1:i*k) = input(i)
         end do
     end if

 end function strjoin_opt

 pure function strjoin_opt_len(input, dotrim)
     character(len=*), intent(in)  :: input(:)
     logical, intent(in)           :: dotrim
     integer                       :: strjoin_opt_len
     integer                       :: i
     if (.not. dotrim) then
         strjoin_opt_len = size(input) * len(input(1))
     else
         strjoin_opt_len = sum([(len_trim(input(i)), i=1, size(input))])
     end if
 end function strjoin_opt_len

end module string
