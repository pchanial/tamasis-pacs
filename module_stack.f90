! A stack implementation
module module_stack

    implicit none

    type stack_int
        integer                      :: nelements = 0
        type(stackcell_int), pointer :: head => null()
    contains
        procedure :: pop
        procedure :: push
        procedure :: is_empty
        procedure :: to_array
        procedure :: print
    end type stack_int

    type stackcell_int
        integer                      :: value
        type(stackcell_int), pointer :: next => null()
    end type stackcell_int


contains


    function pop(this)

        class(stack_int), intent(inout) :: this
        integer                         :: pop
        type(stackcell_int), pointer    :: current

        if (this%is_empty()) then
            stop "Stack is empty."
        end if
        current => this%head
        pop = current%value
        deallocate(this%head)
        this%head => current%next
        this%nelements = this%nelements - 1

    end function pop


    !-------------------------------------------------------------------------------


    subroutine push(this, value)

        class(stack_int), intent(inout) :: this
        integer, intent(in)             :: value
        integer                         :: istat
        type(stackcell_int), pointer    :: current

        allocate(current, stat=istat)
        if (istat /= 0) stop "Error allocating new cell in stack."
        current%value = value
        if (.not. this%is_empty()) then
            current%next => this%head
        end if
        this%head => current
        this%nelements = this%nelements + 1

    end subroutine push


    !-------------------------------------------------------------------------------


    function is_empty(this)

        class(stack_int), intent(in) :: this
        logical                  :: is_empty

        is_empty = this%nelements == 0

    end function


    !-------------------------------------------------------------------------------


    subroutine print(this)

        class(stack_int), intent(in) :: this
        integer :: i
        type(stackcell_int), pointer :: current

        if (this%is_empty()) then
           write (*,*) "Stack is empty"
           return
        end if

        current => this%head
        i = 1
        do while (.true.)
           write (*,*) i, current%value
           i = i + 1
           if (.not. associated(current%next)) exit
           current => current%next
        end do

    end subroutine print


    !-------------------------------------------------------------------------------


    subroutine to_array(this, array)

        class(stack_int), intent(in)      :: this
        integer, intent(out), allocatable :: array(:)
        type(stackcell_int), pointer      :: current
        integer                           :: i

        allocate(array(this%nelements))

        current => this%head
        do i = 1, this%nelements
           array(i) = current%value
           current => current%next
        end do

     end subroutine to_array


end module module_stack
