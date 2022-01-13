module m_block

    integer, parameter :: SIZE_OF_BLOCK_ARRAY = 6

    type t_Block
        integer :: start(3)
        integer :: end(3)
        integer :: sizes(3)
    contains
        procedure :: size => block_size
        procedure :: overlapped => block_overlapped
        procedure :: to_array => block_to_array
        procedure :: from_array => block_from_array
    end type

    private
    public SIZE_OF_BLOCK_ARRAY
    public t_Block
    public new_block

contains

    function new_block(start, end) result(obj)
        integer, intent(in) :: start(3)
        integer, intent(in) :: end(3)
        type(t_Block) :: obj

        obj%start(1:3) = start(1:3)
        obj%end(1:3) = end(1:3)
        obj%sizes(1:3) = max(obj%end(1:3) - obj%start(1:3) + 1, 0)
    end function

    function block_size(self) result(size)
        class(t_Block), intent(in) :: self
        integer :: size

        size = product(self%sizes(1:3))
    end function

    function block_overlapped(self, other) result(overlapped)
        class(t_Block), intent(in) :: self
        class(t_Block), intent(in) :: other
        type(t_Block) :: overlapped

        integer :: start(3), end(3)

        if (self%size() == 0 .or. other%size() == 0) then
            overlapped = new_block([0, 0, 0], [-1, -1, -1])
            return
        end if 

        start(1:3) = max(self%start(1:3), other%start(1:3))
        end(1:3) = min(self%end(1:3), other%end(1:3))

        overlapped = new_block(start, end)
    end function

    subroutine block_to_array(self, array)
        class(t_Block), intent(in) :: self
        integer, intent(out) :: array(SIZE_OF_BLOCK_ARRAY)

        array(1:3) = self%start(1:3)
        array(4:6) = self%end(1:3)
    end subroutine

    subroutine block_from_array(self, array)
        class(t_Block), intent(inout) :: self
        integer, intent(in) :: array(SIZE_OF_BLOCK_ARRAY)

        self%start(1:3) = array(1:3)
        self%end(1:3) = array(4:6)
        self%sizes(1:3) = max(self%end(1:3) - self%start(1:3) + 1, 0)
    end subroutine

end module
