
module m_block_list
    use m_list, only: t_List, init_list
    use m_block, only: t_Block
    implicit none

    type, extends(t_List) :: t_BlockList
        type(t_Block), allocatable :: buffer(:)
        type(t_Block), allocatable :: tmp_buffer(:)

    contains
        procedure :: allocate_tmp_buffer => blockList_allocate_tmp_buffer
        procedure :: copy_to_tmp_buffer => blockList_copy_to_tmp_buffer
        procedure :: switch_to_tmp_buffer => blockList_switch_to_tmp_buffer
        procedure :: destroy => blockList_destroy
        procedure :: append => blockList_append
        procedure :: get => blockList_get
    end type

    private
    public t_BlockList
    public new_BlockList

contains

    function new_BlockList(max_size, growth_factor) result(obj)
        integer, intent(in), optional :: max_size
        double precision, intent(in), optional :: growth_factor
        type(t_BlockList) :: obj

        call init_list(obj, max_size, growth_factor)
    end function

    subroutine blockList_allocate_tmp_buffer(self, n)
        class(t_BlockList), intent(inout) :: self
        integer, intent(in) :: n

        allocate (self%tmp_buffer(n))
    end subroutine

    subroutine blockList_copy_to_tmp_buffer(self)
        class(t_BlockList), intent(inout) :: self

        integer :: i

        do i = 1, self%current_size
            self%tmp_buffer(i) = self%buffer(i)
        end do

        deallocate (self%buffer)
    end subroutine

    subroutine blockList_switch_to_tmp_buffer(self)
        class(t_BlockList), intent(inout) :: self

        allocate(self%buffer, source=self%tmp_buffer)
        deallocate (self%tmp_buffer)
    end subroutine

    subroutine blockList_append(self, blk)
        class(t_BlockList), intent(inout) :: self
        type(t_Block), intent(in) :: blk

        if (self%current_size == self%max_size) then
            call self%extent_size()
        end if

        self%current_size = self%current_size + 1
        self%buffer(self%current_size) = blk
    end subroutine

    subroutine blockList_destroy(self)
        class(t_BlockList), intent(inout) :: self

        deallocate (self%buffer)
    end subroutine

    function blockList_get(self, i) result(ret)
        class(t_BlockList), intent(in) :: self
        integer, intent(in) :: i
        type(t_Block) :: ret

        ret = self%buffer(i)
    end function

end module
