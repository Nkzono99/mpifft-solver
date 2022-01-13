
module m_comm_setting_list
    use m_list, only: t_List, init_list
    use m_comm_setting, only: t_CommSetting
    implicit none

    type, extends(t_List) :: t_CommSettingList
        type(t_CommSetting), allocatable :: buffer(:)
        type(t_CommSetting), allocatable :: tmp_buffer(:)

    contains
        procedure :: allocate_tmp_buffer => commSettingList_allocate_tmp_buffer
        procedure :: copy_to_tmp_buffer => commSettingList_copy_to_tmp_buffer
        procedure :: switch_to_tmp_buffer => commSettingList_switch_to_tmp_buffer
        procedure :: destroy => commSettingList_destroy
        procedure :: append => commSettingList_append
        procedure :: get => commSettingList_get
    end type

    private
    public t_CommSettingList
    public new_CommSettingList

contains

    function new_CommSettingList(max_size, growth_factor) result(obj)
        integer, intent(in), optional :: max_size
        double precision, intent(in), optional :: growth_factor
        type(t_CommSettingList) :: obj

        call init_list(obj, max_size, growth_factor)
    end function

    subroutine commSettingList_allocate_tmp_buffer(self, n)
        class(t_CommSettingList), intent(inout) :: self
        integer, intent(in) :: n

        allocate (self%tmp_buffer(n))
    end subroutine

    subroutine commSettingList_copy_to_tmp_buffer(self)
        class(t_CommSettingList), intent(inout) :: self

        integer :: i

        do i = 1, self%current_size
            self%tmp_buffer(i) = self%buffer(i)
        end do

        deallocate (self%buffer)
    end subroutine

    subroutine commSettingList_switch_to_tmp_buffer(self)
        class(t_CommSettingList), intent(inout) :: self

        allocate (self%buffer, source=self%tmp_buffer)
        deallocate (self%tmp_buffer)
    end subroutine

    subroutine commSettingList_append(self, setting)
        class(t_CommSettingList), intent(inout) :: self
        type(t_CommSetting), intent(in) :: setting

        if (self%current_size == self%max_size) then
            call self%extent_size()
        end if

        self%current_size = self%current_size + 1
        self%buffer(self%current_size) = setting
    end subroutine

    subroutine commSettingList_destroy(self)
        class(t_CommSettingList), intent(inout) :: self

        type(t_CommSetting) :: setting
        integer :: i

        do i = 1, self%current_size
            setting = self%get(i)
            call setting%destroy
        end do

        deallocate (self%buffer)
    end subroutine

    function commSettingList_get(self, i) result(ret)
        class(t_CommSettingList), intent(in) :: self
        integer, intent(in) :: i
        type(t_CommSetting) :: ret

        ret = self%buffer(i)
    end function

end module
