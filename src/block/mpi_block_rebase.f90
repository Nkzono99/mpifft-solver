module m_mpi_block_rebase
    use m_block
    use m_block_list
    use m_comm_setting
    use m_comm_setting_list
    use mpi
    use m_get_default, only: get_default
    implicit none

    type t_MPIBlockRebasor
        type(t_Block) :: local_block
        type(t_Block) :: require_block

        type(t_CommSettingList), private :: send_settings
        type(t_CommSettingList), private :: recv_settings

        integer(kind=kind(MPI_COMM_WORLD)), private :: comm
    contains
        procedure :: run => mpiBlockRebasor_run
        procedure :: destroy => mpiBlockRebasor_destroy
    end type

    private
    public t_MPIBlockRebasor
    public new_MPIBlockRebasor

contains

    function new_MPIBlockRebasor(blocks, require_blocks, pids, ipid, comm) result(obj)
        type(t_BlockList), intent(in) :: blocks
        type(t_BlockList), intent(in) :: require_blocks
        !> Array of process ID (= rank = 0, 1, 2, ..., nproc-1).
        !> It should be the same as the order of blocks and require_blocks.
        integer, intent(in) :: pids(:)
        !> Index to identify the current process ID from pids(1:len(pids)).
        integer, intent(in) :: ipid
        !> MPI Communicator.
        integer(kind=kind(MPI_COMM_WORLD)), intent(in) :: comm
        type(t_MPIBlockRebasor) :: obj

        type(t_Block) :: overlapped
        type(t_CommSetting) :: setting

        integer :: i

        obj%local_block = blocks%get(ipid)
        obj%require_block = require_blocks%get(ipid)

        obj%send_settings = new_CommSettingList()
        obj%recv_settings = new_CommSettingList()

        !> Send settings
        do i = 1, require_blocks%current_size
            overlapped = obj%local_block%overlapped(require_blocks%get(i))

            if (overlapped%size() == 0) then
                cycle
            end if

            setting = new_CommSetting(obj%local_block, overlapped, pids(i), comm)
            call obj%send_settings%append(setting)
        end do

        ! Recv settings
        do i = 1, blocks%current_size
            overlapped = obj%require_block%overlapped(blocks%get(i))
            if (overlapped%size() == 0) then
                cycle
            end if

            setting = new_CommSetting(obj%require_block, overlapped, pids(i), comm)
            call obj%recv_settings%append(setting)
        end do

        obj%comm = comm
    end function

    subroutine mpiBlockRebasor_run(self, send_data, recv_data, tag)
        class(t_MPIBlockRebasor), intent(in) :: self
        double precision, intent(in) :: send_data(self%local_block%start(1):self%local_block%end(1), &
                                                  self%local_block%start(2):self%local_block%end(2), &
                                                  self%local_block%start(3):self%local_block%end(3))
        double precision, intent(inout) :: recv_data(self%require_block%start(1):self%require_block%end(1), &
                                                     self%require_block%start(2):self%require_block%end(2), &
                                                     self%require_block%start(3):self%require_block%end(3))
        integer, intent(in), optional :: tag

        type(t_CommSetting) :: setting
        integer, allocatable :: send_requests(:)
        integer, allocatable :: recv_requests(:)

        integer :: i
        integer :: ierr
        integer, allocatable :: send_status(:, :), recv_status(:, :)

        allocate (send_requests(self%send_settings%current_size))
        allocate (recv_requests(self%recv_settings%current_size))

        do i = 1, self%send_settings%current_size
            setting = self%send_settings%get(i)
            call setting%isend(send_data(:, :, :), get_default(tag, 0), send_requests(i))
        end do

        do i = 1, self%recv_settings%current_size
            setting = self%recv_settings%get(i)
            call setting%irecv(recv_data(:, :, :), get_default(tag, 0), recv_requests(i))
        end do

        allocate(recv_status(MPI_STATUS_SIZE, self%recv_settings%current_size))
        call MPI_Waitall(self%recv_settings%current_size, recv_requests(:), recv_status(:, :), ierr)

        allocate(send_status(MPI_STATUS_SIZE, self%send_settings%current_size))
        call MPI_Waitall(self%send_settings%current_size, send_requests(:), send_status(:, :), ierr)
    end subroutine

    subroutine mpiBlockRebasor_destroy(self)
        class(t_MPIBlockRebasor), intent(inout) :: self

        call self%send_settings%destroy
        call self%recv_settings%destroy
    end subroutine

end module
