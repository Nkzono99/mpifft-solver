module m_comm_setting
    use m_block
    use mpi
    implicit none

    type t_CommSetting
        type(t_Block) :: local_block
        type(t_Block) :: comm_block
        integer :: mpi_type
        integer :: pid
        integer(kind(MPI_COMM_WORLD)) :: comm
    contains
        procedure :: isend => commSetting_isend
        procedure :: irecv => commSetting_irecv
        procedure :: destroy => commSetting_destroy
    end type

contains

    function new_CommSetting(local_block, comm_block, pid, comm) result(obj)
        type(t_Block), intent(in) :: local_block
        type(t_Block), intent(in) :: comm_block
        integer, intent(in) :: pid
        integer(kind(MPI_COMM_WORLD)), intent(in) :: comm
        type(t_CommSetting) :: obj

        integer :: sizes(3)
        integer :: subsizes(3)
        integer :: starts(3)

        integer :: ierr

        obj%local_block = local_block
        obj%comm_block = comm_block
        obj%pid = pid
        obj%comm = comm

        if (obj%comm_block%size() == 0) then
            return
        end if

        sizes(:) = local_block%sizes(:)
        subsizes(:) = comm_block%sizes(:)
        starts(:) = comm_block%start(:) - local_block%start(:)

        call MPI_Type_create_subarray(3, sizes, subsizes, starts, MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, obj%mpi_type, ierr)
        call MPI_Type_commit(obj%mpi_type, ierr)
    end function

    subroutine commSetting_isend(self, senddata, tag, request)
        class(t_CommSetting), intent(inout) :: self
        double precision, intent(in) :: senddata(self%local_block%start(1):self%local_block%end(1), &
                                                 self%local_block%start(2):self%local_block%end(2), &
                                                 self%local_block%start(3):self%local_block%end(3))
        integer, intent(in) :: tag
        integer, intent(out) :: request

        integer :: start(3)
        integer :: ierr

        if (self%comm_block%size() == 0) then
            return
        end if

        start(:) = self%local_block%start(:)
        call MPI_Isend(senddata(start(1), start(2), start(3)), 1, self%mpi_type, self%pid, tag, self%comm, request, ierr)
    end subroutine

    subroutine commSetting_irecv(self, recvdata, tag, request)
        class(t_CommSetting), intent(inout) :: self
        integer, intent(in) :: tag
        double precision, intent(inout) :: recvdata(self%local_block%start(1):self%local_block%end(1), &
                                                 self%local_block%start(2):self%local_block%end(2), &
                                                 self%local_block%start(3):self%local_block%end(3))
        integer, intent(out) :: request

        integer :: start(3)
        integer :: ierr

        if (self%comm_block%size() == 0) then
            return
        end if

        start(:) = self%local_block%start(:)
        call MPI_Irecv(recvdata(start(1), start(2), start(3)), 1, self%mpi_type, self%pid, tag, self%comm, request, ierr)
    end subroutine

    subroutine commSetting_destroy(self)
        class(t_CommSetting), intent(inout) :: self

        integer :: ierr

        call MPI_Type_free(self%mpi_type, ierr)
    end subroutine

end module
