module m_mpifftw3_factory
    use m_block, only: t_Block, new_block, SIZE_OF_BLOCK_ARRAY
    use m_block_list, only: t_BlockList, new_BlockList
    use m_mpi_block_rebase, only: t_MPIBlockRebasor, new_MPIBlockRebasor
    use m_fft_wrapper, only: t_FFTExecutor1d
    use m_fftw3_wrapper_1d, only: t_FFTW3Executor1d, new_FFTW3Executor1d
    use m_mpifft3_wrapper, only: t_MPI_FFTExecutor3d, new_MPI_FFTExecutor3d
    use m_get_default, only: get_default
    use mpi
    implicit none

    type :: t_MPIFFTW3_Factory
        type(t_Block), private :: global_block
        type(t_MPIBlockRebasor), private :: s2x
        type(t_MPIBlockRebasor), private :: x2y
        type(t_MPIBlockRebasor), private :: y2z
        type(t_MPIBlockRebasor), private :: z2s
    contains
        procedure :: create => mpifftw3_Factory_create
    end type

    private
    public t_MPIFFTW3_Factory
    public new_MPIFFTW3_Factory
    public new_MPIFFTW3_Factory_with_blocks

contains

    function new_MPIFFTW3_Factory(local_block, global_block, myrank, nproc, comm, tag) result(obj)
        type(t_Block), intent(in) :: local_block
        type(t_Block), intent(in) :: global_block
        integer, intent(in) :: myrank
        integer, intent(in) :: nproc
        integer(kind(MPI_COMM_WORLD)), intent(in) :: comm
        integer, intent(in), optional :: tag
        type(t_MPIFFTW3_Factory) :: obj

        integer :: send_data(SIZE_OF_BLOCK_ARRAY)
        integer, allocatable :: recv_datas(:, :)
        integer, allocatable :: send_requests(:)
        integer, allocatable :: recv_requests(:)
        integer, allocatable :: send_status(:, :), recv_status(:, :)
        integer :: ierr

        integer :: i, ip

        type(t_Block) :: blk
        type(t_BlockList) :: blocks

        allocate (recv_datas(SIZE_OF_BLOCK_ARRAY, nproc))

        allocate (send_requests(nproc))
        allocate (recv_requests(nproc))

        call local_block%to_array(send_data(:))

        do i = 1, nproc
            ip = i - 1
            call MPI_Isend(send_data(1), SIZE_OF_BLOCK_ARRAY, MPI_INTEGER, ip, get_default(tag, 0), &
                           comm, send_requests(i), ierr)
            call MPI_Irecv(recv_datas(1, i), SIZE_OF_BLOCK_ARRAY, MPI_INTEGER, ip, get_default(tag, 0), &
                           comm, recv_requests(i), ierr)
        end do

        allocate (recv_status(MPI_STATUS_SIZE, nproc))
        call MPI_Waitall(nproc, recv_requests(:), recv_status(:, :), ierr)

        allocate (send_status(MPI_STATUS_SIZE, nproc))
        call MPI_Waitall(nproc, send_requests(:), send_status(:, :), ierr)

        blocks = new_BlockList()
        do i = 1, nproc
            call blk%from_array(recv_datas(:, i))
            call blocks%append(blk)
        end do

        obj = new_MPIFFTW3_Factory_with_blocks(blocks, global_block, myrank, comm)
    end function

    function new_MPIFFTW3_Factory_with_blocks(blocks, global_block, myrank, comm) result(obj)
        type(t_BlockList), intent(in) :: blocks
        type(t_Block), intent(in) :: global_block
        integer, intent(in) :: myrank
        integer(kind(MPI_COMM_WORLD)), intent(in) :: comm
        type(t_MPIFFTW3_Factory) :: obj

        integer :: nproc
        integer, allocatable :: pids(:)
        type(t_BlockList) :: block4xffts, block4yffts, block4zffts

        integer :: i

        obj%global_block = global_block

        nproc = blocks%current_size
        allocate (pids(nproc))
        do i = 1, nproc
            pids(i) = i - 1
        end do

        block4xffts = create_blocks_for_fft(1, global_block, nproc)
        block4yffts = create_blocks_for_fft(2, global_block, nproc)
        block4zffts = create_blocks_for_fft(3, global_block, nproc)

        obj%s2x = new_MPIBlockRebasor(blocks, block4xffts, pids, myrank + 1, comm)
        obj%x2y = new_MPIBlockRebasor(block4xffts, block4yffts, pids, myrank + 1, comm)
        obj%y2z = new_MPIBlockRebasor(block4yffts, block4zffts, pids, myrank + 1, comm)
        obj%z2s = new_MPIBlockRebasor(block4zffts, blocks, pids, myrank + 1, comm)
    end function

    function mpifftw3_Factory_create(self, boundary_types) result(mpifft3d)
        class(t_MPIFFTW3_Factory), target, intent(in) :: self
        integer, intent(in) :: boundary_types(3)
        type(t_MPI_FFTExecutor3d) :: mpifft3d

        type(t_FFTW3Executor1d), pointer :: xfft, yfft, zfft

        block
            integer :: nx, ny, nz

            nx = self%global_block%sizes(1)
            ny = self%global_block%sizes(2)
            nz = self%global_block%sizes(3)

            allocate (xfft)
            allocate (yfft)
            allocate (zfft)
            xfft = new_FFTW3Executor1d(nx, boundary_types(1))
            yfft = new_FFTW3Executor1d(ny, boundary_types(2))
            zfft = new_FFTW3Executor1d(nz, boundary_types(3))
        end block

        block
            class(t_FFTExecutor1d), pointer :: pxfft, pyfft, pzfft
            class(t_MPIBlockRebasor), pointer :: s2x, x2y, y2z, z2s

            pxfft => xfft
            pyfft => yfft
            pzfft => zfft
            s2x => self%s2x
            x2y => self%x2y
            y2z => self%y2z
            z2s => self%z2s
            mpifft3d = new_MPI_FFTExecutor3d(pxfft, pyfft, pzfft, s2x, x2y, y2z, z2s)
        end block
    end function

    function create_blocks_for_fft(axis, global_block, nproc) result(blocks)
        integer, intent(in):: axis
        type(t_Block), intent(in) :: global_block
        integer, intent(in)::nproc
        type(t_BlockList) :: blocks

        type(t_Block) :: blk

        integer :: nx, ny, nz
        integer :: ip
        integer :: npx, npy, npz

        nx = global_block%sizes(1)
        ny = global_block%sizes(2)
        nz = global_block%sizes(3)

        select case (axis)
        case (1)
            npx = 1
            npy = max(1, min(ny, nproc - min(nz, nproc)))
            npz = min(nz, nproc)
        case (2)
            npx = max(1, min(nx, nproc - min(nz, nproc)))
            npy = 1
            npz = min(nz, nproc)
        case (3)
            npx = max(1, min(nx, nproc - min(ny, nproc)))
            npy = min(ny, nproc)
            npz = 1
        end select

        blocks = new_BlockList()
        do ip = 0, nproc - 1
            if (npx*npy*npz < ip + 1) then
                blk = new_Block(global_block%start(:) - 1, global_block%start(:) - 2)
                call blocks%append(blk)
                cycle
            end if

            block
                integer :: ipx, ipy, ipz
                integer :: start(3), end(3)

                ! Divide a sequence of numbers from 0 to nx-1 (= 0, 1, 2, ..., nx-1) into npx pieces (similarly for y and z).
                ipx = mod(ip, npx)
                ipy = mod(int(ip/npx), npy)
                ipz = int(ip/(npx*npy))

                start(:) = [int(nx/npx)*ipx, int(ny/npy)*ipy, int(nz/npz)*ipz]
                end(:) = [int(nx/npx)*(ipx + 1) - 1, &
                          int(ny/npy)*(ipy + 1) - 1, &
                          int(nz/npz)*(ipz + 1) - 1]

                if (ipx == npx - 1) then
                    end(1) = nx - 1
                end if
                if (ipy == npy - 1) then
                    end(2) = ny - 1
                end if
                if (ipz == npz - 1) then
                    end(3) = nz - 1
                end if

                ! Shift by global_block offset.
                blk = new_Block(start(:) + global_block%start(:), end(:) + global_block%start(:))
            end block

            call blocks%append(blk)
        end do
    end function

end module
