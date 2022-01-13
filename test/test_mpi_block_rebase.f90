! This test program have to be run by mpiexec.exe.
program test_block
    use m_block
    use m_block_list
    use m_mpi_block_rebase
    use mpi
    implicit none

    integer :: rank
    integer :: psize
    integer :: ierr

    integer(kind(MPI_COMM_WORLD)) :: comm

    comm = MPI_COMM_WORLD

    call MPI_Init(ierr)
    call MPI_Comm_rank(comm, rank, ierr)
    call MPI_Comm_size(comm, psize, ierr)

    call test_block_comm_run(4, 1, 1, 2, 1, 2)

    call MPI_Finalize(ierr)

contains

    subroutine test_block_comm_run(npx1, npy1, npz1, npx2, npy2, npz2)
        integer, intent(in) :: npx1, npy1, npz1
        integer, intent(in) :: npx2, npy2, npz2
        type(t_MPIBlockRebasor) :: bcomm

        type(t_Block) :: blk
        type(t_Block) :: local_block, require_block
        type(t_BlockList) :: blocks, require_blocks

        integer :: ix, iy, iz
        integer :: ip, ipx, ipy, ipz
        double precision, allocatable :: send_data(:, :, :), recv_data(:, :, :)
        integer :: start(3), end(3)

        integer :: nx, ny, nz
        double precision :: diff
        integer, allocatable :: pids(:)

        integer :: npids

        nx = 32
        ny = 32
        nz = 32

        allocate(pids(psize))

        blocks = new_BlockList()
        require_blocks = new_BlockList()

        npids = 0
        do ip = 0, psize - 1
            ipx = mod(ip, npx1)
            ipy = mod(int(ip/npx1), npy1)
            ipz = int(ip/(npx1*npy1))

            start(:) = [nx/npx1*ipx, ny/npy1*ipy, nz/npz1*ipz]
            end(:) = [nx/npx1*(ipx + 1) - 1, ny/npy1*(ipy + 1) - 1, nz/npz1*(ipz + 1) - 1]
            blk = new_Block(start, end)

            npids = npids + 1
            pids(npids) = ip
            call blocks%append(blk)
        end do

        do ip = 0, psize - 1
            ipx = mod(ip, npx2)
            ipy = mod(int(ip/npx2), npy2)
            ipz = int(ip/(npx2*npy2))

            start(:) = [nx/npx2*ipx, ny/npy2*ipy, nz/npz2*ipz]
            end(:) = [nx/npx2*(ipx + 1) - 1, ny/npy2*(ipy + 1) - 1, nz/npz2*(ipz + 1) - 1]
            blk = new_Block(start, end)

            call require_blocks%append(blk)
        end do

        local_block = blocks%get(rank+1)
        require_block = require_blocks%get(rank+1)

        allocate (send_data(local_block%start(1):local_block%end(1), &
                            local_block%start(2):local_block%end(2), &
                            local_block%start(3):local_block%end(3)))

        allocate (recv_data(require_block%start(1):require_block%end(1), &
                            require_block%start(2):require_block%end(2), &
                            require_block%start(3):require_block%end(3)))

        bcomm = new_MPIBlockRebasor(blocks, require_blocks, pids, rank+1, comm)

        do iz = local_block%start(3), local_block%end(3)
            do iy = local_block%start(2), local_block%end(2)
                do ix = local_block%start(1), local_block%end(1)
                    send_data(ix, iy, iz) = ix + iy*nx + iz*nx*ny
                end do
            end do
        end do

        recv_data = 0
        call bcomm%run(send_data, recv_data)

        diff = 0
        do iz = require_block%start(3), require_block%end(3)
            do iy = require_block%start(2), require_block%end(2)
                do ix = require_block%start(1), require_block%end(1)
                    diff = diff + abs(recv_data(ix, iy, iz) - (ix + iy*nx + iz*nx*ny))
                end do
            end do
        end do

        call MPI_Barrier(comm, ierr)

        print *, 'Diff =', diff, sum(send_data), sum(recv_data), rank
    end subroutine

end program
