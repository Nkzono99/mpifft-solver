! This test program have to be run by mpiexec.exe.
program test_mpifft
    use mpifft
    use mpi
    implicit none

    integer :: nx, ny, nz
    integer :: npx, npy, npz

    integer(kind(MPI_COMM_WORLD)) :: comm
    integer :: myrank
    integer :: psize
    integer :: ierr

    type(t_MPIFFTW3_Factory) :: fft_factory
    type(t_Block) :: local_block
    type(t_Block) :: global_block
    type(t_MPI_FFTExecutor3d), pointer :: fft, fft2

    double precision, allocatable :: original(:, :, :)
    double precision, allocatable :: forwarded(:, :, :)
    double precision, allocatable :: backwarded(:, :, :)

    integer :: ix, iy, iz
    double precision :: diff

    double precision :: stime
    double precision :: etime

    comm = MPI_COMM_WORLD

    call MPI_Init(ierr)
    call MPI_Comm_rank(comm, myrank, ierr)
    call MPI_Comm_size(comm, psize, ierr)

    nx = 128
    ny = 64
    nz = 32

    npx = 2
    npy = 2
    npz = 1

    global_block = new_block([1, 1, 1], [nx, ny, nz])
    ! local_block settings.
    if (npx*npy*npz < myrank + 1) then
        local_block = new_Block(global_block%start(:) - 1, global_block%start(:) - 2)
    else
        block
            integer :: ip
            integer :: ipx, ipy, ipz
            integer :: start(3), end(3)

            ip = myrank

            ! Divide a sequence of numbers from 0 to nx-1 (= 0, 1, 2, ..., nx-1) into npx pieces (similarly for y and z).
            ipx = mod(ip, npx)
            ipy = mod(int(ip/npx), npy)
            ipz = int(ip/(npx*npy))

            start(:) = [int(nx/npx)*ipx, int(ny/npy)*ipy, int(nz/npz)*ipz]
            end(:) = [nx/npx*(ipx + 1) - 1, &
                      ny/npy*(ipy + 1) - 1, &
                      nz/npz*(ipz + 1) - 1]

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
            local_block = new_Block(start(:) + global_block%start(:), end(:) + global_block%start(:))
        end block
    end if

    allocate (original(local_block%start(1):local_block%end(1), &
                       local_block%start(2):local_block%end(2), &
                       local_block%start(3):local_block%end(3)))
    allocate (forwarded(local_block%start(1):local_block%end(1), &
                        local_block%start(2):local_block%end(2), &
                        local_block%start(3):local_block%end(3)))
    allocate (backwarded(local_block%start(1):local_block%end(1), &
                         local_block%start(2):local_block%end(2), &
                         local_block%start(3):local_block%end(3)))

    forwarded(:, :, :) = 0.0d0
    backwarded(:, :, :) = 0.0d0

    do concurrent(ix=local_block%start(1):local_block%end(1), &
                  iy=local_block%start(2):local_block%end(2), &
                  iz=local_block%start(3):local_block%end(3))
        original(ix, iy, iz) = sin((2d0*3.14d0*ix)/nx)*sin((1.2d0*2d0*3.14d0*iy)/ny)*sin((-0.6d0*2d0*3.14d0*iz)/nz)
    end do

    fft_factory = new_MPIFFTW3_Factory(local_block, global_block, myrank, psize, comm)
    fft = fft_factory%create([BoundaryType_Periodic, BoundaryType_Dirichlet, BoundaryType_Neumann])

    call fft%forward(original(:, :, :), forwarded(:, :, :))
    call fft%backward(forwarded(:, :, :), backwarded(:, :, :))

    diff = sum(abs(backwarded(:, :, :) - original(:, :, :)))/sum(abs(original))
    print *, "Diff =", diff, sum(abs(backwarded(:, :, :) - original(:, :, :))), sum(abs(original)), myrank

    call MPI_Finalize(ierr)

end program
