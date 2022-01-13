module m_mpifft3_wrapper
    use m_block, only: t_Block, new_block
    use m_block_list, only: t_BlockList, new_BlockList
    use m_mpi_block_rebase, only: t_MPIBlockRebasor
    use m_fft_wrapper
    implicit none

    type, extends(t_FFTExecutor3d) :: t_MPI_FFTExecutor3d
        class(t_FFTExecutor1d), pointer :: fft_executor_x
        class(t_FFTExecutor1d), pointer :: fft_executor_y
        class(t_FFTExecutor1d), pointer :: fft_executor_z

        double precision, allocatable, private :: buf4xfft(:, :, :)
        double precision, allocatable, private :: buf4yfft(:, :, :)
        double precision, allocatable, private :: buf4zfft(:, :, :)
        double precision, allocatable, private :: buf4xffted(:, :, :)
        double precision, allocatable, private :: buf4yffted(:, :, :)
        double precision, allocatable, private :: buf4zffted(:, :, :)

        class(t_MPIBlockRebasor), pointer, private :: subdomain2xfft
        class(t_MPIBlockRebasor), pointer, private :: xfft2yfft
        class(t_MPIBlockRebasor), pointer, private :: yfft2zfft
        class(t_MPIBlockRebasor), pointer, private :: zfft2subdomain

    contains
        procedure :: forward => mpi_fftExecutor3d_forward
        procedure :: backward => mpi_fftExecutor3d_backward
        procedure :: destroy => mpi_fftExecutor3d_destroy
    end type

contains

    function new_MPI_FFTExecutor3d(fft_executor_x, &
                                   fft_executor_y, &
                                   fft_executor_z, &
                                   subdomain2xfft, &
                                   xfft2yfft, &
                                   yfft2zfft, &
                                   zfft2subdomain) result(obj)
        class(t_FFTExecutor1d), pointer, intent(in) :: fft_executor_x
        class(t_FFTExecutor1d), pointer, intent(in) :: fft_executor_y
        class(t_FFTExecutor1d), pointer, intent(in) :: fft_executor_z

        class(t_MPIBlockRebasor), pointer, intent(in) :: subdomain2xfft
        class(t_MPIBlockRebasor), pointer, intent(in) :: xfft2yfft
        class(t_MPIBlockRebasor), pointer, intent(in) :: yfft2zfft
        class(t_MPIBlockRebasor), pointer, intent(in) :: zfft2subdomain
        type(t_MPI_FFTExecutor3d) :: obj

        integer :: nx, ny, nz
        integer :: boundary_types(3)

        obj%fft_executor_x => fft_executor_x
        obj%fft_executor_y => fft_executor_y
        obj%fft_executor_z => fft_executor_z
        obj%subdomain2xfft => subdomain2xfft
        obj%xfft2yfft => xfft2yfft
        obj%yfft2zfft => yfft2zfft
        obj%zfft2subdomain => zfft2subdomain

        nx = fft_executor_x%n
        ny = fft_executor_y%n
        nz = fft_executor_z%n

        boundary_types(1:3) = [fft_executor_x%boundary_type, &
                               fft_executor_y%boundary_type, &
                               fft_executor_z%boundary_type]

        allocate (obj%buf4xfft(nx, &
                               obj%subdomain2xfft%require_block%sizes(2), &
                               obj%subdomain2xfft%require_block%sizes(3)))

        allocate (obj%buf4yfft(obj%xfft2yfft%require_block%sizes(1), &
                               ny, &
                               obj%xfft2yfft%require_block%sizes(3)))

        allocate (obj%buf4zfft(obj%yfft2zfft%require_block%sizes(1), &
                               obj%yfft2zfft%require_block%sizes(2), &
                               nz))

        allocate (obj%buf4xffted, mold=obj%buf4xfft)
        allocate (obj%buf4yffted, mold=obj%buf4yfft)
        allocate (obj%buf4zffted, mold=obj%buf4zfft)

        call init_FFTExecutor3d(obj, nx, ny, nz, boundary_types(1:3))
    end function

    subroutine mpi_fftExecutor3d_forward(self, in, out)
        class(t_MPI_FFTExecutor3d), intent(inout) :: self
        double precision, intent(in) :: in(:, :, :)
        double precision, intent(inout) :: out(:, :, :)

        integer :: ix, iy, iz

        call self%subdomain2xfft%run(in(:, :, :), self%buf4xfft(:, :, :))
        do iz = 1, self%subdomain2xfft%require_block%sizes(3)
            do iy = 1, self%subdomain2xfft%require_block%sizes(2)
                call self%fft_executor_x%forward(self%buf4xfft(:, iy, iz), self%buf4xffted(:, iy, iz))
            end do
        end do

        call self%xfft2yfft%run(self%buf4xffted(:, :, :), self%buf4yfft(:, :, :))
        do iz = 1, self%xfft2yfft%require_block%sizes(3)
            do ix = 1, self%xfft2yfft%require_block%sizes(1)
                call self%fft_executor_y%forward(self%buf4yfft(ix, :, iz), self%buf4yffted(ix, :, iz))
            end do
        end do

        call self%yfft2zfft%run(self%buf4yffted(:, :, :), self%buf4zfft(:, :, :))
        do iy = 1, self%yfft2zfft%require_block%sizes(2)
            do ix = 1, self%yfft2zfft%require_block%sizes(1)
                call self%fft_executor_z%forward(self%buf4zfft(ix, iy, :), self%buf4zffted(ix, iy, :))
            end do
        end do

        call self%zfft2subdomain%run(self%buf4zffted(:, :, :), out(:, :, :))
    end subroutine

    subroutine mpi_fftExecutor3d_backward(self, in, out)
        class(t_MPI_FFTExecutor3d), intent(inout) :: self
        double precision, intent(in) :: in(:, :, :)
        double precision, intent(inout) :: out(:, :, :)

        integer :: ix, iy, iz

        call self%subdomain2xfft%run(in(:, :, :), self%buf4xfft(:, :, :))
        do iz = 1, self%subdomain2xfft%require_block%sizes(3)
            do iy = 1, self%subdomain2xfft%require_block%sizes(2)
                call self%fft_executor_x%backward(self%buf4xfft(:, iy, iz), self%buf4xffted(:, iy, iz))
            end do
        end do

        call self%xfft2yfft%run(self%buf4xffted(:, :, :), self%buf4yfft(:, :, :))
        do iz = 1, self%xfft2yfft%require_block%sizes(3)
            do ix = 1, self%xfft2yfft%require_block%sizes(1)
                call self%fft_executor_y%backward(self%buf4yfft(ix, :, iz), self%buf4yffted(ix, :, iz))
            end do
        end do

        call self%yfft2zfft%run(self%buf4yffted(:, :, :), self%buf4zfft(:, :, :))
        do iy = 1, self%yfft2zfft%require_block%sizes(2)
            do ix = 1, self%yfft2zfft%require_block%sizes(1)
                call self%fft_executor_z%backward(self%buf4zfft(ix, iy, :), self%buf4zffted(ix, iy, :))
            end do
        end do

        call self%zfft2subdomain%run(self%buf4zffted(:, :, :), out(:, :, :))
    end subroutine

    subroutine mpi_fftExecutor3d_destroy(self)
        class(t_MPI_FFTExecutor3d), intent(inout) :: self

        call self%fft_executor_x%destroy
        call self%fft_executor_y%destroy
        call self%fft_executor_z%destroy
    end subroutine

end module
