module m_fft_wrapper
    implicit none

    type, abstract :: t_FFTExecutor1d
        !> Length of array.
        integer :: n
        !> Boundary type.
        integer :: boundary_type
    contains
        procedure(fftExecutor1d_transform), deferred :: forward
        procedure(fftExecutor1d_transform), deferred :: backward
        procedure(fftExecutor1d_destroy), deferred :: destroy
    end type

    interface
        subroutine fftExecutor1d_transform(self, in, out)
            import t_FFTExecutor1d
            class(t_FFTExecutor1d), intent(inout) :: self
            double precision, intent(in) :: in(:)
            double precision, intent(inout) :: out(:)
        end subroutine

        subroutine fftExecutor1d_destroy(self)
            import t_FFTExecutor1d
            class(t_FFTExecutor1d), intent(inout) :: self
        end subroutine
    end interface

    type, abstract :: t_FFTExecutor3d
        !> Length of array.
        integer :: nx, ny, nz
        !> Boundary type.
        integer :: boundary_types(3)
    contains
        procedure(fftExecutor3d_transform), deferred :: forward
        procedure(fftExecutor3d_transform), deferred :: backward
        procedure(fftExecutor3d_destroy), deferred :: destroy
    end type

    interface
        subroutine fftExecutor3d_transform(self, in, out)
            import t_FFTExecutor3d
            class(t_FFTExecutor3d), intent(inout) :: self
            double precision, intent(in) :: in(:, :, :)
            double precision, intent(inout) :: out(:, :, :)
        end subroutine

        subroutine fftExecutor3d_destroy(self)
            import t_FFTExecutor3d
            class(t_FFTExecutor3d), intent(inout) :: self
        end subroutine
    end interface

    private

    public t_FFTExecutor1d
    public init_FFTExecutor1d

    public t_FFTExecutor3d
    public init_FFTExecutor3d

contains

    subroutine init_FFTExecutor1d(self, n, boundary_type)
        class(t_FFTExecutor1d), intent(inout) :: self
        integer, intent(in) :: n
        integer :: boundary_type

        self%n = n
        self%boundary_type = boundary_type
    end subroutine

    subroutine init_FFTExecutor3d(self, nx, ny, nz, boundary_types)
        class(t_FFTExecutor3d), intent(inout) :: self
        integer, intent(in) :: nx
        integer, intent(in) :: ny
        integer, intent(in) :: nz
        integer :: boundary_types(3)

        self%nx = nx
        self%ny = ny
        self%nz = nz
        self%boundary_types(:) = boundary_types(:)
    end subroutine

end module
