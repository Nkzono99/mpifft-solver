module m_poisson_solver
    use mpifft
    use m_constants, only: PI
    use m_get_default, only: get_default
    implicit none

    !> 3d poisson equation solver.
    !>
    !> Poisson equation:
    !>     ∂^2p/∂^2 + ∂^2p/∂y^2 + ∂^2p/∂z^2 = f(x, y, z)
    type :: t_PoissonSolver3d
        double precision, private :: dx
        double precision, private :: dy
        double precision, private :: dz
        class(t_MPI_FFTExecutor3d), allocatable, private :: fft3d
        double precision, allocatable, private :: modified_wave_number(:, :, :)
    contains
        procedure :: solve => poissonSolver3d_solve
    end type

contains

    function new_PoissonSolver3d(local_block, global_block, fft3d, dx, dy, dz) result(obj)
        type(t_Block), intent(in) :: local_block
        type(t_Block), intent(in) :: global_block
        class(t_MPI_FFTExecutor3d), intent(in) :: fft3d
        double precision, intent(in) :: dx
        double precision, intent(in), optional :: dy
        double precision, intent(in), optional :: dz
        type(t_PoissonSolver3d) :: obj

        integer :: start(3), end(3)
        integer :: kx, ky, kz

        obj%fft3d = fft3d
        obj%dx = dx
        obj%dy = get_default(dy, dx)
        obj%dz = get_default(dz, dx)

        allocate (obj%modified_wave_number(local_block%sizes(1), &
                                           local_block%sizes(2), &
                                           local_block%sizes(3)))

        start(:) = local_block%start(:) - global_block%start(:)
        end(:) = local_block%end(:) - global_block%start(:)
        do concurrent(kz=start(3):end(3), ky=start(2):end(2), kx=start(1):end(1))
            block
                double precision :: wx, wy, wz

                integer :: ix, iy, iz
                double precision :: wn

                wx = calc_wave_number(kx, fft3d%nx, fft3d%boundary_types(1))
                wy = calc_wave_number(ky, fft3d%ny, fft3d%boundary_types(2))
                wz = calc_wave_number(kz, fft3d%nz, fft3d%boundary_types(3))

                ix = kx + global_block%start(1) - local_block%start(1) + 1
                iy = ky + global_block%start(2) - local_block%start(2) + 1
                iz = kz + global_block%start(3) - local_block%start(3) + 1
                wn = wx/(obj%dx*obj%dx) + wy/(obj%dy*obj%dy) + wz/(obj%dz*obj%dz)

                obj%modified_wave_number(ix, iy, iz) = wn
            end block
        end do
    end function

    function calc_wave_number(k, n, boundary_type) result(wn)
        integer, intent(in) :: k
        integer, intent(in) :: n
        integer, intent(in) :: boundary_type
        double precision :: wn

        select case (boundary_type)
        case (BoundaryType_Periodic)
            if (k <= int(n/2)) then
                wn = 2.0d0*sin(PI*k/dble(n))
            else
                wn = 2.0d0*sin(PI*(k - int(n/2))/dble(n))
            end if

        case (BoundaryType_Dirichlet)
            wn = 2.0d0*(cos(PI*(k + 1)/dble(n + 1)) - 1.0d0)

        case (BoundaryType_Neumann)
            wn = 2.0d0*(cos(PI*k/dble(n)) - 1.0d0)

        case (BoundaryType_Dirichlet_Neumann)
            wn = 2.0d0*(cos(PI*(k + 0.5d0)/dble(n + 1)) - 1.0d0)

        case (BoundaryType_Neumann_Dirichlet)
            wn = 2.0d0*(cos(PI*(k + 0.5d0)/dble(n + 1)) - 1.0d0)

        end select
    end function

    subroutine poissonSolver3d_solve(self, f, p)
        class(t_PoissonSolver3d), intent(inout) :: self
        double precision, intent(in) :: f(:, :, :)
        double precision, intent(out) :: p(:, :, :)

        double precision, allocatable :: fk(:, :, :)
        double precision, allocatable :: pk(:, :, :)

        allocate (fk, mold=f)
        allocate (pk, mold=p)

        call self%fft3d%forward(f(:, :, :), fk(:, :, :))

        pk(:, :, :) = f(:, :, :)/self%modified_wave_number(:, :, :)

        call self%fft3d%backward(pk(:, :, :), p(:, :, :))
    end subroutine
end module
