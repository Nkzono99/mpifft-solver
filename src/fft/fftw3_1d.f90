module m_fftw3_wrapper_1d
    use m_boundary_type
    use m_fft_wrapper
    use fftw3
    implicit none

    type, extends(t_FFTExecutor1d) :: t_FFTW3Executor1d
        type(C_PTR), private :: forward_plan
        type(C_PTR), private :: backward_plan

        real(C_DOUBLE), pointer, private :: in(:)
        real(C_DOUBLE), pointer, private :: out(:)
        type(C_PTR), private :: p_in
        type(C_PTR), private :: p_out

        double precision, private :: normalizer
    contains
        procedure :: forward => fftw3Executor1d_forward
        procedure :: backward => fftw3Executor1d_backward
        procedure :: destroy => fftw3Executor1d_destroy
    end type

    private
    public t_FFTW3Executor1d
    public new_FFTW3Executor1d

contains

    subroutine set_fft_type(n, boundary_type, forward_fft_type, backward_fft_type, normalizer)
        integer(C_INT), intent(in) :: n
        integer, intent(in) :: boundary_type
        integer(C_FFTW_R2R_KIND), intent(out) :: forward_fft_type
        integer(C_FFTW_R2R_KIND), intent(out) :: backward_fft_type
        double precision, intent(out) :: normalizer

        select case (boundary_type)
        case (BoundaryType_Periodic)
            forward_fft_type = FFTW_R2HC
            backward_fft_type = FFTW_HC2R
            normalizer = n

        case (BoundaryType_Dirichlet)
            forward_fft_type = FFTW_RODFT00
            backward_fft_type = FFTW_RODFT00
            normalizer = 2.0d0*(n + 1.0d0)

        case (BoundaryType_Neumann)
            forward_fft_type = FFTW_REDFT00
            backward_fft_type = FFTW_REDFT00
            normalizer = 2.0d0*(n - 1.0d0)

        case (BoundaryType_Dirichlet_Neumann)
            forward_fft_type = FFTW_RODFT01
            backward_fft_type = FFTW_RODFT10
            normalizer = 2.0d0*n

        case (BoundaryType_Neumann_Dirichlet)
            forward_fft_type = FFTW_REDFT01
            backward_fft_type = FFTW_REDFT10
            normalizer = 2.0d0*n
        end select
    end subroutine

    function new_FFTW3Executor1d(n, boundary_type) result(obj)
        integer(C_INT), intent(in) :: n
        integer, intent(in) :: boundary_type
        type(t_FFTW3Executor1d) :: obj

        integer(C_FFTW_R2R_KIND) :: forward_fft_type, backward_fft_type

        call init_FFTExecutor1d(obj, n, boundary_type)

        call set_fft_type(n, boundary_type, forward_fft_type, backward_fft_type, obj%normalizer)

        obj%p_in = fftw_alloc_real(int(n, C_SIZE_T))
        call c_f_pointer(obj%p_in, obj%in, [n])

        obj%p_out = fftw_alloc_real(int(n, C_SIZE_T))
        call c_f_pointer(obj%p_out, obj%out, [n])

        obj%forward_plan = fftw_plan_r2r_1d(obj%n, obj%in, obj%out, forward_fft_type, FFTW_ESTIMATE)
        obj%backward_plan = fftw_plan_r2r_1d(obj%n, obj%in, obj%out, backward_fft_type, FFTW_ESTIMATE)
    end function

    subroutine fftw3Executor1d_forward(self, in, out)
        class(t_FFTW3Executor1d), intent(inout) :: self
        double precision, intent(in) :: in(:)
        double precision, intent(inout) :: out(:)

        self%in(1:self%n) = in(1:self%n)
        call dfftw_execute_r2r(self%forward_plan, self%in, self%out)
        out(1:self%n) = self%out(1:self%n)
    end subroutine

    subroutine fftw3Executor1d_backward(self, in, out)
        class(t_FFTW3Executor1d), intent(inout) :: self
        double precision, intent(in) :: in(:)
        double precision, intent(inout) :: out(:)

        self%in(1:self%n) = in(1:self%n)
        call dfftw_execute_r2r(self%backward_plan, self%in, self%out)
        out(1:self%n) = self%out(1:self%n)/self%normalizer
    end subroutine

    subroutine fftw3Executor1d_destroy(self)
        class(t_FFTW3Executor1d), intent(inout) :: self

        call fftw_free(self%p_in)
        call fftw_free(self%p_out)
        call fftw_destroy_plan(self%forward_plan)
        call fftw_destroy_plan(self%backward_plan)
    end subroutine

end module
