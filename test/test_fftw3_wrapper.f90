program test_fftw3_wrapper
    use, intrinsic :: iso_c_binding
    use m_boundary_type
    use m_fft_wrapper
    use m_fftw3_wrapper_1d
    implicit none

    print *, "test_fftw31d periodic"
    call test_fftw31d(63, BoundaryType_Periodic)

    print *, "test_fftw31d dirichlet"
    call test_fftw31d(63, BoundaryType_Dirichlet)

    print *, "test_fftw31d neumman"
    call test_fftw31d(63, BoundaryType_Neumann)

    print *, "test_fftw31d dirichlet neumman"
    call test_fftw31d(63, BoundaryType_Dirichlet_Neumann)

    print *, "test_fftw31d neumman dirichlet"
    call test_fftw31d(63, BoundaryType_Neumann_Dirichlet)

contains

    subroutine test_fftw31d(n, boundary_type)
        integer, intent(in):: n
        integer, intent(in) :: boundary_type

        double precision, allocatable :: original(:), forwarded(:), backwarded(:)
        integer :: i

        type(t_FFTW3Executor1d) :: fft
        double precision :: diff

        allocate (original(n))
        allocate (forwarded(n))
        allocate (backwarded(n))

        do i = 1, n
            original(i) = sin(1.0d0*i/n*2*3.14)
        end do

        fft%boundary_type = boundary_type

        fft = new_FFTW3Executor1d(n, boundary_type)

        call fft%forward(original, forwarded)
        call fft%backward(forwarded, backwarded)

        diff = sum(abs(backwarded(:) - original(:))) / size(original)

        call fft%destroy

        print *, "Diff =", diff
    end subroutine

end program
