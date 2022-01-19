module turbulence 
    use iso_fortran_env
    use mpi 
    use fftwfftc
    implicit none
    real(real32),parameter,private::pi=3.141592653589793_real32

    public::generate_wavenumbers
    public::generate_isotropic
    private::vonkarman_2D
    private::vonkarman_2Dscalar
    public::spectral_parameters
    private::rpm 
    private::rpm_scalar
    private::transpose2D
    private::root_mean_square
    private::random_uniform
    private::random_stduniform

    contains

    function generate_wavenumbers(x, vp, t_type)
        !***********************************************************************
        ! Purpose:                                                             *
        !       Given the domain of a function g(x), obtain the wavenumber     *
        !       bins corresponding to x. This function can generate the        *
        !       bins for 9 different types of transforms                       *
        !                                                                      *
        ! Input: ___________________________________________________________   *
        !       |___Variable____|_______________Description_________________|  *
        !       |       x       |   Domain of some function g(x)            |  *
        !       |      vp       |   Phase speed (c0 for time) (optional)    |  *
        !       |  t_type       |   Type of transform (FFT, DST1, etc.)     |  *
        !       |_______________|___________________________________________|  *
        !                                                                      *
        ! Output:___________________________________________________________   *
        !       |___Variable____|_______________Description_________________|  *
        !       | fourier_space |   Fourier domain (frequency vector)       |  *
        !       |_______________|___________________________________________|  *
        !***********************************************************************
        real(real32), intent(in)    :: x(:)
        real(real32)                :: generate_wavenumbers(size(x,dim=1))
        real(real32), optional      :: vp 
        character(len=*), optional  :: t_type
        ! Local variables ******************************************************
        real(real32)                :: dx
        integer(int32)              :: i, j, n, m
        real(real32)                :: c                                       ! Phase speed 
        character(len=180)          :: transform_type                          ! Type of transform
        !***********************************************************************

        if ( present(vp) ) then
            c = vp
        else
            c = 1.0
        end if 

        if ( present(t_type) ) then
            transform_type = trim(t_type)
        else
            transform_type = "FFT"
        end if

        n = size(x,dim=1)
        dx = x(2) - x(1)

        if ( trim(transform_type) .eq. "FFT" ) then 
            if ( modulo(n,2) .eq. 0 ) then
                i = 0
                do j = -n/2,n/2 - 1
                    i = i+1
                    generate_wavenumbers(i) = (2.0*pi)/(n*dx*c) * (1.0*j)
                end do
                generate_wavenumbers = cshift(generate_wavenumbers,n/2)
            else
                i = 0
                do j = -(n-1)/2,(n-1)/2
                    i = i+1
                    generate_wavenumbers(i) = (2.0*pi)/(n*dx*c) * (1.0*j)
                end do
                generate_wavenumbers = cshift(generate_wavenumbers,(n-1)/2)
            end if
        else if ( trim(transform_type) .eq. "DCT1" ) then
            i = 0
            m = 2*(n-1)
            do j = 0,m/2
                i = i + 1
                generate_wavenumbers(i) = (2.0*pi)/(m*dx*c) * (1.0*j)
            end do
        else if ( trim(transform_type) .eq. "DCT2" ) then
            i = 0
            m = 2*n
            do j = 0,m/2 - 1
                i = i + 1
                generate_wavenumbers(i) = (2.0*pi)/(m*dx*c) * (1.0*j)
            end do
        else if ( trim(transform_type) .eq. "DCT3" ) then
            i = 0
            m = 2*n
            do j = 0,m/2 - 1
                i = i + 1
                generate_wavenumbers(i) = (2.0*pi)/(m*dx*c) * (1.0*j + 0.5)
            end do
        else if ( trim(transform_type) .eq. "DCT4" ) then
            i = 0
            m = 2*n
            do j = 0,m/2 - 1
                i = i + 1
                generate_wavenumbers(i) = (2.0*pi)/(m*dx*c) * (1.0*j + 0.5)
            end do
        else if ( trim(transform_type) .eq. "DST1" ) then
            i = 0
            m = 2*(n+1)
            do j = 1,m/2 - 1
                i = i + 1
                generate_wavenumbers(i) = (2.0*pi)/(m*dx*c) * (1.0*j)
            end do
        else if ( trim(transform_type) .eq. "DST2" ) then
            i = 0
            m = 2*n
            do j = 1,m/2
                i = i + 1
                generate_wavenumbers(i) = (2.0*pi)/(m*dx*c) * (1.0*j)
            end do
        else if ( trim(transform_type) .eq. "DST3" ) then
            i = 0
            m = 2*n
            do j = 0,m/2 - 1
                i = i + 1
                generate_wavenumbers(i) = (2.0*pi)/(m*dx*c) * (1.0*j + 0.5)
            end do
        else if ( trim(transform_type) .eq. "DST4" ) then
            i = 0
            m = 2*n
            do j = 0,m/2 - 1
                i = i + 1
                generate_wavenumbers(i) = (2.0*pi)/(m*dx*c) * (1.0*j + 0.5)
            end do
        else if ( trim(transform_type) .eq. "mixed" ) then
            i = 0
            m = 2*(n+1)
            do j = 1,m/2 - 1
                i = i + 1
                generate_wavenumbers(i) = (2.0*pi)/(m*dx*c) * (1.0*j)
            end do
        else
            if ( modulo(n,2) .eq. 0 ) then
                i = 0
                do j = -n/2,n/2 - 1
                    i = i+1
                    generate_wavenumbers(i) = (2.0*pi)/(n*dx*c) * (1.0*j)
                end do
                generate_wavenumbers = cshift(generate_wavenumbers,n/2)
            else
                i = 0
                do j = -(n-1)/2,(n-1)/2
                    i = i+1
                    generate_wavenumbers(i) = (2.0*pi)/(n*dx*c) * (1.0*j)
                end do
                generate_wavenumbers = cshift(generate_wavenumbers,(n-1)/2)
            end if
        end if

    end function generate_wavenumbers

    ! Including density fluctuations and speed of sound flucuations
    subroutine generate_isotropic(u1, u2, temperature, q, kx1, kx2, sigmaV, LV, sigmaT, LT, sigmaQ, LQ, mode, comm)
        real(real32), intent(out)   :: u1(:,:), u2(:,:), temperature(:,:), q(:,:)
        real(real32), intent(in)    :: kx1(:), kx2(:), sigmaV, LV, sigmaT, LT, sigmaQ, LQ
        character(len=*), optional  :: mode
        integer(int32), optional    :: comm 
        ! Local 
        character(len=180)          :: turb_mode
        real(real32)                :: ell_V, ell_T, ell_Q
        integer(int32)              :: nproc, proc, mpi_err

        if ( present(comm) ) then
            call MPI_COMM_SIZE(comm, nproc, mpi_err)
            call MPI_COMM_RANK(comm, proc, mpi_err)
        else 
            proc = 0
            nproc = 1
        end if

        if ( present(mode) ) then
            turb_mode = trim(mode)
        else
            turb_mode = "all"
        end if

        ! Convert Integral length scales to von Karman length scales (per Frehlich 2000)
        ell_V = (1.0/0.746834)*LV 
        ell_T = (1.0/0.746834)*LT 
        ell_Q = (1.0/0.746834)*LQ

        if ( turb_mode .eq. "all" ) then
            if ( present(comm) ) then
                call rpm(u1, kx1, kx2, sigmaV, ell_V, 1, comm)
                call rpm(u2, kx1, kx2, sigmaV, ell_V, 2, comm)
                call rpm_scalar(temperature, kx1, kx2, sigmaT, ell_T, comm)
                call rpm_scalar(q, kx1, kx2, sigmaQ, ell_Q, comm)
            else 
                call rpm(u1, kx1, kx2, sigmaV, ell_V, 1)
                call rpm(u2, kx1, kx2, sigmaV, ell_V, 2)
                call rpm_scalar(temperature, kx1, kx2, sigmaT, ell_T)
                call rpm_scalar(q, kx1, kx2, sigmaQ, ell_Q)
            end if
        else if ( turb_mode .eq. "velocity" ) then
            if ( present(comm) ) then
                call rpm(u1, kx1, kx2, sigmaV, ell_V, 1, comm)
                call rpm(u2, kx1, kx2, sigmaV, ell_V, 2, comm)
            else 
                call rpm(u1, kx1, kx2, sigmaV, ell_V, 1)
                call rpm(u2, kx1, kx2, sigmaV, ell_V, 2)
            end if
        else if ( turb_mode .eq. "temperature" ) then
            if ( present(comm) ) then
                call rpm_scalar(temperature, kx1, kx2, sigmaT, ell_T, comm)
            else 
                call rpm_scalar(temperature, kx1, kx2, sigmaT, ell_T)
            end if
        else if ( turb_mode .eq. "velocityandtemperature" ) then
            if ( present(comm) ) then
                call rpm(u1, kx1, kx2, sigmaV, ell_V, 1, comm)
                call rpm(u2, kx1, kx2, sigmaV, ell_V, 2, comm)
                call rpm_scalar(temperature, kx1, kx2, sigmaT, ell_T, comm)
            else 
                call rpm(u1, kx1, kx2, sigmaV, ell_V, 1)
                call rpm(u2, kx1, kx2, sigmaV, ell_V, 2)
                call rpm_scalar(temperature, kx1, kx2, sigmaT, ell_T)
            end if
        else if ( turb_mode .eq. "velocityandhumidity" ) then
            if ( present(comm) ) then
                call rpm(u1, kx1, kx2, sigmaV, ell_V, 1, comm)
                call rpm(u2, kx1, kx2, sigmaV, ell_V, 2, comm)
                call rpm_scalar(q, kx1, kx2, sigmaQ, ell_Q, comm)
            else 
                call rpm(u1, kx1, kx2, sigmaV, ell_V, 1)
                call rpm(u2, kx1, kx2, sigmaV, ell_V, 2)
                call rpm_scalar(q, kx1, kx2, sigmaQ, ell_Q)
            end if
        else if ( turb_mode .eq. "humidity" ) then
            if ( present(comm) ) then
                call rpm_scalar(q, kx1, kx2, sigmaQ, ell_Q, comm)
            else 
                call rpm_scalar(q, kx1, kx2, sigmaQ, ell_Q)
            end if
        else if ( turb_mode .eq. "temperatureandhumidity" ) then
            if ( present(comm) ) then
                call rpm_scalar(temperature, kx1, kx2, sigmaT, ell_T, comm)
                call rpm_scalar(q, kx1, kx2, sigmaQ, ell_Q, comm)
            else 
                call rpm_scalar(temperature, kx1, kx2, sigmaT, ell_T)
                call rpm_scalar(q, kx1, kx2, sigmaQ, ell_Q)
            end if
        else
            if ( present(comm) ) then
                call rpm(u1, kx1, kx2, sigmaV, ell_V, 1, comm)
                call rpm(u2, kx1, kx2, sigmaV, ell_V, 2, comm)
                call rpm_scalar(temperature, kx1, kx2, sigmaT, ell_T, comm)
                call rpm_scalar(q, kx1, kx2, sigmaQ, ell_Q, comm)
            else 
                call rpm(u1, kx1, kx2, sigmaV, ell_V, 1)
                call rpm(u2, kx1, kx2, sigmaV, ell_V, 2)
                call rpm_scalar(temperature, kx1, kx2, sigmaT, ell_T)
                call rpm_scalar(q, kx1, kx2, sigmaQ, ell_Q)
            end if
        end if

    end subroutine generate_isotropic

    ! Private
    subroutine vonkarman_2D(phi, kx1, kx2, L, index1, index2, comm)
        real(real32), intent(out)   :: phi(:,:)
        real(real32), intent(in)    :: kx1(:), kx2(:), L 
        integer(int32), optional    :: index1, index2 
        integer(int32), optional    :: comm
        ! Local 
        integer(int32)              :: direction_1, direction_2, nproc, proc, mpi_err, n, m, i, j, global_ind
        real(real32)                :: kh

        if ( present(index1) ) then
            direction_1 = index1 
        else
            direction_1 = 1 
        end if

        if ( present(index2) ) then
            direction_2 = index2 
        else
            direction_2 = 1 
        end if

        if ( present(comm) ) then
            call MPI_COMM_SIZE(comm, nproc, mpi_err)
            call MPI_COMM_RANK(comm, proc, mpi_err)
        else 
            proc = 0
            nproc = 1
        end if

        n = size(phi,dim=1)
        m = size(phi,dim=2)

        select case (direction_1)
            case (1) 
                select case (direction_2)
                    case (1) 
                        ! Phi_11 
                        do j = 1,m 
                            global_ind = j + proc*m 
                            do i = 1,n 
                                kh = sqrt( kx1(i)**2 + kx2(global_ind)**2 )
                                phi(i,j) = ( (L**2)/(6.0*pi*( (1.0 + (kh**2 * L**2))**(4.0/3.0) ) ) )*(1.0 + &
                                    (8.0/3.0)*( (kx2(global_ind)**2 * L**2)/(1.0 + (kh**2 * L**2) ) ))
                                ! phi(i,j) = ((L**2)/(3.0*pi*((1.0 + ( kh**2 * L**2))**(7.0/3.0))))*( (4.0/3.0)*(kx1(i)**2 * L**2) &
                                !     + ((1.0 + (kh**2 * L**2))/2.0) )
                            end do
                        end do
                    case (2) 
                        ! Phi_12
                        do j = 1,m 
                            global_ind = j + proc*m 
                            do i = 1,n 
                                kh = sqrt( kx1(i)**2 + kx2(global_ind)**2 )
                                phi(i,j) = (-4.0 * kx1(i) * kx2(global_ind) * L**4)/(9.0*pi*((1.0 + (kh**2 * L**2))**(7.0/3.0)))
                            end do
                        end do
                end select
            case (2)
                select case (direction_2)
                    case (1) 
                        ! Phi_21
                        do j = 1,m 
                            global_ind = j + proc*m 
                            do i = 1,n 
                                kh = sqrt( kx1(i)**2 + kx2(global_ind)**2 )
                                phi(i,j) = (-4.0 * kx1(i) * kx2(global_ind) * L**4)/(9.0*pi*((1.0 + (kh**2 * L**2))**(7.0/3.0)))
                            end do
                        end do
                    case (2) 
                        ! Phi_22
                        do j = 1,m 
                            global_ind = j + proc*m 
                            do i = 1,n 
                                kh = sqrt( kx1(i)**2 + kx2(global_ind)**2 )
                                phi(i,j) = ( (L**2)/(6.0*pi*( (1.0 + (kh**2 * L**2))**(4.0/3.0) ) ) )*(1.0 + &
                                    (8.0/3.0)*( (kx1(i)**2 * L**2)/(1.0 + (kh**2 * L**2) ) ))
                                ! phi(i,j) = ((L**2)/(3.0*pi*((1.0 + ( kh**2 * L**2))**(7.0/3.0))))*( &
                                !     (4.0/3.0)*(kx2(global_ind)**2 * L**2) + ((1.0 + (kh**2 * L**2))/2.0) )
                            end do
                        end do
                end select
        end select

    end subroutine vonkarman_2D

    subroutine vonkarman_2Dscalar(phi, kx1, kx2, L, comm)
        real(real32), intent(out)   :: phi(:,:)
        real(real32), intent(in)    :: kx1(:), kx2(:), L 
        integer(int32), optional    :: comm
        ! Local 
        integer(int32)              :: nproc, proc, mpi_err, n, m, i, j, global_ind
        real(real32)                :: kh

        if ( present(comm) ) then
            call MPI_COMM_SIZE(comm, nproc, mpi_err)
            call MPI_COMM_RANK(comm, proc, mpi_err)
        else 
            proc = 0
            nproc = 1
        end if

        n = size(phi,dim=1)
        m = size(phi,dim=2)

        do j = 1,m 
            global_ind = j + proc*m 
            do i = 1,n 
                kh = sqrt( kx1(i)**2 + kx2(global_ind)**2 )
                phi(i,j) = ( (L**2)/(6.0*pi*( (1.0 + (kh**2 * L**2))**(4.0/3.0) ) ) )*(1.0 + &
                    (8.0/3.0)*( (kh**2 * L**2)/(1.0 + (kh**2 * L**2) ) ))
            end do
        end do

    end subroutine vonkarman_2Dscalar 

    subroutine spectral_parameters(sigVs, LVs, sigVb, LVb, sigT, LT, z, ustar, wstar, Tstar, L_Obukhov, zi)
        real(real32), intent(out)   :: sigVs(:), LVs(:), sigVb(:), LVb(:), sigT(:), LT(:)
        real(real32), intent(in)    :: ustar, wstar, z(:), zi, Tstar, L_obukhov 
        integer(int32)              :: i, n

        n = size(z,dim=1)
        do i = 1,n
            if ( z(i) .gt. 200.0 ) then
                ! Above 200 m, shear and temperature scales should be kept the same as they are at z = 200 m
                sigVs(i) = sqrt(3.0)*ustar
                LVs(i)   = 1.8*200.0
                sigT(i)  = sqrt( ( 4.0*(Tstar**2) )/( (1.0 - 10.0*(200.0/L_obukhov))**(2.0/3.0) ) )
                LT(i)    = 2.0*200.0*( (1.0 - 7.0*(200.0/L_obukhov) )/(1.0 - 10.0*(200.0/L_obukhov))  )
            else
                sigVs(i) = sqrt(3.0)*ustar
                LVs(i)   = 1.8*z(i) 
                sigT(i)  = sqrt( ( 4.0*(Tstar**2) )/( (1.0 - 10.0*(z(i)/L_obukhov))**(2.0/3.0) ) )
                LT(i)    = 2.0*z(i)*( (1.0 - 7.0*(z(i)/L_obukhov) )/(1.0 - 10.0*(z(i)/L_obukhov))  )
            end if
            sigVb(i) = sqrt(0.35)*wstar
            LVb(i)   = 0.23*zi 
        end do

    end subroutine spectral_parameters

    subroutine rpm(u, kx1, kx2, sigma, L, velocity_index, comm)
        real(real32), intent(out)   :: u(:,:) 
        real(real32), intent(in)    :: kx1(:), kx2(:), sigma, L 
        integer(int32), optional    :: velocity_index, comm 
        ! Local 
        integer(int32)              :: n, m, i, j, nproc, proc, mpi_err, direction
        real(real32)                :: g1(1), g2(1), mean_proc, std, std_proc, mean
        real(real32), allocatable   :: phi(:,:), usquare(:,:)
        complex(real32), allocatable    :: s(:,:), st(:,:)
        complex(real32)             :: random_phase

        if ( present(comm) ) then
            call MPI_COMM_SIZE(comm, nproc, mpi_err)
            call MPI_COMM_RANK(comm, proc, mpi_err)
        else 
            proc = 0
            nproc = 1
        end if

        if ( present(velocity_index) ) then
            direction = velocity_index 
        else
            direction = 1 
        end if

        n = size(u,dim=1)
        m = size(u,dim=2)

        ! allocate 
        allocate ( phi(n,m), s(n,m), st(m*nproc,n/nproc), usquare(n,m) )

        ! Von Karman vector model 
        if ( present(comm) ) then
            call vonkarman_2D(phi, kx1, kx2, L, direction, direction, comm)
        else
            call vonkarman_2D(phi, kx1, kx2, L, direction, direction)
        end if

        ! Random phase method for homogeneous isotropic turbulence 
        call random_seed()
        do j = 1,m 
            do i = 1,n 
                call random_uniform(-1.73205,1.73205,g1)
                call random_uniform(-1.73205,1.73205,g2)
                random_phase = cmplx(g1(1), g2(1))
                s(i,j) = sigma*sqrt(phi(i,j))*random_phase*(1.0*n*(m*nproc)) ! Last two multiplications are scaling factors for FFT
            end do
        end do

        ! Fast Fourier transform to velocity
        do j = 1,m 
            call fftw_fft_1d(s(:,j), inverse=.true.)
        end do

        ! transpose
        if ( present(comm) ) then
            call transpose2D(st, s, comm)
        else
            call transpose2D(st, s)
        end if

        ! Fast Fourier transform
        do i = 1,n/nproc 
            call fftw_fft_1d(st(:,i), inverse=.true.)
        end do

        ! transpose
        if ( present(comm) ) then
            call transpose2D(s, st, comm)
        else
            call transpose2D(s, st)
        end if

        ! Velocity 
        u = real(s)

        ! Set to correct variance
        if ( present(comm) ) then
            if ( sigma .eq. 0.0 ) then
                ! Skip
                continue
            else
                mean_proc = sum(u)/(1.0*n*m*nproc)
                call MPI_Allreduce(mean_proc, mean, 1, MPI_REAL, MPI_SUM, comm, mpi_err)
                ! Find squared differences
                do j = 1,m 
                    do i = 1,n 
                        usquare(i,j) = (u(i,j) - mean)**2
                    end do
                end do
                std_proc = sum(usquare)/(1.0*n*m*nproc)
                call MPI_Allreduce(std_proc, std, 1, MPI_REAL, MPI_SUM, comm, mpi_err)
                std = sqrt(std) 

                u = (sigma/std)*u
            end if
        else
            if ( sigma .eq. 0.0 ) then
                ! Skip
                continue
            else
                mean = sum(u)/(1.0*n*m*nproc)
                do j = 1,m 
                    do i = 1,n 
                        usquare(i,j) = (u(i,j) - mean)**2
                    end do
                end do
                std = sum(usquare)/(1.0*n*m*nproc)
                std = sqrt(std) 

                u = (sigma/std)*u
            end if
        end if

        ! deallocate 
        deallocate ( phi, s, st, usquare )

    end subroutine rpm
        
    subroutine rpm_scalar(u, kx1, kx2, sigma, L, comm)
        real(real32), intent(out)   :: u(:,:) 
        real(real32), intent(in)    :: kx1(:), kx2(:), sigma, L 
        integer(int32), optional    :: comm 
        ! Local 
        integer(int32)              :: n, m, i, j, nproc, proc, mpi_err
        real(real32)                :: g1(1), g2(1), mean_proc, std, std_proc, mean
        real(real32), allocatable   :: phi(:,:), usquare(:,:)
        complex(real32), allocatable    :: s(:,:), st(:,:)
        complex(real32)             :: random_phase

        if ( present(comm) ) then
            call MPI_COMM_SIZE(comm, nproc, mpi_err)
            call MPI_COMM_RANK(comm, proc, mpi_err)
        else 
            proc = 0
            nproc = 1
        end if

        n = size(u,dim=1)
        m = size(u,dim=2)

        ! allocate 
        allocate ( phi(n,m), s(n,m), st(m*nproc,n/nproc), usquare(n,m) )

        ! Von Karman vector model 
        if ( present(comm) ) then
            call vonkarman_2Dscalar(phi, kx1, kx2, L, comm)
        else
            call vonkarman_2Dscalar(phi, kx1, kx2, L)
        end if

        ! Random phase method for homogeneous isotropic turbulence 
        call random_seed()
        do j = 1,m 
            do i = 1,n 
                call random_uniform(-1.73205,1.73205,g1)
                call random_uniform(-1.73205,1.73205,g2)
                random_phase = cmplx(g1(1), g2(1))
                s(i,j) = sigma*sqrt(phi(i,j))*random_phase*(1.0*n*(m*nproc)) ! Last two multiplications are scaling factors for FFT
            end do
        end do

        ! Fast Fourier transform to velocity
        do j = 1,m 
            call fftw_fft_1d(s(:,j), inverse=.true.)
        end do

        ! transpose
        if ( present(comm) ) then
            call transpose2D(st, s, comm)
        else
            call transpose2D(st, s)
        end if

        ! Fast Fourier transform
        do i = 1,n/nproc 
            call fftw_fft_1d(st(:,i), inverse=.true.)
        end do

        ! transpose
        if ( present(comm) ) then
            call transpose2D(s, st, comm)
        else
            call transpose2D(s, st)
        end if

        ! scalar
        u = real(s)

        ! Set to correct variance
        if ( present(comm) ) then
            if ( sigma .eq. 0.0 ) then
                ! Skip
                continue
            else
                mean_proc = sum(u)/(1.0*n*m*nproc)
                call MPI_Allreduce(mean_proc, mean, 1, MPI_REAL, MPI_SUM, comm, mpi_err)
                ! Find squared differences
                do j = 1,m 
                    do i = 1,n 
                        usquare(i,j) = (u(i,j) - mean)**2
                    end do
                end do
                std_proc = sum(usquare)/(1.0*n*m*nproc)
                call MPI_Allreduce(std_proc, std, 1, MPI_REAL, MPI_SUM, comm, mpi_err)
                std = sqrt(std) 

                u = (sigma/std)*u
            end if
        else
            if ( sigma .eq. 0.0 ) then
                ! Skip
                continue
            else
                mean = sum(u)/(1.0*n*m*nproc)
                do j = 1,m 
                    do i = 1,n 
                        usquare(i,j) = (u(i,j) - mean)**2
                    end do
                end do
                std = sum(usquare)/(1.0*n*m*nproc)
                std = sqrt(std) 

                u = (sigma/std)*u
            end if
        end if

        ! deallocate 
        deallocate ( phi, s, st, usquare )

    end subroutine rpm_scalar

    subroutine transpose2D(out, in, comm)
        complex(real32), intent(out)    :: out(:,:) 
        complex(real32), intent(in)     :: in(:,:)
        integer(int32), optional        :: comm 
        ! Local 
        integer(int32)                  :: i, j, n, m, nproc, proc, mpi_err
        integer(int32), allocatable     :: counts_send(:), disp_send(:), counts_recv(:)
        integer(int32), allocatable     :: disp_recv(:), types_send(:), types_recv(:)
        integer(int32)                  :: type_send, transpose_case

        if ( present(comm) ) then
            call MPI_COMM_SIZE(comm, nproc, mpi_err)
            call MPI_COMM_RANK(comm, proc, mpi_err)
            transpose_case = 1
        else 
            proc = 0
            nproc = 1
            transpose_case = 2
        end if

        n = size(in,dim=1)
        m = size(in,dim=2)

        select case (transpose_case)
            case (1)
                allocate ( counts_send(nproc), disp_send(nproc), counts_recv(nproc) )
                allocate ( disp_recv(nproc), types_send(nproc), types_recv(nproc) )

                ! Create datatype that has appropriate strides
                call MPI_Type_vector(m, 1, n, MPI_COMPLEX, type_send, mpi_err)
                call MPI_Type_commit(type_send, mpi_err)

                do j = 1,n/nproc
                    ! Setup transpose
                    do i = 1,nproc 
                        counts_send(i) = 1 
                        disp_send(i) = (i-1)*(n/nproc)*int(sizeof(in(1,1)), kind=int32) + (j-1)*int(sizeof(in(1,1)), kind=int32)
                        types_send(i) = type_send
                        counts_recv(i) = m
                        disp_recv(i) = (i-1)*m*int(sizeof(in(1,1)), kind=int32) + (j-1)*(m*nproc)*int(sizeof(in(1,1)), kind=int32)
                        types_recv(i) = MPI_COMPLEX
                    end do

                    ! Perform all to all communication
                    call MPI_Alltoallw(in, counts_send, disp_send, types_send, out, counts_recv, disp_recv, types_recv, &
                        comm, mpi_err)
                end do
                
                ! Cleanup
                deallocate ( counts_send, disp_send, counts_recv )
                deallocate ( disp_recv, types_send, types_recv )
            case (2)
                out = transpose(in)
        end select

    end subroutine transpose2D

    function root_mean_square(f)
        real(real32), intent(in)    :: f(:)
        real(real32)                :: root_mean_square 
        ! local 
        integer(int32)              :: n

        n = size(f,dim=1)
        root_mean_square = sqrt( real(1.0/n, kind=real32)*sum(f**2) )
    end function root_mean_square

    subroutine random_uniform(a,b,x)
        ! Return a uniform random number x, where x varies
        ! from a < u <= b
        real(real32), intent(in)    :: a,b
        real(real32), intent(out)   :: x(:)
        real(real32), allocatable   :: u(:)
        allocate ( u( size(x,dim=1) ) )
        call random_stduniform(u)
        x = (b-a)*u + a
    end subroutine random_uniform

    subroutine random_stduniform(u)
        ! Return a uniform random number u, where u varies
        ! from 0 < u <= 1
        real(real32), intent(out)   :: u(:)
        real(real32), allocatable   :: r(:)
        allocate ( r(size(u,dim=1)) )
        call random_number(r)
        u = 1 - r
    end subroutine random_stduniform

end module turbulence