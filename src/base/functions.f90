module functions
    use mpi 
    use hdf5
    use mkl_dfti
    use special_functions
    implicit none
    integer, parameter          :: sp = selected_real_kind(6, 37)
    integer, parameter          :: dp = selected_real_kind(15, 307)
    integer, parameter          :: isp = selected_int_kind(6)
    integer, parameter          :: idp = selected_int_kind(10)
    real(kind=sp), parameter    :: pi = 3.141592653_sp
    complex(kind=sp), parameter :: img = (0.0_sp,1.0_sp)


    contains

    subroutine gaussian_boxmuller(y, u1, u2, a, s)
        real, intent(out)   :: y
        real, intent(in)    :: u1, u2, a, s
        real                :: X

        x = sqrt( -2.0*log(u1) )*sin(2.0*pi*u2)
        y = s*x + a

    end subroutine gaussian_boxmuller

    subroutine vonKarman_velocity(E, k, L, sigma)
        real, intent(out)   :: E
        real, intent(in)    :: k
        real, intent(in)    :: L, sigma
        real                :: a, b, c

        a = 1.12878702991
        b = 2.67893853470
        c = (55.0*a)/(36.0*(pi**(1.5))*b)

        E = c*(sigma**2 * k**2 * L**5)/( ( 1.0 + (k**2 * L**2) )**(17.0/6.0))

    end subroutine vonKarman_velocity

    subroutine vonKarman_temperature(phi, k, L, sigma)
      implicit none
      real, intent(out) :: phi
      real, intent(in) :: k
      real, intent(in) :: L, sigma
      real :: a, b
  
      a = 0.94065585825
      b = 2.67893853471 
  
      phi = (a/(pi**(3.0/2.0) * b))*( (sigma**2 * L**3)/((1.0 &
        + (k**2 * L**2) )**(11.0/6.0)) )
  
    end subroutine vonKarman_temperature

    subroutine velocitySpectrum(phi, ki, sigma, L, n, m)
        real, intent(out)   :: phi
        real, intent(in)    :: ki(:)
        real, intent(in)    :: sigma, L
        integer, intent(in) :: n, m
        real                :: E, k


        k = sqrt( ki(1)**2 + ki(2)**2 + ki(3)**2 )
        call vonKarman_velocity(E, k, L, sigma)
        phi = E*( (k**2)*kronecker(n, m) - ki(n)*ki(m) )

    end subroutine velocitySpectrum

    function kronecker(n,m)
        integer :: n, m
        integer :: kronecker

        if ( n .eq. m ) then
          kronecker = 1
        else
          kronecker = 0
        end if

    end function kronecker

    ! subroutine fft(w, array_size, dim, inverse)
    !     ! ========================================================================
    !     ! Note: The forward Fourier transform in Frelich has the scaling factor
    !     !
    !     ! Note: This subroutine follows the same convention
    !     ! ========================================================================

    !     ! Complex to complex fft
    !     complex, intent(inout)          :: w(:)
    !     integer, intent(in)             :: array_size(:)
    !     integer, intent(in)             :: dim
    !     logical, intent(in)             :: inverse
    !     type(DFTI_DESCRIPTOR), pointer  :: My_Desc_Handle
    !     integer(isp)                    :: status = 0_isp, ignored_status

    !     ! Perform a complex to complex transform
    !     status = DftiCreateDescriptor( My_Desc_Handle, DFTI_SINGLE,&
    !       DFTI_COMPLEX, dim, array_size )
    !     if (status .ne. 0) goto 999
    !     status = DftiSetValue(My_Desc_Handle, DFTI_FORWARD_SCALE, &
    !       1.0/( 1.0*product( array_size ) ) )
    !     if (status .ne. 0) goto 999
    !     status = DftiCommitDescriptor( My_Desc_Handle )
    !     if (status .ne. 0) goto 999
    !     if ( .not. inverse ) then
    !       status = DftiComputeForward( My_Desc_Handle, w )
    !       if (status .ne. 0) goto 999
    !     else
    !       status = DftiComputeBackward( My_Desc_Handle, w )
    !       if (status .ne. 0) goto 999
    !     end if
    !     100 continue
    !     ignored_status = DftiFreeDescriptor(My_Desc_Handle)
    !     goto 200

    !     999 continue
    !     print '(" Error, status = ",I0)', status

    !     goto 100

    !     200 continue

    ! end subroutine fft

    ! subroutine construct_wavenumbers(kx, ky, kz, x, y, z)
    !     real(sp), intent(out)       :: kx(:), ky(:), kz(:)
    !     real(sp), intent(in)        :: x(:), y(:), z(:)
    !     real(sp)                    :: fxs, fys, fzs
    !     real(sp)                    :: dfx, dfy, dfz
    !     real(sp), allocatable       :: fx(:), fy(:), fz(:)
    !     integer                     :: i

    !     allocate ( fx(size(x,dim=1)), fy(size(y,dim=1)), fz(size(z,dim=1)) )

    !     ! Sampling frequencies ===================================================
    !     fxs = 1.0*size(x)/( 1.0_sp*abs( x(size(x)) - x(1) ) )
    !     fys = 1.0*size(y)/( 1.0_sp*abs( y(size(y)) - y(1) ) )
    !     fzs = 1.0*size(z)/( 1.0_sp*abs( z(size(z)) - z(1) ) )
    !     ! ========================================================================
    !     ! Bin Widths =============================================================
    !     dfx = fxs/(size(x)*1.0_sp)
    !     dfy = fys/(size(y)*1.0_sp)
    !     dfz = fzs/(size(z)*1.0_sp)
    !     ! ========================================================================
    !     ! Construct frequencies ==================================================
    !     do i = 1,size(x)
    !       if (i .le. int(size(x)/2.0) + 1) then
    !         fx(i) = dfx*real(i-1) ! dft*real(i-1)
    !       else
    !         fx(i) = fx(2*(int(size(x)/2.0) + 1) - i)
    !       end if
    !     end do
    !     do i = 1,size(y)
    !       if (i .le. int(size(y)/2.0) + 1) then
    !         fy(i) = dfy*real(i-1) ! dfy*real(i-1)
    !       else
    !         fy(i) = fy(2*(int(size(y)/2.0) + 1) - i)
    !       end if
    !     end do
    !     do i = 1,size(z)
    !       if (i .le. int(size(z)/2.0) + 1) then
    !         fz(i) = dfz*real(i-1) ! dfz*real(i-1)
    !       else
    !         fz(i) = fz(2*(int(size(z)/2.0) + 1) - i)
    !       end if
    !     end do
    !     ! ========================================================================

    !     ! Wavenumbers ============================================================
    !     kx = 2.0*pi*fx
    !     ky = 2.0*pi*fy
    !     kz = 2.0*pi*fz
    !     ! ========================================================================

    !     deallocate ( fx, fy, fz )

    ! end subroutine construct_wavenumbers

    ! subroutine construct_w(w1, w2, w3, x, y, z, kx, ky, kz, mean, sigma, L,type)
    !     complex(sp), intent(out)      :: w1(:,:,:), w2(:,:,:), w3(:,:,:)
    !     real(sp), intent(in)          :: kx(:), ky(:), kz(:)
    !     real(sp), intent(in)          :: x(:), y(:), z(:)
    !     real(sp), intent(in)          :: mean, sigma, L
    !     character(len=*), intent(in)  :: type
    !     real(sp)                      :: a1(1), b1(1), a2(1), b2(1), a(1), b(1)
    !     real(sp)                      :: h11, h12, h22, h13, h23, h33
    !     real(sp)                      :: phi11, phi12, phi22, phi13, phi23, phi33
    !     complex(sp)                   :: N1, N2, N3
    !     real(sp)                      :: dkx, dky, dkz, G, Gd, r
    !     real(dp)                      :: vj1, vj2, vy1, vy2, vi1, vi2, vk1, vk2
    !     integer(isp)                  :: i, j, k
    !     complex(sp), allocatable      :: B11(:,:,:), B22(:,:,:), B33(:,:,:)

    !     dkx   = kx(2) - kx(1)
    !     dky   = ky(2) - ky(1)
    !     dkz   = kz(2) - kz(1)

    !     allocate ( B11(size(w1,dim=1),size(w1,dim=2),size(w1,dim=3)) )
    !     allocate ( B22(size(w2,dim=1),size(w2,dim=2),size(w2,dim=3)) )
    !     allocate ( B33(size(w3,dim=1),size(w3,dim=2),size(w3,dim=3)) )

    !     if ( type .eq. "gaussian" ) then
    !       ! Loop through wavenumbers with gaussian phase correlation =============
    !       call random_seed()
    !       do k = 1,size(w1,dim=3)
    !         do j = 1,size(w1,dim=2)
    !           do i = 1,size(w1,dim=1)
    !             ! Random vectors
    !             call random_number(a1)
    !             call random_number(b1)
    !             call gaussian_boxmuller(a, a1, b1, 0.0, 1.0)
    !             call random_number(a2)
    !             call random_number(b2)
    !             call gaussian_boxmuller(b, a2, b2, 0.0, 1.0)
    !             N1 = a(1) + img*b(1)
    !             call random_number(a1)
    !             call random_number(b1)
    !             call gaussian_boxmuller(a, a1, b1, 0.0, 1.0)
    !             call random_number(a2)
    !             call random_number(b2)
    !             call gaussian_boxmuller(b, a2, b2, 0.0, 1.0)
    !             N2 = a(1) + img*b(1)
    !             call random_number(a1)
    !             call random_number(b1)
    !             call gaussian_boxmuller(a, a1, b1, 0.0, 1.0)
    !             call random_number(a2)
    !             call random_number(b2)
    !             call gaussian_boxmuller(b, a2, b2, 0.0, 1.0)
    !             N3 = a(1) + img*b(1)
    !             ! Velocity spectrum
    !             call velocitySpectrum(phi11, kx(i), ky(j), kz(k), sigma, L, 1, 1)
    !             h11 = sqrt( phi11*dkx*dky*dkz )
    !             w1(i,j,k)  = cmplx(h11,0.0)*N1
    !             call velocitySpectrum(phi12, kx(i), ky(j), kz(k), sigma, L, 1, 2)
    !             h12 = phi12*sqrt( dkx*dky*dkz )/sqrt(phi11)
    !             call velocitySpectrum(phi22, kx(i), ky(j), kz(k), sigma, L, 2, 2)
    !             h22 = sqrt( phi22*dkx*dky*dkz  - h12**2 )
    !             w2(i,j,k)  = (cmplx(h12,0.0)*N1 + cmplx(h22,0.0)*N2)
    !             call velocitySpectrum(phi13, kx(i), ky(j), kz(k), sigma, L, 1, 3)
    !             h13 = phi13*sqrt( dkx*dky*dkz )/sqrt(phi11)
    !             call velocitySpectrum(phi23, kx(i), ky(j), kz(k), sigma, L, 2, 3)
    !             h23 = ( phi23*dkx*dky*dkz - h12*h13 )/sqrt(phi22)
    !             call velocitySpectrum(phi33, kx(i), ky(j), kz(k), sigma, L, 3, 3)
    !             h33 = sqrt( phi33*dkx*dky*dkz - (h13**2) - (h23**2) )
    !             w3(i,j,k)  = (cmplx(h13,0.0)*N1 + cmplx(h23,0.0)*N2 &
    !               + cmplx(h33,0.0)*N3)
    !           end do
    !         end do
    !       end do
    !       ! ======================================================================
    !     else if ( type .eq. "vonkarman" ) then
    !       ! Model Correlation function
    !       do k = 1,size(w1,dim=3)
    !         do j = 1,size(w1,dim=2)
    !           do i = 1,size(w1,dim=1)
    !             r = sqrt( x(i)**2 + y(j)**2 + z(k)**2 )
    !             call ajyik ( dble(r/L), vj1, vj2, vy1, vy2, vi1, vi2, vk1, vk2 )
    !             G = 0.592548*(r**(1.0/3.0))*vk1
    !             Gd = 0.29627426*(r**(4.0/3.0))*vk2
    !             B11(i,j,k) = sigma**2 * G
    !             B22(i,j,k) = sigma**2 * ( G - Gd)
    !             B33(i,j,k) = sigma**2 * ( G - Gd)
    !           end do
    !         end do
    !       end do

    !       call fft( B11(:,1,1), [size(B11,dim=1), size(B11,dim=2), &
    !         size(B11,dim=3)], 3, .false. )
    !       call fft( B22(:,1,1), [size(B22,dim=1), size(B22,dim=2), &
    !         size(B22,dim=3)], 3, .false. )
    !       call fft( B33(:,1,1), [size(B33,dim=1), size(B33,dim=2), &
    !         size(B33,dim=3)], 3, .false. )

    !       w1 = B11
    !       w2 = B22
    !       w3 = B33

    !     else if ( type .eq. "uncorrelated" ) then
    !       ! Loop through wavenumbers with gaussian phase correlation =============
    !       call random_seed()
    !       do k = 1,size(w1,dim=3)
    !         do j = 1,size(w1,dim=2)
    !           do i = 1,size(w1,dim=1)
    !             ! Random vectors
    !             call random_number(a)
    !             call random_number(b)
    !             N1 = a(1) + img*b(1)
    !             call random_number(a)
    !             call random_number(b)
    !             N2 = a(1) + img*b(1)
    !             call random_number(a)
    !             call random_number(b)
    !             N3 = a(1) + img*b(1)
    !             ! Velocity spectrum
    !             call velocitySpectrum(phi11, kx(i), ky(j), kz(k), sigma, L, 1, 1)
    !             h11 = sqrt( phi11*dkx*dky*dkz )
    !             w1(i,j,k)  = cmplx(h11,0.0)*N1
    !             call velocitySpectrum(phi12, kx(i), ky(j), kz(k), sigma, L, 1, 2)
    !             h12 = phi12*sqrt( dkx*dky*dkz )/sqrt(phi11)
    !             call velocitySpectrum(phi22, kx(i), ky(j), kz(k), sigma, L, 2, 2)
    !             h22 = sqrt( phi22*dkx*dky*dkz  - h12**2 )
    !             w2(i,j,k)  = (cmplx(h12,0.0)*N1 + cmplx(h22,0.0)*N2)
    !             call velocitySpectrum(phi13, kx(i), ky(j), kz(k), sigma, L, 1, 3)
    !             h13 = phi13*sqrt( dkx*dky*dkz )/sqrt(phi11)
    !             call velocitySpectrum(phi23, kx(i), ky(j), kz(k), sigma, L, 2, 3)
    !             h23 = ( phi23*dkx*dky*dkz - h12*h13 )/sqrt(phi22)
    !             call velocitySpectrum(phi33, kx(i), ky(j), kz(k), sigma, L, 3, 3)
    !             h33 = sqrt( phi33*dkx*dky*dkz - (h13**2) - (h23**2) )
    !             w3(i,j,k)  = (cmplx(h13,0.0)*N1 + cmplx(h23,0.0)*N2 &
    !               + cmplx(h33,0.0)*N3)
    !           end do
    !         end do
    !       end do
    !       ! ======================================================================
    !     else

    !     end if

    !     deallocate ( B11, B22, B33 )

    ! end subroutine construct_w

    subroutine most(u, v, T, h, z, L, z0, zr, Tr, ustar, Tstar, Gammad, theta,&
        qr, q_star)
        !**********************************************************************
        ! Monin-Obukhov Similarity Theory                                     *
        !   Universal functions chosen according to D. K. Wilson (2001)       *
        !   https://link.springer.com/article/10.1023/A:1018718707419         *
        !                                                                     *
        ! Inputs: ___________________________________________________________ *
        !        |__Variable__|_______________Details________________________|*
        !        |     z      | Altitude                                     |*
        !        |     L      | Obukhov length scale                         |*
        !        |    z0      | Roughness height                             |*
        !        |    zr      | Reference height                             |*
        !        |    Tr      | Reference temperature                        |*
        !        |   ustar    | Shear velocity                               |*
        !        |   Tstar    |                                              |*
        !        |  Gammad    |                                              |*
        !        |   theta    | Direction angle ( in radians )               |*
        !        |     qr     | Reference humidity                           |*
        !        |   qstar    |                                              !*
        !        |____________|______________________________________________|*
        !                                                                     *
        ! Output: ___________________________________________________________ *
        !        |__Variable__|_______________Details________________________|*
        !        |     u      | Long. Velocity as a function of altitude     |*
        !        |     v      | Transverse velocity as a function of altitude|*
        !        |     T      | Temeperature as a function of altitude       |*
        !        |     h      | Humidity                                     |*
        !        |____________|______________________________________________|*
        !**********************************************************************
        implicit none
        real(sp), intent(out)   :: u(:), v(:), T(:), h(:)
        real(sp), intent(in)    :: z(:), L, z0, zr, Tr, ustar, Tstar, Gammad
        real(sp), intent(in)    :: theta, qr, q_star
        ! Local Variables *****************************************************
        integer(isp)            :: i 
        real(sp), parameter     :: kappa    = 0.4
        real(sp), parameter     :: Pt       = 0.95 
        real(sp), parameter     :: gamma_m  = 3.6 
        real(sp), parameter     :: gamma_h  = 7.9
        real(sp)                :: psi_u, psi_T, Vel(size(u,dim=1))
        !**********************************************************************

        do i = 1,size(z,dim=1)
            psi_u = ( 1.0 + sqrt(1.0 + gamma_m*(abs(z(i)/L)**(2.0/3.0)) ) & 
                )/( 1.0 + sqrt( 1.0 + gamma_m*(abs(z0/L)**(2.0/3.0)) ) )
            psi_T = ( 1.0 + sqrt(1.0 + gamma_h*(abs(z(i)/L)**(2.0/3.0)) ) & 
                )/( 1.0 + sqrt( 1.0 + gamma_h*(abs(zr/L)**(2.0/3.0)) ) )
            Vel(i) = (ustar/kappa)*( log(z(i)/z0) - 3.0*log(psi_u) )
            T(i) = Tr - (z(i) - zr)*Gammad + (Pt*Tstar/kappa)*(log(z(i)/zr) &
                - 3.0*log( psi_T ) )
            u(i) = Vel(i)*cos(theta) 
            v(i) = Vel(i)*sin(theta)
            h(i) = qr + (Pt*q_star/kappa)*( log(z(i)/zr) - 3.0*log( psi_T ) )
        end do

    end subroutine most 

    function fourier_space(x)
        !**********************************************************************
        ! Purpose:                                                            *
        !       Given the domain of a function g(x), obtain the frequencies   *
        !       and wavenumbers (Fourier domain) corresponding to x.          *
        !                                                                     *
        ! Input: ___________________________________________________________  *
        !       |___Variable____|_______________Description_________________| *
        !       |       x       |   Domain of some function g(x)            | *
        !       |_______________|___________________________________________| *
        !                                                                     *
        ! Output:___________________________________________________________  *
        !       |___Variable____|_______________Description_________________| *
        !       | fourier_space |   Fourier domain (frequency vector)       | *
        !       |_______________|___________________________________________| *
        !**********************************************************************
        implicit none
        real(sp), intent(in)    :: x(:)
        real(sp)                :: fourier_space(size(x,dim=1))
        ! Local variables *****************************************************
        real(sp)                :: f(size(x,dim=1)), fs, df
        integer(isp)            :: i
        !**********************************************************************

        fs = 1.0*size(x)/( 1.0*abs( x(size(x,dim=1)) - x(1) ) ) ! Sampling freq
        df = fs/( 1.0*size(x,dim=1) )                           ! Bin width

        ! Construct frequency array *******************************************
        do i = 1,size(x,dim=1)
            if (i .le. int(size(x,dim=1)/2.0) + 1) then
              f(i) = df*real(i-1)
            else
              f(i) = f( 2*(int(size(x,dim=1)/2.0) + 1) - i )
            end if
        end do
        !**********************************************************************

        fourier_space = f

    end function fourier_space

    subroutine construct_w(w1, w2, w3, wT, kx, ky, kz, sigma, L, sigmaT, LT, &
        comm)
        !**************************************************************************
        !   Purpose:                                                              *
        !           Construct Fourier velocities.                                 *
        !                                                                         *
        !   Method:                                                               *  
        !           The method used to construct the turbulent field is a spectral*
        !           domain algorithm by Frelich et al. 2000.                      *
        !                                                                         *
        ! Link: https://doi.org/10.1175/1520-0450(2001)040<0246:SOTDTV>2.0.CO;2   *
        !                                                                         *
        !                                                                         *
        !   Inputs: ____________________________________________________________  *
        !          |____Variables___|_______________Description_________________| *
        !          |        x       |   Array of x values                       | *
        !          |        y       |   Array of y values                       | *
        !          |        z       |   Array of z values                       | *
        !          |    sigma       |   Variance of the velocity (meters/second)| *
        !          |        L       |   Integral length scale (meters)          | *
        !          |     comm       |   MPI_COMM_WORLD                          | *
        !          |________________|___________________________________________| *
        !                                                                         *
        !                                                                         *
        !   Outputs:____________________________________________________________  *
        !          |____Variables___|_______________Description_________________| *
        !          |       ux       |   Turbulent u velocity                    | *
        !          |       uy       |   Turbulent v velocity                    | *
        !          |       uz       |   Turbulent w velocity                    | *
        !          |________________|___________________________________________| *
        !                                                                         *
        !**************************************************************************
        implicit none
        ! I/O variables ***********************************************************
        complex(sp), intent(out)      :: w1(:,:,:), w2(:,:,:), w3(:,:,:), wT(:,:,:)
        real(sp), intent(in)          :: kx(:), ky(:), kz(:)
        real(sp), intent(in)          :: sigma(:), L(:)
        real(sp), intent(in)          :: sigmaT(:), LT(:)
        integer(isp), intent(in)      :: comm
        ! Local variables *********************************************************
        real(sp)                      :: a1, b1, a2, b2, a, b
        real(sp)                      :: h11, h12, h22, h13, h23, h33, hT
        real(sp)                      :: phi11, phi12, phi22, phi13, phi23, phi33
        real(sp)                      :: phiT
        complex(sp)                   :: N1, N2, N3, NT
        integer(isp)                  :: i, j, k, temp3
        real(sp)                      :: dkx, dky, dkz
        integer(isp)                  :: dim1(3), dim2(3), dim3(3), dimT(3)
        ! MPI variables ***********************************************************
        integer(isp)                  :: proc, nproc, mpi_err
        !**************************************************************************
    
        ! Get processor information ***********************************************
        call MPI_COMM_SIZE(comm, nproc, mpi_err)
        call MPI_COMM_RANK(comm, proc, mpi_err)
        !**************************************************************************
    
        dkx   = kx(2) - kx(1)
        dky   = ky(2) - ky(1)
        dkz   = kz(2) - kz(1)
    
        dim1 = [size(w1,dim=1),size(w1,dim=2),size(w1,dim=3)*nproc]
        dim2 = [size(w2,dim=1),size(w2,dim=2),size(w2,dim=3)*nproc]
        dim3 = [size(w3,dim=1),size(w3,dim=2),size(w3,dim=3)*nproc]
        dimT = [size(wT,dim=1),size(wT,dim=2),size(wT,dim=3)*nproc]
    
        call random_seed()
        do k = 1,size(w1,dim=3)
          temp3 = k + proc*size(w1,dim=3)
          do j = 1,size(w1,dim=2)
            do i = 1,size(w1,dim=1)
              ! Random vectors
              call random_number(a1)
              do while ( a1 .eq. 0 )
                call random_number(a1)
              end do
              call random_number(b1)
              do while ( b1 .eq. 0 )
                call random_number(b1)
              end do
              call gaussian_boxmuller(a, a1, b1, 0.0_sp, 1.0_sp)
              call random_number(a2)
              do while ( a2 .eq. 0 )
                call random_number(a2)
              end do
              call random_number(b2)
              do while ( b1 .eq. 0 )
                call random_number(b1)
              end do
              call gaussian_boxmuller(b, a2, b2, 0.0_sp, 1.0_sp)
              N1 = a + img*b
              call random_number(a1)
              do while ( a1 .eq. 0 )
                call random_number(a1)
              end do
              call random_number(b1)
              do while ( b1 .eq. 0 )
                call random_number(b1)
              end do
              call gaussian_boxmuller(a, a1, b1, 0.0_sp, 1.0_sp)
              call random_number(a2)
              do while ( a2 .eq. 0 )
                call random_number(a2)
              end do
              call random_number(b2)
              do while ( b1 .eq. 0 )
                call random_number(b1)
              end do
              call gaussian_boxmuller(b, a2, b2, 0.0_sp, 1.0_sp)
              N2 = a + img*b
              call random_number(a1)
              do while ( a1 .eq. 0 )
                call random_number(a1)
              end do
              call random_number(b1)
              do while ( b1 .eq. 0 )
                call random_number(b1)
              end do
              call gaussian_boxmuller(a, a1, b1, 0.0_sp, 1.0_sp)
              call random_number(a2)
              do while ( a2 .eq. 0 )
                call random_number(a2)
              end do
              call random_number(b2)
              do while ( b1 .eq. 0 )
                call random_number(b1)
              end do
              call gaussian_boxmuller(b, a2, b2, 0.0_sp, 1.0_sp)
              N3 = a + img*b
              call random_number(a1)
              do while ( a1 .eq. 0 )
                call random_number(a1)
              end do
              call random_number(b1)
              do while ( b1 .eq. 0 )
                call random_number(b1)
              end do
              call gaussian_boxmuller(a, a1, b1, 0.0_sp, 1.0_sp)
              call random_number(a2)
              do while ( a2 .eq. 0 )
                call random_number(a2)
              end do
              call random_number(b2)
              do while ( b1 .eq. 0 )
                call random_number(b1)
              end do
              call gaussian_boxmuller(b, a2, b2, 0.0_sp, 1.0_sp)
              NT = a + img*b
              ! Velocity spectrum
              call velocitySpectrum(phi11, [kx(i), ky(j), kz(temp3)] , &
                sigma(temp3),L(temp3), 1, 1)
              h11 = sqrt( phi11*dkx*dky*dkz )
              w1(i,j,k)  = cmplx(h11*real(N1),h11*imag(N1))*real(product(dim1))
              call velocitySpectrum(phi12, [kx(i), ky(j), kz(temp3)], &
                sigma(temp3), L(temp3), 1, 2)
              h12 = phi12*sqrt( dkx*dky*dkz )/sqrt(phi11)
              call velocitySpectrum(phi22, [kx(i), ky(j), kz(temp3)], &
                sigma(temp3), L(temp3), 2, 2)
              h22 = sqrt( phi22*dkx*dky*dkz  - h12**2 )
              w2(i,j,k)  = (cmplx(h12,0.0)*N1 + &
                cmplx(h22,0.0)*N2)*real(product(dim2))
              call velocitySpectrum(phi13, [kx(i), ky(j), kz(temp3)], &
                sigma(temp3), L(temp3), 1, 3)
              h13 = phi13*sqrt( dkx*dky*dkz )/sqrt(phi11)
              call velocitySpectrum(phi23, [kx(i), ky(j), kz(temp3)], &
                sigma(temp3), L(temp3), 2, 3)
              h23 = ( phi23*dkx*dky*dkz - h12*h13 )/sqrt(phi22)
              call velocitySpectrum(phi33, [kx(i), ky(j), kz(temp3)], &
                sigma(temp3), L(temp3), 3, 3)
              h33 = sqrt( phi33*dkx*dky*dkz - (h13**2) - (h23**2) )
              w3(i,j,k)  = (cmplx(h13,0.0)*N1 + cmplx(h23,0.0)*N2 &
                + cmplx(h33,0.0)*N3)*real(product(dim3))
              call vonKarman_temperature(phiT, &
                sqrt(kx(i)**2 + ky(j)**2 + kz(temp3)**2), LT(temp3), sigmaT(temp3))
              hT = sqrt( phiT*dkx*dky*dkz )
              wT(i,j,k)  = cmplx(hT,0.0)*NT*real(product(dimT))
              if ( isnan(real(w1(i,j,k))) ) then
                write (*,*) real(w1(i,j,k)), kx(i), ky(j), kz(temp3), proc, N1
              end if
            end do
          end do
        end do

    end subroutine construct_w

    subroutine mkl_fft(f, dims, fft_dim, inverse, err_code)
        !**********************************************************************
        ! Purpose: 
        !           Perform Fast Fourier Transform using Intel Math Kernel
        !           library.                                                        
        !                                                                        
        ! Inputs:    _________________________________________________________  
        !           |___Variable___|_______________Purpose____________________|  
        !           |      f       |  Complex potential array                 |
        !           |     dims     |  bounds of processor data on local scale |
        !           |    fft_dim   |  Dimension of fft to take                |
        !           |    inverse   |  Logical arguement (if .true. take       | 
        !           |              |  inverse fft). Default is .false.        |
        !           |     comm     |  MPI_COMM_WORLD                          |
        !           |______________|__________________________________________|  
        !          
        ! Outputs:   _________________________________________________________
        !           |___Variable____|_______________Purpose___________________|  
        !           |      f        |  Complex potential array (in-place)     |
        !           |   err_code    |  Error code                             |
        !           |     time      |  Wall time to step through routine.     |
        !           |_______________|_________________________________________|  
        !                                                                     
        !**********************************************************************
        use mkl_dfti
        implicit none
        complex(sp), intent(inout)          :: f(:)
        integer(sp), intent(in)             :: dims(:)
        integer(sp), intent(in)             :: fft_dim
        logical, intent(in), optional       :: inverse
        integer(isp), intent(out), optional :: err_code
        ! Local variables *****************************************************
        integer(isp)                    :: i
        type(DFTI_DESCRIPTOR), pointer  :: My_Desc_Handle
        integer(isp)                    :: status = 0_isp, ignored_status, num_fft
        logical                         :: backward
        integer(isp), allocatable       :: L(:)

        allocate ( L(fft_dim) )

        if ( present(err_code) ) then
            err_code = 0_isp
        end if

        if ( present(inverse) ) then
            backward = inverse
        else
            backward = .false.
        end if

        do i = 1, fft_dim
            L(i) = dims(i)
        end do

        num_fft = product( dims )/product( L )

        ! Perform a real to complex transform
        status = DftiCreateDescriptor( My_Desc_Handle, DFTI_SINGLE,&
            DFTI_COMPLEX, fft_dim, L )
        if (status .ne. 0) goto 999

        ! Perform grid%ys*im_b%zs_l one dimensional transforms in t dimension
        status = DftiSetValue(My_Desc_Handle, DFTI_NUMBER_OF_TRANSFORMS, num_fft)
        if (status .ne. 0) goto 999
        status = DftiSetValue(My_Desc_Handle, DFTI_BACKWARD_SCALE, 1.0/( 1.0*product( L ) ) )
        if (status .ne. 0) goto 999
        status = DftiSetValue(My_Desc_Handle, DFTI_INPUT_DISTANCE, product(L) )
        if (status .ne. 0) goto 999
        status = DftiCommitDescriptor( My_Desc_Handle )
        if (status .ne. 0) goto 999
        if ( .not. backward ) then
            status = DftiComputeForward( My_Desc_Handle, f )
            if (status .ne. 0) goto 999
        else
            status = DftiComputeBackward( My_Desc_Handle, f )
            if (status .ne. 0) goto 999
        end if
        100 continue
        ignored_status = DftiFreeDescriptor(My_Desc_Handle)
        goto 200

        999 continue
        print '(" Error, status = ",I0)', status
        if ( present(err_code) ) then
            err_code = 1_isp
        end if
        goto 100

        200 continue
      
    end subroutine mkl_fft

    subroutine transpose(out, temp, in, dout, dtemp, din, &
        index_switch, comm, err, time)
        !**********************************************************************
        ! Purpose: 
        !           Perform forward and backwards transposition of a mutli-
        !           dimensional array that is distributed across nproc
        !           processors in a slab decomposition. That is, the last
        !           dimension of the array has a local length on each
        !           processor l = L/nproc. This subroutine can transpose arrays
        !           up to 3 dimensions. For more dimensions, an in-place 
        !           transposition routine would have to be written, this would
        !           inevitably take longer as there would be two write calls.                                                          
        !                                                                        
        ! Inputs:    _________________________________________________________  
        !           |___Variable___|_______________Purpose____________________|  
        !           |      out     |   The transpose array.                   |
        !           |     temp     |   Temporary array for storage.           |  
        !           |       in     |   Input array                            | 
        !           |     dout     |   Dimensions of the transpose array.     |
        !           |    dtemp     |   Dims of temporary array for storage.   |  
        !           |      din     |   Dimensions of the input array          |
        !           | index_switch |   Order in which to re-arrange indices   |
        !           |              |   (i.e. 123 -> 231, index_switch=231)    |
        !           |     comm     |   MPI communicator                       | 
        !           |______________|__________________________________________|  
        !          
        ! Outputs:   _________________________________________________________  
        !           |___Variable___|_______________Purpose____________________|  
        !           |      out     |   The transpose array.                   |
        !           |     temp     |   Temporary array for storage.           |  
        !           |       in     |   Input array                            |  
        !           |     time     |   Time of transposition                  |
        !           |______________|__________________________________________|  
        !                                                                     
        !**********************************************************************
        use mpi
        implicit none
        complex(sp), intent(inout)  :: out(:)
        complex(sp), intent(inout)  :: temp(:)
        complex(sp), intent(inout)  :: in(:)
        integer(isp), intent(in)    :: dout(:)
        integer(isp), intent(in)    :: dtemp(:)
        integer(isp), intent(in)    :: din(:)
        integer(isp), intent(in)    :: index_switch
        integer(isp), intent(in)    :: comm
        integer(isp), intent(out), optional :: err
        real(dp), intent(out), optional     :: time
        ! Local variables *****************************************************
        integer(isp)    :: i, j                                ! Loop variables
        integer(isp)    :: nproc, proc, mpi_err                ! MPI variables
        real(dp)        :: t1, t2                              ! Time vars
        integer(isp), allocatable :: sendcounts(:,:), sdisps(:,:), &
            recvcounts(:,:), rdisps(:,:), sendtype(:), recvtype(:)
        integer(isp)    :: type_send
        !**********************************************************************
        
        ! Obtain processor info ***********************************************
        call MPI_COMM_SIZE(comm, nproc, mpi_err)
        call MPI_COMM_RANK(comm, proc, mpi_err)
        !**********************************************************************

        ! Wait for all processes to reach this point, since this is a global
        ! operation we must ensure that all processors are "up to date."
        call MPI_BARRIER(comm, mpi_err)

        ! Options *************************************************************
        if ( present(time) ) then
            t1 = MPI_WTIME()
        end if

        if ( present(err) ) then
            err = 0_isp
        end if
        !**********************************************************************

        ! Allocate the MPI send and receive info ******************************
        allocate ( sendcounts(nproc, dout(3)), sdisps(nproc,dout(3)) )
        allocate ( recvcounts(nproc,dout(3)), rdisps(nproc,dout(3)) )
        allocate ( sendtype(nproc), recvtype(nproc) )
        !**********************************************************************

        !**********************************************************************
        if ( index_switch .eq. 231 ) then
            ! Forward transpose ***********************************************

            ! Create datatype *************************************************
            call MPI_Type_vector(din(2)*din(3), 1, din(1), &
                MPI_COMPLEX, type_send, mpi_err)
            call MPI_Type_commit(type_send, mpi_err)
            !******************************************************************

            ! Counts and displacements ****************************************
            do j = 1,dout(3)
                do i = 1,nproc
                    sendcounts(i,j) = 1
                    sdisps(i,j)     = ((i-1)*dout(3) + &
                        (j-1))*int(sizeof(in(1)),4)
                    sendtype(i)     = type_send
                    recvcounts(i,j) = din(2)*din(3)
                    rdisps(i,j)     = ((i-1)*din(2)*din(3) + &
                        (j-1)*dout(1)*dout(2))*int(sizeof(out(1)),4)
                    recvtype(i)     = MPI_COMPLEX
                end do
            end do
            !******************************************************************

            ! All to All communication ****************************************
            do j = 1,dout(3)
                call MPI_ALLTOALLW(in, sendcounts(:,j), sdisps(:,j), sendtype,&
                    out, recvcounts(:,j), rdisps(:,j), recvtype, comm, mpi_err)
            end do
            !******************************************************************

            ! Free the type ***************************************************
            call MPI_Type_free(type_send, mpi_err)
            !******************************************************************

            !******************************************************************
        else if ( index_switch .eq. 312) then
            ! Backward transpose **********************************************

            ! Note: It is actually quicker to compute two forward transposes.
            !       It also scales much better with processes.

            ! Forward transpose to temp ***************************************

            ! Create datatype *************************************************
            call MPI_Type_vector(din(2)*din(3), 1, din(1), &
                MPI_COMPLEX, type_send, mpi_err)
            call MPI_Type_commit(type_send, mpi_err)
            !******************************************************************

            ! Counts and displacements ****************************************
            do j = 1,dtemp(3)
                do i = 1,nproc
                    sendcounts(i,j) = 1
                    sdisps(i,j)     = ((i-1)*dtemp(3) + &
                        (j-1))*int(sizeof(in(1)),4)
                    sendtype(i)     = type_send
                    recvcounts(i,j) = din(2)*din(3)
                    rdisps(i,j)     = ((i-1)*din(2)*din(3) + &
                        (j-1)*dtemp(1)*dtemp(2))*int(sizeof(temp(1)),4)
                    recvtype(i)     = MPI_COMPLEX
                end do
            end do
            !******************************************************************

            ! All to All communication ****************************************
            do j = 1,dout(3)
                call MPI_ALLTOALLW(in, sendcounts(:,j), sdisps(:,j), sendtype,&
                    temp, recvcounts(:,j), rdisps(:,j), recvtype, comm, mpi_err)
            end do
            !******************************************************************

            ! Free the type ***************************************************
            call MPI_Type_free(type_send, mpi_err)
            !******************************************************************

            !******************************************************************

            ! Forward transpose ***********************************************

            ! Create datatype *************************************************
            call MPI_Type_vector(dtemp(2)*dtemp(3), 1, dtemp(1), &
                MPI_COMPLEX, type_send, mpi_err)
            call MPI_Type_commit(type_send, mpi_err)
            !******************************************************************

            ! Counts and displacements ****************************************
            do j = 1,dout(3)
                do i = 1,nproc
                    sendcounts(i,j) = 1
                    sdisps(i,j)     = ((i-1)*dout(3) + &
                        (j-1))*int(sizeof(temp(1)),4)
                    sendtype(i)     = type_send
                    recvcounts(i,j) = dtemp(2)*dtemp(3)
                    rdisps(i,j)     = ((i-1)*dtemp(2)*dtemp(3) + &
                        (j-1)*dout(1)*dout(2))*int(sizeof(out(1)),4)
                    recvtype(i)     = MPI_COMPLEX
                end do
            end do
            !******************************************************************

            ! All to All communication ****************************************
            do j = 1,dout(3)
                call MPI_ALLTOALLW(temp, sendcounts(:,j),sdisps(:,j),sendtype,&
                    out, recvcounts(:,j), rdisps(:,j), recvtype, comm, mpi_err)
            end do
            !******************************************************************

            ! Free the type ***************************************************
            call MPI_Type_free(type_send, mpi_err)
            !******************************************************************

            !******************************************************************

            !******************************************************************
        else if ( index_switch .eq. 132 ) then
            ! Not offered yet
            continue
        else if ( index_switch .eq. 213 ) then
            ! Not offered yet
            continue
        else if ( index_switch .eq. 321 ) then
            ! Not offered yet
            continue
        else if ( index_switch .eq. 21 ) then
            ! Not offered yet
            continue
        else
            ! Error
            continue
        end if
        !**********************************************************************

        ! Deallocate everything ***********************************************
        deallocate ( sendcounts, sdisps, recvcounts )
        deallocate ( rdisps, sendtype, recvtype )
        !**********************************************************************

        call MPI_BARRIER(comm, mpi_err)

        if ( present(time) ) then
            t2 = MPI_WTIME()
            time = t2 - t1
        end if

    end subroutine transpose

    subroutine create_h5_f(filename, comm, groupnames)
      character(len=*), intent(in)         :: filename  ! file name
      character(len=*), dimension(:), optional, intent(in) :: groupnames ! dataset name
      integer(isp), intent(in)                :: comm
      integer                               :: mpi_err, error
      integer                               :: i
      integer(HID_T)                        :: file_id   ! File Identifier
      integer(HID_T)                        :: group_id   ! Dataset Identifier
      integer(HID_T)                        :: plist_id

      call h5open_f(error)

      ! Set up parallel i/o access
      call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
      call h5pset_fapl_mpio_f(plist_id, comm, MPI_INFO_NULL, error)

      ! Collectively create file
      call h5fcreate_f(trim(filename), H5F_ACC_TRUNC_F, file_id, error, &
        access_prp = plist_id)
      call h5pclose_f(plist_id, error)

      if ( present(groupnames) ) then
          do i = 1,size(groupnames)
              call h5gcreate_f(file_id, trim(groupnames(i)), group_id, error)
              call h5gclose_f(group_id, error)
          end do
      end if

      call h5fclose_f(file_id, error)

      call h5close_f(error)

      call MPI_BARRIER(comm, mpi_err)

    end subroutine create_h5_f

    subroutine create_h5_d(filename, dsetname, dim, comm, f, chunked, &
        dim_chunk, offset, extendable, loff)
        character(len=*), intent(in)           :: filename  ! file name
        character(len=*), intent(in)           :: dsetname    ! dataset name
        integer(isp), intent(in)             :: dim(:)      ! total data dimensions
        integer(isp), intent(in)             :: comm        ! MPI communicator
        real(sp),  optional, intent(in)         :: f(:)        ! data buffer
        logical, optional, intent(in)          :: chunked     ! .true. = dataset is chunked
        integer(isp), optional, intent(in)   :: dim_chunk(:)   ! Local chunk size
        integer(isp), optional, intent(in)   :: offset(:)      ! Global position of chunk
        logical, optional, intent(in)          :: extendable  ! .true. = dataset is chunked
        integer(isp), intent(in), optional   :: loff(:)
        integer                            :: mpi_err, error
        integer                            :: errcode, i
        integer(HID_T)                     :: file_id     ! File Identifier
        integer(HID_T)                     :: dset_id     ! Dataset Identifier
        integer(HID_T)                     :: dataspace   ! Dataspace ID
        integer(HID_T)                     :: memspace
        integer(HID_T)                     :: plist_id
        integer(HSIZE_T), allocatable      :: dims(:), stride(:)
        integer(HSIZE_T), allocatable      :: maxdims(:), dimt(:), off_set(:), offset_out(:), loc_off(:)
        logical, dimension(2)              :: cases

        allocate ( dims(size(dim)), stride(size(dim)), off_set(size(dim)), &
            maxdims(size(dim)), dimt(size(dim)) , offset_out(size(dim)), loc_off(size(dim)))

        if ( present(loff) ) then
            loc_off = loff
        else
            loc_off = 0
        end if

        if ( present(chunked) ) then
            cases(1) = chunked
        else
            cases(1) = .false.
        end if

        if ( present(extendable) ) then
            cases(2) = extendable
        else
            cases(2) = .false.
        end if

        call h5open_f(error)

        ! Set up parallel i/o access
        call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
        call h5pset_fapl_mpio_f(plist_id, comm, MPI_INFO_NULL, error)

        ! Collectively open file
        call h5fopen_f(trim(filename), H5F_ACC_RDWR_F, file_id, error, access_prp = plist_id)
        call h5pclose_f(plist_id, error)

        if ( cases(1) .and. cases(2) ) then
            ! Chunked and Extendable

            ! Dimensions
            do i = 1, size(dim)
                dimt(i)     = dim(i)
                dims(i)     = dim_chunk(i)    ! Dimensions of chunks
                stride(i)   = 1               ! Continuous data
                off_set(i)  = offset(i)       ! Zero indexed
                offset_out(i) = loc_off(i)
                maxdims(i)  = H5S_UNLIMITED_F
            end do

            call h5screate_simple_f(size(dimt), dimt, dataspace, error, maxdims)
            call h5screate_simple_f(size(dims), dimt, memspace, error)

            ! Create chunked dataset
            call h5pcreate_f(H5P_DATASET_CREATE_F, plist_id, error)
            call h5pset_chunk_f(plist_id, size(dims), dims, error)   ! Only works if blmax = lmax
            call h5pset_deflate_f(plist_id, 6, error) ! Compression level
            call h5dcreate_f(file_id, dsetname, H5T_NATIVE_REAL, dataspace, &
                dset_id, error, plist_id)

            if ( present(f) ) then
                ! Select hyperslab in file
                call h5dget_space_f(dset_id, dataspace, error)
                call h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, &
                    off_set, dims, error)
                call h5dget_space_f(dset_id, memspace, error)
                call h5sselect_hyperslab_f(memspace, H5S_SELECT_SET_F, &
                    offset_out, dims, error)

                ! Create property list for collective dataset write
                call h5pclose_f(plist_id, error)
                call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
                call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)

                ! Write the dataset collectively
                call h5dwrite_f(dset_id, H5T_NATIVE_REAL, f, dimt, error, &
                    file_space_id = dataspace, mem_space_id = memspace, &
                    xfer_prp = plist_id)
            end if
            call h5sclose_f(memspace, error)
        else if ( cases(1) .and. (.not. cases(2) ) ) then
            ! Chunked

            ! Dimensions
            do i = 1, size(dim)
                dimt(i)     = dim(i)
                dims(i)     = dim_chunk(i)    ! Dimensions of chunks
                stride(i)   = 1               ! Continuous data
                off_set(i)  = offset(i)       ! Zero indexed
                maxdims(i)  = H5S_UNLIMITED_F
            end do

            call h5screate_simple_f(size(dimt), dimt, dataspace, error)
            call h5screate_simple_f(size(dims), dims, memspace, error)

            ! Create chunked dataset
            call h5pcreate_f(H5P_DATASET_CREATE_F, plist_id, error)
            call h5pset_chunk_f(plist_id, size(dims), dims, error)   ! Only works if blmax = lmax
            call h5pset_deflate_f(plist_id, 6, error) ! Compression level
            call h5dcreate_f(file_id, dsetname, H5T_NATIVE_REAL, dataspace, &
                dset_id, error, plist_id)

            if ( present(f) ) then
                ! Select hyperslab in file
                call h5dget_space_f(dset_id, dataspace, error)
                call h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, &
                    off_set, dims, error)

                ! Create property list for collective dataset write
                call h5pclose_f(plist_id, error)
                call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
                call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)

                ! Write the dataset collectively
                call h5dwrite_f(dset_id, H5T_NATIVE_REAL, f, dimt, error, &
                    file_space_id = dataspace, mem_space_id = memspace, &
                    xfer_prp = plist_id)

            end if
            call h5sclose_f(memspace, error)

        else if ( (.not. cases(1) ) .and. cases(2) ) then
            print *, "Error in subroutine create_h5_d(): Cannot extend dataset without chunking"
            call h5fclose_f(file_id, error)

            call h5close_f(error)
            deallocate ( dims, stride, off_set, maxdims, dimt )
            call MPI_Abort(comm, errcode, mpi_err)
        else
            ! Normal

            ! Dimensions
            do i = 1, size(dim)
                dimt(i)      = dim(i)
                dims(i)      = dim(i)    ! Dimensions
                stride(i)    = 1               ! Continuous data
                off_set(i)   = 0       ! Zero indexed
                maxdims(i)   = H5S_UNLIMITED_F
            end do

            ! Create dataset
            call h5screate_simple_f(size(dims), dims, dataspace, error)
            call h5pcreate_f(H5P_DATASET_CREATE_F, plist_id, error)
            call h5dcreate_f(file_id, dsetname, H5T_NATIVE_REAL, dataspace, &
                dset_id, error, plist_id)

            if ( present(f) ) then
                call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
                call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)

                ! Write the datasets collectively
                call h5dwrite_f(dset_id, H5T_NATIVE_REAL, f, dims, error, &
                  xfer_prp = plist_id)
            end if
        end if

        call h5sclose_f(dataspace, error)
        call h5pclose_f(plist_id, error)
        call h5dclose_f(dset_id, error)
        call h5fclose_f(file_id, error)

        call h5close_f(error)

        deallocate ( dims, stride, off_set, maxdims, dimt )

        call MPI_BARRIER(comm, mpi_err)

    end subroutine create_h5_d

    subroutine write_h5_d(filename, dsetname, dim, comm, f, chunked, &
        dim_chunk, offset, extendable, extension)
        character(len=*), intent(in)           :: filename  ! file name
        character(len=*), intent(in)           :: dsetname    ! dataset name
        integer(isp), intent(in)             :: dim(:)      ! total data dimensions
        integer(isp), intent(in)             :: comm        ! MPI communicator
        real(sp), intent(in)               :: f(:)        ! data buffer
        logical, optional, intent(in)          :: chunked     ! .true. = dataset is chunked
        integer(isp), optional, intent(in)   :: dim_chunk(:)   ! Local chunk size
        integer(isp), optional, intent(in)   :: offset(:)      ! Global position of chunk
        logical, optional, intent(in)          :: extendable  ! .true. = dataset is chunked
        integer(isp), optional, intent(in)   :: extension(:)   ! dimension to extend the dataset
        integer                            :: mpi_err, error
        integer                            :: errcode, i
        integer(HID_T)                     :: file_id     ! File Identifier
        integer(HID_T)                     :: dset_id     ! Dataset Identifier
        integer(HID_T)                     :: dataspace   ! Dataspace ID
        integer(HID_T)                     :: memspace
        integer(HID_T)                     :: plist_id
        integer(HSIZE_T), allocatable      :: dims(:), stride(:)
        integer(HSIZE_T), allocatable      :: maxdims(:), dimt(:), off_set(:)
        logical, dimension(2)              :: cases

        allocate ( dims(size(dim)), stride(size(dim)), off_set(size(dim)), &
            maxdims(size(dim)), dimt(size(dim)) )



        if ( present(chunked) ) then
            cases(1) = chunked
        else
            cases(1) = .false.
        end if

        if ( present(extendable) ) then
            cases(2) = extendable
        else
            cases(2) = .false.
        end if

        call h5open_f(error)

        ! Set up parallel i/o access
        call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
        call h5pset_fapl_mpio_f(plist_id, comm, MPI_INFO_NULL, error)

        ! Collectively open file
        call h5fopen_f(trim(filename), H5F_ACC_RDWR_F, file_id, error, access_prp = plist_id)
        call h5pclose_f(plist_id, error)

        if ( cases(1) .and. cases(2) ) then
            ! Chunked and Extendable

            ! Dimensions
            do i = 1, size(dim)
                dimt(i)     = dim(i)
                dims(i)     = dim_chunk(i)    ! Dimensions of chunks
                stride(i)   = 1               ! Continuous data
                off_set(i)  = offset(i)       ! Zero indexed
                maxdims(i)  = extension(i)    !
            end do

            ! Open dataset and get dataset properties
            call h5dopen_f(file_id, dsetname, dset_id, error)
            call h5screate_simple_f(size(dims), dims, memspace, error)

            ! Extend the dataset
            call h5dset_extent_f(dset_id, maxdims, error)

            ! Select hyperslab in file
            call h5dget_space_f(dset_id, dataspace, error)
            call h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, &
                off_set, dims, error)

            ! Create property list for collective dataset write
            call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
            call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)

            ! Write the dataset collectively
            call h5dwrite_f(dset_id, H5T_NATIVE_REAL, f, dimt, error, &
                file_space_id = dataspace, mem_space_id = memspace, &
                xfer_prp = plist_id)
            call h5sclose_f(memspace, error)

        else if ( cases(1) .and. (.not. cases(2) ) ) then
            ! Chunked

            ! Dimensions
            do i = 1, size(dim)
                dimt(i)     = dim(i)
                dims(i)     = dim_chunk(i)    ! Dimensions of chunks
                stride(i)   = 1               ! Continuous data
                off_set(i)  = offset(i)       ! Zero indexed
                maxdims(i)  = 0               !
            end do

            ! Open dataset and get dataset properties
            call h5dopen_f(file_id, dsetname, dset_id, error)
            call h5screate_simple_f(size(dims), dims, memspace, error)

            ! Select hyperslab in file
            call h5dget_space_f(dset_id, dataspace, error)
            call h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, &
                off_set, dims, error)

            ! Create property list for collective dataset write
            call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
            call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)

            ! Write the dataset collectively
            call h5dwrite_f(dset_id, H5T_NATIVE_REAL, f, dimt, error, &
                file_space_id = dataspace, mem_space_id = memspace, &
                xfer_prp = plist_id)
            call h5sclose_f(memspace, error)

        else if ( (.not. cases(1) ) .and. cases(2) ) then
            print *, "Error in subroutine create_h5_d(): Cannot extend dataset without chunking"
            call h5fclose_f(file_id, error)

            call h5close_f(error)
            deallocate ( dims, stride, off_set, maxdims, dimt )
            call MPI_Abort(comm, errcode, mpi_err)
        else
            ! Normal

            ! Dimensions
            do i = 1, size(dim)
                dimt(i)      = dim(i)
                dims(i)      = dim(i)    ! Dimensions
                stride(i)    = 1               ! Continuous data
                off_set(i)   = 0       ! Zero indexed
                maxdims(i)   = 0
            end do

            ! Open dataset
            call h5dopen_f(file_id, dsetname, dset_id, error)

            call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
            call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)

            ! Write the datasets collectively
            call h5dwrite_f(dset_id, H5T_NATIVE_REAL, f, dims, error, xfer_prp = plist_id)

        end if

        call h5sclose_f(dataspace, error)
        call h5pclose_f(plist_id, error)
        call h5dclose_f(dset_id, error)
        call h5fclose_f(file_id, error)

        call h5close_f(error)

        deallocate ( dims, stride, off_set, maxdims, dimt )

        call MPI_BARRIER(comm, mpi_err)

    end subroutine write_h5_d

end module functions