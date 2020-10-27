module functions
    use mkl_dfti
    use special_functions
    implicit none
    integer, parameter       :: sp = selected_real_kind(6, 37)
    integer, parameter       :: dp = selected_real_kind(15, 307)
    integer, parameter       :: isp = selected_int_kind(6)
    integer, parameter       :: idp = selected_int_kind(10)
    real(kind=dp), parameter :: pi = 3.141592653589793_dp
    complex(kind=sp), parameter:: img = (0.0_sp,1.0_sp)


  contains

    subroutine gaussian_boxmuller(y, u1, u2, a, s)
      real, intent(out) :: y(:)
      real, intent(in)  :: u1(:), u2(:), a, s
      real :: x(size(u1,dim=1))
      integer :: i

      do i = 1,size(u1,dim=1)
        x(i) = sqrt( -2.0*log(u1(i)))*sin(2.0*pi*u2(i))
        y(i)  = s*x(i) + a
      end do

    end subroutine gaussian_boxmuller

    subroutine vonKarman_velocity(E, k, L, sigma)
      real, intent(out) :: E
      real, intent(in) :: k
      real, intent(in) :: L, sigma
      real :: a, b

      a = 1.12878702991
      b = 2.67893853470

      E = ( 55.0*a/(9.0*(sqrt(pi))*b ) )*( &
        (sigma**2 * k**4 * L**5)/( ( 1.0 + (k**2 * L**2) )**(17.0/6.0)))

    end subroutine vonKarman_velocity

    subroutine velocitySpectrum(phi, kx, ky, kz, sigma, L, n, m)
      real, intent(out) :: phi
      real, intent(in)  :: kx, ky, kz
      real, intent(in)  :: sigma, L
      integer, intent(in) :: n, m
      real                :: E, k


      k = sqrt( kx**2 + ky**2 + kz**2 )
      call vonKarman_velocity(E, k, L, sigma)
      if ( k .eq. 0.0 ) then
        phi = 0.0
      else if ( (n .eq. 1) .and. (m .eq. 1) ) then
        phi = (3.0/(4.0*pi*(k**4)) )*E*(kronecker(n,m)*(k**2) - &
          kx*kx)
      else if ( (n .eq. 1) .and. (m .eq. 2) ) then
        phi = (3.0/(4.0*pi*(k**4)) )*E*(kronecker(n,m)*(k**2) - &
          kx*ky)
      else if ( (n .eq. 1) .and. (m .eq. 3) ) then
        phi = (3.0/(4.0*pi*(k**4)) )*E*(kronecker(n,m)*(k**2) - &
          kx*kz)
      else if ( (n .eq. 2) .and. (m .eq. 1) ) then
        phi = (3.0/(4.0*pi*(k**4)) )*E*(kronecker(n,m)*(k**2) - &
          ky*kx)
      else if ( (n .eq. 2) .and. (m .eq. 2) ) then
        phi = (3.0/(4.0*pi*(k**4)) )*E*(kronecker(n,m)*(k**2) - &
          ky*ky)
      else if ( (n .eq. 2) .and. (m .eq. 3) ) then
        phi = (3.0/(4.0*pi*(k**4)) )*E*(kronecker(n,m)*(k**2) - &
          ky*kz)
      else if ( (n .eq. 3) .and. (m .eq. 1) ) then
        phi = (3.0/(4.0*pi*(k**4)) )*E*(kronecker(n,m)*(k**2) - &
          kz*kx)
      else if ( (n .eq. 3) .and. (m .eq. 2) ) then
        phi = (3.0/(4.0*pi*(k**4)) )*E*(kronecker(n,m)*(k**2) - &
          kz*ky)
      else if ( (n .eq. 3) .and. (m .eq. 3) ) then
        phi = (3.0/(4.0*pi*(k**4)) )*E*(kronecker(n,m)*(k**2) - &
          kz*kz )
      end if

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

    subroutine fft(w, array_size, dim, inverse)
      ! ========================================================================
      ! Note: The forward Fourier transform in Frelich has the scaling factor
      !
      ! Note: This subroutine follows the same convention
      ! ========================================================================

      ! Complex to complex fft
      complex, intent(inout)          :: w(:)
      integer, intent(in)             :: array_size(:)
      integer, intent(in)             :: dim
      logical, intent(in)             :: inverse
      type(DFTI_DESCRIPTOR), pointer  :: My_Desc_Handle
      integer(isp)                    :: status = 0_isp, ignored_status

      ! Perform a complex to complex transform
      status = DftiCreateDescriptor( My_Desc_Handle, DFTI_SINGLE,&
        DFTI_COMPLEX, dim, array_size )
      if (status .ne. 0) goto 999
      status = DftiSetValue(My_Desc_Handle, DFTI_FORWARD_SCALE, &
        1.0/( 1.0*product( array_size ) ) )
      if (status .ne. 0) goto 999
      status = DftiCommitDescriptor( My_Desc_Handle )
      if (status .ne. 0) goto 999
      if ( .not. inverse ) then
        status = DftiComputeForward( My_Desc_Handle, w )
        if (status .ne. 0) goto 999
      else
        status = DftiComputeBackward( My_Desc_Handle, w )
        if (status .ne. 0) goto 999
      end if
      100 continue
      ignored_status = DftiFreeDescriptor(My_Desc_Handle)
      goto 200

      999 continue
      print '(" Error, status = ",I0)', status

      goto 100

      200 continue

    end subroutine fft

    subroutine construct_wavenumbers(kx, ky, kz, x, y, z)
      real(sp), intent(out)     :: kx(:), ky(:), kz(:)
      real(sp), intent(in)      :: x(:), y(:), z(:)
      real(sp)                  :: fxs, fys, fzs
      real(sp)                  :: dfx, dfy, dfz
      real(sp), allocatable     :: fx(:), fy(:), fz(:)
      integer :: i, j, k

      allocate ( fx(size(x,dim=1)), fy(size(y,dim=1)), fz(size(z,dim=1)) )

      ! Sampling frequencies ===================================================
      fxs = 1.0*size(x)/( 1.0_sp*abs( x(size(x)) - x(1) ) )
      fys = 1.0*size(y)/( 1.0_sp*abs( y(size(y)) - y(1) ) )
      fzs = 1.0*size(z)/( 1.0_sp*abs( z(size(z)) - z(1) ) )
      ! ========================================================================
      ! Bin Widths =============================================================
      dfx = fxs/(size(x)*1.0_sp)
      dfy = fys/(size(y)*1.0_sp)
      dfz = fzs/(size(z)*1.0_sp)
      ! ========================================================================
      ! Construct frequencies ==================================================
      do i = 1,size(x)
        if (i .le. int(size(x)/2.0) + 1) then
          fx(i) = dfx*real(i-1) ! dft*real(i-1)
        else
          fx(i) = fx(2*(int(size(x)/2.0) + 1) - i)
        end if
      end do
      do i = 1,size(y)
        if (i .le. int(size(y)/2.0) + 1) then
          fy(i) = dfy*real(i-1) ! dfy*real(i-1)
        else
          fy(i) = fy(2*(int(size(y)/2.0) + 1) - i)
        end if
      end do
      do i = 1,size(z)
        if (i .le. int(size(z)/2.0) + 1) then
          fz(i) = dfz*real(i-1) ! dfz*real(i-1)
        else
          fz(i) = fz(2*(int(size(z)/2.0) + 1) - i)
        end if
      end do
      ! ========================================================================

      ! Wavenumbers ============================================================
      kx = 2.0*pi*fx
      ky = 2.0*pi*fy
      kz = 2.0*pi*fz
      ! ========================================================================

      deallocate ( fx, fy, fz )

    end subroutine construct_wavenumbers

    subroutine construct_w(w1, w2, w3, x, y, z, kx, ky, kz, mean, sigma, L,type)
      complex(sp), intent(out)      :: w1(:,:,:), w2(:,:,:), w3(:,:,:)
      real(sp), intent(in)          :: kx(:), ky(:), kz(:)
      real(sp), intent(in)          :: x(:), y(:), z(:)
      real(sp), intent(in)          :: mean, sigma, L
      character(len=*), intent(in)  :: type
      real(sp)                      :: a1(1), b1(1), a2(1), b2(1), a(1), b(1)
      real(sp)                      :: h11, h12, h22, h13, h23, h33
      real(sp)                      :: phi11, phi12, phi22, phi13, phi23, phi33
      complex(sp)                   :: N1, N2, N3
      real(sp)                      :: dkx, dky, dkz, G, Gd, r
      real(dp)                      :: vj1, vj2, vy1, vy2, vi1, vi2, vk1, vk2
      integer(isp)                  :: i, j, k
      complex(sp), allocatable      :: B11(:,:,:), B22(:,:,:), B33(:,:,:)

      dkx   = kx(2) - kx(1)
      dky   = ky(2) - ky(1)
      dkz   = kz(2) - kz(1)

      allocate ( B11(size(w1,dim=1),size(w1,dim=2),size(w1,dim=3)) )
      allocate ( B22(size(w2,dim=1),size(w2,dim=2),size(w2,dim=3)) )
      allocate ( B33(size(w3,dim=1),size(w3,dim=2),size(w3,dim=3)) )

      if ( type .eq. "gaussian" ) then
        ! Loop through wavenumbers with gaussian phase correlation =============
        call random_seed()
        do k = 1,size(w1,dim=3)
          do j = 1,size(w1,dim=2)
            do i = 1,size(w1,dim=1)
              ! Random vectors
              call random_number(a1)
              call random_number(b1)
              call gaussian_boxmuller(a, a1, b1, 0.0, 1.0)
              call random_number(a2)
              call random_number(b2)
              call gaussian_boxmuller(b, a2, b2, 0.0, 1.0)
              N1 = a(1) + img*b(1)
              call random_number(a1)
              call random_number(b1)
              call gaussian_boxmuller(a, a1, b1, 0.0, 1.0)
              call random_number(a2)
              call random_number(b2)
              call gaussian_boxmuller(b, a2, b2, 0.0, 1.0)
              N2 = a(1) + img*b(1)
              call random_number(a1)
              call random_number(b1)
              call gaussian_boxmuller(a, a1, b1, 0.0, 1.0)
              call random_number(a2)
              call random_number(b2)
              call gaussian_boxmuller(b, a2, b2, 0.0, 1.0)
              N3 = a(1) + img*b(1)
              ! Velocity spectrum
              call velocitySpectrum(phi11, kx(i), ky(j), kz(k), sigma, L, 1, 1)
              h11 = sqrt( phi11*dkx*dky*dkz )
              w1(i,j,k)  = cmplx(h11,0.0)*N1
              call velocitySpectrum(phi12, kx(i), ky(j), kz(k), sigma, L, 1, 2)
              h12 = phi12*sqrt( dkx*dky*dkz )/sqrt(phi11)
              call velocitySpectrum(phi22, kx(i), ky(j), kz(k), sigma, L, 2, 2)
              h22 = sqrt( phi22*dkx*dky*dkz  - h12**2 )
              w2(i,j,k)  = (cmplx(h12,0.0)*N1 + cmplx(h22,0.0)*N2)
              call velocitySpectrum(phi13, kx(i), ky(j), kz(k), sigma, L, 1, 3)
              h13 = phi13*sqrt( dkx*dky*dkz )/sqrt(phi11)
              call velocitySpectrum(phi23, kx(i), ky(j), kz(k), sigma, L, 2, 3)
              h23 = ( phi23*dkx*dky*dkz - h12*h13 )/sqrt(phi22)
              call velocitySpectrum(phi33, kx(i), ky(j), kz(k), sigma, L, 3, 3)
              h33 = sqrt( phi33*dkx*dky*dkz - (h13**2) - (h23**2) )
              w3(i,j,k)  = (cmplx(h13,0.0)*N1 + cmplx(h23,0.0)*N2 &
                + cmplx(h33,0.0)*N3)
            end do
          end do
        end do
        ! ======================================================================
      else if ( type .eq. "vonkarman" ) then
        ! Model Correlation function
        do k = 1,size(w1,dim=3)
          do j = 1,size(w1,dim=2)
            do i = 1,size(w1,dim=1)
              r = sqrt( x(i)**2 + y(j)**2 + z(k)**2 )
              call ajyik ( dble(r/L), vj1, vj2, vy1, vy2, vi1, vi2, vk1, vk2 )
              G = 0.592548*(r**(1.0/3.0))*vk1
              Gd = 0.29627426*(r**(4.0/3.0))*vk2
              B11(i,j,k) = sigma**2 * G
              B22(i,j,k) = sigma**2 * ( G - Gd)
              B33(i,j,k) = sigma**2 * ( G - Gd)
            end do
          end do
        end do

        call fft( B11(:,1,1), [size(B11,dim=1), size(B11,dim=2), &
          size(B11,dim=3)], 3, .false. )
        call fft( B22(:,1,1), [size(B22,dim=1), size(B22,dim=2), &
          size(B22,dim=3)], 3, .false. )
        call fft( B33(:,1,1), [size(B33,dim=1), size(B33,dim=2), &
          size(B33,dim=3)], 3, .false. )

        w1 = B11
        w2 = B22
        w3 = B33

      else if ( type .eq. "uncorrelated" ) then
        ! Loop through wavenumbers with gaussian phase correlation =============
        call random_seed()
        do k = 1,size(w1,dim=3)
          do j = 1,size(w1,dim=2)
            do i = 1,size(w1,dim=1)
              ! Random vectors
              call random_number(a)
              call random_number(b)
              N1 = a(1) + img*b(1)
              call random_number(a)
              call random_number(b)
              N2 = a(1) + img*b(1)
              call random_number(a)
              call random_number(b)
              N3 = a(1) + img*b(1)
              ! Velocity spectrum
              call velocitySpectrum(phi11, kx(i), ky(j), kz(k), sigma, L, 1, 1)
              h11 = sqrt( phi11*dkx*dky*dkz )
              w1(i,j,k)  = cmplx(h11,0.0)*N1
              call velocitySpectrum(phi12, kx(i), ky(j), kz(k), sigma, L, 1, 2)
              h12 = phi12*sqrt( dkx*dky*dkz )/sqrt(phi11)
              call velocitySpectrum(phi22, kx(i), ky(j), kz(k), sigma, L, 2, 2)
              h22 = sqrt( phi22*dkx*dky*dkz  - h12**2 )
              w2(i,j,k)  = (cmplx(h12,0.0)*N1 + cmplx(h22,0.0)*N2)
              call velocitySpectrum(phi13, kx(i), ky(j), kz(k), sigma, L, 1, 3)
              h13 = phi13*sqrt( dkx*dky*dkz )/sqrt(phi11)
              call velocitySpectrum(phi23, kx(i), ky(j), kz(k), sigma, L, 2, 3)
              h23 = ( phi23*dkx*dky*dkz - h12*h13 )/sqrt(phi22)
              call velocitySpectrum(phi33, kx(i), ky(j), kz(k), sigma, L, 3, 3)
              h33 = sqrt( phi33*dkx*dky*dkz - (h13**2) - (h23**2) )
              w3(i,j,k)  = (cmplx(h13,0.0)*N1 + cmplx(h23,0.0)*N2 &
                + cmplx(h33,0.0)*N3)
            end do
          end do
        end do
        ! ======================================================================
      else

      end if

      deallocate ( B11, B22, B33 )

    end subroutine construct_w



end module functions
