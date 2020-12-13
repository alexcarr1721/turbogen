! Test of Random Generation Method

program main
    use hdf5
    use functions
    implicit none
    ! General variables ==========================================================
    integer(isp)              :: i, j, k, ii, jj, kk
    real(sp)                  :: a1(1), b1(1), a2(1), b2(1), a(1), b(1)
    character(len=80)         :: trash
    ! ============================================================================
    ! Frequency space variables ==================================================
    real(sp), parameter       :: k_max = 2.0_sp*pi
    real(sp), allocatable     :: kx(:), ky(:), kz(:)
    integer(isp)              :: domain_size, modes
    ! ============================================================================
    ! RGM method variables =======================================================
    complex(sp), allocatable  :: w1(:,:,:), w2(:,:,:), w3(:,:,:)
    real(sp)                  :: h11, h12, h22, h13, h23, h33, phi11, phi12, phi22
    real(sp)                  :: phi13, phi23, phi33
    complex(sp)               :: N1, N2, N3
    ! ============================================================================
    ! Field variables ============================================================
    complex(sp), allocatable  :: u1(:,:,:), u2(:,:,:), u3(:,:,:)
    real(sp), allocatable     :: ux(:,:,:), uy(:,:,:), uz(:,:,:)
    real(sp), allocatable     :: x(:), y(:), z(:)
    real(sp)                  :: box_size
    ! ============================================================================
    ! Mean Flow variables ********************************************************
    real(sp), allocatable     :: Vx(:), Vy(:), Vz(:), T(:), sigV(:), sigT(:)
    real(sp), allocatable     :: LV(:), LT(:), h(:)
    real(sp)                  :: z0, zr, ustar, Tstar, Tr, Gammad, g, cp, Qh, Ts
    real(sp)                  :: Qe, rho0, L0, qr, L_v, q_star, wstar, zi
    ! Measurement variables ======================================================
    real(sp)                  :: L, mean, sigma
    ! ============================================================================
    ! HDF5 variables =============================================================
    CHARACTER(LEN=12), PARAMETER :: filename = "turbo.h5"
    INTEGER, PARAMETER  :: rank = 3   ! Rank of the data set

    INTEGER(hid_t) :: file_id, dataset_id, dataspace_id ! Identifiers
    INTEGER(hid_t) :: dataset_id1, dataspace_id1 ! Identifiers
    INTEGER(hid_t) :: plist_id ! Property list identifier

    INTEGER :: herror, chunk
    INTEGER(hsize_t), DIMENSION(1:rank) :: dims ! dimensions of data
    INTEGER(hsize_t), DIMENSION(1:rank) :: cdims ! sizes of chunked data
    INTEGER(HSIZE_T), DIMENSION(1)      :: data_dims ! dimensions of data buffers
    ! ============================================================================

    ! Read input file *********************************************************
    open(10, file="input.inp")
    read(10,*) trash, ustar 
    read(10,*) trash, wstar
    read(10,*) trash, z0 
    read(10,*) trash, zr 
    read(10,*) trash, zi
    read(10,*) trash, Tr
    read(10,*) trash, g 
    read(10,*) trash, rho0
    read(10,*) trash, cp
    read(10,*) trash, Qh 
    read(10,*) trash, Ts 
    read(10,*) trash, Qe 
    read(10,*) trash, qr 
    read(10,*) trash, L_v
    close(10)
    !**************************************************************************

    ! Set parameters **********************************************************
    Tstar = -1.0*Qh/(rho0*cp*ustar)
    L0    = -1.0*(ustar**3)*Ts*rho0*cp/(g*0.4*Qh)
    Gammad = g/cp 
    q_star = -1.0*Qe/(rho0*L_v*ustar)
    !**************************************************************************

    ! From measurements ==========================================================
    mean  = 0.0     ! m/s
    sigma = 22.2   ! m/s
    L     = 2.40  ! m
    ! ============================================================================

    ! Initialize wavenumber and spatial arrays ===================================
    box_size    = 30.0*L
    domain_size = 512
    chunk       = domain_size/16
    allocate ( kx(domain_size), ky(domain_size), kz(domain_size) )
    allocate ( ux(domain_size,domain_size,domain_size) )
    allocate ( uy(domain_size,domain_size,domain_size) )
    allocate ( uz(domain_size,domain_size,domain_size) )
    allocate ( u1(domain_size,domain_size,domain_size) )
    allocate ( u2(domain_size,domain_size,domain_size) )
    allocate ( u3(domain_size,domain_size,domain_size) )
    allocate ( w1(domain_size,domain_size,domain_size) )
    allocate ( w2(domain_size,domain_size,domain_size) )
    allocate ( w3(domain_size,domain_size,domain_size) )
    allocate ( x(domain_size), y(domain_size), z(domain_size) )
    allocate ( Vx(domain_size), Vy(domain_size), Vz(domain_size) )
    allocate ( sigV(domain_size), sigT(domain_size), LV(domain_size))
    allocate ( LT(domain_size), h(domain_size), T(domain_size) )
    do i = 1,domain_size
        x(i) = (i-1)*(1.0/(domain_size - 1.0) )*box_size
        y(i) = (i-1)*(1.0/(domain_size - 1.0) )*box_size
        z(i) = (i-1)*(1.0/(domain_size - 1.0) )*box_size
    end do
    do i = 1,size(z,dim=1)
        sigT(i) = sqrt( (Tstar**2)*( 4.0/( (1.0 + 10.0*(-1.0*z(i)/L0) &
            )**(2.0/3.0) ) ) )
        LT(i)   = 2.0*z(i)*( ( 1.0 + 7.0*(-1.0*z(i)/L0) )/( 1.0 + &
            10.0*(-1.0*z(i)/L0) ) )
        sigV(i) = sqrt( 3.0*(ustar**2) + 0.35*(wstar**2) )
        LV(i)   = 1.8*z(i) + 0.23*zi 
    end do

    ! Model mean flow *********************************************************
    call most(Vx, Vy, T, h, z, L0, z0, zr, Tr, ustar, Tstar, Gammad, 0.0, qr,&
        q_star)
    Vz = 0.0
    do i = 1,size(Vx,dim=1)
        write(*,*) sigT(i), LT(i), sigV(i), LV(i)
    end do
    !**************************************************************************

    call construct_wavenumbers(kx, ky, kz, x, y, z)

    ! ============================================================================

    ! Construct w ================================================================
    call construct_w(w1, w2, w3, x, y, z, kx, ky, kz, mean, sigma, L, "uncorrelated")
    ! ============================================================================

    ! Inverse FFT ================================================================
    call fft( w1(:,1,1), [domain_size, domain_size, domain_size], 3, .true. )
    call fft( w2(:,1,1), [domain_size, domain_size, domain_size], 3, .true. )
    call fft( w3(:,1,1), [domain_size, domain_size, domain_size], 3, .true. )
    ux = real(w1)
    uy = real(w2)
    uz = real(w3)
    ! ============================================================================

    ! HDF5 storage ===============================================================
    call h5open_f(herror)
    !
    ! Create a file
    call h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, herror)
    !
    ! Create dataset "Compressed Data" in the group using absolute name.
    dims(1:3) = [domain_size, domain_size, domain_size]
    call h5screate_simple_f(rank, dims, dataspace_id, herror)
    call h5pcreate_f(H5P_DATASET_CREATE_F, plist_id, herror)
    !
    ! Dataset must be chunked for compression
    cdims(1:3) = [chunk, chunk, chunk]
    CALL h5pset_chunk_f(plist_id, 3, cdims, herror)

    ! Set ZLIB / DEFLATE Compression using compression level 6.
    ! To use SZIP Compression comment out these lines.
    CALL h5pset_deflate_f(plist_id, 6, herror)

    ! Create data set
    CALL h5dcreate_f(file_id, "ux", H5T_NATIVE_REAL, dataspace_id, &
        dataset_id, herror, dcpl_id=plist_id)

    CALL h5dwrite_f(dataset_id, H5T_NATIVE_REAL, ux, dims, herror)

    ! Create data set
    CALL h5dcreate_f(file_id, "uy", H5T_NATIVE_REAL, dataspace_id, &
        dataset_id, herror, dcpl_id=plist_id)

    CALL h5dwrite_f(dataset_id, H5T_NATIVE_REAL, uy, dims, herror)

    ! Create data set
    CALL h5dcreate_f(file_id, "uz", H5T_NATIVE_REAL, dataspace_id, &
        dataset_id, herror, dcpl_id=plist_id)

    CALL h5dwrite_f(dataset_id, H5T_NATIVE_REAL, uz, dims, herror)

    ! Create data set
    data_dims = [domain_size]
    call h5screate_simple_f(1, data_dims, dataspace_id1, herror)
    CALL h5dcreate_f(file_id, "x", H5T_NATIVE_REAL, dataspace_id1, &
        dataset_id1, herror)

    CALL h5dwrite_f(dataset_id1, H5T_NATIVE_REAL, x, data_dims, herror)

    ! Create data set
    CALL h5dcreate_f(file_id, "y", H5T_NATIVE_REAL, dataspace_id1, &
        dataset_id1, herror)

    CALL h5dwrite_f(dataset_id1, H5T_NATIVE_REAL, y, data_dims, herror)

    ! Create data set
    CALL h5dcreate_f(file_id, "z", H5T_NATIVE_REAL, dataspace_id1, &
        dataset_id1, herror)

    CALL h5dwrite_f(dataset_id1, H5T_NATIVE_REAL, z, data_dims, herror)

    ! Close resources
    CALL h5sclose_f(dataspace_id, herror)
    CALL h5pclose_f(plist_id, herror)
    CALL h5dclose_f(dataset_id, herror)
    CALL h5sclose_f(dataspace_id1, herror)
    CALL h5dclose_f(dataset_id1, herror)
    CALL h5fclose_f(file_id, herror)
    ! ============================================================================

    deallocate ( kx, ky, kz, ux, uy, uz, u1, u2, u3, w1, w2, w3, x, y, z )


end program main
