! Simulation of a 3D turbulent velocity and temperature field *****************
!                                                                             *
!   Purpose: This code simulates a 3D turbulent field using the method        *
!            of Frelich et al. (2000) [Journal of Applied Meteorology].       *
!            In its current form, the code simulates an incompressible        *
!            isotropic field whose TKE spectrum is goverened by the           *
!            von Karman model spectrum. Future work will focus on             *
!            generalizing the code to simulate fields that are inhomogenous   *
!            in the vertical direction, and accept TKE spectra as inputs.     *
!                                                                             *
!   Author: Alexander N. Carr   ( alexcarr1721@gmail.com )                    *
!           Ph.D. Student, Theoretical Fluid Dynamics and Turbulence Group    *
!           University of Florida                                             *
!           Department of Mechanical and Aerospace Engineering                *
!                                                                             *
!******************************************************************************

program main
    use mpi
    use hdf5
    use functions
    implicit none
    ! General variables *******************************************************
    integer(isp)                :: i, proc, nproc, mpi_err, err
    character(len=80)           :: trash, groupname
    character(len=8), parameter :: filename = "turbo.h5"
    logical                     :: file_exists
    !**************************************************************************
    ! Frequency space variables ***********************************************
    real(sp), allocatable     :: k1(:), k2(:), k3(:), f1(:), f2(:), f3(:)
    !**************************************************************************
    ! RGM method variables ****************************************************
    complex(sp), allocatable  :: w1(:,:,:), w2(:,:,:), w3(:,:,:), wT(:,:,:)
    complex(sp), allocatable  :: temp1(:,:,:), temp2(:,:,:), temp3(:,:,:)
    !**************************************************************************
    ! Field variables *********************************************************
    complex(sp), allocatable  :: u1(:,:,:), u2(:,:,:), u3(:,:,:), T1(:,:,:)
    real(sp), allocatable     :: ux(:,:,:), uy(:,:,:), uz(:,:,:), Temp(:,:,:)
    real(sp), allocatable     :: x(:), y(:), z(:)
    real(sp)                  :: xmin, ymin, zmin, xmax, ymax, zmax
    integer(isp)              :: xsize, ysize, zsize
    !**************************************************************************
    ! Mean Flow variables *****************************************************
    real(sp), allocatable     :: Vx(:), Vy(:), Vz(:), T(:), sigV(:), sigT(:)
    real(sp), allocatable     :: LV(:), LT(:), h(:)
    real(sp)                  :: z0, zr, ustar, Tstar, Tr, Gammad, g, cp, Qh
    real(sp)                  :: Qe, rho0, L0, qr, L_v, q_star, wstar, zi, Ts
    !**************************************************************************

    ! Get processor information ***********************************************
    call MPI_INIT(mpi_err)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, nproc, mpi_err)
    call MPI_COMM_RANK(MPI_COMM_WORLD, proc, mpi_err)
    !**************************************************************************

    ! err = 0
    ! if ( proc .eq. 0 ) then
    !     if(command_argument_count() .ne. 1)THEN
    !         write (*,*) 'error: command line argument required, stopping.'
    !         err = 1
    !     end if
    ! end if
    ! call MPI_Bcast(err, 1, MPI_INTEGER, 0, MPI_COMM_WORLD)
    ! write (*,*) err
    ! if ( err .ne. 0 ) goto 999

    ! call get_command_argument(1,groupname)

    if ( proc .eq. 0 ) then
        open(unit = 10, file="number.inp")
        read(10,*) groupname
        groupname = trim(groupname)
    end if
    call MPI_Bcast(groupname, len(groupname), MPI_CHARACTER, 0, &
        MPI_COMM_WORLD, mpi_err)

    ! Read input file *********************************************************
    if ( proc .eq. 0 ) then
        open(unit = 10, file="input.inp")
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
    end if
    call MPI_Bcast(ustar, 1, MPI_REAL, 0, MPI_COMM_WORLD, mpi_err)
    call MPI_Bcast(wstar, 1, MPI_REAL, 0, MPI_COMM_WORLD, mpi_err)
    call MPI_Bcast(z0, 1, MPI_REAL, 0, MPI_COMM_WORLD, mpi_err)
    call MPI_Bcast(zr, 1, MPI_REAL, 0, MPI_COMM_WORLD, mpi_err)
    call MPI_Bcast(zi, 1, MPI_REAL, 0, MPI_COMM_WORLD, mpi_err)
    call MPI_Bcast(Tr, 1, MPI_REAL, 0, MPI_COMM_WORLD, mpi_err)
    call MPI_Bcast(g, 1, MPI_REAL, 0, MPI_COMM_WORLD, mpi_err)
    call MPI_Bcast(rho0, 1, MPI_REAL, 0, MPI_COMM_WORLD, mpi_err)
    call MPI_Bcast(cp, 1, MPI_REAL, 0, MPI_COMM_WORLD, mpi_err)
    call MPI_Bcast(Qh, 1, MPI_REAL, 0, MPI_COMM_WORLD, mpi_err)
    call MPI_Bcast(Ts, 1, MPI_REAL, 0, MPI_COMM_WORLD, mpi_err)
    call MPI_Bcast(Qe, 1, MPI_REAL, 0, MPI_COMM_WORLD, mpi_err)
    call MPI_Bcast(qr, 1, MPI_REAL, 0, MPI_COMM_WORLD, mpi_err)
    call MPI_Bcast(L_v, 1, MPI_REAL, 0, MPI_COMM_WORLD, mpi_err)
    call MPI_Barrier(MPI_COMM_WORLD, mpi_err)
    !**************************************************************************

    ! Set parameters **********************************************************
    Tstar = -1.0*Qh/(rho0*cp*ustar)
    L0    = -1.0*(ustar**3)*Ts*rho0*cp/(g*0.4*Qh)
    Gammad = g/cp 
    q_star = -1.0*Qe/(rho0*L_v*ustar)
    !**************************************************************************

    ! Construct Domain ********************************************************
    if ( proc .eq. 0 ) then
        open (unit = 10, file = "domain.inp" )
        read (10,*) trash, trash, trash
        read (10,*) xmin, ymin, zmin
        read (10,*) xmax, ymax, zmax
        close(10)
        open(10, file="grid.inp")
        read(10,*) trash, xsize
        read(10,*) trash, ysize 
        read(10,*) trash, zsize 
        close(10)
    end if
    call MPI_Bcast(xmin, 1, MPI_REAL, 0, MPI_COMM_WORLD, mpi_err)
    call MPI_Bcast(ymin, 1, MPI_REAL, 0, MPI_COMM_WORLD, mpi_err)
    call MPI_Bcast(zmin, 1, MPI_REAL, 0, MPI_COMM_WORLD, mpi_err)
    call MPI_Bcast(xmax, 1, MPI_REAL, 0, MPI_COMM_WORLD, mpi_err)
    call MPI_Bcast(ymax, 1, MPI_REAL, 0, MPI_COMM_WORLD, mpi_err)
    call MPI_Bcast(zmax, 1, MPI_REAL, 0, MPI_COMM_WORLD, mpi_err)
    call MPI_Bcast(xsize, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpi_err)
    call MPI_Bcast(ysize, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpi_err)
    call MPI_Bcast(zsize, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpi_err)
    call MPI_Barrier(MPI_COMM_WORLD, mpi_err)
    allocate ( x(xsize), y(ysize), z(zsize) )
    x = (/ ( ( xmin + (i-1)*( (xmax - xmin)/(size(x,dim=1) - 1) ) &
        ),i=1,size(x,dim=1) ) /)
    y = (/ ( ( ymin + (i-1)*( (ymax - ymin)/(size(y,dim=1) - 1) ) &
        ),i=1,size(y,dim=1) ) /)
    z = (/ ( ( zmin + (i-1)*( (zmax - zmin)/(size(z,dim=1) - 1) ) &
        ),i=1,size(z,dim=1) ) /)
    !**************************************************************************

    ! Initialize wavenumber and spatial arrays ********************************
    allocate ( k1(xsize), k2(ysize), k3(zsize) )
    allocate ( f1(xsize), f2(ysize), f3(zsize) )
    allocate ( ux(xsize,ysize,zsize/nproc) )
    allocate ( uy(xsize,ysize,zsize/nproc) )
    allocate ( uz(xsize,ysize,zsize/nproc) )
    allocate ( Temp(xsize,ysize,zsize/nproc) )
    allocate ( u1(xsize,ysize,zsize/nproc) )
    allocate ( u2(xsize,ysize,zsize/nproc) )
    allocate ( u3(xsize,ysize,zsize/nproc) )
    allocate ( T1(xsize,ysize,zsize/nproc) )
    allocate ( w1(size(ux,dim=3)*nproc,size(ux,dim=1),size(ux,dim=2)/nproc) )
    allocate ( w2(size(uy,dim=3)*nproc,size(uy,dim=1),size(uy,dim=2)/nproc) )
    allocate ( w3(size(uz,dim=3)*nproc,size(uz,dim=1),size(uz,dim=2)/nproc) )
    allocate ( wT(size(T1,dim=3)*nproc,size(T1,dim=1),size(T1,dim=2)/nproc) )
    allocate ( Vx(xsize), Vy(ysize), Vz(zsize) )
    allocate ( sigV(zsize), sigT(zsize), LV(zsize))
    allocate ( LT(zsize), h(zsize), T(zsize) )
    allocate ( temp1(2,2,2))
    allocate ( temp2(2,2,2))
    allocate ( temp3(2,2,2))
    do i = 1,size(z,dim=1)
        sigT(i) = 2.0 !sqrt( (Tstar**2)*( 4.0/( (1.0 + 10.0*(-1.0*z(i)/L0) &
            !)**(2.0/3.0) ) ) )
        LT(i)   = 200.0!2.0*z(i)*( ( 1.0 + 7.0*(-1.0*z(i)/L0) )/( 1.0 + &
            !10.0*(-1.0*z(i)/L0) ) )
        sigV(i) = 5.0!sqrt( 3.0*(ustar**2) + 0.35*(wstar**2) )
        LV(i)   = 300.0!1.8*z(i) + 0.23*zi 
    end do
    !**************************************************************************

    ! Model mean flow *********************************************************
    call most(Vx, Vy, T, h, z, L0, z0, zr, Tr, ustar, Tstar, Gammad, 0.0, qr,&
        q_star)
    Vz = 0.0
    !**************************************************************************

    ! Construct wavenumbers ***************************************************
    f1 = fourier_space(x)
    f2 = fourier_space(y)
    f3 = fourier_space(z)
    k1 = 2.0*pi*f1
    k2 = 2.0*pi*f2 
    k3 = 2.0*pi*f3
    !**************************************************************************

    ! Construct w *************************************************************
    call construct_w(w1, w2, w3, wT, k1, k2, k3, sigV, LV, sigT, LT, "ii", &
        MPI_COMM_WORLD)
    !**************************************************************************

    ! Compute iFFT in z ( 1st dimension of w ) ********************************
    call mkl_fft(w1(:,1,1), [size(ux,dim=3)*nproc,size(ux,dim=1),&
        size(ux,dim=2)/nproc], 1, inverse=.true. )
    call mkl_fft(w2(:,1,1), [size(uy,dim=3)*nproc,size(uy,dim=1),&
        size(uy,dim=2)/nproc], 1, inverse=.true. )
    call mkl_fft(w3(:,1,1), [size(uz,dim=3)*nproc,size(uz,dim=1),&
        size(uz,dim=2)/nproc], 1, inverse=.true. )
    call mkl_fft(wT(:,1,1), [size(T1,dim=3)*nproc,size(T1,dim=1),&
        size(T1,dim=2)/nproc], 1, inverse=.true. )
    !**************************************************************************

    ! Transpose out of place **************************************************
    call transpose(u1(:,1,1), temp1(:,1,1), w1(:,1,1), [size(ux,dim=1),&
        size(ux,dim=2),size(ux,dim=3)], [2,2,2], [size(ux,dim=3)*nproc,&
        size(ux,dim=1),size(ux,dim=2)/nproc], 231, MPI_COMM_WORLD )
    call transpose(u2(:,1,1), temp2(:,1,1), w2(:,1,1), [size(uy,dim=1),&
        size(uy,dim=2),size(uy,dim=3)], [2,2,2], [size(uy,dim=3)*nproc,&
        size(uy,dim=1),size(uy,dim=2)/nproc], 231, MPI_COMM_WORLD  )
    call transpose(u3(:,1,1), temp3(:,1,1), w3(:,1,1), [size(uz,dim=1),&
        size(uz,dim=2),size(uz,dim=3)], [2,2,2], [size(uz,dim=3)*nproc,&
        size(uz,dim=1),size(uz,dim=2)/nproc], 231, MPI_COMM_WORLD  )
    call transpose(T1(:,1,1), temp3(:,1,1), wT(:,1,1), [size(T1,dim=1),&
        size(T1,dim=2),size(T1,dim=3)], [2,2,2], [size(T1,dim=3)*nproc,&
        size(T1,dim=1),size(T1,dim=2)/nproc], 231, MPI_COMM_WORLD  )
    !**************************************************************************

    ! Compute iFFT in x, y ****************************************************
    call mkl_fft(u1(:,1,1), [size(ux,dim=1),size(ux,dim=2),size(ux,dim=3)],&
        2, inverse=.true. )
    call mkl_fft(u2(:,1,1), [size(uy,dim=1),size(uy,dim=2),size(uy,dim=3)],&
        2, inverse=.true. )
    call mkl_fft(u3(:,1,1), [size(ux,dim=1),size(ux,dim=2),size(ux,dim=3)],&
        2, inverse=.true. )
    call mkl_fft(T1(:,1,1), [size(T1,dim=1),size(T1,dim=2),size(T1,dim=3)],&
        2, inverse=.true. )
    !**************************************************************************

    ux = real(u1)
    uy = real(u2)
    uz = real(u3)
    Temp = real(T1)

    ! Write to file ***********************************************************
    inquire(file="turbo.h5",exist=file_exists)
    if ( file_exists ) then
        call create_h5_f(filename, MPI_COMM_WORLD )
        call create_h5_d(filename, trim(groupname)//"/"//"ux", [size(ux,dim=1),&
            size(uy,dim=2), size(uz,dim=3)*nproc], MPI_COMM_WORLD, ux(:,1,1),&
            chunked=.true., dim_chunk=shape(ux), &
            offset=[0, 0, proc*size(ux,dim=3)])
        call create_h5_d(filename, trim(groupname)//"/"//"uy", [size(uy,dim=1),&
            size(uy,dim=2), size(uy,dim=3)*nproc], MPI_COMM_WORLD, uy(:,1,1),&
            chunked=.true., dim_chunk=shape(uy), &
            offset=[0, 0, proc*size(uy,dim=3)])
        call create_h5_d(filename, trim(groupname)//"/"//"uz", [size(uz,dim=1),&
            size(uz,dim=2), size(uz,dim=3)*nproc], MPI_COMM_WORLD, uz(:,1,1),&
            chunked=.true., dim_chunk=shape(uz), &
            offset=[0, 0, proc*size(uz,dim=3)])
        call create_h5_d(filename, trim(groupname)//"/"//"T", [size(T1,dim=1),&
            size(T1,dim=2), size(T1,dim=3)*nproc], MPI_COMM_WORLD, Temp(:,1,1),&
            chunked=.true., dim_chunk=shape(Temp), &
            offset=[0, 0,proc*size(T1,dim=3)])
    else
        call create_h5_d(filename, trim(groupname)//"/"//"ux", [size(ux,dim=1),&
            size(uy,dim=2), size(uz,dim=3)*nproc], MPI_COMM_WORLD, ux(:,1,1),&
            chunked=.true., dim_chunk=shape(ux), &
            offset=[0, 0, proc*size(ux,dim=3)])
        call create_h5_d(filename, trim(groupname)//"/"//"uy", [size(uy,dim=1),&
            size(uy,dim=2), size(uy,dim=3)*nproc], MPI_COMM_WORLD, uy(:,1,1),&
            chunked=.true., dim_chunk=shape(uy), &
            offset=[0, 0, proc*size(uy,dim=3)])
        call create_h5_d(filename, trim(groupname)//"/"//"uz", [size(uz,dim=1),&
            size(uz,dim=2), size(uz,dim=3)*nproc], MPI_COMM_WORLD, uz(:,1,1),&
            chunked=.true., dim_chunk=shape(uz), &
            offset=[0, 0, proc*size(uz,dim=3)])
        call create_h5_d(filename, trim(groupname)//"/"//"T", [size(T1,dim=1),&
            size(T1,dim=2), size(T1,dim=3)*nproc], MPI_COMM_WORLD, Temp(:,1,1),&
            chunked=.true., dim_chunk=shape(Temp), &
            offset=[0, 0,proc*size(T1,dim=3)])
    end if
    !**************************************************************************


    deallocate ( k1, k2, k3, f1, f2, f3, ux, uy, uz, Temp, u1, u2, u3, T1, &
        w1, w2, w3, wT, x, y, z )
    call MPI_FINALIZE(mpi_err)

end program main