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
    integer(isp)              :: i, j, k, proc, nproc, mpi_err, count, stat
    integer(isp)              :: zglobal
    character(len=80)         :: trash, groupname, filename, case, mkdirCmd
    !**************************************************************************
    ! Frequency space variables ***********************************************
    real(sp), allocatable     :: k1(:), k2(:), k3(:), f1(:), f2(:), f3(:)
    !**************************************************************************
    ! RGM method variables ****************************************************
    complex(sp), allocatable  :: w1(:,:,:), w2(:,:,:), w3(:,:,:), wT(:,:,:)
    ! complex(sp), allocatable  :: wTz(:,:,:)
    complex(sp), allocatable  :: temp1(:,:,:), temp2(:,:,:), temp3(:,:,:)
    !**************************************************************************
    ! Field variables *********************************************************
    complex(sp), allocatable  :: u1(:,:,:), u2(:,:,:), u3(:,:,:), T1(:,:,:)
    ! complex(sp), allocatable  :: T1y(:,:,:), T1x(:,:,:)!, T1z(:,:,:)
    ! real(sp), allocatable     :: Tempy(:,:,:), Tempx(:,:,:), Tempz(:,:,:)
    real(sp), allocatable     :: ux(:,:,:), uy(:,:,:), uz(:,:,:), Temp(:,:,:)
    real(sp), allocatable     :: cturb(:,:,:), rhoturb(:,:,:)!, rhoturbx(:,:,:)
    ! real(sp), allocatable     :: rhoturby(:,:,:), rhoturbz(:,:,:)
    real(sp), allocatable     :: x(:), y(:), z(:), tau(:)
    real(sp)                  :: xmin, ymin, zmin, xmax, ymax, zmax
    integer(isp)              :: xsize, ysize, zsize, tsize
    real(sp)                  :: mic_loc(3,1000)
    !**************************************************************************
    ! Mean Flow variables *****************************************************
    real(sp)                  :: Ltemp, LTtemp, sigtemp, sigTtemp
    real(sp)                  :: sigV, sigT, LV, LT
    real(sp)                  :: c0, rho0, p0, T0, h0
    real(sp)                  :: psat, pvap, pdry
    real(sp), allocatable     :: meanzeros(:,:), hbar(:), pbar(:), Tbar(:)
    real(sp), allocatable     :: rhobar(:), cbar(:), p(:,:,:), pressure(:)
    !**************************************************************************

    ! Get processor information ***********************************************
    call MPI_INIT(mpi_err)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, nproc, mpi_err)
    call MPI_COMM_RANK(MPI_COMM_WORLD, proc, mpi_err)
    !**************************************************************************

    if (command_argument_count() .ne. 1) then
        groupname = "1"
    else
        call get_command_argument(1,groupname)
    end if
    call MPI_Barrier(MPI_COMM_WORLD, mpi_err)

    ! Read input file *********************************************************
    count = 0
    if ( proc .eq. 0 ) then
        open(unit = 10, file="input.inp")
        read(10,*) trash
        read(10,*) trash, case
        read(10,*) trash
        read(10,*) trash, xsize 
        read(10,*) trash, ysize
        read(10,*) trash, zsize
        read(10,*) trash, tsize
        read(10,*) trash
        read(10,*) trash, xmin
        read(10,*) trash, xmax
        read(10,*) trash, ymin
        read(10,*) trash, ymax
        read(10,*) trash, zmin
        read(10,*) trash, zmax
        read(10,*) trash
        read(10,*) trash, sigtemp
        read(10,*) trash, Ltemp 
        read(10,*) trash, sigTtemp 
        read(10,*) trash, LTtemp
        read(10,*) trash
        read(10,*) trash
        do while (.not. IS_IOSTAT_END(stat))
            count = count + 1
            read (10,*,IOSTAT=stat) mic_loc(1,count), mic_loc(2,count), &
                mic_loc(3,count)
        end do
        close(10)
    end if
    call MPI_Bcast(case, 80, MPI_CHARACTER, 0, MPI_COMM_WORLD, mpi_err)
    case = trim(case)
    call MPI_Bcast(xmin, 1, MPI_REAL, 0, MPI_COMM_WORLD, mpi_err)
    call MPI_Bcast(ymin, 1, MPI_REAL, 0, MPI_COMM_WORLD, mpi_err)
    call MPI_Bcast(zmin, 1, MPI_REAL, 0, MPI_COMM_WORLD, mpi_err)
    call MPI_Bcast(xmax, 1, MPI_REAL, 0, MPI_COMM_WORLD, mpi_err)
    call MPI_Bcast(ymax, 1, MPI_REAL, 0, MPI_COMM_WORLD, mpi_err)
    call MPI_Bcast(zmax, 1, MPI_REAL, 0, MPI_COMM_WORLD, mpi_err)
    call MPI_Bcast(xsize, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpi_err)
    call MPI_Bcast(ysize, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpi_err)
    call MPI_Bcast(zsize, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpi_err)
    call MPI_Bcast(tsize, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpi_err)
    call MPI_Bcast(sigtemp, 1, MPI_REAL, 0, MPI_COMM_WORLD, mpi_err)
    call MPI_Bcast(Ltemp, 1, MPI_REAL, 0, MPI_COMM_WORLD, mpi_err)
    call MPI_Bcast(sigTtemp, 1, MPI_REAL, 0, MPI_COMM_WORLD, mpi_err)
    call MPI_Bcast(LTtemp, 1, MPI_REAL, 0, MPI_COMM_WORLD, mpi_err)
    call MPI_Bcast(count , 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpi_err)
    call MPI_Bcast(mic_loc, count*3, MPI_REAL, 0, MPI_COMM_WORLD, mpi_err)
    call MPI_Barrier(MPI_COMM_WORLD, mpi_err)
    !**************************************************************************

    ! Set parameters **********************************************************
    sigT = sigTtemp
    LT   = LTtemp*(1.0/0.746834)
    sigV = sigtemp
    LV   = Ltemp*(1.0/0.746834)
    !**************************************************************************

    ! Construct Domain ********************************************************
    allocate ( x(xsize), y(ysize), z(zsize), tau(tsize) )
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
    allocate ( ux(ysize,zsize,xsize/nproc) )
    allocate ( uy(ysize,zsize,xsize/nproc) )
    allocate ( uz(ysize,zsize,xsize/nproc) )
    allocate ( Temp(ysize,zsize,xsize/nproc) )
    allocate ( u1(ysize,zsize,xsize/nproc) )
    allocate ( u2(ysize,zsize,xsize/nproc) )
    allocate ( u3(ysize,zsize,xsize/nproc) )
    allocate ( T1(ysize,zsize,xsize/nproc) )
    allocate ( w1(size(ux,dim=3)*nproc,size(ux,dim=1),size(ux,dim=2)/nproc) )
    allocate ( w2(size(uy,dim=3)*nproc,size(uy,dim=1),size(uy,dim=2)/nproc) )
    allocate ( w3(size(uz,dim=3)*nproc,size(uz,dim=1),size(uz,dim=2)/nproc) )
    allocate ( wT(size(T1,dim=3)*nproc,size(T1,dim=1),size(T1,dim=2)/nproc) )
    ! allocate ( wTz(size(T1,dim=3)*nproc,size(T1,dim=1),size(T1,dim=2)/nproc) )
    ! allocate ( T1z(xsize,ysize,zsize/nproc) )
    ! allocate ( T1y(xsize,ysize,zsize/nproc) )
    ! allocate ( T1x(xsize,ysize,zsize/nproc) )
    ! allocate ( Tempz(xsize,ysize,zsize/nproc) )
    ! allocate ( Tempy(xsize,ysize,zsize/nproc) )
    ! allocate ( Tempx(xsize,ysize,zsize/nproc) )
    allocate ( temp1(2,2,2))
    allocate ( temp2(2,2,2))
    allocate ( temp3(2,2,2))
    allocate ( meanzeros(xsize, ysize), hbar(xsize), pbar(xsize), Tbar(xsize) )
    allocate ( cbar(xsize), rhobar(xsize), p(tsize, ysize, zsize/nproc) )
    allocate ( pressure(tsize), cturb(ysize, zsize, xsize/nproc) )
    allocate ( rhoturb(ysize, zsize, xsize/nproc) )
    ! allocate ( rhoturbx(xsize, ysize, zsize/nproc) )
    ! allocate ( rhoturby(xsize, ysize, zsize/nproc) )
    ! allocate ( rhoturbz(xsize, ysize, zsize/nproc) )
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
    call construct_w(w1, w2, w3, wT, k1, k2, k3, sigV, LV, sigT, LT, "ii", &  ! y is the last index of w, that is why k2 is considered kz here
        MPI_COMM_WORLD)
    !**************************************************************************

    ! do k = 1,size(w1,dim=3)
    !     do j = 1,size(w1,dim=2)
    !         do i = 1,size(w1,dim=1)
    !             wTz(i,j,k) = img*k3(i)*wT(i,j,k)
    !         end do
    !     end do
    ! end do

    ! Compute iFFT in x ( 1st dimension of w ) ********************************
    call mkl_fft(w1(:,1,1), [size(ux,dim=3)*nproc,size(ux,dim=1),&
        size(ux,dim=2)/nproc], 1, inverse=.true. )
    call mkl_fft(w2(:,1,1), [size(uy,dim=3)*nproc,size(uy,dim=1),&
        size(uy,dim=2)/nproc], 1, inverse=.true. )
    call mkl_fft(w3(:,1,1), [size(uz,dim=3)*nproc,size(uz,dim=1),&
        size(uz,dim=2)/nproc], 1, inverse=.true. )
    call mkl_fft(wT(:,1,1), [size(T1,dim=3)*nproc,size(T1,dim=1),&
        size(T1,dim=2)/nproc], 1, inverse=.true. )
    ! call mkl_fft(wTz(:,1,1), [size(T1,dim=3)*nproc,size(T1,dim=1),&
    !     size(T1,dim=2)/nproc], 1, inverse=.true. )
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
    ! call transpose(T1z(:,1,1), temp3(:,1,1), wTz(:,1,1), [size(T1z,dim=1),&
    !     size(T1z,dim=2),size(T1z,dim=3)], [2,2,2], [size(T1z,dim=3)*nproc,&
    !     size(T1z,dim=1),size(T1z,dim=2)/nproc], 231, MPI_COMM_WORLD  )
    !**************************************************************************

    ! do k = 1,size(T1y,dim=3)
    !     do j = 1,size(T1y,dim=2)
    !         do i = 1,size(T1y,dim=1)
    !             T1y(i,j,k) = img*k2(j)*T1(i,j,k)
    !             T1x(i,j,k) = img*k1(i)*T1(i,j,k)
    !         end do 
    !     end do 
    ! end do

    ! Compute iFFT in z, y ****************************************************
    call mkl_fft(u1(:,1,1), [size(ux,dim=1),size(ux,dim=2),size(ux,dim=3)],&
        2, inverse=.true. )
    call mkl_fft(u2(:,1,1), [size(uy,dim=1),size(uy,dim=2),size(uy,dim=3)],&
        2, inverse=.true. )
    call mkl_fft(u3(:,1,1), [size(ux,dim=1),size(ux,dim=2),size(ux,dim=3)],&
        2, inverse=.true. )
    call mkl_fft(T1(:,1,1), [size(T1,dim=1),size(T1,dim=2),size(T1,dim=3)],&
        2, inverse=.true. )
    ! call mkl_fft(T1z(:,1,1), [size(T1z,dim=1),size(T1z,dim=2),size(T1z,dim=3)],&
    !     2, inverse=.true. )
    ! call mkl_fft(T1y(:,1,1), [size(T1y,dim=1),size(T1y,dim=2),size(T1y,dim=3)],&
    !     2, inverse=.true. )
    ! call mkl_fft(T1x(:,1,1), [size(T1x,dim=1),size(T1x,dim=2),size(T1x,dim=3)],&
    !     2, inverse=.true. )
    !**************************************************************************

    ux = real(u1)
    uy = real(u2)
    uz = real(u3)
    Temp = real(T1)
    ! Tempx = real(T1x)
    ! Tempy = real(T1y)
    ! Tempz = real(T1z)

    ! Set mean values
    c0 = 343.0 ! m/s
    rho0 = 1.225 
    T0 = 288.15
    p0 = 101325.0
    h0 = 0
    pbar = p0
    Tbar = T0
    hbar = h0
    meanzeros = 0.0

    ! Turbulent fields of density and speed of sound **************************
    ! Mean flow
    do i = 1,size(x,dim=1)
        psat = 6.1078_sp * 10.0_sp**(7.5_sp*(Tbar(i))/( &
            (Tbar(i) ) + 237.3_sp))
        pvap = psat*hbar(i)
        pdry = p0 - pvap 
        rhobar(i) = pdry/(287.058*( Tbar(i) ) ) + &
            pvap/(461.495*( Tbar(i) ) )
        cbar(i) = sqrt( 1.4*p0/rhobar(i) )
    end do
    do k = 1,size(ux,dim=3)
        zglobal = k + proc*size(ux,dim=3)
        do j = 1,size(ux,dim=2)
            do i = 1,size(ux,dim=1)
                ! Total field *************************************************
                psat = 6.1078_sp * 10.0_sp**(7.5_sp*(Tbar(zglobal) + &
                    Temp(i,j,k))/( (Tbar(zglobal) + Temp(i,j,k) ) + 237.3_sp))
                pvap = psat*hbar(zglobal)
                pdry = p0 - pvap 
                rhoturb(i,j,k) = pdry/(287.058*( Tbar(zglobal) + Temp(i,j,k) ) ) + &
                    pvap/(461.495*( Tbar(zglobal) + Temp(i,j,k) ) )
                cturb(i,j,k) = sqrt( 1.4*p0/rhoturb(i,j,k) )
                ! Subtract mean flow from turbulent variables to get true *****
                ! fluctuations ************************************************
                rhoturb(i,j,k) = rhoturb(i,j,k) - rhobar(zglobal)
                cturb(i,j,k)   = cturb(i,j,k) - cbar(zglobal)
                ! rhoturbx(i,j,k) = (-1.0/( Tbar(i) + Temp(i,j,k) ))*( &
                !     (pdry/287.058) + (pvap/461.495) )*Tempx(i,j,k)
                ! rhoturby(i,j,k) = (-1.0/( Tbar(i) + Temp(i,j,k) ))*( &
                !     (pdry/287.058) + (pvap/461.495) )*Tempy(i,j,k)
                ! rhoturbz(i,j,k) = (-1.0/( Tbar(i) + Temp(i,j,k) ))*( &
                !     (pdry/287.058) + (pvap/461.495) )*Tempz(i,j,k)
            end do 
        end do
    end do
    !**************************************************************************

    ! Pressure input **********************************************************
    call read_waveform("waveform.inp", pressure, tau, MPI_COMM_WORLD)
    do j = 1,zsize/nproc
        do i = 1,ysize 
            p(:,i,j) = pressure(:)
        end do 
    end do

    ! ! Write to text file
    ! open(unit=10, file="F18_waveform.dat")
    ! write (10,*) "t (s)   ", "pressure (pa)   ", "sampling frequency (Hz) = ", size(tau,1)/( tau(size(tau,1)) - tau(1) )
    ! do i = 1,size(tau,dim=1)
    !     write (10,*) tau(i), p(i,1,1)
    ! end do
    ! close(10)
    !**************************************************************************

    ! Write to file ***********************************************************
    mkdirCmd = 'mkdir -p '//trim(case)
    call system( mkdirCmd )
    filename = trim(case)//"/"//"field"//trim(groupname)//".h5"
    call create_h5_f(filename, MPI_COMM_WORLD, &
        (/"pressure   ","grid       ","turb       ","mean       ",&
        "microphones","routines   "/) )
    call create_h5_d(filename, "pressure/initial", [size(p,dim=1),&
        size(p,dim=2), size(p,dim=3)*nproc], MPI_COMM_WORLD, p(:,1,1),&
        chunked=.true., dim_chunk=shape(p), &
        offset=[0, 0, proc*size(p,dim=3)])
    call create_h5_d(filename, "turb/ux1", [size(ux,dim=1),&
        size(uy,dim=2), size(uz,dim=3)*nproc], MPI_COMM_WORLD, ux(:,1,1),&
        chunked=.true., dim_chunk=shape(ux), &
        offset=[0, 0, proc*size(ux,dim=3)])
    call create_h5_d(filename, "turb/ux2", [size(uy,dim=1),&
        size(uy,dim=2), size(uy,dim=3)*nproc], MPI_COMM_WORLD, uy(:,1,1),&
        chunked=.true., dim_chunk=shape(uy), &
        offset=[0, 0, proc*size(uy,dim=3)])
    call create_h5_d(filename, "turb/ux3", [size(uz,dim=1),&
        size(uz,dim=2), size(uz,dim=3)*nproc], MPI_COMM_WORLD, uz(:,1,1),&
        chunked=.true., dim_chunk=shape(uz), &
        offset=[0, 0, proc*size(uz,dim=3)])
    call create_h5_d(filename, "turb/temperature", [size(T1,dim=1),&
        size(T1,dim=2), size(T1,dim=3)*nproc], MPI_COMM_WORLD, Temp(:,1,1),&
        chunked=.true., dim_chunk=shape(Temp), &
        offset=[0, 0,proc*size(T1,dim=3)])
    call create_h5_d(filename, "turb/c", [size(cturb,dim=1),&
        size(cturb,dim=2), size(cturb,dim=3)*nproc], MPI_COMM_WORLD, &
        cturb(:,1,1), chunked=.true., dim_chunk=shape(cturb), &
        offset=[0, 0,proc*size(cturb,dim=3)])
    call create_h5_d(filename, "turb/rho", [size(rhoturb,dim=1),&
        size(rhoturb,dim=2), size(rhoturb,dim=3)*nproc], MPI_COMM_WORLD, &
        rhoturb(:,1,1), chunked=.true., dim_chunk=shape(rhoturb), &
        offset=[0, 0,proc*size(rhoturb,dim=3)])
    ! call create_h5_d(filename, "turb/drhodx1", [size(rhoturbx,dim=1),&
    !     size(rhoturbx,dim=2), size(rhoturbx,dim=3)*nproc], MPI_COMM_WORLD, &
    !     rhoturbx(:,1,1), chunked=.true., dim_chunk=shape(rhoturbx), &
    !     offset=[0, 0,proc*size(rhoturbx,dim=3)])
    ! call create_h5_d(filename, "turb/drhodx2", [size(rhoturby,dim=1),&
    !     size(rhoturby,dim=2), size(rhoturby,dim=3)*nproc], MPI_COMM_WORLD, &
    !     rhoturby(:,1,1), chunked=.true., dim_chunk=shape(rhoturby), &
    !     offset=[0, 0,proc*size(rhoturby,dim=3)])
    ! call create_h5_d(filename, "turb/drhodx3", [size(rhoturbz,dim=1),&
    !     size(rhoturbz,dim=2), size(rhoturbz,dim=3)*nproc], MPI_COMM_WORLD, &
    !     rhoturbz(:,1,1), chunked=.true., dim_chunk=shape(rhoturbz), &
    !     offset=[0, 0,proc*size(rhoturbz,dim=3)])
    call create_h5_d(filename, "mean/sigV", [1], &
        MPI_COMM_WORLD,[sigV])
    call create_h5_d(filename, "mean/sigT", [1], &
        MPI_COMM_WORLD,[sigT])
    call create_h5_d(filename, "mean/LV", [1], MPI_COMM_WORLD, [LV])
    call create_h5_d(filename, "mean/LT", [1], MPI_COMM_WORLD, [LT])
    call create_h5_d(filename, "mean/Vx1", [size(meanzeros,dim=1), &
        size(meanzeros,dim=2)], MPI_COMM_WORLD, meanzeros(:,1))
    call create_h5_d(filename, "mean/Vx2", [size(meanzeros,dim=1), &
        size(meanzeros,dim=2)], MPI_COMM_WORLD,  meanzeros(:,1))
    call create_h5_d(filename, "mean/Vx3", [size(meanzeros,dim=1), &
        size(meanzeros,dim=2)], MPI_COMM_WORLD, meanzeros(:,1))
    call create_h5_d(filename, "mean/Tbar", [size(Tbar,dim=1)], MPI_COMM_WORLD, Tbar)
    call create_h5_d(filename, "mean/hbar", [size(hbar,dim=1)], MPI_COMM_WORLD, hbar)
    call create_h5_d(filename, "mean/cbar", [size(cbar,dim=1)], &
        MPI_COMM_WORLD, cbar)
    call create_h5_d(filename, "mean/rhobar", [size(rhobar,dim=1)], &
        MPI_COMM_WORLD, rhobar)
    call create_h5_d(filename, "mean/pbar", [size(pbar,dim=1)], &
        MPI_COMM_WORLD, pbar)
    call create_h5_d(filename, "mean/c0", [ 1 ], &
        MPI_COMM_WORLD, [c0] )
    call create_h5_d(filename, "mean/rho0", [ 1 ], &
        MPI_COMM_WORLD, [rho0] )
    call create_h5_d(filename, "grid/x1", [size(x,dim=1)], MPI_COMM_WORLD, x)
    call create_h5_d(filename, "grid/x2", [size(y,dim=1)], MPI_COMM_WORLD, y)
    call create_h5_d(filename, "grid/x3", [size(z,dim=1)], MPI_COMM_WORLD, z)
    call create_h5_d(filename, "grid/t", [size(tau,dim=1)], MPI_COMM_WORLD, &
        tau)
    call create_h5_d(filename,"routines/nonlinear",(/1/), &
        MPI_COMM_WORLD, [1.0] )
    call create_h5_d(filename,"routines/diffraction",(/1/), &
        MPI_COMM_WORLD, [1.0] )
    call create_h5_d(filename,"routines/absorption",(/1/), &
        MPI_COMM_WORLD, [1.0] )
    call create_h5_d(filename,"routines/phase",(/1/), &
        MPI_COMM_WORLD, [1.0] )
    call create_h5_d(filename,"routines/coupling",(/1/), &
        MPI_COMM_WORLD, [1.0] )
    call create_h5_d(filename,"routines/boundaryconditions",(/1/), &
        MPI_COMM_WORLD, [0.0] )
    call create_h5_d(filename,"microphones/location",(/3, count-1/), &
        MPI_COMM_WORLD, mic_loc(:,1))
    call create_h5_d(filename,"propagation_direction",(/1/), MPI_COMM_WORLD,&
        (/0.0/) )
    !**************************************************************************


    deallocate ( k1, k2, k3, f1, f2, f3, ux, uy, uz, Temp, u1, u2, u3, T1, &
        w1, w2, w3, wT, x, y, z )
    call MPI_FINALIZE(mpi_err)

end program main