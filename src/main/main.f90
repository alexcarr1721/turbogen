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
    integer(isp)                :: i, proc, nproc, mpi_err
    character(len=80)           :: trash, groupname, filename, case, mkdirCmd
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
    real(sp)                  :: Ltemp, LTtemp, sigtemp, sigTtemp
    real(sp)                  :: sigV, sigT, LV, LT
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
    if ( proc .eq. 0 ) then
        open(unit = 10, file="input.inp")
        read(10,*) trash
        read(10,*) trash, case
        read(10,*) trash
        read(10,*) trash, xsize 
        read(10,*) trash, ysize
        read(10,*) trash, zsize
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
    call MPI_Bcast(sigtemp, 1, MPI_REAL, 0, MPI_COMM_WORLD, mpi_err)
    call MPI_Bcast(Ltemp, 1, MPI_REAL, 0, MPI_COMM_WORLD, mpi_err)
    call MPI_Bcast(sigTtemp, 1, MPI_REAL, 0, MPI_COMM_WORLD, mpi_err)
    call MPI_Bcast(LTtemp, 1, MPI_REAL, 0, MPI_COMM_WORLD, mpi_err)
    call MPI_Barrier(MPI_COMM_WORLD, mpi_err)
    !**************************************************************************

    ! Set parameters **********************************************************
    sigT = sigTtemp
    LT   = LTtemp*(1.0/0.746834)
    sigV = sigtemp
    LV   = Ltemp*(1.0/0.746834)
    !**************************************************************************

    ! Construct Domain ********************************************************
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
    allocate ( temp1(2,2,2))
    allocate ( temp2(2,2,2))
    allocate ( temp3(2,2,2))
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
    mkdirCmd = 'mkdir -p '//trim(case)
    call system( mkdirCmd )
    filename = trim(case)//"/"//"field"//trim(groupname)//".h5"
    call create_h5_f(filename, MPI_COMM_WORLD )
    call create_h5_d(filename, "ux", [size(ux,dim=1),&
        size(uy,dim=2), size(uz,dim=3)*nproc], MPI_COMM_WORLD, ux(:,1,1),&
        chunked=.true., dim_chunk=shape(ux), &
        offset=[0, 0, proc*size(ux,dim=3)])
    call create_h5_d(filename, "uy", [size(uy,dim=1),&
        size(uy,dim=2), size(uy,dim=3)*nproc], MPI_COMM_WORLD, uy(:,1,1),&
        chunked=.true., dim_chunk=shape(uy), &
        offset=[0, 0, proc*size(uy,dim=3)])
    call create_h5_d(filename, "uz", [size(uz,dim=1),&
        size(uz,dim=2), size(uz,dim=3)*nproc], MPI_COMM_WORLD, uz(:,1,1),&
        chunked=.true., dim_chunk=shape(uz), &
        offset=[0, 0, proc*size(uz,dim=3)])
    call create_h5_d(filename, "T", [size(T1,dim=1),&
        size(T1,dim=2), size(T1,dim=3)*nproc], MPI_COMM_WORLD, Temp(:,1,1),&
        chunked=.true., dim_chunk=shape(Temp), &
        offset=[0, 0,proc*size(T1,dim=3)])
    call create_h5_d(filename, "sigV", [1], MPI_COMM_WORLD,[sigV])
    call create_h5_d(filename, "sigT", [1], MPI_COMM_WORLD,[sigT])
    call create_h5_d(filename, "LV", [1], MPI_COMM_WORLD, [LV] )
    call create_h5_d(filename, "LT", [1], MPI_COMM_WORLD, [LT] )
    call create_h5_d(filename, "x", [size(x,dim=1)], MPI_COMM_WORLD, x)
    call create_h5_d(filename, "y", [size(x,dim=1)], MPI_COMM_WORLD, y)
    call create_h5_d(filename, "z", [size(x,dim=1)], MPI_COMM_WORLD, z)
    !**************************************************************************


    deallocate ( k1, k2, k3, f1, f2, f3, ux, uy, uz, Temp, u1, u2, u3, T1, &
        w1, w2, w3, wT, x, y, z )
    call MPI_FINALIZE(mpi_err)

end program main