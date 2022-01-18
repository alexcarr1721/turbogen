!******************************************************************************
!       Generation of turbulent velocity, temperature, and humidity           *
!       fields using the method of random phase outlined by Ostashev          *
!       and Wilson (2015).                                                    *
!                                                                             *
!                                                                             *
!   Author: Alexander N. Carr   ( alexcarr1721@gmail.com )                    *
!           Ph.D. Candidate, Theoretical Fluid Dynamics and Turbulence Group  *
!           University of Florida                                             *
!           Department of Mechanical and Aerospace Engineering                *
!******************************************************************************

program main
    use iso_fortran_env
    use mpi 
    use sbabl_hdf5
    use turbulence
    implicit none 
    ! Program variables *******************************************************
    integer(int32)      :: i, j, k, global_ind, proc, nproc, mpi_err
    integer(int32)      :: dimensions
    character(len=180)  :: trash, fieldname, filename
    character(len=180)  :: directoryname, mkdirCmd
    !**************************************************************************
    ! Grid ********************************************************************
    real(real32), allocatable   :: x1(:), x2(:), x3(:)
    real(real32)                :: x1min, x2min, x3min, x1max, x2max, x3max
    integer(int32)              :: x1size, x2size, x3size
    real(real32), allocatable   :: k1(:), k2(:), k3(:)
    !**************************************************************************
    ! Velocities, temperature, and humidity ***********************************
    real(real32), allocatable   :: u1_2d(:,:), u2_2d(:,:)
    real(real32), allocatable   :: temperature_2d(:,:), q_2d(:,:)
    real(real32), allocatable   :: u1_3d(:,:,:), u2_3d(:,:,:), u3_3d(:,:,:)
    real(real32), allocatable   :: temperature_3d(:,:,:), q_3d(:,:,:)
    !**************************************************************************
    ! Isotropic turbulence parameters *****************************************
    real(real32)                :: sigV, sigT, sigQ, LV, LT, LQ
    !**************************************************************************

    ! Get processor information ***********************************************
    call MPI_INIT(mpi_err)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, nproc, mpi_err)
    call MPI_COMM_RANK(MPI_COMM_WORLD, proc, mpi_err)
    !**************************************************************************

    ! Get field number information from command line **************************
    if (command_argument_count() .ne. 1) then
        fieldname = "1"
    else
        call get_command_argument(1,fieldname)
    end if
    call MPI_Barrier(MPI_COMM_WORLD, mpi_err)
    filename = "field"
    !**************************************************************************

    ! Read input file *********************************************************
    if ( proc .eq. 0 ) then
        open(unit = 10, file="input.inp")
        read(10,*) trash
        read(10,*) trash, directoryname
        read(10,*) trash
        read(10,*) trash, x1size 
        read(10,*) trash, x2size
        read(10,*) trash, x3size
        read(10,*) trash
        read(10,*) trash, x1min
        read(10,*) trash, x1max
        read(10,*) trash, x2min
        read(10,*) trash, x2max
        read(10,*) trash, x3min
        read(10,*) trash, x3max
        read(10,*) trash
        read(10,*) trash, sigV
        read(10,*) trash, LV 
        read(10,*) trash, sigT 
        read(10,*) trash, LT
        read(10,*) trash, sigQ 
        read(10,*) trash, LQ
        read(10,*) trash 
        read(10,*) trash, dimensions
        close(10)
    end if
    call MPI_Bcast(directoryname, 180, MPI_CHARACTER, 0, MPI_COMM_WORLD, mpi_err)
    directoryname = trim(directoryname)
    call MPI_Bcast(dimensions, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpi_err)
    call MPI_Bcast(x1min, 1, MPI_REAL, 0, MPI_COMM_WORLD, mpi_err)
    call MPI_Bcast(x2min, 1, MPI_REAL, 0, MPI_COMM_WORLD, mpi_err)
    call MPI_Bcast(x3min, 1, MPI_REAL, 0, MPI_COMM_WORLD, mpi_err)
    call MPI_Bcast(x1max, 1, MPI_REAL, 0, MPI_COMM_WORLD, mpi_err)
    call MPI_Bcast(x2max, 1, MPI_REAL, 0, MPI_COMM_WORLD, mpi_err)
    call MPI_Bcast(x3max, 1, MPI_REAL, 0, MPI_COMM_WORLD, mpi_err)
    call MPI_Bcast(x1size, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpi_err)
    call MPI_Bcast(x2size, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpi_err)
    call MPI_Bcast(x3size, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpi_err)
    call MPI_Bcast(sigV, 1, MPI_REAL, 0, MPI_COMM_WORLD, mpi_err)
    call MPI_Bcast(LV, 1, MPI_REAL, 0, MPI_COMM_WORLD, mpi_err)
    call MPI_Bcast(sigT, 1, MPI_REAL, 0, MPI_COMM_WORLD, mpi_err)
    call MPI_Bcast(LT, 1, MPI_REAL, 0, MPI_COMM_WORLD, mpi_err)
    call MPI_Bcast(sigQ, 1, MPI_REAL, 0, MPI_COMM_WORLD, mpi_err)
    call MPI_Bcast(LQ, 1, MPI_REAL, 0, MPI_COMM_WORLD, mpi_err)
    call MPI_Barrier(MPI_COMM_WORLD, mpi_err)
    !**************************************************************************

    if ( dimensions .eq. 2 ) then
        goto 100
    else if ( dimensions .eq. 3 ) then
        goto 200
    else
        goto 100
    end if

    100 continue

    ! Construct domain ********************************************************
    allocate ( x1(x1size), x2(x2size) )
    x1 = (/ ( ( x1min + (i-1)*( (x1max - x1min)/(size(x1,dim=1) - 1) ) ),i=1,size(x1,dim=1) ) /)
    x2 = (/ ( ( x2min + (i-1)*( (x2max - x2min)/(size(x2,dim=1) - 1) ) ),i=1,size(x2,dim=1) ) /)
    !**************************************************************************

    ! Construct wavenumbers ***************************************************
    allocate ( k1(x1size), k2(x2size) )
    k1 = generate_wavenumbers(x1)
    k2 = generate_wavenumbers(x2)
    !**************************************************************************

    ! Generate isotropic ******************************************************
    allocate ( u1_2d(x1size, x2size/nproc), u2_2d(x1size, x2size/nproc) )
    allocate ( temperature_2d(x1size, x2size/nproc), q_2d(x1size, x2size/nproc) )
    call generate_isotropic(u1_2d, u2_2d, temperature_2d, q_2d, k1, k2, sigV, LV, sigT, LT, sigQ, LQ)
    !**************************************************************************

    ! Write to file ***********************************************************
    mkdirCmd = 'mkdir -p '//trim(directoryname)
    call system( mkdirCmd )
    filename = trim(directoryname)//"/"//trim(filename)//trim(fieldname)//".h5"
    ! Create file
    call create_h5_f(filename, MPI_COMM_WORLD)
    call create_h5_d(filename, "u1", [x1size, x2size], MPI_COMM_WORLD, u1_2d(:,1), &
        chunked=.true., dim_chunk=[x1size, x2size/nproc], offset=[0, proc*(x2size/nproc)] )
    call create_h5_d(filename, "u2", [x1size, x2size], MPI_COMM_WORLD, u2_2d(:,1), &
        chunked=.true., dim_chunk=[x1size, x2size/nproc], offset=[0, proc*(x2size/nproc)] )
    call create_h5_d(filename, "temperature", [x1size, x2size], MPI_COMM_WORLD, temperature_2d(:,1), &
        chunked=.true., dim_chunk=[x1size, x2size/nproc], offset=[0, proc*(x2size/nproc)] )
    call create_h5_d(filename, "q", [x1size, x2size], MPI_COMM_WORLD, q_2d(:,1), &
        chunked=.true., dim_chunk=[x1size, x2size/nproc], offset=[0, proc*(x2size/nproc)] )
    call create_h5_d(filename, "x1", [x1size], MPI_COMM_WORLD, x1)
    call create_h5_d(filename, "x2", [x2size], MPI_COMM_WORLD, x2)
    !**************************************************************************

    deallocate ( x1, x2, k1, k2, u1_2d, u2_2d, temperature_2d, q_2d )

    200 continue

    call MPI_FINALIZE(mpi_err)

end program main
