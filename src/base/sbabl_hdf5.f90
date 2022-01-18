module sbabl_hdf5
    use iso_fortran_env
    use mpi
    use hdf5
    implicit none


    contains

    subroutine read_h5_rank(filename, dsetname, rank)
      character(len=*), intent(in)          :: filename  ! file name
      character(len=*), intent(in)          :: dsetname  ! dataset name
      integer(int32), intent(inout)         :: rank      ! rank of array
      ! HDF5 type
      integer(HID_T) :: file_id
      integer(HID_T) :: dset_id
      integer(HID_T) :: dspace_id
      ! error
      integer        :: error

      call h5open_f(error)
      call h5fopen_f(trim(filename), H5F_ACC_RDONLY_F, file_id, error)
      call h5dopen_f(file_id, trim(dsetname), dset_id, error)
      call h5dget_space_f(dset_id, dspace_id, error)
      call h5sget_simple_extent_ndims_f(dspace_id, rank, error)
      call h5sclose_f(dspace_id, error)
      call h5dclose_f(dset_id, error)
      call h5fclose_f(file_id, error)
      call h5close_f(error)

    end subroutine read_h5_rank

    subroutine read_h5_dims(filename, dsetname, dims, maxdims)
      character(len=*), intent(in)          :: filename  ! file name
      character(len=*), intent(in)          :: dsetname  ! dataset name
      integer(int32), intent(inout)         :: dims(:)   ! dims of array
      integer(int32), intent(inout)         :: maxdims(:)! maxdims of array
      ! HDF5 type
      integer(HID_T)   :: file_id
      integer(HID_T)   :: dset_id
      integer(HID_T)   :: dspace_id
      integer(HSIZE_T) :: dimread(size(dims)), maxdimread(size(maxdims))
      ! error
      integer(int32)   :: error

      call h5open_f(error)
      call h5fopen_f(trim(filename), H5F_ACC_RDONLY_F, file_id, error)
      call h5dopen_f(file_id, trim(dsetname), dset_id, error)
      call h5dget_space_f(dset_id, dspace_id, error)
      call h5sget_simple_extent_dims_f(dspace_id, dimread, maxdimread, error)
      dims(:) = Int(dimread(:),int32)
      maxdims(:) = Int(maxdimread(:),int32)
      call h5sclose_f(dspace_id, error)
      call h5dclose_f(dset_id, error)
      call h5fclose_f(file_id, error)
      call h5close_f(error)

    end subroutine read_h5_dims

    subroutine read_h5(f, filename, dsetname, dim, gmin, comm, loff)
        real(real32), dimension(:), intent(inout) :: f         ! Data buffer
        character(len=*), intent(in)         :: filename  ! file name
        character(len=*), intent(in)         :: dsetname  ! dataset name
        integer(int32), dimension(:), intent(in) :: dim       ! dimensions of dataset to be read
        integer(int32), dimension(:), intent(in) :: gmin      ! global array dimensions (offset)
        integer(int32), dimension(:), intent(in), optional :: loff ! Offset in local memory
        integer                               :: mpi_err, error, i, comm
        integer(HID_T)                        :: file_id   ! File Identifier
        integer(HID_T)                        :: dset_id   ! Dataset Identifier
        integer(HID_T)                        :: dataspace ! Dataspace ID
        integer(HID_T)                        :: memspace
        integer(HSIZE_T), allocatable         :: dims(:), stride(:), offset(:)
        integer(HSIZE_T), allocatable         :: offset_out(:), loc_off(:)

        allocate ( dims(size(dim)), stride(size(dim)), offset(size(dim)), &
            offset_out(size(dim)) , loc_off(size(dim)) )

        if ( present(loff) ) then
            loc_off = loff
        else
            loc_off = 0
        end if

        ! Set dimensions
        do i = 1, size(dim)
            dims(i) = dim(i)           ! Dimensions of slab to read
            stride(i) = 1              ! Continuous data
            offset(i) = gmin(i) - 1    ! Zero indexed
            offset_out(i) = loc_off(i) ! No offset in memspace memory
        end do

        call h5open_f(error)

        call h5fopen_f(trim(filename), H5F_ACC_RDONLY_F, file_id, &
            error)

        ! Read data (Works for chunked datasets as well) =======================
        call h5dopen_f(file_id, trim(dsetname), dset_id, error)       ! Open dataset
        call h5dget_space_f(dset_id, dataspace, error)          ! Get identifier
        call h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, &
            offset, dims, error)                              ! Select hyperslab
        call h5screate_simple_f(size(dim), dims, memspace, error)! memspace
        call h5sselect_hyperslab_f(memspace, H5S_SELECT_SET_F, &
            offset_out, dims, error)              ! Select hyperslab in memory
        call h5dread_f(dset_id, H5T_NATIVE_REAL, f, dims, error, &
            memspace, dataspace) ! Read hyperslab on file to hyperslab in memory
        ! ======================================================================

        ! Close objects
        call h5dclose_f(dset_id, error)
        call h5sclose_f(dataspace, error)
        call h5sclose_f(memspace, error)
        call h5fclose_f(file_id, error)

        call h5close_f(error)

        deallocate ( dims, stride, offset, offset_out )

        call MPI_BARRIER(comm, mpi_err)

    end subroutine read_h5

    subroutine read_h5_int(f, filename, dsetname, dim, gmin, comm, loff)
        integer(int32), dimension(:), intent(inout) :: f         ! Data buffer
        character(len=*), intent(in)         :: filename  ! file name
        character(len=*), intent(in)         :: dsetname  ! dataset name
        integer(int32), dimension(:), intent(in) :: dim       ! dimensions of dataset to be read
        integer(int32), dimension(:), intent(in) :: gmin      ! global array dimensions (offset)
        integer(int32), dimension(:), intent(in), optional :: loff ! Offset in local memory
        integer                               :: mpi_err, error, i, comm
        integer(HID_T)                        :: file_id   ! File Identifier
        integer(HID_T)                        :: dset_id   ! Dataset Identifier
        integer(HID_T)                        :: dataspace ! Dataspace ID
        integer(HID_T)                        :: memspace
        integer(HSIZE_T), allocatable         :: dims(:), stride(:), offset(:)
        integer(HSIZE_T), allocatable         :: offset_out(:), loc_off(:)

        allocate ( dims(size(dim)), stride(size(dim)), offset(size(dim)), &
            offset_out(size(dim)) , loc_off(size(dim)) )

        if ( present(loff) ) then
            loc_off = loff
        else
            loc_off = 0
        end if

        ! Set dimensions
        do i = 1, size(dim)
            dims(i) = dim(i)           ! Dimensions of slab to read
            stride(i) = 1              ! Continuous data
            offset(i) = gmin(i) - 1    ! Zero indexed
            offset_out(i) = loc_off(i) ! No offset in memspace memory
        end do

        call h5open_f(error)

        call h5fopen_f(trim(filename), H5F_ACC_RDONLY_F, file_id, &
            error)

        ! Read data (Works for chunked datasets as well) =======================
        call h5dopen_f(file_id, trim(dsetname), dset_id, error)       ! Open dataset
        call h5dget_space_f(dset_id, dataspace, error)          ! Get identifier
        call h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, &
            offset, dims, error)                              ! Select hyperslab
        call h5screate_simple_f(size(dim), dims, memspace, error)! memspace
        call h5sselect_hyperslab_f(memspace, H5S_SELECT_SET_F, &
            offset_out, dims, error)              ! Select hyperslab in memory
        call h5dread_f(dset_id, H5T_NATIVE_INTEGER, f, dims, error, &
            memspace, dataspace) ! Read hyperslab on file to hyperslab in memory
        ! ======================================================================

        ! Close objects
        call h5dclose_f(dset_id, error)
        call h5sclose_f(dataspace, error)
        call h5sclose_f(memspace, error)
        call h5fclose_f(file_id, error)

        call h5close_f(error)

        deallocate ( dims, stride, offset, offset_out )

        call MPI_BARRIER(comm, mpi_err)

    end subroutine read_h5_int

    subroutine create_h5_f(filename, comm, groupnames)
        character(len=*), intent(in)         :: filename  ! file name
        character(len=*), dimension(:), optional, intent(in) :: groupnames ! dataset name
        integer(4), intent(in)                :: comm
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
        integer(int32), intent(in)             :: dim(:)      ! total data dimensions
        integer(int32), intent(in)             :: comm        ! MPI communicator
        real(real32),  optional, intent(in)         :: f(:)        ! data buffer
        logical, optional, intent(in)          :: chunked     ! .true. = dataset is chunked
        integer(int32), optional, intent(in)   :: dim_chunk(:)   ! Local chunk size
        integer(int32), optional, intent(in)   :: offset(:)      ! Global position of chunk
        logical, optional, intent(in)          :: extendable  ! .true. = dataset is chunked
        integer(int32), intent(in), optional   :: loff(:)
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

    subroutine create_h5_d_int(filename, dsetname, dim, comm, f, chunked, &
        dim_chunk, offset, extendable, loff)
        character(len=*), intent(in)           :: filename  ! file name
        character(len=*), intent(in)           :: dsetname    ! dataset name
        integer(int32), intent(in)             :: dim(:)      ! total data dimensions
        integer(int32), intent(in)             :: comm        ! MPI communicator
        integer(int32),  optional, intent(in)         :: f(:)        ! data buffer
        logical, optional, intent(in)          :: chunked     ! .true. = dataset is chunked
        integer(int32), optional, intent(in)   :: dim_chunk(:)   ! Local chunk size
        integer(int32), optional, intent(in)   :: offset(:)      ! Global position of chunk
        logical, optional, intent(in)          :: extendable  ! .true. = dataset is chunked
        integer(int32), intent(in), optional   :: loff(:)
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
            call h5dcreate_f(file_id, dsetname, H5T_NATIVE_INTEGER, dataspace, &
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
                call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, f, dimt, error, &
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
            call h5dcreate_f(file_id, dsetname, H5T_NATIVE_INTEGER, dataspace, &
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
                call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, f, dimt, error, &
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
            call h5dcreate_f(file_id, dsetname, H5T_NATIVE_INTEGER, dataspace, &
                dset_id, error, plist_id)

            if ( present(f) ) then
                call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
                call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)

                ! Write the datasets collectively
                call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, f, dims, error, &
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

    end subroutine create_h5_d_int

    subroutine write_h5_d(filename, dsetname, dim, comm, f, chunked, &
        dim_chunk, offset, extendable, extension)
        character(len=*), intent(in)           :: filename  ! file name
        character(len=*), intent(in)           :: dsetname    ! dataset name
        integer(int32), intent(in)             :: dim(:)      ! total data dimensions
        integer(int32), intent(in)             :: comm        ! MPI communicator
        real(real32), intent(in)               :: f(:)        ! data buffer
        logical, optional, intent(in)          :: chunked     ! .true. = dataset is chunked
        integer(int32), optional, intent(in)   :: dim_chunk(:)   ! Local chunk size
        integer(int32), optional, intent(in)   :: offset(:)      ! Global position of chunk
        logical, optional, intent(in)          :: extendable  ! .true. = dataset is chunked
        integer(int32), optional, intent(in)   :: extension(:)   ! dimension to extend the dataset
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


end module sbabl_hdf5
