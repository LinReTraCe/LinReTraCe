program hdf5test
  use hdf5
  implicit none

  write(*,*) 'begin HDF5 tests'
  call hdf5_interface()
  call hdf5_file()
  call hdf5_logical() ! implicitly tests for version of HDF5 library
  call hdf5_complex()
  write(*,*) 'end HDF5 Tests'

  contains

  subroutine hdf5_interface()
    implicit none
    integer :: hdf_err
    write(*,'(A)', advance='no') '    opening HDF5 interface:'
    call h5open_f(hdf_err)
    if (hdf_err) then
      write(*,*) 'failed'
    else
      write(*,*) 'success'
    endif
  end subroutine

  subroutine hdf5_file()
    implicit none
    integer(hid_t) ifile
    integer :: hdf_err

    write(*,'(A)', advance='no') '    creating HDF5 file:'
    call h5fcreate_f('test.hdf5', h5f_acc_trunc_f, ifile, hdf_err)
    if (hdf_err) then
      write(*,*) 'failed'
    else
      write(*,*) 'success'
    endif

    write(*,'(A)', advance='no') '    opening HDF5 file:'
    call h5fopen_f('test.hdf5', h5f_acc_rdwr_f, ifile, hdf_err)
    if (hdf_err) then
      write(*,*) 'failed'
    else
      write(*,*) 'success'
    endif

    write(*,'(A)', advance='no') '    closing HDF5 file:'
    call h5fclose_f(ifile, hdf_err)
    if (hdf_err) then
      write(*,*) 'failed'
    else
      write(*,*) 'success'
    endif

  end subroutine

  subroutine hdf5_complex()
    implicit none
    integer                        :: hdf_err
    integer(size_t), parameter     :: zero = 0
    integer(hid_t)                 :: complex_id_r_dp, complex_id_i_dp, complex_id_dp
    integer(hid_t)                 :: plist_id, dspace_id, dset_id, grp_id
    integer(hid_t)                 :: ifile
    integer(size_t)                :: compound_size, type_sized
    complex(8)                     :: darray = (0,2.d0)
    complex(8), parameter          :: ci = (0.d0, 1.d0)
    real(8)                        :: tmp_r, tmp_i
    integer(8), dimension(0)       :: dims
    integer(hsize_t), dimension(0) :: dset_dims

    write(*,'(A)', advance='yes') '    complex datatype tests.'
    write(*,'(A)', advance='no')  '      creating complex structure:'
    ! double precision
    call h5pcreate_f(h5p_dataset_xfer_f, plist_id, hdf_err)
    call h5pset_preserve_f(plist_id, .true., hdf_err)
    ! compound
    call h5tget_size_f(h5t_native_double, type_sized, hdf_err)
    compound_size = 2*type_sized
    call h5tcreate_f(h5t_compound_f, compound_size, complex_id_dp, hdf_err)
    call h5tinsert_f(complex_id_dp, "r", zero, h5t_native_double, hdf_err)
    call h5tinsert_f(complex_id_dp, "i", type_sized, h5t_native_double, hdf_err)
    ! separate
    call h5tcreate_f(h5t_compound_f, type_sized, complex_id_r_dp, hdf_err)
    call h5tinsert_f(complex_id_r_dp, "r", zero, h5t_native_double, hdf_err)
    call h5tcreate_f(h5t_compound_f, type_sized, complex_id_i_dp, hdf_err)
    call h5tinsert_f(complex_id_i_dp, "i", zero, h5t_native_double, hdf_err)
    call h5pclose_f(plist_id, hdf_err)
    if (hdf_err) then
      write(*,*) 'failed'
    else
      write(*,*) 'success'
    endif

    write(*,'(A)', advance='no') '      opening file:'
    call h5fopen_f('test.hdf5', h5f_acc_rdwr_f, ifile, hdf_err)
    if (hdf_err) then
      write(*,*) 'failed'
    else
      write(*,*) 'success'
    endif

    write(*,'(A)', advance='no') '      writing complex to file:'
    ! writing logical 0D array
    call h5screate_simple_f(0, dims, dspace_id, hdf_err)
    grp_id = ifile ! write into root
    call h5dcreate_f(grp_id, 'complex', complex_id_dp, dspace_id, dset_id, hdf_err)
    call h5dwrite_f(dset_id, complex_id_r_dp, real(darray), dims, hdf_err)
    call h5dwrite_f(dset_id, complex_id_i_dp, aimag(darray), dims, hdf_err)
    call h5dclose_f(dset_id, hdf_err)
    call h5sclose_f(dspace_id, hdf_err)
    if (hdf_err) then
      write(*,*) 'failed'
    else
      write(*,*) 'success'
    endif


    write(*,'(A)', advance='no') '      reading logical from file:'
    ! writing logical 0D array
    darray = 123.d0
    call h5dopen_f(ifile, 'complex', dset_id, hdf_err)
    call h5dread_f(dset_id, complex_id_r_dp, tmp_r, dset_dims, hdf_err)
    call h5dread_f(dset_id, complex_id_i_dp, tmp_i, dset_dims, hdf_err)
    call h5dclose_f(dset_id, hdf_err)

    darray = tmp_r + ci*tmp_i
    if (real(darray) == 0.d0 .and. aimag(darray) == 2.d0) then
      write(*,*) 'success'
    else
      write(*,*) 'failed'
    endif

    write(*,'(A)', advance='no') '      closing file:'
    call h5fclose_f(ifile, hdf_err)
    if (hdf_err) then
      write(*,*) 'failed'
    else
      write(*,*) 'success'
    endif

  end subroutine

  subroutine hdf5_logical()
    implicit none
    integer, parameter             :: logical_size = 1
    integer(hid_t)                 :: grp_id, logical_id, dset_id, ifile, dspace_id
    integer, parameter             :: zero = 0
    integer, parameter             :: one = 1
    integer(size_t)                :: type_sized
    integer                        :: hdf_err
    integer(8), dimension(0)       :: dims
    integer(hsize_t), dimension(0) :: dset_dims

    integer(logical_size)    :: darray = 1

    write(*,'(A)', advance='yes') '    logical datatype tests.'
    write(*,'(A)', advance='no')  '      creating logical enumeration:'
    ! create logical datatype : enumeration 0: FALSE; 1: TRUE
    call h5tget_size_f(h5t_native_b8, type_sized, hdf_err) ! h5py uses numpy 1-byte bool
    call h5tcreate_f(h5t_enum_f, type_sized, logical_id, hdf_err)
    call h5tenum_insert_f(logical_id, "FALSE", zero, hdf_err)
    call h5tenum_insert_f(logical_id, "TRUE", one, hdf_err)
    if (hdf_err) then
      write(*,*) 'failed'
    else
      write(*,*) 'success'
    endif

    write(*,'(A)', advance='no') '      opening file:'
    call h5fopen_f('test.hdf5', h5f_acc_rdwr_f, ifile, hdf_err)
    if (hdf_err) then
      write(*,*) 'failed'
    else
      write(*,*) 'success'
    endif

    write(*,'(A)', advance='no') '      writing logical to file:'
    ! writing logical 0D array
    call h5screate_simple_f(0, dims, dspace_id, hdf_err)
    grp_id = ifile ! write into root
    call h5dcreate_f(grp_id, 'logical', logical_id, dspace_id, dset_id, hdf_err)
    call h5dwrite_f(dset_id, logical_id, darray, dims, hdf_err)
    call h5dclose_f(dset_id, hdf_err)
    call h5sclose_f(dspace_id, hdf_err)
    if (hdf_err) then
      write(*,*) 'failed'
    else
      write(*,*) 'success'
    endif


    write(*,'(A)', advance='no') '      reading logical from file:'
    ! writing logical 0D array
    darray = 123
    call h5dopen_f(ifile, 'logical', dset_id, hdf_err)
    call h5dread_f(dset_id, logical_id, darray, dset_dims, hdf_err)
    call h5dclose_f(dset_id, hdf_err)
    if (darray == 1) then
      write(*,*) 'success'
    else
      write(*,*) 'failed'
    endif

    write(*,'(A)', advance='no') '      closing file:'
    call h5fclose_f(ifile, hdf_err)
    if (hdf_err) then
      write(*,*) 'failed'
    else
      write(*,*) 'success'
    endif

  end subroutine

end program
