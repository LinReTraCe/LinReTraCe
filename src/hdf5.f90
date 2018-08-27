module Mhdf5
  use hdf5
  use Mparams
  implicit none

  ! for an easy access to common hdf_error and special complex data types
  integer         :: hdf_err
  integer(hid_t)  :: compound_id ! complex compound dtype
  integer(hid_t)  :: type_r_id, type_i_id

  ! may want to extend this to higher dimensions
  interface h5_loaddata
    module procedure h5_loaddata_1d_complex, h5_loaddata_2d_complex, h5_loaddata_3d_complex, &
                     h5_loaddata_1d_real,    h5_loaddata_2d_real,    h5_loaddata_3d_real
  end interface h5_loaddata

  interface h5_flushdata
    module procedure h5_flushdata_1d_complex, h5_flushdata_2d_complex, h5_flushdata_3d_complex, &
                     h5_flushdata_1d_real,    h5_flushdata_2d_real,    h5_flushdata_3d_real
  end interface h5_flushdata

  contains

  subroutine h5_init()
    call h5open_f(hdf_err)
    call h5_create_complex_datatype()
  end subroutine h5_init

  subroutine h5_finalize()
    call h5close_f(hdf_err)
  end subroutine h5_finalize

  subroutine h5_create_file(fname)
    character(len=*), intent(in) :: fname

    integer(hid_t) :: file_id
    ! truncate -> if it already exists, erase all data
    ! other possibility: h5f_acc_excl_f -> fail if already exists
    call h5fcreate_f(trim(adjustl(fname)), h5f_acc_trunc_f, file_id, hdf_err)
  end subroutine h5_create_file

  ! this code part is from ADGA
  ! github.com/abinitiodga/adga
  ! GPLv3
  subroutine h5_create_complex_datatype()
    integer(size_t), parameter :: zero = 0
    integer(hid_t)  :: plist_id
    integer(size_t) :: compound_size, type_sized

    call h5pcreate_f(h5p_dataset_xfer_f, plist_id, hdf_err)
    call h5pset_preserve_f(plist_id, .true., hdf_err)

    ! create compound datatype for complex arrays:
    call h5tget_size_f(h5t_native_double, type_sized, hdf_err)
    compound_size = 2*type_sized
    call h5tcreate_f(h5t_compound_f, compound_size, compound_id, hdf_err)
    call h5tinsert_f(compound_id, "r", zero, h5t_native_double, hdf_err)
    call h5tinsert_f(compound_id, "i", type_sized, h5t_native_double, hdf_err)

    !complex type to write real and imaginary individually:
    call h5tcreate_f(h5t_compound_f, type_sized, type_r_id, hdf_err)
    call h5tinsert_f(type_r_id, "r", zero, h5t_native_double, hdf_err)
    call h5tcreate_f(h5t_compound_f, type_sized, type_i_id, hdf_err)
    call h5tinsert_f(type_i_id, "i", zero, h5t_native_double, hdf_err)
  end subroutine h5_create_complex_datatype

  subroutine h5_loaddata_1d_real(fname, dset, darray)
    character(len=*), intent(in)       :: fname
    character(len=*), intent(in)       :: dset
    real(8), allocatable, dimension(:) :: darray

    integer(hid_t) :: file_id, dset_id, dset_space_id
    integer(hsize_t), dimension(1) :: dset_dims, dset_maxdims

    call h5fopen_f(trim(adjustl(fname)), h5f_acc_rdonly_f, file_id, hdf_err)
    call h5dopen_f(file_id, trim(adjustl(dset)), dset_id, hdf_err)
    call h5dget_space_f(dset_id, dset_space_id, hdf_err)
    call h5sget_simple_extent_dims_f(dset_space_id, dset_dims, dset_maxdims, hdf_err)
    allocate(darray(dset_dims(1)))
    call h5dread_f(dset_id, h5t_native_double, darray, dset_dims, hdf_err)
    call h5dclose_f(dset_id, hdf_err)
    call h5fclose_f(file_id, hdf_err)
  end subroutine h5_loaddata_1d_real

  subroutine h5_loaddata_1d_complex(fname, dset, darray)
    character(len=*), intent(in)          :: fname
    character(len=*), intent(in)          :: dset
    complex(8), allocatable, dimension(:) :: darray

    real(8), allocatable, dimension(:)    :: tmp_i, tmp_r
    integer(hid_t) :: file_id, dset_id, dset_space_id
    integer(hsize_t), dimension(1) :: dset_dims, dset_maxdims

    call h5fopen_f(trim(adjustl(fname)), h5f_acc_rdonly_f, file_id, hdf_err)
    call h5dopen_f(file_id, trim(adjustl(dset)), dset_id, hdf_err)
    call h5dget_space_f(dset_id, dset_space_id, hdf_err)
    call h5sget_simple_extent_dims_f(dset_space_id, dset_dims, dset_maxdims, hdf_err)
    allocate(darray(dset_dims(1)))
    allocate(tmp_i(dset_dims(1)),tmp_r(dset_dims(1)))
    call h5dread_f(dset_id, type_r_id, tmp_r, dset_dims, hdf_err)
    call h5dread_f(dset_id, type_i_id, tmp_i, dset_dims, hdf_err)
    call h5dclose_f(dset_id, hdf_err)
    call h5fclose_f(file_id, hdf_err)
    darray = tmp_r + ci* tmp_i
    deallocate(tmp_i, tmp_r)
  end subroutine h5_loaddata_1d_complex

  subroutine h5_loaddata_2d_real(fname, dset, darray)
    character(len=*), intent(in)         :: fname
    character(len=*), intent(in)         :: dset
    real(8), allocatable, dimension(:,:) :: darray

    integer(hid_t) :: file_id, dset_id, dset_space_id
    integer(hsize_t), dimension(2) :: dset_dims, dset_maxdims

    call h5fopen_f(trim(adjustl(fname)), h5f_acc_rdonly_f, file_id, hdf_err)
    call h5dopen_f(file_id, trim(adjustl(dset)), dset_id, hdf_err)
    call h5dget_space_f(dset_id, dset_space_id, hdf_err)
    call h5sget_simple_extent_dims_f(dset_space_id, dset_dims, dset_maxdims, hdf_err)
    allocate(darray(dset_dims(1), dset_dims(2)))
    call h5dread_f(dset_id, h5t_native_double, darray, dset_dims, hdf_err)
    call h5dclose_f(dset_id, hdf_err)
    call h5fclose_f(file_id, hdf_err)
  end subroutine h5_loaddata_2d_real

  subroutine h5_loaddata_2d_complex(fname, dset, darray)
    character(len=*), intent(in)            :: fname
    character(len=*), intent(in)            :: dset
    complex(8), allocatable, dimension(:,:) :: darray

    real(8), allocatable, dimension(:,:)    :: tmp_i, tmp_r
    integer(hid_t) :: file_id, dset_id, dset_space_id
    integer(hsize_t), dimension(2) :: dset_dims, dset_maxdims

    call h5fopen_f(trim(adjustl(fname)), h5f_acc_rdonly_f, file_id, hdf_err)
    call h5dopen_f(file_id, trim(adjustl(dset)), dset_id, hdf_err)
    call h5dget_space_f(dset_id, dset_space_id, hdf_err)
    call h5sget_simple_extent_dims_f(dset_space_id, dset_dims, dset_maxdims, hdf_err)
    allocate(darray(dset_dims(1), dset_dims(2)))
    allocate(tmp_r(dset_dims(1), dset_dims(2)),tmp_i(dset_dims(1), dset_dims(2)))
    call h5dread_f(dset_id, type_r_id, tmp_r, dset_dims, hdf_err)
    call h5dread_f(dset_id, type_i_id, tmp_i, dset_dims, hdf_err)
    call h5dclose_f(dset_id, hdf_err)
    call h5fclose_f(file_id, hdf_err)
    darray = tmp_r + ci* tmp_i
    deallocate(tmp_i, tmp_r)
  end subroutine h5_loaddata_2d_complex

  subroutine h5_loaddata_3d_real(fname, dset, darray)
    character(len=*), intent(in)           :: fname
    character(len=*), intent(in)           :: dset
    real(8), allocatable, dimension(:,:,:) :: darray

    integer(hid_t) :: file_id, dset_id, dset_space_id
    integer(hsize_t), dimension(3) :: dset_dims, dset_maxdims

    call h5fopen_f(trim(adjustl(fname)), h5f_acc_rdonly_f, file_id, hdf_err)
    call h5dopen_f(file_id, trim(adjustl(dset)), dset_id, hdf_err)
    call h5dget_space_f(dset_id, dset_space_id, hdf_err)
    call h5sget_simple_extent_dims_f(dset_space_id, dset_dims, dset_maxdims, hdf_err)
    allocate(darray(dset_dims(1), dset_dims(2), dset_dims(3)))
    call h5dread_f(dset_id, h5t_native_double, darray, dset_dims, hdf_err)
    call h5dclose_f(dset_id, hdf_err)
    call h5fclose_f(file_id, hdf_err)
  end subroutine h5_loaddata_3d_real

  subroutine h5_loaddata_3d_complex(fname, dset, darray)
    character(len=*), intent(in)              :: fname
    character(len=*), intent(in)              :: dset
    complex(8), allocatable, dimension(:,:,:) :: darray

    real(8), allocatable, dimension(:,:,:)    :: tmp_i, tmp_r
    integer(hid_t)                            :: file_id, dset_id, dset_space_id
    integer(hsize_t), dimension(3)            :: dset_dims, dset_maxdims

    call h5fopen_f(trim(adjustl(fname)), h5f_acc_rdonly_f, file_id, hdf_err)
    call h5dopen_f(file_id, trim(adjustl(dset)), dset_id, hdf_err)
    call h5dget_space_f(dset_id, dset_space_id, hdf_err)
    call h5sget_simple_extent_dims_f(dset_space_id, dset_dims, dset_maxdims, hdf_err)
    allocate(darray(dset_dims(1), dset_dims(2), dset_dims(3)))
    allocate(tmp_r(dset_dims(1), dset_dims(2), dset_dims(3)),tmp_i(dset_dims(1), dset_dims(2), dset_dims(3)))
    call h5dread_f(dset_id, type_r_id, tmp_r, dset_dims, hdf_err)
    call h5dread_f(dset_id, type_i_id, tmp_i, dset_dims, hdf_err)
    call h5dclose_f(dset_id, hdf_err)
    call h5fclose_f(file_id, hdf_err)
    darray = tmp_r + ci* tmp_i
    deallocate(tmp_i, tmp_r)
  end subroutine h5_loaddata_3d_complex

  subroutine h5_flushdata_1d_real(fname, dset, darray, group)
    character(len=*), intent(in)           :: fname
    character(len=*), intent(in)           :: dset
    real(8), dimension(:)                  :: darray
    character(len=*), intent(in), optional :: group

    integer(hid_t)                         :: file_id, grp_id, dset_id, dspace_id
    integer(hsize_t), dimension(1)         :: dset_dims, dset_maxdims
    integer, parameter                     :: rank = 1
    integer(8), dimension(:), allocatable  :: dims

    allocate(dims(rank))
    dims(1) = size(darray)

    call h5fopen_f(trim(adjustl(fname)),h5f_acc_rdwr_f,file_id,hdf_err)
    call h5screate_f(H5S_SIMPLE_F, dspace_id, hdf_err)
    call h5sset_extent_simple_f(dspace_id, rank, dims, dims, hdf_err)
    if (present(group)) then
      call h5gopen_f(file_id, trim(adjustl(group)), grp_id, hdf_err)
      call h5dcreate_f(grp_id, trim(adjustl(dset)), h5t_native_double, dspace_id, dset_id, hdf_err)
    else
      call h5dcreate_f(file_id,trim(adjustl(dset)), h5t_native_double, dspace_id, dset_id, hdf_err)
    endif
    call h5dwrite_f(dset_id, h5t_native_double, darray, dims, hdf_err)
    call h5sclose_f(dspace_id, hdf_err)
    call h5dclose_f(dset_id, hdf_err)
    if (present(group)) then
      call h5gclose_f(grp_id, hdf_err)
    endif
    call h5fclose_f(file_id, hdf_err)

  end subroutine h5_flushdata_1d_real

  subroutine h5_flushdata_1d_complex(fname, dset, darray, group)
    character(len=*), intent(in)           :: fname
    character(len=*), intent(in)           :: dset
    complex(8), dimension(:)               :: darray
    character(len=*), intent(in), optional :: group

    integer(hid_t)                         :: file_id, grp_id, dset_id, dspace_id
    integer(hsize_t), dimension(1)         :: dset_dims, dset_maxdims
    integer, parameter                     :: rank = 1
    integer(8), dimension(:), allocatable  :: dims

    allocate(dims(rank))
    dims(1) = size(darray)

    call h5fopen_f(trim(adjustl(fname)),h5f_acc_rdwr_f,file_id,hdf_err)
    call h5screate_f(H5S_SIMPLE_F, dspace_id, hdf_err)
    call h5sset_extent_simple_f(dspace_id, rank, dims, dims, hdf_err)
    if (present(group)) then
      call h5gopen_f(file_id, trim(adjustl(group)), grp_id, hdf_err)
      call h5dcreate_f(grp_id, trim(adjustl(dset)), h5t_native_double, dspace_id, dset_id, hdf_err)
    else
      call h5dcreate_f(file_id,trim(adjustl(dset)), h5t_native_double, dspace_id, dset_id, hdf_err)
    endif
    call h5dwrite_f(dset_id, type_r_id, real(darray), dims, hdf_err)
    call h5dwrite_f(dset_id, type_i_id, aimag(darray), dims, hdf_err)
    call h5sclose_f(dspace_id, hdf_err)
    call h5dclose_f(dset_id, hdf_err)
    if (present(group)) then
      call h5gclose_f(grp_id, hdf_err)
    endif
    call h5fclose_f(file_id, hdf_err)

  end subroutine h5_flushdata_1d_complex

  subroutine h5_flushdata_2d_real(fname, dset, darray, group)
    character(len=*), intent(in)           :: fname
    character(len=*), intent(in)           :: dset
    real(8), dimension(:,:)                :: darray
    character(len=*), intent(in), optional :: group

    integer(hid_t)                         :: file_id, grp_id, dset_id, dspace_id
    integer(hsize_t), dimension(2)         :: dset_dims, dset_maxdims
    integer, parameter                     :: rank = 2
    integer(8), dimension(:), allocatable  :: dims

    allocate(dims(rank))
    dims(1) = size(darray,1)
    dims(2) = size(darray,2)

    call h5fopen_f(trim(adjustl(fname)),h5f_acc_rdwr_f,file_id,hdf_err)
    call h5screate_f(H5S_SIMPLE_F, dspace_id, hdf_err)
    call h5sset_extent_simple_f(dspace_id, rank, dims, dims, hdf_err)
    if (present(group)) then
      call h5gopen_f(file_id, trim(adjustl(group)), grp_id, hdf_err)
      call h5dcreate_f(grp_id, trim(adjustl(dset)), h5t_native_double, dspace_id, dset_id, hdf_err)
    else
      call h5dcreate_f(file_id,trim(adjustl(dset)), h5t_native_double, dspace_id, dset_id, hdf_err)
    endif
    call h5dwrite_f(dset_id, h5t_native_double, darray, dims, hdf_err)
    call h5sclose_f(dspace_id, hdf_err)
    call h5dclose_f(dset_id, hdf_err)
    if (present(group)) then
      call h5gclose_f(grp_id, hdf_err)
    endif
    call h5fclose_f(file_id, hdf_err)

  end subroutine h5_flushdata_2d_real

  subroutine h5_flushdata_2d_complex(fname, dset, darray, group)
    character(len=*), intent(in)           :: fname
    character(len=*), intent(in)           :: dset
    complex(8), dimension(:,:)             :: darray
    character(len=*), intent(in), optional :: group

    integer(hid_t)                         :: file_id, grp_id, dset_id, dspace_id
    integer(hsize_t), dimension(2)         :: dset_dims, dset_maxdims
    integer, parameter                     :: rank = 2
    integer(8), dimension(:), allocatable  :: dims

    allocate(dims(rank))
    dims(1) = size(darray,1)
    dims(2) = size(darray,2)

    call h5fopen_f(trim(adjustl(fname)),h5f_acc_rdwr_f,file_id,hdf_err)
    call h5screate_f(H5S_SIMPLE_F, dspace_id, hdf_err)
    call h5sset_extent_simple_f(dspace_id, rank, dims, dims, hdf_err)
    if (present(group)) then
      call h5gopen_f(file_id, trim(adjustl(group)), grp_id, hdf_err)
      call h5dcreate_f(grp_id, trim(adjustl(dset)), h5t_native_double, dspace_id, dset_id, hdf_err)
    else
      call h5dcreate_f(file_id,trim(adjustl(dset)), h5t_native_double, dspace_id, dset_id, hdf_err)
    endif
    call h5dwrite_f(dset_id, type_r_id, real(darray), dims, hdf_err)
    call h5dwrite_f(dset_id, type_i_id, aimag(darray), dims, hdf_err)
    call h5sclose_f(dspace_id, hdf_err)
    call h5dclose_f(dset_id, hdf_err)
    if (present(group)) then
      call h5gclose_f(grp_id, hdf_err)
    endif
    call h5fclose_f(file_id, hdf_err)

  end subroutine h5_flushdata_2d_complex

  subroutine h5_flushdata_3d_real(fname, dset, darray, group)
    character(len=*), intent(in)           :: fname
    character(len=*), intent(in)           :: dset
    real(8), dimension(:,:,:)              :: darray
    character(len=*), intent(in), optional :: group

    integer(hid_t)                         :: file_id, grp_id, dset_id, dspace_id
    integer(hsize_t), dimension(3)         :: dset_dims, dset_maxdims
    integer                                :: rank = 3
    integer(8), dimension(:), allocatable  :: dims

    allocate(dims(rank))
    dims(1) = size(darray,1)
    dims(2) = size(darray,2)
    dims(3) = size(darray,3)

    call h5fopen_f(trim(adjustl(fname)),h5f_acc_rdwr_f,file_id,hdf_err)
    call h5screate_f(H5S_SIMPLE_F, dspace_id, hdf_err)
    call h5sset_extent_simple_f(dspace_id, rank, dims, dims, hdf_err)
    if (present(group)) then
      call h5gopen_f(file_id, trim(adjustl(group)), grp_id, hdf_err)
      call h5dcreate_f(grp_id, trim(adjustl(dset)), h5t_native_double, dspace_id, dset_id, hdf_err)
    else
      call h5dcreate_f(file_id,trim(adjustl(dset)), h5t_native_double, dspace_id, dset_id, hdf_err)
    endif
    call h5dwrite_f(dset_id, h5t_native_double, darray, dims, hdf_err)
    call h5sclose_f(dspace_id, hdf_err)
    call h5dclose_f(dset_id, hdf_err)
    if (present(group)) then
      call h5gclose_f(grp_id, hdf_err)
    endif
    call h5fclose_f(file_id, hdf_err)

  end subroutine h5_flushdata_3d_real

  subroutine h5_flushdata_3d_complex(fname, dset, darray, group)
    character(len=*), intent(in)           :: fname
    character(len=*), intent(in)           :: dset
    complex(8), dimension(:,:,:)           :: darray
    character(len=*), intent(in), optional :: group

    integer(hid_t)                         :: file_id, grp_id, dset_id, dspace_id
    integer(hsize_t), dimension(3)         :: dset_dims, dset_maxdims
    integer                                :: rank = 3
    integer(8), dimension(:), allocatable  :: dims

    allocate(dims(rank))
    dims(1) = size(darray,1)
    dims(2) = size(darray,2)
    dims(3) = size(darray,3)

    call h5fopen_f(trim(adjustl(fname)),h5f_acc_rdwr_f,file_id,hdf_err)
    call h5screate_f(H5S_SIMPLE_F, dspace_id, hdf_err)
    call h5sset_extent_simple_f(dspace_id, rank, dims, dims, hdf_err)
    if (present(group)) then
      call h5gopen_f(file_id, trim(adjustl(group)), grp_id, hdf_err)
      call h5dcreate_f(grp_id, trim(adjustl(dset)), h5t_native_double, dspace_id, dset_id, hdf_err)
    else
      call h5dcreate_f(file_id,trim(adjustl(dset)), h5t_native_double, dspace_id, dset_id, hdf_err)
    endif
    call h5dwrite_f(dset_id, type_r_id, real(darray), dims, hdf_err)
    call h5dwrite_f(dset_id, type_r_id, aimag(darray), dims, hdf_err)
    call h5sclose_f(dspace_id, hdf_err)
    call h5dclose_f(dset_id, hdf_err)
    if (present(group)) then
      call h5gclose_f(grp_id, hdf_err)
    endif
    call h5fclose_f(file_id, hdf_err)

  end subroutine h5_flushdata_3d_complex

  subroutine h5_create_group(fname, gname, parentgroup)
    character(len=*), intent(in)           :: fname
    character(len=*), intent(in)           :: gname
    character(len=*), intent(in), optional :: parentgroup

    integer(hid_t) :: file_id
    integer(hid_t) :: grp_id_parent, grp_id

    call h5fopen_f(trim(adjustl(fname)), h5f_acc_rdwr_f, file_id, hdf_err)
    if (present(parentgroup) .and. (trim(adjustl(parentgroup)) .ne. "")) then
      call h5gopen_f(file_id, trim(adjustl(parentgroup)), grp_id_parent, hdf_err)
      call h5gcreate_f(grp_id_parent, trim(adjustl(gname)), grp_id, hdf_err) ! on parent
    else
      call h5gcreate_f(file_id, trim(adjustl(gname)), grp_id, hdf_err) ! on root
    endif
    call h5fclose_f(file_id, hdf_err)
  end subroutine h5_create_group

end module Mhdf5
