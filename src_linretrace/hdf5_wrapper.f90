module hdf5_wrapper
  use hdf5
  implicit none

  integer         :: hdf_err
  integer(hid_t)  :: complex_id_dp, complex_id_sp
  integer(hid_t)  :: complex_id_r_dp, complex_id_i_dp, complex_id_r_sp, complex_id_i_sp

  interface hdf5_read_data
    module procedure hdf5_read_data_0d_logical, &
                     hdf5_read_data_1d_logical, &
                     hdf5_read_data_2d_logical, &
                     hdf5_read_data_3d_logical, &
                     hdf5_read_data_4d_logical, &
                     hdf5_read_data_5d_logical, &
                     hdf5_read_data_6d_logical, &
                     hdf5_read_data_7d_logical, &
                     hdf5_read_data_0d_int4,    &
                     hdf5_read_data_1d_int4,    &
                     hdf5_read_data_2d_int4,    &
                     hdf5_read_data_3d_int4,    &
                     hdf5_read_data_4d_int4,    &
                     hdf5_read_data_5d_int4,    &
                     hdf5_read_data_6d_int4,    &
                     hdf5_read_data_7d_int4,    &
                     hdf5_read_data_0d_real4,   &
                     hdf5_read_data_1d_real4,   &
                     hdf5_read_data_2d_real4,   &
                     hdf5_read_data_3d_real4,   &
                     hdf5_read_data_4d_real4,   &
                     hdf5_read_data_5d_real4,   &
                     hdf5_read_data_6d_real4,   &
                     hdf5_read_data_7d_real4,   &
                     hdf5_read_data_0d_complex4,&
                     hdf5_read_data_1d_complex4,&
                     hdf5_read_data_2d_complex4,&
                     hdf5_read_data_3d_complex4,&
                     hdf5_read_data_4d_complex4,&
                     hdf5_read_data_5d_complex4,&
                     hdf5_read_data_6d_complex4,&
                     hdf5_read_data_7d_complex4,&
                     hdf5_read_data_0d_real8,   &
                     hdf5_read_data_1d_real8,   &
                     hdf5_read_data_2d_real8,   &
                     hdf5_read_data_3d_real8,   &
                     hdf5_read_data_4d_real8,   &
                     hdf5_read_data_5d_real8,   &
                     hdf5_read_data_6d_real8,   &
                     hdf5_read_data_7d_real8,   &
                     hdf5_read_data_0d_complex8,&
                     hdf5_read_data_1d_complex8,&
                     hdf5_read_data_2d_complex8,&
                     hdf5_read_data_3d_complex8,&
                     hdf5_read_data_4d_complex8,&
                     hdf5_read_data_5d_complex8,&
                     hdf5_read_data_6d_complex8,&
                     hdf5_read_data_7d_complex8
  end interface hdf5_read_data

  interface hdf5_write_data
    module procedure hdf5_write_data_0d_logical, &
                     hdf5_write_data_1d_logical, &
                     hdf5_write_data_2d_logical, &
                     hdf5_write_data_3d_logical, &
                     hdf5_write_data_4d_logical, &
                     hdf5_write_data_5d_logical, &
                     hdf5_write_data_6d_logical, &
                     hdf5_write_data_7d_logical, &
                     hdf5_write_data_0d_int4,    &
                     hdf5_write_data_1d_int4,    &
                     hdf5_write_data_2d_int4,    &
                     hdf5_write_data_3d_int4,    &
                     hdf5_write_data_4d_int4,    &
                     hdf5_write_data_5d_int4,    &
                     hdf5_write_data_6d_int4,    &
                     hdf5_write_data_7d_int4,    &
                     hdf5_write_data_0d_real4,   &
                     hdf5_write_data_1d_real4,   &
                     hdf5_write_data_2d_real4,   &
                     hdf5_write_data_3d_real4,   &
                     hdf5_write_data_4d_real4,   &
                     hdf5_write_data_5d_real4,   &
                     hdf5_write_data_6d_real4,   &
                     hdf5_write_data_7d_real4,   &
                     hdf5_write_data_0d_complex4,&
                     hdf5_write_data_1d_complex4,&
                     hdf5_write_data_2d_complex4,&
                     hdf5_write_data_3d_complex4,&
                     hdf5_write_data_4d_complex4,&
                     hdf5_write_data_5d_complex4,&
                     hdf5_write_data_6d_complex4,&
                     hdf5_write_data_7d_complex4,&
                     hdf5_write_data_0d_real8,   &
                     hdf5_write_data_1d_real8,   &
                     hdf5_write_data_2d_real8,   &
                     hdf5_write_data_3d_real8,   &
                     hdf5_write_data_4d_real8,   &
                     hdf5_write_data_5d_real8,   &
                     hdf5_write_data_6d_real8,   &
                     hdf5_write_data_7d_real8,   &
                     hdf5_write_data_0d_complex8,&
                     hdf5_write_data_1d_complex8,&
                     hdf5_write_data_2d_complex8,&
                     hdf5_write_data_3d_complex8,&
                     hdf5_write_data_4d_complex8,&
                     hdf5_write_data_5d_complex8,&
                     hdf5_write_data_6d_complex8,&
                     hdf5_write_data_7d_complex8
  end interface hdf5_write_data

  interface hdf5_write_attribute
    module procedure hdf5_write_attribute_logical, &
                     hdf5_write_attribute_int4,    &
                     hdf5_write_attribute_real4,   &
                     hdf5_write_attribute_real8,   &
                     hdf5_write_attribute_complex4,&
                     hdf5_write_attribute_complex8,&
                     hdf5_write_attribute_string
  end interface hdf5_write_attribute

  interface hdf5_read_attribute
    module procedure hdf5_read_attribute_logical, &
                     hdf5_read_attribute_int4,    &
                     hdf5_read_attribute_real4,   &
                     hdf5_read_attribute_real8,   &
                     hdf5_read_attribute_complex4,&
                     hdf5_read_attribute_complex8,&
                     hdf5_read_attribute_string
  end interface hdf5_read_attribute

  contains

  subroutine hdf5_init()
    call h5open_f(hdf_err)
    call hdf5_create_complex_datatype()
  end subroutine hdf5_init

  subroutine hdf5_finalize()
    call h5close_f(hdf_err)
  end subroutine hdf5_finalize

  subroutine hdf5_create_file(fname)
    character(len=*), intent(in) :: fname

    integer(hid_t)               :: ifile
    ! truncate -> if it already exists, erase all data
    ! other possibility: h5f_acc_excl_f -> fail if already exists
    call h5fcreate_f(trim(adjustl(fname)), h5f_acc_trunc_f, ifile, hdf_err)
    ! we only create here
    call h5fclose_f(ifile, hdf_err)
  end subroutine hdf5_create_file

  subroutine hdf5_open_file(fname, ifile, rdonly)
    character(len=*), intent(in)  :: fname
    integer(hid_t), intent(out)    :: ifile
    logical, optional :: rdonly

    if (present(rdonly)) then
      if (rdonly) then
        call h5fopen_f(trim(adjustl(fname)), h5f_acc_rdonly_f, ifile, hdf_err)
      else
        call h5fopen_f(trim(adjustl(fname)), h5f_acc_rdwr_f, ifile, hdf_err)
      endif
    else
      call h5fopen_f(trim(adjustl(fname)), h5f_acc_rdwr_f, ifile, hdf_err)
    endif
  end subroutine hdf5_open_file

  subroutine hdf5_close_file(ifile)
    integer(hid_t), intent(in)    :: ifile
    call h5fclose_f(ifile, hdf_err)
  end subroutine hdf5_close_file

  subroutine hdf5_create_group(ifile, gname, iforce_error)
    integer(hid_t), intent(in)             :: ifile
    character(len=*), intent(in)           :: gname
    logical, intent(in), optional          :: iforce_error

    integer(hid_t) :: file_id
    integer(hid_t) :: grp_parent_id, grp_id

    integer :: grpcnt
    integer :: pos,i
    logical :: force_error
    character(len=50), allocatable :: ngroup(:)
    character(len=150) :: strim

    ! adjusting string and removing possible leading and trailing '/'
    strim = adjustl(trim(gname))
    if (index(strim, '/', back=.true.) .eq. len_trim(strim)) then
      strim = strim(:(len_trim(strim)-1))
    endif
    if (index(strim, '/') .eq. 1) then
      strim = strim(2:)
    endif
    ! now we have a string with looks like 'group1/group2/group3'

    if (len_trim(strim) .eq. 0) return

    allocate(ngroup(30))
    grpcnt = 1
    do
      pos = index(strim, '/')
      if (pos .eq. 0) then
        ngroup(grpcnt) = strim(:len_trim(strim)) ! until the end
        exit
      else
        ngroup(grpcnt) = strim(:(pos-1)) ! without the '/'
        strim = strim((pos+1):) ! truncating the left
        grpcnt = grpcnt + 1
      endif
    enddo

    ! if we call create_group from the outside we call errors if we stumble upon
    ! existing groups
    if (present(iforce_error)) then
      if (iforce_error) then
        force_error = .true.
      else
        force_error = .false.
      endif
    else
      force_error = .true.
    endif
    ! now we have the number of groups (grpcnt)
    ! the according group names are stored in the array ngroup(1...grpcnt)
    call h5eset_auto_f(0,hdf_err) ! deactivate error printing
    grp_parent_id = ifile ! for the first entry into the file

    ! the last group creating is outside of the error supression.
    if (force_error) then
      do i=1,grpcnt-1
        ! check for existance
        ! this check causes all the unecessary error messages
        ! I tried doing it with h5oexists_by_name_f, but there one
        ! cannot specifically check for group ... so
        ! I will stick to this 'hack'
        call h5gopen_f(grp_parent_id, trim(ngroup(i)), grp_id, hdf_err)
        if (hdf_err .ne. 0) then ! failed
          call h5gcreate_f(grp_parent_id, trim(ngroup(i)), grp_id, hdf_err)
        endif
        if (grp_parent_id .ne. ifile) then
          call h5gclose_f(grp_parent_id, hdf_err)
        endif
        grp_parent_id = grp_id
      enddo
      call h5eclear_f(hdf_err)      ! clear the error
      call h5eset_auto_f(1,hdf_err) ! acitvate it again
      ! this will cause the error message if the group already exists
      call h5gcreate_f(grp_parent_id, trim(ngroup(grpcnt)), grp_id, hdf_err)
      if (grp_parent_id .ne. ifile) then
        call h5gclose_f(grp_parent_id, hdf_err)
      endif
    ! here everything is inside the error supression.
    else
      do i=1,grpcnt
        call h5gopen_f(grp_parent_id, trim(ngroup(i)), grp_id, hdf_err)
        if (hdf_err .ne. 0) then ! failed
          call h5gcreate_f(grp_parent_id, trim(ngroup(i)), grp_id, hdf_err)
        endif
        if (grp_parent_id .ne. ifile) then
          call h5gclose_f(grp_parent_id, hdf_err)
        endif
        grp_parent_id = grp_id
      enddo
      call h5eclear_f(hdf_err)      ! clear the error
      call h5eset_auto_f(1,hdf_err) ! acitvate it again
    endif
    call h5gclose_f(grp_id, hdf_err)
    deallocate(ngroup)

  end subroutine hdf5_create_group

  ! this code part is from ADGA
  ! github.com/abinitiodga/adga
  ! GPLv3
  subroutine hdf5_create_complex_datatype()
    integer(size_t), parameter :: zero = 0
    integer(hid_t)  :: plist_id
    integer(size_t) :: compound_size, type_sized

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

    ! single precision
    call h5pcreate_f(h5p_dataset_xfer_f, plist_id, hdf_err)
    call h5pset_preserve_f(plist_id, .true., hdf_err)
    ! compound
    call h5tget_size_f(h5t_native_real, type_sized, hdf_err)
    compound_size = 2*type_sized
    call h5tcreate_f(h5t_compound_f, compound_size, complex_id_sp, hdf_err)
    call h5tinsert_f(complex_id_sp, "r", zero, h5t_native_real, hdf_err)
    call h5tinsert_f(complex_id_sp, "i", type_sized, h5t_native_real, hdf_err)
    ! separate
    call h5tcreate_f(h5t_compound_f, type_sized, complex_id_r_sp, hdf_err)
    call h5tinsert_f(complex_id_r_sp, "r", zero, h5t_native_real, hdf_err)
    call h5tcreate_f(h5t_compound_f, type_sized, complex_id_i_sp, hdf_err)
    call h5tinsert_f(complex_id_i_sp, "i", zero, h5t_native_real, hdf_err)
    call h5pclose_f(plist_id, hdf_err)
  end subroutine hdf5_create_complex_datatype

  integer function hdf5_get_dimensions(ifile, dset)
    integer(hid_t), intent(in)   :: ifile
    character(len=*), intent(in) :: dset

    integer :: rank
    integer(hid_t) :: dset_id, dset_space_id

    call h5dopen_f(ifile, trim(adjustl(dset)), dset_id, hdf_err)
    call h5dget_space_f(dset_id, dset_space_id, hdf_err)
    call h5sget_simple_extent_ndims_f(dset_space_id, rank, hdf_err)
    call h5sclose_f(dset_space_id, hdf_err)
    call h5dclose_f(dset_id, hdf_err)

    hdf5_get_dimensions = rank
  end function hdf5_get_dimensions

  subroutine hdf5_get_shape(ifile, dset, dshape)
    integer(hid_t), intent(in)   :: ifile
    character(len=*), intent(in) :: dset
    integer, allocatable         :: dshape(:)

    integer(hsize_t), allocatable :: dshape5(:)
    integer(hsize_t), allocatable :: dmaxshape5(:)
    integer :: rank
    integer(hid_t) :: dset_id, dset_space_id

    rank = hdf5_get_dimensions(ifile, dset)
    allocate(dshape(rank), dshape5(rank), dmaxshape5(rank))

    call h5dopen_f(ifile, trim(adjustl(dset)), dset_id, hdf_err)
    call h5dget_space_f(dset_id, dset_space_id, hdf_err)
    call h5sget_simple_extent_dims_f(dset_space_id, dshape5, dmaxshape5, hdf_err)
    call h5sclose_f(dset_space_id, hdf_err)
    call h5dclose_f(dset_id, hdf_err)

    dshape = int(dshape5)
    deallocate(dshape5, dmaxshape5)
  end subroutine hdf5_get_shape

  subroutine hdf5_read_data_0d_logical(ifile, dset, darray)
    integer(hid_t), intent(in)   :: ifile
    character(len=*), intent(in) :: dset
    logical                      :: darray

    integer(4)                   :: darrayi

    call hdf5_read_data_0d_int4(ifile, dset, darrayi)
    if (darrayi == 0) then
      darray = .false.
    else
      darray = .true.
    endif
  end subroutine hdf5_read_data_0d_logical

  subroutine hdf5_read_data_1d_logical(ifile, dset, darray)
    integer(hid_t), intent(in)         :: ifile
    character(len=*), intent(in)       :: dset
    logical, allocatable, dimension(:) :: darray

    integer(4), allocatable            :: darrayi(:)
    integer(8), dimension(1)           :: dims
    integer                            :: a

    call hdf5_read_data_1d_int4(ifile, dset, darrayi)

    dims(1) = size(darrayi,1)
    allocate(darray(dims(1)))
    do a=1,dims(1)
      if (darrayi(a) == 0) then
        darray(a) = .false.
      else
        darray(a) = .true.
      endif
    enddo
    deallocate(darrayi)
  end subroutine hdf5_read_data_1d_logical

  subroutine hdf5_read_data_2d_logical(ifile, dset, darray)
    integer(hid_t), intent(in)           :: ifile
    character(len=*), intent(in)         :: dset
    logical, allocatable, dimension(:,:) :: darray

    integer(4), allocatable              :: darrayi(:,:)
    integer(8), dimension(2)             :: dims
    integer                              :: a,b

    call hdf5_read_data_2d_int4(ifile, dset, darrayi)

    dims(1) = size(darrayi,1)
    dims(2) = size(darrayi,2)
    allocate(darray(dims(1),dims(2)))
    do b=1,dims(2)
    do a=1,dims(1)
      if (darrayi(a,b) == 0) then
        darray(a,b) = .false.
      else
        darray(a,b) = .true.
      endif
    enddo
    enddo
    deallocate(darrayi)
  end subroutine hdf5_read_data_2d_logical

  subroutine hdf5_read_data_3d_logical(ifile, dset, darray)
    integer(hid_t), intent(in)             :: ifile
    character(len=*), intent(in)           :: dset
    logical, allocatable, dimension(:,:,:) :: darray

    integer(4), allocatable                :: darrayi(:,:,:)
    integer(8), dimension(3)               :: dims
    integer                                :: a,b,c

    call hdf5_read_data_3d_int4(ifile, dset, darrayi)

    dims(1) = size(darrayi,1)
    dims(2) = size(darrayi,2)
    dims(3) = size(darrayi,3)
    allocate(darray(dims(1),dims(2),dims(3)))
    do c=1,dims(3)
    do b=1,dims(2)
    do a=1,dims(1)
      if (darrayi(a,b,c) == 0) then
        darray(a,b,c) = .false.
      else
        darray(a,b,c) = .true.
      endif
    enddo
    enddo
    enddo
    deallocate(darrayi)
  end subroutine hdf5_read_data_3d_logical

  subroutine hdf5_read_data_4d_logical(ifile, dset, darray)
    integer(hid_t), intent(in)               :: ifile
    character(len=*), intent(in)             :: dset
    logical, allocatable, dimension(:,:,:,:) :: darray

    integer(4), allocatable                  :: darrayi(:,:,:,:)
    integer(8), dimension(4)                 :: dims
    integer                                  :: a,b,c,d

    call hdf5_read_data_4d_int4(ifile, dset, darrayi)

    dims(1) = size(darrayi,1)
    dims(2) = size(darrayi,2)
    dims(3) = size(darrayi,3)
    dims(4) = size(darrayi,4)
    allocate(darray(dims(1),dims(2),dims(3),dims(4)))
    do d=1,dims(4)
    do c=1,dims(3)
    do b=1,dims(2)
    do a=1,dims(1)
      if (darrayi(a,b,c,d) == 0) then
        darray(a,b,c,d) = .false.
      else
        darray(a,b,c,d) = .true.
      endif
    enddo
    enddo
    enddo
    enddo
    deallocate(darrayi)
  end subroutine hdf5_read_data_4d_logical

  subroutine hdf5_read_data_5d_logical(ifile, dset, darray)
    integer(hid_t), intent(in)                 :: ifile
    character(len=*), intent(in)               :: dset
    logical, allocatable, dimension(:,:,:,:,:) :: darray

    integer(4), allocatable                    :: darrayi(:,:,:,:,:)
    integer(8), dimension(5)                   :: dims
    integer                                    :: a,b,c,d,e

    call hdf5_read_data_5d_int4(ifile, dset, darrayi)

    dims(1) = size(darrayi,1)
    dims(2) = size(darrayi,2)
    dims(3) = size(darrayi,3)
    dims(4) = size(darrayi,4)
    dims(5) = size(darrayi,5)
    allocate(darray(dims(1),dims(2),dims(3),dims(4),dims(5)))
    do e=1,dims(5)
    do d=1,dims(4)
    do c=1,dims(3)
    do b=1,dims(2)
    do a=1,dims(1)
      if (darrayi(a,b,c,d,e) == 0) then
        darray(a,b,c,d,e) = .false.
      else
        darray(a,b,c,d,e) = .true.
      endif
    enddo
    enddo
    enddo
    enddo
    enddo
    deallocate(darrayi)
  end subroutine hdf5_read_data_5d_logical

  subroutine hdf5_read_data_6d_logical(ifile, dset, darray)
    integer(hid_t), intent(in)                   :: ifile
    character(len=*), intent(in)                 :: dset
    logical, allocatable, dimension(:,:,:,:,:,:) :: darray

    integer(4), allocatable                      :: darrayi(:,:,:,:,:,:)
    integer(8), dimension(6)                     :: dims
    integer                                      :: a,b,c,d,e,f

    call hdf5_read_data_6d_int4(ifile, dset, darrayi)

    dims(1) = size(darrayi,1)
    dims(2) = size(darrayi,2)
    dims(3) = size(darrayi,3)
    dims(4) = size(darrayi,4)
    dims(5) = size(darrayi,5)
    dims(6) = size(darrayi,6)
    allocate(darray(dims(1),dims(2),dims(3),dims(4),dims(5),dims(6)))
    do f=1,dims(6)
    do e=1,dims(5)
    do d=1,dims(4)
    do c=1,dims(3)
    do b=1,dims(2)
    do a=1,dims(1)
      if (darrayi(a,b,c,d,e,f) == 0) then
        darray(a,b,c,d,e,f) = .false.
      else
        darray(a,b,c,d,e,f) = .true.
      endif
    enddo
    enddo
    enddo
    enddo
    enddo
    enddo
    deallocate(darrayi)
  end subroutine hdf5_read_data_6d_logical

  subroutine hdf5_read_data_7d_logical(ifile, dset, darray)
    integer(hid_t), intent(in)                     :: ifile
    character(len=*), intent(in)                   :: dset
    logical, allocatable, dimension(:,:,:,:,:,:,:) :: darray

    integer(4), allocatable                        :: darrayi(:,:,:,:,:,:,:)
    integer(8), dimension(7)                       :: dims
    integer                                        :: a,b,c,d,e,f,g

    call hdf5_read_data_7d_int4(ifile, dset, darrayi)

    dims(1) = size(darrayi,1)
    dims(2) = size(darrayi,2)
    dims(3) = size(darrayi,3)
    dims(4) = size(darrayi,4)
    dims(5) = size(darrayi,5)
    dims(6) = size(darrayi,6)
    dims(7) = size(darrayi,7)
    allocate(darray(dims(1),dims(2),dims(3),dims(4),dims(5),dims(6),dims(7)))
    do g=1,dims(7)
    do f=1,dims(6)
    do e=1,dims(5)
    do d=1,dims(4)
    do c=1,dims(3)
    do b=1,dims(2)
    do a=1,dims(1)
      if (darrayi(a,b,c,d,e,f,g) == 0) then
        darray(a,b,c,d,e,f,g) = .false.
      else
        darray(a,b,c,d,e,f,g) = .true.
      endif
    enddo
    enddo
    enddo
    enddo
    enddo
    enddo
    enddo
    deallocate(darrayi)
  end subroutine hdf5_read_data_7d_logical

  subroutine hdf5_read_data_0d_int4(ifile, dset, darray)
    integer(hid_t), intent(in)   :: ifile
    character(len=*), intent(in) :: dset
    integer(4)                   :: darray

    integer(hid_t) :: dset_id
    integer(hsize_t), dimension(0) :: dset_dims

    call h5dopen_f(ifile, trim(adjustl(dset)), dset_id, hdf_err)
    call h5dread_f(dset_id, h5t_native_integer, darray, dset_dims, hdf_err)
    call h5dclose_f(dset_id, hdf_err)
  end subroutine hdf5_read_data_0d_int4

  subroutine hdf5_read_data_1d_int4(ifile, dset, darray)
    integer(hid_t), intent(in)            :: ifile
    character(len=*), intent(in)          :: dset
    integer(4), allocatable, dimension(:) :: darray

    integer(hid_t) :: dset_id, dset_space_id
    integer(hsize_t), dimension(1) :: dset_dims, dset_maxdims

    call h5dopen_f(ifile, trim(adjustl(dset)), dset_id, hdf_err)
    call h5dget_space_f(dset_id, dset_space_id, hdf_err)
    call h5sget_simple_extent_dims_f(dset_space_id, dset_dims, dset_maxdims, hdf_err)
    allocate(darray(dset_dims(1)))
    call h5dread_f(dset_id, h5t_native_integer, darray, dset_dims, hdf_err)
    call h5sclose_f(dset_space_id, hdf_err)
    call h5dclose_f(dset_id, hdf_err)
  end subroutine hdf5_read_data_1d_int4

  subroutine hdf5_read_data_2d_int4(ifile, dset, darray)
    integer(hid_t), intent(in)              :: ifile
    character(len=*), intent(in)            :: dset
    integer(4), allocatable, dimension(:,:) :: darray

    integer(hid_t) :: dset_id, dset_space_id
    integer(hsize_t), dimension(2) :: dset_dims, dset_maxdims

    call h5dopen_f(ifile, trim(adjustl(dset)), dset_id, hdf_err)
    call h5dget_space_f(dset_id, dset_space_id, hdf_err)
    call h5sget_simple_extent_dims_f(dset_space_id, dset_dims, dset_maxdims, hdf_err)
    allocate(darray(dset_dims(1), dset_dims(2)))
    call h5dread_f(dset_id, h5t_native_integer, darray, dset_dims, hdf_err)
    call h5sclose_f(dset_space_id, hdf_err)
    call h5dclose_f(dset_id, hdf_err)
  end subroutine hdf5_read_data_2d_int4

  subroutine hdf5_read_data_3d_int4(ifile, dset, darray)
    integer(hid_t), intent(in)                :: ifile
    character(len=*), intent(in)              :: dset
    integer(4), allocatable, dimension(:,:,:) :: darray

    integer(hid_t) :: dset_id, dset_space_id
    integer(hsize_t), dimension(3) :: dset_dims, dset_maxdims

    call h5dopen_f(ifile, trim(adjustl(dset)), dset_id, hdf_err)
    call h5dget_space_f(dset_id, dset_space_id, hdf_err)
    call h5sget_simple_extent_dims_f(dset_space_id, dset_dims, dset_maxdims, hdf_err)
    allocate(darray(dset_dims(1), dset_dims(2), dset_dims(3)))
    call h5dread_f(dset_id, h5t_native_integer, darray, dset_dims, hdf_err)
    call h5sclose_f(dset_space_id, hdf_err)
    call h5dclose_f(dset_id, hdf_err)
  end subroutine hdf5_read_data_3d_int4

  subroutine hdf5_read_data_4d_int4(ifile, dset, darray)
    integer(hid_t), intent(in)                  :: ifile
    character(len=*), intent(in)                :: dset
    integer(4), allocatable, dimension(:,:,:,:) :: darray

    integer(hid_t) :: dset_id, dset_space_id
    integer(hsize_t), dimension(4) :: dset_dims, dset_maxdims

    call h5dopen_f(ifile, trim(adjustl(dset)), dset_id, hdf_err)
    call h5dget_space_f(dset_id, dset_space_id, hdf_err)
    call h5sget_simple_extent_dims_f(dset_space_id, dset_dims, dset_maxdims, hdf_err)
    allocate(darray(dset_dims(1), dset_dims(2), dset_dims(3), dset_dims(4)))
    call h5dread_f(dset_id, h5t_native_integer, darray, dset_dims, hdf_err)
    call h5sclose_f(dset_space_id, hdf_err)
    call h5dclose_f(dset_id, hdf_err)
  end subroutine hdf5_read_data_4d_int4

  subroutine hdf5_read_data_5d_int4(ifile, dset, darray)
    integer(hid_t), intent(in)                    :: ifile
    character(len=*), intent(in)                  :: dset
    integer(4), allocatable, dimension(:,:,:,:,:) :: darray

    integer(hid_t) :: dset_id, dset_space_id
    integer(hsize_t), dimension(5) :: dset_dims, dset_maxdims

    call h5dopen_f(ifile, trim(adjustl(dset)), dset_id, hdf_err)
    call h5dget_space_f(dset_id, dset_space_id, hdf_err)
    call h5sget_simple_extent_dims_f(dset_space_id, dset_dims, dset_maxdims, hdf_err)
    allocate(darray(dset_dims(1), dset_dims(2), dset_dims(3), dset_dims(4), dset_dims(5)))
    call h5dread_f(dset_id, h5t_native_integer, darray, dset_dims, hdf_err)
    call h5sclose_f(dset_space_id, hdf_err)
    call h5dclose_f(dset_id, hdf_err)
  end subroutine hdf5_read_data_5d_int4

  subroutine hdf5_read_data_6d_int4(ifile, dset, darray)
    integer(hid_t), intent(in)                      :: ifile
    character(len=*), intent(in)                    :: dset
    integer(4), allocatable, dimension(:,:,:,:,:,:) :: darray

    integer(hid_t) :: dset_id, dset_space_id
    integer(hsize_t), dimension(6) :: dset_dims, dset_maxdims

    call h5dopen_f(ifile, trim(adjustl(dset)), dset_id, hdf_err)
    call h5dget_space_f(dset_id, dset_space_id, hdf_err)
    call h5sget_simple_extent_dims_f(dset_space_id, dset_dims, dset_maxdims, hdf_err)
    allocate(darray(dset_dims(1), dset_dims(2), dset_dims(3), dset_dims(4), dset_dims(5), dset_dims(6)))
    call h5dread_f(dset_id, h5t_native_integer, darray, dset_dims, hdf_err)
    call h5sclose_f(dset_space_id, hdf_err)
    call h5dclose_f(dset_id, hdf_err)
  end subroutine hdf5_read_data_6d_int4

  subroutine hdf5_read_data_7d_int4(ifile, dset, darray)
    integer(hid_t), intent(in)                        :: ifile
    character(len=*), intent(in)                      :: dset
    integer(4), allocatable, dimension(:,:,:,:,:,:,:) :: darray

    integer(hid_t) :: dset_id, dset_space_id
    integer(hsize_t), dimension(7) :: dset_dims, dset_maxdims

    call h5dopen_f(ifile, trim(adjustl(dset)), dset_id, hdf_err)
    call h5dget_space_f(dset_id, dset_space_id, hdf_err)
    call h5sget_simple_extent_dims_f(dset_space_id, dset_dims, dset_maxdims, hdf_err)
    allocate(darray(dset_dims(1), dset_dims(2), dset_dims(3), dset_dims(4), dset_dims(5), dset_dims(6), dset_dims(7)))
    call h5dread_f(dset_id, h5t_native_integer, darray, dset_dims, hdf_err)
    call h5sclose_f(dset_space_id, hdf_err)
    call h5dclose_f(dset_id, hdf_err)
  end subroutine hdf5_read_data_7d_int4

  subroutine hdf5_read_data_0d_real4(ifile, dset, darray)
    integer(hid_t), intent(in)   :: ifile
    character(len=*), intent(in) :: dset
    real(4)                      :: darray

    integer(hid_t) :: dset_id
    integer(hsize_t), dimension(0) :: dset_dims

    call h5dopen_f(ifile, trim(adjustl(dset)), dset_id, hdf_err)
    call h5dread_f(dset_id, h5t_native_real, darray, dset_dims, hdf_err)
    call h5dclose_f(dset_id, hdf_err)
  end subroutine hdf5_read_data_0d_real4

  subroutine hdf5_read_data_1d_real4(ifile, dset, darray)
    integer(hid_t), intent(in)            :: ifile
    character(len=*), intent(in)          :: dset
    real(4), allocatable, dimension(:) :: darray

    integer(hid_t) :: dset_id, dset_space_id
    integer(hsize_t), dimension(1) :: dset_dims, dset_maxdims

    call h5dopen_f(ifile, trim(adjustl(dset)), dset_id, hdf_err)
    call h5dget_space_f(dset_id, dset_space_id, hdf_err)
    call h5sget_simple_extent_dims_f(dset_space_id, dset_dims, dset_maxdims, hdf_err)
    allocate(darray(dset_dims(1)))
    call h5dread_f(dset_id, h5t_native_real, darray, dset_dims, hdf_err)
    call h5sclose_f(dset_space_id, hdf_err)
    call h5dclose_f(dset_id, hdf_err)
  end subroutine hdf5_read_data_1d_real4

  subroutine hdf5_read_data_2d_real4(ifile, dset, darray)
    integer(hid_t), intent(in)           :: ifile
    character(len=*), intent(in)         :: dset
    real(4), allocatable, dimension(:,:) :: darray

    integer(hid_t) :: dset_id, dset_space_id
    integer(hsize_t), dimension(2) :: dset_dims, dset_maxdims

    call h5dopen_f(ifile, trim(adjustl(dset)), dset_id, hdf_err)
    call h5dget_space_f(dset_id, dset_space_id, hdf_err)
    call h5sget_simple_extent_dims_f(dset_space_id, dset_dims, dset_maxdims, hdf_err)
    allocate(darray(dset_dims(1), dset_dims(2)))
    call h5dread_f(dset_id, h5t_native_real, darray, dset_dims, hdf_err)
    call h5sclose_f(dset_space_id, hdf_err)
    call h5dclose_f(dset_id, hdf_err)
  end subroutine hdf5_read_data_2d_real4

  subroutine hdf5_read_data_3d_real4(ifile, dset, darray)
    integer(hid_t), intent(in)                :: ifile
    character(len=*), intent(in)              :: dset
    real(4), allocatable, dimension(:,:,:) :: darray

    integer(hid_t) :: dset_id, dset_space_id
    integer(hsize_t), dimension(3) :: dset_dims, dset_maxdims

    call h5dopen_f(ifile, trim(adjustl(dset)), dset_id, hdf_err)
    call h5dget_space_f(dset_id, dset_space_id, hdf_err)
    call h5sget_simple_extent_dims_f(dset_space_id, dset_dims, dset_maxdims, hdf_err)
    allocate(darray(dset_dims(1), dset_dims(2), dset_dims(3)))
    call h5dread_f(dset_id, h5t_native_real, darray, dset_dims, hdf_err)
    call h5sclose_f(dset_space_id, hdf_err)
    call h5dclose_f(dset_id, hdf_err)
  end subroutine hdf5_read_data_3d_real4

  subroutine hdf5_read_data_4d_real4(ifile, dset, darray)
    integer(hid_t), intent(in)               :: ifile
    character(len=*), intent(in)             :: dset
    real(4), allocatable, dimension(:,:,:,:) :: darray

    integer(hid_t) :: dset_id, dset_space_id
    integer(hsize_t), dimension(4) :: dset_dims, dset_maxdims

    call h5dopen_f(ifile, trim(adjustl(dset)), dset_id, hdf_err)
    call h5dget_space_f(dset_id, dset_space_id, hdf_err)
    call h5sget_simple_extent_dims_f(dset_space_id, dset_dims, dset_maxdims, hdf_err)
    allocate(darray(dset_dims(1), dset_dims(2), dset_dims(3), dset_dims(4)))
    call h5dread_f(dset_id, h5t_native_real, darray, dset_dims, hdf_err)
    call h5sclose_f(dset_space_id, hdf_err)
    call h5dclose_f(dset_id, hdf_err)
  end subroutine hdf5_read_data_4d_real4

  subroutine hdf5_read_data_5d_real4(ifile, dset, darray)
    integer(hid_t), intent(in)                 :: ifile
    character(len=*), intent(in)               :: dset
    real(4), allocatable, dimension(:,:,:,:,:) :: darray

    integer(hid_t) :: dset_id, dset_space_id
    integer(hsize_t), dimension(5) :: dset_dims, dset_maxdims

    call h5dopen_f(ifile, trim(adjustl(dset)), dset_id, hdf_err)
    call h5dget_space_f(dset_id, dset_space_id, hdf_err)
    call h5sget_simple_extent_dims_f(dset_space_id, dset_dims, dset_maxdims, hdf_err)
    allocate(darray(dset_dims(1), dset_dims(2), dset_dims(3), dset_dims(4), dset_dims(5)))
    call h5dread_f(dset_id, h5t_native_real, darray, dset_dims, hdf_err)
    call h5sclose_f(dset_space_id, hdf_err)
    call h5dclose_f(dset_id, hdf_err)
  end subroutine hdf5_read_data_5d_real4

  subroutine hdf5_read_data_6d_real4(ifile, dset, darray)
    integer(hid_t), intent(in)                   :: ifile
    character(len=*), intent(in)                 :: dset
    real(4), allocatable, dimension(:,:,:,:,:,:) :: darray

    integer(hid_t) :: dset_id, dset_space_id
    integer(hsize_t), dimension(6) :: dset_dims, dset_maxdims

    call h5dopen_f(ifile, trim(adjustl(dset)), dset_id, hdf_err)
    call h5dget_space_f(dset_id, dset_space_id, hdf_err)
    call h5sget_simple_extent_dims_f(dset_space_id, dset_dims, dset_maxdims, hdf_err)
    allocate(darray(dset_dims(1), dset_dims(2), dset_dims(3), dset_dims(4), dset_dims(5), dset_dims(6)))
    call h5dread_f(dset_id, h5t_native_real, darray, dset_dims, hdf_err)
    call h5sclose_f(dset_space_id, hdf_err)
    call h5dclose_f(dset_id, hdf_err)
  end subroutine hdf5_read_data_6d_real4

  subroutine hdf5_read_data_7d_real4(ifile, dset, darray)
    integer(hid_t), intent(in)                     :: ifile
    character(len=*), intent(in)                   :: dset
    real(4), allocatable, dimension(:,:,:,:,:,:,:) :: darray

    integer(hid_t) :: dset_id, dset_space_id
    integer(hsize_t), dimension(7) :: dset_dims, dset_maxdims

    call h5dopen_f(ifile, trim(adjustl(dset)), dset_id, hdf_err)
    call h5dget_space_f(dset_id, dset_space_id, hdf_err)
    call h5sget_simple_extent_dims_f(dset_space_id, dset_dims, dset_maxdims, hdf_err)
    allocate(darray(dset_dims(1), dset_dims(2), dset_dims(3), dset_dims(4), dset_dims(5), dset_dims(6), dset_dims(7)))
    call h5dread_f(dset_id, h5t_native_real, darray, dset_dims, hdf_err)
    call h5sclose_f(dset_space_id, hdf_err)
    call h5dclose_f(dset_id, hdf_err)
  end subroutine hdf5_read_data_7d_real4

  subroutine hdf5_read_data_0d_real8(ifile, dset, darray)
    integer(hid_t), intent(in)   :: ifile
    character(len=*), intent(in) :: dset
    real(8)                      :: darray

    integer(hid_t) :: dset_id
    integer(hsize_t), dimension(0) :: dset_dims

    call h5dopen_f(ifile, trim(adjustl(dset)), dset_id, hdf_err)
    call h5dread_f(dset_id, h5t_native_double, darray, dset_dims, hdf_err)
    call h5dclose_f(dset_id, hdf_err)
  end subroutine hdf5_read_data_0d_real8

  subroutine hdf5_read_data_1d_real8(ifile, dset, darray)
    integer(hid_t), intent(in)            :: ifile
    character(len=*), intent(in)          :: dset
    real(8), allocatable, dimension(:) :: darray

    integer(hid_t) :: dset_id, dset_space_id
    integer(hsize_t), dimension(1) :: dset_dims, dset_maxdims

    call h5dopen_f(ifile, trim(adjustl(dset)), dset_id, hdf_err)
    call h5dget_space_f(dset_id, dset_space_id, hdf_err)
    call h5sget_simple_extent_dims_f(dset_space_id, dset_dims, dset_maxdims, hdf_err)
    allocate(darray(dset_dims(1)))
    call h5dread_f(dset_id, h5t_native_double, darray, dset_dims, hdf_err)
    call h5sclose_f(dset_space_id, hdf_err)
    call h5dclose_f(dset_id, hdf_err)
  end subroutine hdf5_read_data_1d_real8

  subroutine hdf5_read_data_2d_real8(ifile, dset, darray)
    integer(hid_t), intent(in)           :: ifile
    character(len=*), intent(in)         :: dset
    real(8), allocatable, dimension(:,:) :: darray

    integer(hid_t) :: dset_id, dset_space_id
    integer(hsize_t), dimension(2) :: dset_dims, dset_maxdims

    call h5dopen_f(ifile, trim(adjustl(dset)), dset_id, hdf_err)
    call h5dget_space_f(dset_id, dset_space_id, hdf_err)
    call h5sget_simple_extent_dims_f(dset_space_id, dset_dims, dset_maxdims, hdf_err)
    allocate(darray(dset_dims(1), dset_dims(2)))
    call h5dread_f(dset_id, h5t_native_double, darray, dset_dims, hdf_err)
    call h5sclose_f(dset_space_id, hdf_err)
    call h5dclose_f(dset_id, hdf_err)
  end subroutine hdf5_read_data_2d_real8

  subroutine hdf5_read_data_3d_real8(ifile, dset, darray)
    integer(hid_t), intent(in)                :: ifile
    character(len=*), intent(in)              :: dset
    real(8), allocatable, dimension(:,:,:) :: darray

    integer(hid_t) :: dset_id, dset_space_id
    integer(hsize_t), dimension(3) :: dset_dims, dset_maxdims

    call h5dopen_f(ifile, trim(adjustl(dset)), dset_id, hdf_err)
    call h5dget_space_f(dset_id, dset_space_id, hdf_err)
    call h5sget_simple_extent_dims_f(dset_space_id, dset_dims, dset_maxdims, hdf_err)
    allocate(darray(dset_dims(1), dset_dims(2), dset_dims(3)))
    call h5dread_f(dset_id, h5t_native_double, darray, dset_dims, hdf_err)
    call h5sclose_f(dset_space_id, hdf_err)
    call h5dclose_f(dset_id, hdf_err)
  end subroutine hdf5_read_data_3d_real8

  subroutine hdf5_read_data_4d_real8(ifile, dset, darray)
    integer(hid_t), intent(in)               :: ifile
    character(len=*), intent(in)             :: dset
    real(8), allocatable, dimension(:,:,:,:) :: darray

    integer(hid_t) :: dset_id, dset_space_id
    integer(hsize_t), dimension(4) :: dset_dims, dset_maxdims

    call h5dopen_f(ifile, trim(adjustl(dset)), dset_id, hdf_err)
    call h5dget_space_f(dset_id, dset_space_id, hdf_err)
    call h5sget_simple_extent_dims_f(dset_space_id, dset_dims, dset_maxdims, hdf_err)
    allocate(darray(dset_dims(1), dset_dims(2), dset_dims(3), dset_dims(4)))
    call h5dread_f(dset_id, h5t_native_double, darray, dset_dims, hdf_err)
    call h5sclose_f(dset_space_id, hdf_err)
    call h5dclose_f(dset_id, hdf_err)
  end subroutine hdf5_read_data_4d_real8

  subroutine hdf5_read_data_5d_real8(ifile, dset, darray)
    integer(hid_t), intent(in)                 :: ifile
    character(len=*), intent(in)               :: dset
    real(8), allocatable, dimension(:,:,:,:,:) :: darray

    integer(hid_t) :: dset_id, dset_space_id
    integer(hsize_t), dimension(5) :: dset_dims, dset_maxdims

    call h5dopen_f(ifile, trim(adjustl(dset)), dset_id, hdf_err)
    call h5dget_space_f(dset_id, dset_space_id, hdf_err)
    call h5sget_simple_extent_dims_f(dset_space_id, dset_dims, dset_maxdims, hdf_err)
    allocate(darray(dset_dims(1), dset_dims(2), dset_dims(3), dset_dims(4), dset_dims(5)))
    call h5dread_f(dset_id, h5t_native_double, darray, dset_dims, hdf_err)
    call h5sclose_f(dset_space_id, hdf_err)
    call h5dclose_f(dset_id, hdf_err)
  end subroutine hdf5_read_data_5d_real8

  subroutine hdf5_read_data_6d_real8(ifile, dset, darray)
    integer(hid_t), intent(in)                   :: ifile
    character(len=*), intent(in)                 :: dset
    real(8), allocatable, dimension(:,:,:,:,:,:) :: darray

    integer(hid_t) :: dset_id, dset_space_id
    integer(hsize_t), dimension(6) :: dset_dims, dset_maxdims

    call h5dopen_f(ifile, trim(adjustl(dset)), dset_id, hdf_err)
    call h5dget_space_f(dset_id, dset_space_id, hdf_err)
    call h5sget_simple_extent_dims_f(dset_space_id, dset_dims, dset_maxdims, hdf_err)
    allocate(darray(dset_dims(1), dset_dims(2), dset_dims(3), dset_dims(4), dset_dims(5), dset_dims(6)))
    call h5dread_f(dset_id, h5t_native_double, darray, dset_dims, hdf_err)
    call h5sclose_f(dset_space_id, hdf_err)
    call h5dclose_f(dset_id, hdf_err)
  end subroutine hdf5_read_data_6d_real8

  subroutine hdf5_read_data_7d_real8(ifile, dset, darray)
    integer(hid_t), intent(in)                     :: ifile
    character(len=*), intent(in)                   :: dset
    real(8), allocatable, dimension(:,:,:,:,:,:,:) :: darray

    integer(hid_t) :: dset_id, dset_space_id
    integer(hsize_t), dimension(7) :: dset_dims, dset_maxdims

    call h5dopen_f(ifile, trim(adjustl(dset)), dset_id, hdf_err)
    call h5dget_space_f(dset_id, dset_space_id, hdf_err)
    call h5sget_simple_extent_dims_f(dset_space_id, dset_dims, dset_maxdims, hdf_err)
    allocate(darray(dset_dims(1), dset_dims(2), dset_dims(3), dset_dims(4), dset_dims(5), dset_dims(6), dset_dims(7)))
    call h5dread_f(dset_id, h5t_native_double, darray, dset_dims, hdf_err)
    call h5sclose_f(dset_space_id, hdf_err)
    call h5dclose_f(dset_id, hdf_err)
  end subroutine hdf5_read_data_7d_real8

  subroutine hdf5_read_data_0d_complex4(ifile, dset, darray)
    integer(hid_t), intent(in)   :: ifile
    character(len=*), intent(in) :: dset
    complex(4)                   :: darray

    complex(4) :: ci = (0.0, 1.0)
    real(4) :: tmp_r, tmp_i
    integer(hid_t) :: dset_id, dset_space_id
    integer(hsize_t), dimension(0) :: dset_dims

    call h5dopen_f(ifile, trim(adjustl(dset)), dset_id, hdf_err)
    call h5dread_f(dset_id, complex_id_r_sp, tmp_r, dset_dims, hdf_err)
    call h5dread_f(dset_id, complex_id_i_sp, tmp_i, dset_dims, hdf_err)
    call h5dclose_f(dset_id, hdf_err)
    darray = tmp_r + ci*tmp_i
  end subroutine hdf5_read_data_0d_complex4

  subroutine hdf5_read_data_1d_complex4(ifile, dset, darray)
    integer(hid_t), intent(in)            :: ifile
    character(len=*), intent(in)          :: dset
    complex(4), allocatable, dimension(:) :: darray

    complex(4) :: ci = (0.0, 1.0)
    real(4), allocatable, dimension(:)    :: tmp_i, tmp_r
    integer(hid_t) :: dset_id, dset_space_id
    integer(hsize_t), dimension(1) :: dset_dims, dset_maxdims

    call h5dopen_f(ifile, trim(adjustl(dset)), dset_id, hdf_err)
    call h5dget_space_f(dset_id, dset_space_id, hdf_err)
    call h5sget_simple_extent_dims_f(dset_space_id, dset_dims, dset_maxdims, hdf_err)
    allocate(darray(dset_dims(1)))
    allocate(tmp_i(dset_dims(1)))
    allocate(tmp_r(dset_dims(1)))
    call h5dread_f(dset_id, complex_id_r_sp, tmp_r, dset_dims, hdf_err)
    call h5dread_f(dset_id, complex_id_i_sp, tmp_i, dset_dims, hdf_err)
    call h5sclose_f(dset_space_id, hdf_err)
    call h5dclose_f(dset_id, hdf_err)
    darray = tmp_r + ci*tmp_i
    deallocate(tmp_i, tmp_r)
  end subroutine hdf5_read_data_1d_complex4

  subroutine hdf5_read_data_2d_complex4(ifile, dset, darray)
    integer(hid_t), intent(in)              :: ifile
    character(len=*), intent(in)            :: dset
    complex(4), allocatable, dimension(:,:) :: darray

    complex(4) :: ci = (0.0, 1.0)
    real(4), allocatable, dimension(:,:)    :: tmp_i, tmp_r
    integer(hid_t) :: dset_id, dset_space_id
    integer(hsize_t), dimension(2) :: dset_dims, dset_maxdims

    call h5dopen_f(ifile, trim(adjustl(dset)), dset_id, hdf_err)
    call h5dget_space_f(dset_id, dset_space_id, hdf_err)
    call h5sget_simple_extent_dims_f(dset_space_id, dset_dims, dset_maxdims, hdf_err)
    allocate(darray(dset_dims(1), dset_dims(2)))
    allocate(tmp_i(dset_dims(1), dset_dims(2)))
    allocate(tmp_r(dset_dims(1), dset_dims(2)))
    call h5dread_f(dset_id, complex_id_r_sp, tmp_r, dset_dims, hdf_err)
    call h5dread_f(dset_id, complex_id_i_sp, tmp_i, dset_dims, hdf_err)
    call h5sclose_f(dset_space_id, hdf_err)
    call h5dclose_f(dset_id, hdf_err)
    darray = tmp_r + ci*tmp_i
    deallocate(tmp_i, tmp_r)
  end subroutine hdf5_read_data_2d_complex4

  subroutine hdf5_read_data_3d_complex4(ifile, dset, darray)
    integer(hid_t), intent(in)                :: ifile
    character(len=*), intent(in)              :: dset
    complex(4), allocatable, dimension(:,:,:) :: darray

    complex(4) :: ci = (0.0, 1.0)
    real(4), allocatable, dimension(:,:,:)    :: tmp_i, tmp_r
    integer(hid_t) :: dset_id, dset_space_id
    integer(hsize_t), dimension(3) :: dset_dims, dset_maxdims

    call h5dopen_f(ifile, trim(adjustl(dset)), dset_id, hdf_err)
    call h5dget_space_f(dset_id, dset_space_id, hdf_err)
    call h5sget_simple_extent_dims_f(dset_space_id, dset_dims, dset_maxdims, hdf_err)
    allocate(darray(dset_dims(1), dset_dims(2), dset_dims(3)))
    allocate(tmp_i(dset_dims(1), dset_dims(2), dset_dims(3)))
    allocate(tmp_r(dset_dims(1), dset_dims(2), dset_dims(3)))
    call h5dread_f(dset_id, complex_id_r_sp, tmp_r, dset_dims, hdf_err)
    call h5dread_f(dset_id, complex_id_i_sp, tmp_i, dset_dims, hdf_err)
    call h5sclose_f(dset_space_id, hdf_err)
    call h5dclose_f(dset_id, hdf_err)
    darray = tmp_r + ci*tmp_i
    deallocate(tmp_i, tmp_r)
  end subroutine hdf5_read_data_3d_complex4

  subroutine hdf5_read_data_4d_complex4(ifile, dset, darray)
    integer(hid_t), intent(in)                  :: ifile
    character(len=*), intent(in)                :: dset
    complex(4), allocatable, dimension(:,:,:,:) :: darray

    complex(4) :: ci = (0.0, 1.0)
    real(4), allocatable, dimension(:,:,:,:)    :: tmp_i, tmp_r
    integer(hid_t) :: dset_id, dset_space_id
    integer(hsize_t), dimension(4) :: dset_dims, dset_maxdims

    call h5dopen_f(ifile, trim(adjustl(dset)), dset_id, hdf_err)
    call h5dget_space_f(dset_id, dset_space_id, hdf_err)
    call h5sget_simple_extent_dims_f(dset_space_id, dset_dims, dset_maxdims, hdf_err)
    allocate(darray(dset_dims(1), dset_dims(2), dset_dims(3), dset_dims(4)))
    allocate(tmp_i(dset_dims(1), dset_dims(2), dset_dims(3), dset_dims(4)))
    allocate(tmp_r(dset_dims(1), dset_dims(2), dset_dims(3), dset_dims(4)))
    call h5dread_f(dset_id, complex_id_r_sp, tmp_r, dset_dims, hdf_err)
    call h5dread_f(dset_id, complex_id_i_sp, tmp_i, dset_dims, hdf_err)
    call h5sclose_f(dset_space_id, hdf_err)
    call h5dclose_f(dset_id, hdf_err)
    darray = tmp_r + ci*tmp_i
    deallocate(tmp_i, tmp_r)
  end subroutine hdf5_read_data_4d_complex4

  subroutine hdf5_read_data_5d_complex4(ifile, dset, darray)
    integer(hid_t), intent(in)                    :: ifile
    character(len=*), intent(in)                  :: dset
    complex(4), allocatable, dimension(:,:,:,:,:) :: darray

    complex(4) :: ci = (0.0, 1.0)
    real(4), allocatable, dimension(:,:,:,:,:)    :: tmp_i, tmp_r
    integer(hid_t) :: dset_id, dset_space_id
    integer(hsize_t), dimension(5) :: dset_dims, dset_maxdims

    call h5dopen_f(ifile, trim(adjustl(dset)), dset_id, hdf_err)
    call h5dget_space_f(dset_id, dset_space_id, hdf_err)
    call h5sget_simple_extent_dims_f(dset_space_id, dset_dims, dset_maxdims, hdf_err)
    allocate(darray(dset_dims(1), dset_dims(2), dset_dims(3), dset_dims(4), dset_dims(5)))
    allocate(tmp_i(dset_dims(1), dset_dims(2), dset_dims(3), dset_dims(4), dset_dims(5)))
    allocate(tmp_r(dset_dims(1), dset_dims(2), dset_dims(3), dset_dims(4), dset_dims(5)))
    call h5dread_f(dset_id, complex_id_r_sp, tmp_r, dset_dims, hdf_err)
    call h5dread_f(dset_id, complex_id_i_sp, tmp_i, dset_dims, hdf_err)
    call h5sclose_f(dset_space_id, hdf_err)
    call h5dclose_f(dset_id, hdf_err)
    darray = tmp_r + ci*tmp_i
    deallocate(tmp_i, tmp_r)
  end subroutine hdf5_read_data_5d_complex4

  subroutine hdf5_read_data_6d_complex4(ifile, dset, darray)
    integer(hid_t), intent(in)                      :: ifile
    character(len=*), intent(in)                    :: dset
    complex(4), allocatable, dimension(:,:,:,:,:,:) :: darray

    complex(4) :: ci = (0.0, 1.0)
    real(4), allocatable, dimension(:,:,:,:,:,:)    :: tmp_i, tmp_r
    integer(hid_t) :: dset_id, dset_space_id
    integer(hsize_t), dimension(6) :: dset_dims, dset_maxdims

    call h5dopen_f(ifile, trim(adjustl(dset)), dset_id, hdf_err)
    call h5dget_space_f(dset_id, dset_space_id, hdf_err)
    call h5sget_simple_extent_dims_f(dset_space_id, dset_dims, dset_maxdims, hdf_err)
    allocate(darray(dset_dims(1), dset_dims(2), dset_dims(3), dset_dims(4), dset_dims(5), dset_dims(6)))
    allocate(tmp_i(dset_dims(1), dset_dims(2), dset_dims(3), dset_dims(4), dset_dims(5), dset_dims(6)))
    allocate(tmp_r(dset_dims(1), dset_dims(2), dset_dims(3), dset_dims(4), dset_dims(5), dset_dims(6)))
    call h5dread_f(dset_id, complex_id_r_sp, tmp_r, dset_dims, hdf_err)
    call h5dread_f(dset_id, complex_id_i_sp, tmp_i, dset_dims, hdf_err)
    call h5sclose_f(dset_space_id, hdf_err)
    call h5dclose_f(dset_id, hdf_err)
    darray = tmp_r + ci*tmp_i
    deallocate(tmp_i, tmp_r)
  end subroutine hdf5_read_data_6d_complex4

  subroutine hdf5_read_data_7d_complex4(ifile, dset, darray)
    integer(hid_t), intent(in)                        :: ifile
    character(len=*), intent(in)                      :: dset
    complex(4), allocatable, dimension(:,:,:,:,:,:,:) :: darray

    complex(4) :: ci = (0.0, 1.0)
    real(4), allocatable, dimension(:,:,:,:,:,:,:)    :: tmp_i, tmp_r
    integer(hid_t) :: dset_id, dset_space_id
    integer(hsize_t), dimension(7) :: dset_dims, dset_maxdims

    call h5dopen_f(ifile, trim(adjustl(dset)), dset_id, hdf_err)
    call h5dget_space_f(dset_id, dset_space_id, hdf_err)
    call h5sget_simple_extent_dims_f(dset_space_id, dset_dims, dset_maxdims, hdf_err)
    allocate(darray(dset_dims(1), dset_dims(2), dset_dims(3), dset_dims(4), dset_dims(5), dset_dims(6), dset_dims(7)))
    allocate(tmp_i(dset_dims(1), dset_dims(2), dset_dims(3), dset_dims(4), dset_dims(5), dset_dims(6), dset_dims(7)))
    allocate(tmp_r(dset_dims(1), dset_dims(2), dset_dims(3), dset_dims(4), dset_dims(5), dset_dims(6), dset_dims(7)))
    call h5dread_f(dset_id, complex_id_r_sp, tmp_r, dset_dims, hdf_err)
    call h5dread_f(dset_id, complex_id_i_sp, tmp_i, dset_dims, hdf_err)
    call h5sclose_f(dset_space_id, hdf_err)
    call h5dclose_f(dset_id, hdf_err)
    darray = tmp_r + ci*tmp_i
    deallocate(tmp_i, tmp_r)
  end subroutine hdf5_read_data_7d_complex4

  subroutine hdf5_read_data_0d_complex8(ifile, dset, darray)
    integer(hid_t), intent(in)   :: ifile
    character(len=*), intent(in) :: dset
    complex(8)                   :: darray

    complex(8) :: ci = (0.d0, 1.d0)
    real(8) :: tmp_r, tmp_i
    integer(hid_t) :: dset_id, dset_space_id
    integer(hsize_t), dimension(0) :: dset_dims

    call h5dopen_f(ifile, trim(adjustl(dset)), dset_id, hdf_err)
    call h5dread_f(dset_id, complex_id_r_dp, tmp_r, dset_dims, hdf_err)
    call h5dread_f(dset_id, complex_id_i_dp, tmp_i, dset_dims, hdf_err)
    call h5dclose_f(dset_id, hdf_err)
    darray = tmp_r + ci*tmp_i
  end subroutine hdf5_read_data_0d_complex8

  subroutine hdf5_read_data_1d_complex8(ifile, dset, darray)
    integer(hid_t), intent(in)            :: ifile
    character(len=*), intent(in)          :: dset
    complex(8), allocatable, dimension(:) :: darray

    complex(8) :: ci = (0.d0, 1.d0)
    real(8), allocatable, dimension(:)    :: tmp_i, tmp_r
    integer(hid_t) :: dset_id, dset_space_id
    integer(hsize_t), dimension(1) :: dset_dims, dset_maxdims

    call h5dopen_f(ifile, trim(adjustl(dset)), dset_id, hdf_err)
    call h5dget_space_f(dset_id, dset_space_id, hdf_err)
    call h5sget_simple_extent_dims_f(dset_space_id, dset_dims, dset_maxdims, hdf_err)
    allocate(darray(dset_dims(1)))
    allocate(tmp_i(dset_dims(1)))
    allocate(tmp_r(dset_dims(1)))
    call h5dread_f(dset_id, complex_id_r_dp, tmp_r, dset_dims, hdf_err)
    call h5dread_f(dset_id, complex_id_i_dp, tmp_i, dset_dims, hdf_err)
    call h5sclose_f(dset_space_id, hdf_err)
    call h5dclose_f(dset_id, hdf_err)
    darray = tmp_r + ci*tmp_i
    deallocate(tmp_i, tmp_r)
  end subroutine hdf5_read_data_1d_complex8

  subroutine hdf5_read_data_2d_complex8(ifile, dset, darray)
    integer(hid_t), intent(in)              :: ifile
    character(len=*), intent(in)            :: dset
    complex(8), allocatable, dimension(:,:) :: darray

    complex(8) :: ci = (0.d0, 1.d0)
    real(8), allocatable, dimension(:,:)    :: tmp_i, tmp_r
    integer(hid_t) :: dset_id, dset_space_id
    integer(hsize_t), dimension(2) :: dset_dims, dset_maxdims

    call h5dopen_f(ifile, trim(adjustl(dset)), dset_id, hdf_err)
    call h5dget_space_f(dset_id, dset_space_id, hdf_err)
    call h5sget_simple_extent_dims_f(dset_space_id, dset_dims, dset_maxdims, hdf_err)
    allocate(darray(dset_dims(1), dset_dims(2)))
    allocate(tmp_i(dset_dims(1), dset_dims(2)))
    allocate(tmp_r(dset_dims(1), dset_dims(2)))
    call h5dread_f(dset_id, complex_id_r_dp, tmp_r, dset_dims, hdf_err)
    call h5dread_f(dset_id, complex_id_i_dp, tmp_i, dset_dims, hdf_err)
    call h5sclose_f(dset_space_id, hdf_err)
    call h5dclose_f(dset_id, hdf_err)
    darray = tmp_r + ci*tmp_i
    deallocate(tmp_i, tmp_r)
  end subroutine hdf5_read_data_2d_complex8

  subroutine hdf5_read_data_3d_complex8(ifile, dset, darray)
    integer(hid_t), intent(in)                :: ifile
    character(len=*), intent(in)              :: dset
    complex(8), allocatable, dimension(:,:,:) :: darray

    complex(8) :: ci = (0.d0, 1.d0)
    real(8), allocatable, dimension(:,:,:)    :: tmp_i, tmp_r
    integer(hid_t) :: dset_id, dset_space_id
    integer(hsize_t), dimension(3) :: dset_dims, dset_maxdims

    call h5dopen_f(ifile, trim(adjustl(dset)), dset_id, hdf_err)
    call h5dget_space_f(dset_id, dset_space_id, hdf_err)
    call h5sget_simple_extent_dims_f(dset_space_id, dset_dims, dset_maxdims, hdf_err)
    allocate(darray(dset_dims(1), dset_dims(2), dset_dims(3)))
    allocate(tmp_i(dset_dims(1), dset_dims(2), dset_dims(3)))
    allocate(tmp_r(dset_dims(1), dset_dims(2), dset_dims(3)))
    call h5dread_f(dset_id, complex_id_r_dp, tmp_r, dset_dims, hdf_err)
    call h5dread_f(dset_id, complex_id_i_dp, tmp_i, dset_dims, hdf_err)
    call h5sclose_f(dset_space_id, hdf_err)
    call h5dclose_f(dset_id, hdf_err)
    darray = tmp_r + ci*tmp_i
    deallocate(tmp_i, tmp_r)
  end subroutine hdf5_read_data_3d_complex8

  subroutine hdf5_read_data_4d_complex8(ifile, dset, darray)
    integer(hid_t), intent(in)                  :: ifile
    character(len=*), intent(in)                :: dset
    complex(8), allocatable, dimension(:,:,:,:) :: darray

    complex(8) :: ci = (0.d0, 1.d0)
    real(8), allocatable, dimension(:,:,:,:)    :: tmp_i, tmp_r
    integer(hid_t) :: dset_id, dset_space_id
    integer(hsize_t), dimension(4) :: dset_dims, dset_maxdims

    call h5dopen_f(ifile, trim(adjustl(dset)), dset_id, hdf_err)
    call h5dget_space_f(dset_id, dset_space_id, hdf_err)
    call h5sget_simple_extent_dims_f(dset_space_id, dset_dims, dset_maxdims, hdf_err)
    allocate(darray(dset_dims(1), dset_dims(2), dset_dims(3), dset_dims(4)))
    allocate(tmp_i(dset_dims(1), dset_dims(2), dset_dims(3), dset_dims(4)))
    allocate(tmp_r(dset_dims(1), dset_dims(2), dset_dims(3), dset_dims(4)))
    call h5dread_f(dset_id, complex_id_r_dp, tmp_r, dset_dims, hdf_err)
    call h5dread_f(dset_id, complex_id_i_dp, tmp_i, dset_dims, hdf_err)
    call h5sclose_f(dset_space_id, hdf_err)
    call h5dclose_f(dset_id, hdf_err)
    darray = tmp_r + ci*tmp_i
    deallocate(tmp_i, tmp_r)
  end subroutine hdf5_read_data_4d_complex8

  subroutine hdf5_read_data_5d_complex8(ifile, dset, darray)
    integer(hid_t), intent(in)                    :: ifile
    character(len=*), intent(in)                  :: dset
    complex(8), allocatable, dimension(:,:,:,:,:) :: darray

    complex(8) :: ci = (0.d0, 1.d0)
    real(8), allocatable, dimension(:,:,:,:,:)    :: tmp_i, tmp_r
    integer(hid_t) :: dset_id, dset_space_id
    integer(hsize_t), dimension(5) :: dset_dims, dset_maxdims

    call h5dopen_f(ifile, trim(adjustl(dset)), dset_id, hdf_err)
    call h5dget_space_f(dset_id, dset_space_id, hdf_err)
    call h5sget_simple_extent_dims_f(dset_space_id, dset_dims, dset_maxdims, hdf_err)
    allocate(darray(dset_dims(1), dset_dims(2), dset_dims(3), dset_dims(4), dset_dims(5)))
    allocate(tmp_i(dset_dims(1), dset_dims(2), dset_dims(3), dset_dims(4), dset_dims(5)))
    allocate(tmp_r(dset_dims(1), dset_dims(2), dset_dims(3), dset_dims(4), dset_dims(5)))
    call h5dread_f(dset_id, complex_id_r_dp, tmp_r, dset_dims, hdf_err)
    call h5dread_f(dset_id, complex_id_i_dp, tmp_i, dset_dims, hdf_err)
    call h5sclose_f(dset_space_id, hdf_err)
    call h5dclose_f(dset_id, hdf_err)
    darray = tmp_r + ci*tmp_i
    deallocate(tmp_i, tmp_r)
  end subroutine hdf5_read_data_5d_complex8

  subroutine hdf5_read_data_6d_complex8(ifile, dset, darray)
    integer(hid_t), intent(in)                      :: ifile
    character(len=*), intent(in)                    :: dset
    complex(8), allocatable, dimension(:,:,:,:,:,:) :: darray

    complex(8) :: ci = (0.d0, 1.d0)
    real(8), allocatable, dimension(:,:,:,:,:,:)    :: tmp_i, tmp_r
    integer(hid_t) :: dset_id, dset_space_id
    integer(hsize_t), dimension(6) :: dset_dims, dset_maxdims

    call h5dopen_f(ifile, trim(adjustl(dset)), dset_id, hdf_err)
    call h5dget_space_f(dset_id, dset_space_id, hdf_err)
    call h5sget_simple_extent_dims_f(dset_space_id, dset_dims, dset_maxdims, hdf_err)
    allocate(darray(dset_dims(1), dset_dims(2), dset_dims(3), dset_dims(4), dset_dims(5), dset_dims(6)))
    allocate(tmp_i(dset_dims(1), dset_dims(2), dset_dims(3), dset_dims(4), dset_dims(5), dset_dims(6)))
    allocate(tmp_r(dset_dims(1), dset_dims(2), dset_dims(3), dset_dims(4), dset_dims(5), dset_dims(6)))
    call h5dread_f(dset_id, complex_id_r_dp, tmp_r, dset_dims, hdf_err)
    call h5dread_f(dset_id, complex_id_i_dp, tmp_i, dset_dims, hdf_err)
    call h5sclose_f(dset_space_id, hdf_err)
    call h5dclose_f(dset_id, hdf_err)
    darray = tmp_r + ci*tmp_i
    deallocate(tmp_i, tmp_r)
  end subroutine hdf5_read_data_6d_complex8

  subroutine hdf5_read_data_7d_complex8(ifile, dset, darray)
    integer(hid_t), intent(in)                        :: ifile
    character(len=*), intent(in)                      :: dset
    complex(8), allocatable, dimension(:,:,:,:,:,:,:) :: darray

    complex(8) :: ci = (0.d0, 1.d0)
    real(8), allocatable, dimension(:,:,:,:,:,:,:)    :: tmp_i, tmp_r
    integer(hid_t) :: dset_id, dset_space_id
    integer(hsize_t), dimension(7) :: dset_dims, dset_maxdims

    call h5dopen_f(ifile, trim(adjustl(dset)), dset_id, hdf_err)
    call h5dget_space_f(dset_id, dset_space_id, hdf_err)
    call h5sget_simple_extent_dims_f(dset_space_id, dset_dims, dset_maxdims, hdf_err)
    allocate(darray(dset_dims(1), dset_dims(2), dset_dims(3), dset_dims(4), dset_dims(5), dset_dims(6), dset_dims(7)))
    allocate(tmp_i(dset_dims(1), dset_dims(2), dset_dims(3), dset_dims(4), dset_dims(5), dset_dims(6), dset_dims(7)))
    allocate(tmp_r(dset_dims(1), dset_dims(2), dset_dims(3), dset_dims(4), dset_dims(5), dset_dims(6), dset_dims(7)))
    call h5dread_f(dset_id, complex_id_r_dp, tmp_r, dset_dims, hdf_err)
    call h5dread_f(dset_id, complex_id_i_dp, tmp_i, dset_dims, hdf_err)
    call h5sclose_f(dset_space_id, hdf_err)
    call h5dclose_f(dset_id, hdf_err)
    darray = tmp_r + ci*tmp_i
    deallocate(tmp_i, tmp_r)
  end subroutine hdf5_read_data_7d_complex8

  subroutine hdf5_write_data_0d_logical(ifile, dsetfull, darray)
    integer(hid_t), intent(in)   :: ifile
    character(len=*), intent(in) :: dsetfull
    logical                      :: darray

    integer :: darrayi

    if (darray) then
      darrayi = 1
    else
      darrayi = 0
    endif
    call hdf5_write_data_0d_int4(ifile, dsetfull, darrayi)
  end subroutine hdf5_write_data_0d_logical

  subroutine hdf5_write_data_1d_logical(ifile, dsetfull, darray)
    integer(hid_t), intent(in)   :: ifile
    character(len=*), intent(in) :: dsetfull
    logical, dimension(:)        :: darray

    integer(4), allocatable      :: darrayi(:)
    integer(8), dimension(1)     :: dims
    integer                      :: a

    dims(1) = size(darray,1)
    allocate(darrayi(dims(1)))
    do a=1,dims(1)
      if (darray(a)) then
        darrayi(a) = 1
      else
        darrayi(a) = 0
      endif
    enddo
    call hdf5_write_data_1d_int4(ifile, dsetfull, darrayi)
    deallocate(darrayi)
  end subroutine hdf5_write_data_1d_logical

  subroutine hdf5_write_data_2d_logical(ifile, dsetfull, darray)
    integer(hid_t), intent(in)   :: ifile
    character(len=*), intent(in) :: dsetfull
    logical, dimension(:,:)      :: darray

    integer(4), allocatable      :: darrayi(:,:)
    integer(8), dimension(2)     :: dims
    integer                      :: a,b

    dims(1) = size(darray,1)
    dims(2) = size(darray,2)
    allocate(darrayi(dims(1),dims(2)))
    do b=1,dims(2)
    do a=1,dims(1)
      if (darray(a,b)) then
        darrayi(a,b) = 1
      else
        darrayi(a,b) = 0
      endif
    enddo
    enddo
    call hdf5_write_data_2d_int4(ifile, dsetfull, darrayi)
    deallocate(darrayi)
  end subroutine hdf5_write_data_2d_logical

  subroutine hdf5_write_data_3d_logical(ifile, dsetfull, darray)
    integer(hid_t), intent(in)   :: ifile
    character(len=*), intent(in) :: dsetfull
    logical, dimension(:,:,:)    :: darray

    integer(4), allocatable      :: darrayi(:,:,:)
    integer(8), dimension(3)     :: dims
    integer                      :: a,b,c

    dims(1) = size(darray,1)
    dims(2) = size(darray,2)
    dims(3) = size(darray,3)
    allocate(darrayi(dims(1),dims(2),dims(3)))
    do c=1,dims(3)
    do b=1,dims(2)
    do a=1,dims(1)
      if (darray(a,b,c)) then
        darrayi(a,b,c) = 1
      else
        darrayi(a,b,c) = 0
      endif
    enddo
    enddo
    enddo
    call hdf5_write_data_3d_int4(ifile, dsetfull, darrayi)
    deallocate(darrayi)
  end subroutine hdf5_write_data_3d_logical

  subroutine hdf5_write_data_4d_logical(ifile, dsetfull, darray)
    integer(hid_t), intent(in)   :: ifile
    character(len=*), intent(in) :: dsetfull
    logical, dimension(:,:,:,:)  :: darray

    integer(4), allocatable      :: darrayi(:,:,:,:)
    integer(8), dimension(4)     :: dims
    integer                      :: a,b,c,d

    dims(1) = size(darray,1)
    dims(2) = size(darray,2)
    dims(3) = size(darray,3)
    dims(4) = size(darray,4)
    allocate(darrayi(dims(1),dims(2),dims(3),dims(4)))
    do d=1,dims(4)
    do c=1,dims(3)
    do b=1,dims(2)
    do a=1,dims(1)
      if (darray(a,b,c,d)) then
        darrayi(a,b,c,d) = 1
      else
        darrayi(a,b,c,d) = 0
      endif
    enddo
    enddo
    enddo
    enddo
    call hdf5_write_data_4d_int4(ifile, dsetfull, darrayi)
    deallocate(darrayi)
  end subroutine hdf5_write_data_4d_logical

  subroutine hdf5_write_data_5d_logical(ifile, dsetfull, darray)
    integer(hid_t), intent(in)    :: ifile
    character(len=*), intent(in)  :: dsetfull
    logical, dimension(:,:,:,:,:) :: darray

    integer(4), allocatable       :: darrayi(:,:,:,:,:)
    integer(8), dimension(5)      :: dims
    integer                       :: a,b,c,d,e

    dims(1) = size(darray,1)
    dims(2) = size(darray,2)
    dims(3) = size(darray,3)
    dims(4) = size(darray,4)
    dims(5) = size(darray,5)
    allocate(darrayi(dims(1),dims(2),dims(3),dims(4),dims(5)))
    do e=1,dims(5)
    do d=1,dims(4)
    do c=1,dims(3)
    do b=1,dims(2)
    do a=1,dims(1)
      if (darray(a,b,c,d,e)) then
        darrayi(a,b,c,d,e) = 1
      else
        darrayi(a,b,c,d,e) = 0
      endif
    enddo
    enddo
    enddo
    enddo
    enddo
    call hdf5_write_data_5d_int4(ifile, dsetfull, darrayi)
    deallocate(darrayi)
  end subroutine hdf5_write_data_5d_logical

  subroutine hdf5_write_data_6d_logical(ifile, dsetfull, darray)
    integer(hid_t), intent(in)      :: ifile
    character(len=*), intent(in)    :: dsetfull
    logical, dimension(:,:,:,:,:,:) :: darray

    integer(4), allocatable         :: darrayi(:,:,:,:,:,:)
    integer(8), dimension(6)        :: dims
    integer                         :: a,b,c,d,e,f

    dims(1) = size(darray,1)
    dims(2) = size(darray,2)
    dims(3) = size(darray,3)
    dims(4) = size(darray,4)
    dims(5) = size(darray,5)
    dims(6) = size(darray,6)
    allocate(darrayi(dims(1),dims(2),dims(3),dims(4),dims(5),dims(6)))
    do f=1,dims(6)
    do e=1,dims(5)
    do d=1,dims(4)
    do c=1,dims(3)
    do b=1,dims(2)
    do a=1,dims(1)
      if (darray(a,b,c,d,e,f)) then
        darrayi(a,b,c,d,e,f) = 1
      else
        darrayi(a,b,c,d,e,f) = 0
      endif
    enddo
    enddo
    enddo
    enddo
    enddo
    enddo
    call hdf5_write_data_6d_int4(ifile, dsetfull, darrayi)
    deallocate(darrayi)
  end subroutine hdf5_write_data_6d_logical

  subroutine hdf5_write_data_7d_logical(ifile, dsetfull, darray)
    integer(hid_t), intent(in)        :: ifile
    character(len=*), intent(in)      :: dsetfull
    logical, dimension(:,:,:,:,:,:,:) :: darray

    integer(4), allocatable           :: darrayi(:,:,:,:,:,:,:)
    integer(8), dimension(7)          :: dims
    integer                           :: a,b,c,d,e,f,g

    dims(1) = size(darray,1)
    dims(2) = size(darray,2)
    dims(3) = size(darray,3)
    dims(4) = size(darray,4)
    dims(5) = size(darray,5)
    dims(6) = size(darray,6)
    dims(7) = size(darray,7)
    allocate(darrayi(dims(1),dims(2),dims(3),dims(4),dims(5),dims(6),dims(7)))
    do g=1,dims(7)
    do f=1,dims(6)
    do e=1,dims(5)
    do d=1,dims(4)
    do c=1,dims(3)
    do b=1,dims(2)
    do a=1,dims(1)
      if (darray(a,b,c,d,e,f,g)) then
        darrayi(a,b,c,d,e,f,g) = 1
      else
        darrayi(a,b,c,d,e,f,g) = 0
      endif
    enddo
    enddo
    enddo
    enddo
    enddo
    enddo
    enddo
    call hdf5_write_data_7d_int4(ifile, dsetfull, darrayi)
    deallocate(darrayi)
  end subroutine hdf5_write_data_7d_logical

  subroutine hdf5_write_data_0d_int4(ifile, dsetfull, darray)
    integer(hid_t), intent(in)   :: ifile
    character(len=*), intent(in) :: dsetfull
    integer(4)                   :: darray

    integer(hid_t)               :: grp_id, dset_id, dspace_id
    integer, parameter           :: rank = 0
    integer(8), dimension(0)     :: dims

    character(len=150) :: gname, dset

    call hdf5_help_separate_dsetname(dsetfull, gname, dset)
    call hdf5_create_group(ifile, gname, iforce_error=.false.)

    if (trim(adjustl(gname)) == '') then
      grp_id = ifile
    else
      call h5gopen_f(ifile, trim(adjustl(gname)), grp_id, hdf_err)
    endif
    call h5screate_simple_f(0, dims, dspace_id, hdf_err)
    call h5dcreate_f(grp_id, trim(adjustl(dset)), h5t_native_integer, dspace_id, dset_id, hdf_err)
    call h5dwrite_f(dset_id, h5t_native_integer, darray, dims, hdf_err)
    call h5dclose_f(dset_id, hdf_err)
    call h5sclose_f(dspace_id, hdf_err)
    if (grp_id .ne. ifile) then
      call h5gclose_f(grp_id, hdf_err)
    endif
  end subroutine hdf5_write_data_0d_int4

  subroutine hdf5_write_data_1d_int4(ifile, dsetfull, darray)
    integer(hid_t), intent(in)     :: ifile
    character(len=*), intent(in)   :: dsetfull
    integer(4), dimension(:)       :: darray

    integer(hid_t)                 :: grp_id, dset_id, dspace_id
    integer(hsize_t), dimension(1) :: dset_dims, dset_maxdims
    integer, parameter             :: rank = 1
    integer(8), dimension(1)       :: dims

    character(len=150) :: gname, dset

    dims(1) = size(darray)

    call hdf5_help_separate_dsetname(dsetfull, gname, dset)
    call hdf5_create_group(ifile, gname, iforce_error=.false.)

    if (trim(adjustl(gname)) == '') then
      grp_id = ifile
    else
      call h5gopen_f(ifile, trim(adjustl(gname)), grp_id, hdf_err)
    endif
    call h5screate_f(h5s_simple_f, dspace_id, hdf_err)
    call h5sset_extent_simple_f(dspace_id, rank, dims, dims, hdf_err)
    call h5dcreate_f(grp_id, trim(adjustl(dset)), h5t_native_integer, dspace_id, dset_id, hdf_err)
    call h5dwrite_f(dset_id, h5t_native_integer, darray, dims, hdf_err)
    call h5dclose_f(dset_id, hdf_err)
    call h5sclose_f(dspace_id, hdf_err)
    if (grp_id .ne. ifile) then
      call h5gclose_f(grp_id, hdf_err)
    endif
  end subroutine hdf5_write_data_1d_int4

  subroutine hdf5_write_data_2d_int4(ifile, dsetfull, darray)
    integer(hid_t), intent(in)     :: ifile
    character(len=*), intent(in)   :: dsetfull
    integer(4), dimension(:,:)     :: darray

    integer(hid_t)                 :: grp_id, dset_id, dspace_id
    integer(hsize_t), dimension(2) :: dset_dims, dset_maxdims
    integer, parameter             :: rank = 2
    integer(8), dimension(2)       :: dims

    character(len=150) :: gname, dset

    dims(1) = size(darray,1)
    dims(2) = size(darray,2)

    call hdf5_help_separate_dsetname(dsetfull, gname, dset)
    call hdf5_create_group(ifile, gname, iforce_error=.false.)

    if (trim(adjustl(gname)) == '') then
      grp_id = ifile
    else
      call h5gopen_f(ifile, trim(adjustl(gname)), grp_id, hdf_err)
    endif
    call h5screate_f(h5s_simple_f, dspace_id, hdf_err)
    call h5sset_extent_simple_f(dspace_id, rank, dims, dims, hdf_err)
    call h5dcreate_f(grp_id, trim(adjustl(dset)), h5t_native_integer, dspace_id, dset_id, hdf_err)
    call h5dwrite_f(dset_id, h5t_native_integer, darray, dims, hdf_err)
    call h5dclose_f(dset_id, hdf_err)
    call h5sclose_f(dspace_id, hdf_err)
    if (grp_id .ne. ifile) then
      call h5gclose_f(grp_id, hdf_err)
    endif
  end subroutine hdf5_write_data_2d_int4

  subroutine hdf5_write_data_3d_int4(ifile, dsetfull, darray)
    integer(hid_t), intent(in)     :: ifile
    character(len=*), intent(in)   :: dsetfull
    integer(4), dimension(:,:,:)   :: darray

    integer(hid_t)                 :: grp_id, dset_id, dspace_id
    integer(hsize_t), dimension(3) :: dset_dims, dset_maxdims
    integer, parameter             :: rank = 3
    integer(8), dimension(3)       :: dims

    character(len=150) :: gname, dset

    dims(1) = size(darray,1)
    dims(2) = size(darray,2)
    dims(3) = size(darray,3)

    call hdf5_help_separate_dsetname(dsetfull, gname, dset)
    call hdf5_create_group(ifile, gname, iforce_error=.false.)

    if (trim(adjustl(gname)) == '') then
      grp_id = ifile
    else
      call h5gopen_f(ifile, trim(adjustl(gname)), grp_id, hdf_err)
    endif
    call h5screate_f(h5s_simple_f, dspace_id, hdf_err)
    call h5sset_extent_simple_f(dspace_id, rank, dims, dims, hdf_err)
    call h5dcreate_f(grp_id, trim(adjustl(dset)), h5t_native_integer, dspace_id, dset_id, hdf_err)
    call h5dwrite_f(dset_id, h5t_native_integer, darray, dims, hdf_err)
    call h5dclose_f(dset_id, hdf_err)
    call h5sclose_f(dspace_id, hdf_err)
    if (grp_id .ne. ifile) then
      call h5gclose_f(grp_id, hdf_err)
    endif
  end subroutine hdf5_write_data_3d_int4

  subroutine hdf5_write_data_4d_int4(ifile, dsetfull, darray)
    integer(hid_t), intent(in)     :: ifile
    character(len=*), intent(in)   :: dsetfull
    integer(4), dimension(:,:,:,:) :: darray

    integer(hid_t)                 :: grp_id, dset_id, dspace_id
    integer(hsize_t), dimension(4) :: dset_dims, dset_maxdims
    integer, parameter             :: rank = 4
    integer(8), dimension(4)       :: dims

    character(len=150) :: gname, dset

    dims(1) = size(darray,1)
    dims(2) = size(darray,2)
    dims(3) = size(darray,3)
    dims(4) = size(darray,4)

    call hdf5_help_separate_dsetname(dsetfull, gname, dset)
    call hdf5_create_group(ifile, gname, iforce_error=.false.)

    if (trim(adjustl(gname)) == '') then
      grp_id = ifile
    else
      call h5gopen_f(ifile, trim(adjustl(gname)), grp_id, hdf_err)
    endif
    call h5screate_f(h5s_simple_f, dspace_id, hdf_err)
    call h5sset_extent_simple_f(dspace_id, rank, dims, dims, hdf_err)
    call h5dcreate_f(grp_id, trim(adjustl(dset)), h5t_native_integer, dspace_id, dset_id, hdf_err)
    call h5dwrite_f(dset_id, h5t_native_integer, darray, dims, hdf_err)
    call h5dclose_f(dset_id, hdf_err)
    call h5sclose_f(dspace_id, hdf_err)
    if (grp_id .ne. ifile) then
      call h5gclose_f(grp_id, hdf_err)
    endif
  end subroutine hdf5_write_data_4d_int4

  subroutine hdf5_write_data_5d_int4(ifile, dsetfull, darray)
    integer(hid_t), intent(in)       :: ifile
    character(len=*), intent(in)     :: dsetfull
    integer(4), dimension(:,:,:,:,:) :: darray

    integer(hid_t)                   :: grp_id, dset_id, dspace_id
    integer(hsize_t), dimension(5)   :: dset_dims, dset_maxdims
    integer, parameter               :: rank = 5
    integer(8), dimension(5)         :: dims

    character(len=150) :: gname, dset

    dims(1) = size(darray,1)
    dims(2) = size(darray,2)
    dims(3) = size(darray,3)
    dims(4) = size(darray,4)
    dims(5) = size(darray,5)

    call hdf5_help_separate_dsetname(dsetfull, gname, dset)
    call hdf5_create_group(ifile, gname, iforce_error=.false.)

    if (trim(adjustl(gname)) == '') then
      grp_id = ifile
    else
      call h5gopen_f(ifile, trim(adjustl(gname)), grp_id, hdf_err)
    endif
    call h5screate_f(h5s_simple_f, dspace_id, hdf_err)
    call h5sset_extent_simple_f(dspace_id, rank, dims, dims, hdf_err)
    call h5dcreate_f(grp_id, trim(adjustl(dset)), h5t_native_integer, dspace_id, dset_id, hdf_err)
    call h5dwrite_f(dset_id, h5t_native_integer, darray, dims, hdf_err)
    call h5dclose_f(dset_id, hdf_err)
    call h5sclose_f(dspace_id, hdf_err)
    if (grp_id .ne. ifile) then
      call h5gclose_f(grp_id, hdf_err)
    endif
  end subroutine hdf5_write_data_5d_int4

  subroutine hdf5_write_data_6d_int4(ifile, dsetfull, darray)
    integer(hid_t), intent(in)         :: ifile
    character(len=*), intent(in)       :: dsetfull
    integer(4), dimension(:,:,:,:,:,:) :: darray

    integer(hid_t)                     :: grp_id, dset_id, dspace_id
    integer(hsize_t), dimension(6)     :: dset_dims, dset_maxdims
    integer, parameter                 :: rank = 6
    integer(8), dimension(6)           :: dims

    character(len=150) :: gname, dset

    dims(1) = size(darray,1)
    dims(2) = size(darray,2)
    dims(3) = size(darray,3)
    dims(4) = size(darray,4)
    dims(5) = size(darray,5)
    dims(6) = size(darray,6)

    call hdf5_help_separate_dsetname(dsetfull, gname, dset)
    call hdf5_create_group(ifile, gname, iforce_error=.false.)

    if (trim(adjustl(gname)) == '') then
      grp_id = ifile
    else
      call h5gopen_f(ifile, trim(adjustl(gname)), grp_id, hdf_err)
    endif
    call h5screate_f(h5s_simple_f, dspace_id, hdf_err)
    call h5sset_extent_simple_f(dspace_id, rank, dims, dims, hdf_err)
    call h5dcreate_f(grp_id, trim(adjustl(dset)), h5t_native_integer, dspace_id, dset_id, hdf_err)
    call h5dwrite_f(dset_id, h5t_native_integer, darray, dims, hdf_err)
    call h5dclose_f(dset_id, hdf_err)
    call h5sclose_f(dspace_id, hdf_err)
    if (grp_id .ne. ifile) then
      call h5gclose_f(grp_id, hdf_err)
    endif
  end subroutine hdf5_write_data_6d_int4

  subroutine hdf5_write_data_7d_int4(ifile, dsetfull, darray)
    integer(hid_t), intent(in)           :: ifile
    character(len=*), intent(in)         :: dsetfull
    integer(4), dimension(:,:,:,:,:,:,:) :: darray

    integer(hid_t)                       :: grp_id, dset_id, dspace_id
    integer(hsize_t), dimension(7)       :: dset_dims, dset_maxdims
    integer, parameter                   :: rank = 7
    integer(8), dimension(7)             :: dims

    character(len=150) :: gname, dset

    dims(1) = size(darray,1)
    dims(2) = size(darray,2)
    dims(3) = size(darray,3)
    dims(4) = size(darray,4)
    dims(5) = size(darray,5)
    dims(6) = size(darray,6)
    dims(7) = size(darray,7)

    call hdf5_help_separate_dsetname(dsetfull, gname, dset)
    call hdf5_create_group(ifile, gname, iforce_error=.false.)

    if (trim(adjustl(gname)) == '') then
      grp_id = ifile
    else
      call h5gopen_f(ifile, trim(adjustl(gname)), grp_id, hdf_err)
    endif
    call h5screate_f(h5s_simple_f, dspace_id, hdf_err)
    call h5sset_extent_simple_f(dspace_id, rank, dims, dims, hdf_err)
    call h5dcreate_f(grp_id, trim(adjustl(dset)), h5t_native_integer, dspace_id, dset_id, hdf_err)
    call h5dwrite_f(dset_id, h5t_native_integer, darray, dims, hdf_err)
    call h5dclose_f(dset_id, hdf_err)
    call h5sclose_f(dspace_id, hdf_err)
    if (grp_id .ne. ifile) then
      call h5gclose_f(grp_id, hdf_err)
    endif
  end subroutine hdf5_write_data_7d_int4

  subroutine hdf5_write_data_0d_real4(ifile, dsetfull, darray)
    integer(hid_t), intent(in)   :: ifile
    character(len=*), intent(in) :: dsetfull
    real(4)                      :: darray

    integer(hid_t)               :: grp_id, dset_id, dspace_id
    integer, parameter           :: rank = 0
    integer(8), dimension(0)     :: dims

    character(len=150) :: gname, dset

    call hdf5_help_separate_dsetname(dsetfull, gname, dset)
    call hdf5_create_group(ifile, gname, iforce_error=.false.)

    if (trim(adjustl(gname)) == '') then
      grp_id = ifile
    else
      call h5gopen_f(ifile, trim(adjustl(gname)), grp_id, hdf_err)
    endif
    call h5screate_simple_f(0, dims, dspace_id, hdf_err)
    call h5dcreate_f(grp_id, trim(adjustl(dset)), h5t_native_real, dspace_id, dset_id, hdf_err)
    call h5dwrite_f(dset_id, h5t_native_real, darray, dims, hdf_err)
    call h5dclose_f(dset_id, hdf_err)
    call h5sclose_f(dspace_id, hdf_err)
    if (grp_id .ne. ifile) then
      call h5gclose_f(grp_id, hdf_err)
    endif
  end subroutine hdf5_write_data_0d_real4

  subroutine hdf5_write_data_1d_real4(ifile, dsetfull, darray)
    integer(hid_t), intent(in)     :: ifile
    character(len=*), intent(in)   :: dsetfull
    real(4), dimension(:)          :: darray

    integer(hid_t)                 :: grp_id, dset_id, dspace_id
    integer(hsize_t), dimension(1) :: dset_dims, dset_maxdims
    integer, parameter             :: rank = 1
    integer(8), dimension(1)       :: dims

    character(len=150) :: gname, dset

    dims(1) = size(darray)

    call hdf5_help_separate_dsetname(dsetfull, gname, dset)
    call hdf5_create_group(ifile, gname, iforce_error=.false.)

    if (trim(adjustl(gname)) == '') then
      grp_id = ifile
    else
      call h5gopen_f(ifile, trim(adjustl(gname)), grp_id, hdf_err)
    endif
    call h5screate_f(h5s_simple_f, dspace_id, hdf_err)
    call h5sset_extent_simple_f(dspace_id, rank, dims, dims, hdf_err)
    call h5dcreate_f(grp_id, trim(adjustl(dset)), h5t_native_real, dspace_id, dset_id, hdf_err)
    call h5dwrite_f(dset_id, h5t_native_real, darray, dims, hdf_err)
    call h5dclose_f(dset_id, hdf_err)
    call h5sclose_f(dspace_id, hdf_err)
    if (grp_id .ne. ifile) then
      call h5gclose_f(grp_id, hdf_err)
    endif
  end subroutine hdf5_write_data_1d_real4

  subroutine hdf5_write_data_2d_real4(ifile, dsetfull, darray)
    integer(hid_t), intent(in)     :: ifile
    character(len=*), intent(in)   :: dsetfull
    real(4), dimension(:,:)        :: darray

    integer(hid_t)                 :: grp_id, dset_id, dspace_id
    integer(hsize_t), dimension(2) :: dset_dims, dset_maxdims
    integer, parameter             :: rank = 2
    integer(8), dimension(2)       :: dims

    character(len=150) :: gname, dset

    dims(1) = size(darray,1)
    dims(2) = size(darray,2)

    call hdf5_help_separate_dsetname(dsetfull, gname, dset)
    call hdf5_create_group(ifile, gname, iforce_error=.false.)

    if (trim(adjustl(gname)) == '') then
      grp_id = ifile
    else
      call h5gopen_f(ifile, trim(adjustl(gname)), grp_id, hdf_err)
    endif
    call h5screate_f(h5s_simple_f, dspace_id, hdf_err)
    call h5sset_extent_simple_f(dspace_id, rank, dims, dims, hdf_err)
    call h5dcreate_f(grp_id, trim(adjustl(dset)), h5t_native_real, dspace_id, dset_id, hdf_err)
    call h5dwrite_f(dset_id, h5t_native_real, darray, dims, hdf_err)
    call h5dclose_f(dset_id, hdf_err)
    call h5sclose_f(dspace_id, hdf_err)
    if (grp_id .ne. ifile) then
      call h5gclose_f(grp_id, hdf_err)
    endif
  end subroutine hdf5_write_data_2d_real4

  subroutine hdf5_write_data_3d_real4(ifile, dsetfull, darray)
    integer(hid_t), intent(in)     :: ifile
    character(len=*), intent(in)   :: dsetfull
    real(4), dimension(:,:,:)      :: darray

    integer(hid_t)                 :: grp_id, dset_id, dspace_id
    integer(hsize_t), dimension(3) :: dset_dims, dset_maxdims
    integer, parameter             :: rank = 3
    integer(8), dimension(3)       :: dims

    character(len=150) :: gname, dset

    dims(1) = size(darray,1)
    dims(2) = size(darray,2)
    dims(3) = size(darray,3)

    call hdf5_help_separate_dsetname(dsetfull, gname, dset)
    call hdf5_create_group(ifile, gname, iforce_error=.false.)

    if (trim(adjustl(gname)) == '') then
      grp_id = ifile
    else
      call h5gopen_f(ifile, trim(adjustl(gname)), grp_id, hdf_err)
    endif
    call h5screate_f(h5s_simple_f, dspace_id, hdf_err)
    call h5sset_extent_simple_f(dspace_id, rank, dims, dims, hdf_err)
    call h5dcreate_f(grp_id, trim(adjustl(dset)), h5t_native_real, dspace_id, dset_id, hdf_err)
    call h5dwrite_f(dset_id, h5t_native_real, darray, dims, hdf_err)
    call h5dclose_f(dset_id, hdf_err)
    call h5sclose_f(dspace_id, hdf_err)
    if (grp_id .ne. ifile) then
      call h5gclose_f(grp_id, hdf_err)
    endif
  end subroutine hdf5_write_data_3d_real4

  subroutine hdf5_write_data_4d_real4(ifile, dsetfull, darray)
    integer(hid_t), intent(in)     :: ifile
    character(len=*), intent(in)   :: dsetfull
    real(4), dimension(:,:,:,:)    :: darray

    integer(hid_t)                 :: grp_id, dset_id, dspace_id
    integer(hsize_t), dimension(4) :: dset_dims, dset_maxdims
    integer, parameter             :: rank = 4
    integer(8), dimension(4)       :: dims

    character(len=150) :: gname, dset

    dims(1) = size(darray,1)
    dims(2) = size(darray,2)
    dims(3) = size(darray,3)
    dims(4) = size(darray,4)

    call hdf5_help_separate_dsetname(dsetfull, gname, dset)
    call hdf5_create_group(ifile, gname, iforce_error=.false.)

    if (trim(adjustl(gname)) == '') then
      grp_id = ifile
    else
      call h5gopen_f(ifile, trim(adjustl(gname)), grp_id, hdf_err)
    endif
    call h5screate_f(h5s_simple_f, dspace_id, hdf_err)
    call h5sset_extent_simple_f(dspace_id, rank, dims, dims, hdf_err)
    call h5dcreate_f(grp_id, trim(adjustl(dset)), h5t_native_real, dspace_id, dset_id, hdf_err)
    call h5dwrite_f(dset_id, h5t_native_real, darray, dims, hdf_err)
    call h5dclose_f(dset_id, hdf_err)
    call h5sclose_f(dspace_id, hdf_err)
    if (grp_id .ne. ifile) then
      call h5gclose_f(grp_id, hdf_err)
    endif
  end subroutine hdf5_write_data_4d_real4

  subroutine hdf5_write_data_5d_real4(ifile, dsetfull, darray)
    integer(hid_t), intent(in)     :: ifile
    character(len=*), intent(in)   :: dsetfull
    real(4), dimension(:,:,:,:,:)  :: darray

    integer(hid_t)                 :: grp_id, dset_id, dspace_id
    integer(hsize_t), dimension(5) :: dset_dims, dset_maxdims
    integer, parameter             :: rank = 5
    integer(8), dimension(5)       :: dims

    character(len=150) :: gname, dset

    dims(1) = size(darray,1)
    dims(2) = size(darray,2)
    dims(3) = size(darray,3)
    dims(4) = size(darray,4)
    dims(5) = size(darray,5)

    call hdf5_help_separate_dsetname(dsetfull, gname, dset)
    call hdf5_create_group(ifile, gname, iforce_error=.false.)

    if (trim(adjustl(gname)) == '') then
      grp_id = ifile
    else
      call h5gopen_f(ifile, trim(adjustl(gname)), grp_id, hdf_err)
    endif
    call h5screate_f(h5s_simple_f, dspace_id, hdf_err)
    call h5sset_extent_simple_f(dspace_id, rank, dims, dims, hdf_err)
    call h5dcreate_f(grp_id, trim(adjustl(dset)), h5t_native_real, dspace_id, dset_id, hdf_err)
    call h5dwrite_f(dset_id, h5t_native_real, darray, dims, hdf_err)
    call h5dclose_f(dset_id, hdf_err)
    call h5sclose_f(dspace_id, hdf_err)
    if (grp_id .ne. ifile) then
      call h5gclose_f(grp_id, hdf_err)
    endif
  end subroutine hdf5_write_data_5d_real4

  subroutine hdf5_write_data_6d_real4(ifile, dsetfull, darray)
    integer(hid_t), intent(in)      :: ifile
    character(len=*), intent(in)    :: dsetfull
    real(4), dimension(:,:,:,:,:,:) :: darray

    integer(hid_t)                  :: grp_id, dset_id, dspace_id
    integer(hsize_t), dimension(6)  :: dset_dims, dset_maxdims
    integer, parameter              :: rank = 6
    integer(8), dimension(6)        :: dims

    character(len=150) :: gname, dset

    dims(1) = size(darray,1)
    dims(2) = size(darray,2)
    dims(3) = size(darray,3)
    dims(4) = size(darray,4)
    dims(5) = size(darray,5)
    dims(6) = size(darray,6)

    call hdf5_help_separate_dsetname(dsetfull, gname, dset)
    call hdf5_create_group(ifile, gname, iforce_error=.false.)

    if (trim(adjustl(gname)) == '') then
      grp_id = ifile
    else
      call h5gopen_f(ifile, trim(adjustl(gname)), grp_id, hdf_err)
    endif
    call h5screate_f(h5s_simple_f, dspace_id, hdf_err)
    call h5sset_extent_simple_f(dspace_id, rank, dims, dims, hdf_err)
    call h5dcreate_f(grp_id, trim(adjustl(dset)), h5t_native_real, dspace_id, dset_id, hdf_err)
    call h5dwrite_f(dset_id, h5t_native_real, darray, dims, hdf_err)
    call h5dclose_f(dset_id, hdf_err)
    call h5sclose_f(dspace_id, hdf_err)
    if (grp_id .ne. ifile) then
      call h5gclose_f(grp_id, hdf_err)
    endif
  end subroutine hdf5_write_data_6d_real4

  subroutine hdf5_write_data_7d_real4(ifile, dsetfull, darray)
    integer(hid_t), intent(in)        :: ifile
    character(len=*), intent(in)      :: dsetfull
    real(4), dimension(:,:,:,:,:,:,:) :: darray

    integer(hid_t)                    :: grp_id, dset_id, dspace_id
    integer(hsize_t), dimension(7)    :: dset_dims, dset_maxdims
    integer, parameter                :: rank = 7
    integer(8), dimension(7)          :: dims

    character(len=150) :: gname, dset

    dims(1) = size(darray,1)
    dims(2) = size(darray,2)
    dims(3) = size(darray,3)
    dims(4) = size(darray,4)
    dims(5) = size(darray,5)
    dims(6) = size(darray,6)
    dims(7) = size(darray,7)

    call hdf5_help_separate_dsetname(dsetfull, gname, dset)
    call hdf5_create_group(ifile, gname, iforce_error=.false.)

    if (trim(adjustl(gname)) == '') then
      grp_id = ifile
    else
      call h5gopen_f(ifile, trim(adjustl(gname)), grp_id, hdf_err)
    endif
    call h5screate_f(h5s_simple_f, dspace_id, hdf_err)
    call h5sset_extent_simple_f(dspace_id, rank, dims, dims, hdf_err)
    call h5dcreate_f(grp_id, trim(adjustl(dset)), h5t_native_real, dspace_id, dset_id, hdf_err)
    call h5dwrite_f(dset_id, h5t_native_real, darray, dims, hdf_err)
    call h5dclose_f(dset_id, hdf_err)
    call h5sclose_f(dspace_id, hdf_err)
    if (grp_id .ne. ifile) then
      call h5gclose_f(grp_id, hdf_err)
    endif
  end subroutine hdf5_write_data_7d_real4

  subroutine hdf5_write_data_0d_real8(ifile, dsetfull, darray)
    integer(hid_t), intent(in)   :: ifile
    character(len=*), intent(in) :: dsetfull
    real(8)                      :: darray

    integer(hid_t)               :: grp_id, dset_id, dspace_id
    integer, parameter           :: rank = 0
    integer(8), dimension(0)     :: dims

    character(len=150) :: gname, dset

    call hdf5_help_separate_dsetname(dsetfull, gname, dset)
    call hdf5_create_group(ifile, gname, iforce_error=.false.)

    if (trim(adjustl(gname)) == '') then
      grp_id = ifile
    else
      call h5gopen_f(ifile, trim(adjustl(gname)), grp_id, hdf_err)
    endif
    call h5screate_simple_f(0, dims, dspace_id, hdf_err)
    call h5dcreate_f(grp_id, trim(adjustl(dset)), h5t_native_double, dspace_id, dset_id, hdf_err)
    call h5dwrite_f(dset_id, h5t_native_double, darray, dims, hdf_err)
    call h5dclose_f(dset_id, hdf_err)
    call h5sclose_f(dspace_id, hdf_err)
    if (grp_id .ne. ifile) then
      call h5gclose_f(grp_id, hdf_err)
    endif
  end subroutine hdf5_write_data_0d_real8

  subroutine hdf5_write_data_1d_real8(ifile, dsetfull, darray)
    integer(hid_t), intent(in)     :: ifile
    character(len=*), intent(in)   :: dsetfull
    real(8), dimension(:)          :: darray

    integer(hid_t)                 :: grp_id, dset_id, dspace_id
    integer(hsize_t), dimension(1) :: dset_dims, dset_maxdims
    integer, parameter             :: rank = 1
    integer(8), dimension(1)       :: dims

    character(len=150) :: gname, dset

    dims(1) = size(darray)

    call hdf5_help_separate_dsetname(dsetfull, gname, dset)
    call hdf5_create_group(ifile, gname, iforce_error=.false.)

    if (trim(adjustl(gname)) == '') then
      grp_id = ifile
    else
      call h5gopen_f(ifile, trim(adjustl(gname)), grp_id, hdf_err)
    endif
    call h5screate_f(h5s_simple_f, dspace_id, hdf_err)
    call h5sset_extent_simple_f(dspace_id, rank, dims, dims, hdf_err)
    call h5dcreate_f(grp_id, trim(adjustl(dset)), h5t_native_double, dspace_id, dset_id, hdf_err)
    call h5dwrite_f(dset_id, h5t_native_double, darray, dims, hdf_err)
    call h5dclose_f(dset_id, hdf_err)
    call h5sclose_f(dspace_id, hdf_err)
    if (grp_id .ne. ifile) then
      call h5gclose_f(grp_id, hdf_err)
    endif
  end subroutine hdf5_write_data_1d_real8

  subroutine hdf5_write_data_2d_real8(ifile, dsetfull, darray)
    integer(hid_t), intent(in)     :: ifile
    character(len=*), intent(in)   :: dsetfull
    real(8), dimension(:,:)        :: darray

    integer(hid_t)                 :: grp_id, dset_id, dspace_id
    integer(hsize_t), dimension(2) :: dset_dims, dset_maxdims
    integer, parameter             :: rank = 2
    integer(8), dimension(2)       :: dims

    character(len=150) :: gname, dset

    dims(1) = size(darray,1)
    dims(2) = size(darray,2)

    call hdf5_help_separate_dsetname(dsetfull, gname, dset)
    call hdf5_create_group(ifile, gname, iforce_error=.false.)

    if (trim(adjustl(gname)) == '') then
      grp_id = ifile
    else
      call h5gopen_f(ifile, trim(adjustl(gname)), grp_id, hdf_err)
    endif
    call h5screate_f(h5s_simple_f, dspace_id, hdf_err)
    call h5sset_extent_simple_f(dspace_id, rank, dims, dims, hdf_err)
    call h5dcreate_f(grp_id, trim(adjustl(dset)), h5t_native_double, dspace_id, dset_id, hdf_err)
    call h5dwrite_f(dset_id, h5t_native_double, darray, dims, hdf_err)
    call h5dclose_f(dset_id, hdf_err)
    call h5sclose_f(dspace_id, hdf_err)
    if (grp_id .ne. ifile) then
      call h5gclose_f(grp_id, hdf_err)
    endif
  end subroutine hdf5_write_data_2d_real8

  subroutine hdf5_write_data_3d_real8(ifile, dsetfull, darray)
    integer(hid_t), intent(in)     :: ifile
    character(len=*), intent(in)   :: dsetfull
    real(8), dimension(:,:,:)      :: darray

    integer(hid_t)                 :: grp_id, dset_id, dspace_id
    integer(hsize_t), dimension(3) :: dset_dims, dset_maxdims
    integer, parameter             :: rank = 3
    integer(8), dimension(3)       :: dims

    character(len=150) :: gname, dset

    dims(1) = size(darray,1)
    dims(2) = size(darray,2)
    dims(3) = size(darray,3)

    call hdf5_help_separate_dsetname(dsetfull, gname, dset)
    call hdf5_create_group(ifile, gname, iforce_error=.false.)

    if (trim(adjustl(gname)) == '') then
      grp_id = ifile
    else
      call h5gopen_f(ifile, trim(adjustl(gname)), grp_id, hdf_err)
    endif
    call h5screate_f(h5s_simple_f, dspace_id, hdf_err)
    call h5sset_extent_simple_f(dspace_id, rank, dims, dims, hdf_err)
    call h5dcreate_f(grp_id, trim(adjustl(dset)), h5t_native_double, dspace_id, dset_id, hdf_err)
    call h5dwrite_f(dset_id, h5t_native_double, darray, dims, hdf_err)
    call h5dclose_f(dset_id, hdf_err)
    call h5sclose_f(dspace_id, hdf_err)
    if (grp_id .ne. ifile) then
      call h5gclose_f(grp_id, hdf_err)
    endif
  end subroutine hdf5_write_data_3d_real8

  subroutine hdf5_write_data_4d_real8(ifile, dsetfull, darray)
    integer(hid_t), intent(in)     :: ifile
    character(len=*), intent(in)   :: dsetfull
    real(8), dimension(:,:,:,:)    :: darray

    integer(hid_t)                 :: grp_id, dset_id, dspace_id
    integer(hsize_t), dimension(4) :: dset_dims, dset_maxdims
    integer, parameter             :: rank = 4
    integer(8), dimension(4)       :: dims

    character(len=150) :: gname, dset

    dims(1) = size(darray,1)
    dims(2) = size(darray,2)
    dims(3) = size(darray,3)
    dims(4) = size(darray,4)

    call hdf5_help_separate_dsetname(dsetfull, gname, dset)
    call hdf5_create_group(ifile, gname, iforce_error=.false.)

    if (trim(adjustl(gname)) == '') then
      grp_id = ifile
    else
      call h5gopen_f(ifile, trim(adjustl(gname)), grp_id, hdf_err)
    endif
    call h5screate_f(h5s_simple_f, dspace_id, hdf_err)
    call h5sset_extent_simple_f(dspace_id, rank, dims, dims, hdf_err)
    call h5dcreate_f(grp_id, trim(adjustl(dset)), h5t_native_double, dspace_id, dset_id, hdf_err)
    call h5dwrite_f(dset_id, h5t_native_double, darray, dims, hdf_err)
    call h5dclose_f(dset_id, hdf_err)
    call h5sclose_f(dspace_id, hdf_err)
    if (grp_id .ne. ifile) then
      call h5gclose_f(grp_id, hdf_err)
    endif
  end subroutine hdf5_write_data_4d_real8

  subroutine hdf5_write_data_5d_real8(ifile, dsetfull, darray)
    integer(hid_t), intent(in)     :: ifile
    character(len=*), intent(in)   :: dsetfull
    real(8), dimension(:,:,:,:,:)  :: darray

    integer(hid_t)                 :: grp_id, dset_id, dspace_id
    integer(hsize_t), dimension(5) :: dset_dims, dset_maxdims
    integer, parameter             :: rank = 5
    integer(8), dimension(5)       :: dims

    character(len=150) :: gname, dset

    dims(1) = size(darray,1)
    dims(2) = size(darray,2)
    dims(3) = size(darray,3)
    dims(4) = size(darray,4)
    dims(5) = size(darray,5)

    call hdf5_help_separate_dsetname(dsetfull, gname, dset)
    call hdf5_create_group(ifile, gname, iforce_error=.false.)

    if (trim(adjustl(gname)) == '') then
      grp_id = ifile
    else
      call h5gopen_f(ifile, trim(adjustl(gname)), grp_id, hdf_err)
    endif
    call h5screate_f(h5s_simple_f, dspace_id, hdf_err)
    call h5sset_extent_simple_f(dspace_id, rank, dims, dims, hdf_err)
    call h5dcreate_f(grp_id, trim(adjustl(dset)), h5t_native_double, dspace_id, dset_id, hdf_err)
    call h5dwrite_f(dset_id, h5t_native_double, darray, dims, hdf_err)
    call h5dclose_f(dset_id, hdf_err)
    call h5sclose_f(dspace_id, hdf_err)
    if (grp_id .ne. ifile) then
      call h5gclose_f(grp_id, hdf_err)
    endif
  end subroutine hdf5_write_data_5d_real8

  subroutine hdf5_write_data_6d_real8(ifile, dsetfull, darray)
    integer(hid_t), intent(in)      :: ifile
    character(len=*), intent(in)    :: dsetfull
    real(8), dimension(:,:,:,:,:,:) :: darray

    integer(hid_t)                  :: grp_id, dset_id, dspace_id
    integer(hsize_t), dimension(6)  :: dset_dims, dset_maxdims
    integer, parameter              :: rank = 6
    integer(8), dimension(6)        :: dims

    character(len=150) :: gname, dset

    dims(1) = size(darray,1)
    dims(2) = size(darray,2)
    dims(3) = size(darray,3)
    dims(4) = size(darray,4)
    dims(5) = size(darray,5)
    dims(6) = size(darray,6)

    call hdf5_help_separate_dsetname(dsetfull, gname, dset)
    call hdf5_create_group(ifile, gname, iforce_error=.false.)

    if (trim(adjustl(gname)) == '') then
      grp_id = ifile
    else
      call h5gopen_f(ifile, trim(adjustl(gname)), grp_id, hdf_err)
    endif
    call h5screate_f(h5s_simple_f, dspace_id, hdf_err)
    call h5sset_extent_simple_f(dspace_id, rank, dims, dims, hdf_err)
    call h5dcreate_f(grp_id, trim(adjustl(dset)), h5t_native_double, dspace_id, dset_id, hdf_err)
    call h5dwrite_f(dset_id, h5t_native_double, darray, dims, hdf_err)
    call h5dclose_f(dset_id, hdf_err)
    call h5sclose_f(dspace_id, hdf_err)
    if (grp_id .ne. ifile) then
      call h5gclose_f(grp_id, hdf_err)
    endif
  end subroutine hdf5_write_data_6d_real8

  subroutine hdf5_write_data_7d_real8(ifile, dsetfull, darray)
    integer(hid_t), intent(in)        :: ifile
    character(len=*), intent(in)      :: dsetfull
    real(8), dimension(:,:,:,:,:,:,:) :: darray

    integer(hid_t)                    :: grp_id, dset_id, dspace_id
    integer(hsize_t), dimension(7)    :: dset_dims, dset_maxdims
    integer, parameter                :: rank = 7
    integer(8), dimension(7)          :: dims

    character(len=150) :: gname, dset

    dims(1) = size(darray,1)
    dims(2) = size(darray,2)
    dims(3) = size(darray,3)
    dims(4) = size(darray,4)
    dims(5) = size(darray,5)
    dims(6) = size(darray,6)
    dims(7) = size(darray,7)

    call hdf5_help_separate_dsetname(dsetfull, gname, dset)
    call hdf5_create_group(ifile, gname, iforce_error=.false.)

    if (trim(adjustl(gname)) == '') then
      grp_id = ifile
    else
      call h5gopen_f(ifile, trim(adjustl(gname)), grp_id, hdf_err)
    endif
    call h5screate_f(h5s_simple_f, dspace_id, hdf_err)
    call h5sset_extent_simple_f(dspace_id, rank, dims, dims, hdf_err)
    call h5dcreate_f(grp_id, trim(adjustl(dset)), h5t_native_double, dspace_id, dset_id, hdf_err)
    call h5dwrite_f(dset_id, h5t_native_double, darray, dims, hdf_err)
    call h5dclose_f(dset_id, hdf_err)
    call h5sclose_f(dspace_id, hdf_err)
    if (grp_id .ne. ifile) then
      call h5gclose_f(grp_id, hdf_err)
    endif
  end subroutine hdf5_write_data_7d_real8

  subroutine hdf5_write_data_0d_complex4(ifile, dsetfull, darray)
    integer(hid_t), intent(in)   :: ifile
    character(len=*), intent(in) :: dsetfull
    complex(4)                   :: darray

    integer(hid_t)               :: grp_id, dset_id, dspace_id
    integer, parameter           :: rank = 0
    integer(8), dimension(0)     :: dims

    character(len=150) :: gname, dset

    call hdf5_help_separate_dsetname(dsetfull, gname, dset)
    call hdf5_create_group(ifile, gname, iforce_error=.false.)

    if (trim(adjustl(gname)) == '') then
      grp_id = ifile
    else
      call h5gopen_f(ifile, trim(adjustl(gname)), grp_id, hdf_err)
    endif
    call h5screate_simple_f(0, dims, dspace_id, hdf_err)
    call h5dcreate_f(grp_id, trim(adjustl(dset)), complex_id_sp, dspace_id, dset_id, hdf_err)
    call h5dwrite_f(dset_id, complex_id_r_sp, real(darray), dims, hdf_err)
    call h5dwrite_f(dset_id, complex_id_i_sp, aimag(darray), dims, hdf_err)
    call h5dclose_f(dset_id, hdf_err)
    call h5sclose_f(dspace_id, hdf_err)
    if (grp_id .ne. ifile) then
      call h5gclose_f(grp_id, hdf_err)
    endif
  end subroutine hdf5_write_data_0d_complex4

  subroutine hdf5_write_data_1d_complex4(ifile, dsetfull, darray)
    integer(hid_t), intent(in)   :: ifile
    character(len=*), intent(in) :: dsetfull
    complex(4), dimension(:)     :: darray

    integer(hid_t)               :: grp_id, dset_id, dspace_id
    integer, parameter           :: rank = 1
    integer(8), dimension(1)     :: dims

    character(len=150) :: gname, dset

    dims(1) = size(darray,1)

    call hdf5_help_separate_dsetname(dsetfull, gname, dset)
    call hdf5_create_group(ifile, gname, iforce_error=.false.)

    if (trim(adjustl(gname)) == '') then
      grp_id = ifile
    else
      call h5gopen_f(ifile, trim(adjustl(gname)), grp_id, hdf_err)
    endif
    call h5screate_f(h5s_simple_f, dspace_id, hdf_err)
    call h5sset_extent_simple_f(dspace_id, rank, dims, dims, hdf_err)
    call h5dcreate_f(grp_id, trim(adjustl(dset)), complex_id_sp, dspace_id, dset_id, hdf_err)
    call h5dwrite_f(dset_id, complex_id_r_sp, real(darray), dims, hdf_err)
    call h5dwrite_f(dset_id, complex_id_i_sp, aimag(darray), dims, hdf_err)
    call h5dclose_f(dset_id, hdf_err)
    call h5sclose_f(dspace_id, hdf_err)
    if (grp_id .ne. ifile) then
      call h5gclose_f(grp_id, hdf_err)
    endif
  end subroutine hdf5_write_data_1d_complex4

  subroutine hdf5_write_data_2d_complex4(ifile, dsetfull, darray)
    integer(hid_t), intent(in)   :: ifile
    character(len=*), intent(in) :: dsetfull
    complex(4), dimension(:,:)   :: darray

    integer(hid_t)               :: grp_id, dset_id, dspace_id
    integer, parameter           :: rank = 2
    integer(8), dimension(2)     :: dims

    character(len=150) :: gname, dset

    dims(1) = size(darray,1)
    dims(2) = size(darray,2)

    call hdf5_help_separate_dsetname(dsetfull, gname, dset)
    call hdf5_create_group(ifile, gname, iforce_error=.false.)

    if (trim(adjustl(gname)) == '') then
      grp_id = ifile
    else
      call h5gopen_f(ifile, trim(adjustl(gname)), grp_id, hdf_err)
    endif
    call h5screate_f(h5s_simple_f, dspace_id, hdf_err)
    call h5sset_extent_simple_f(dspace_id, rank, dims, dims, hdf_err)
    call h5dcreate_f(grp_id, trim(adjustl(dset)), complex_id_sp, dspace_id, dset_id, hdf_err)
    call h5dwrite_f(dset_id, complex_id_r_sp, real(darray), dims, hdf_err)
    call h5dwrite_f(dset_id, complex_id_i_sp, aimag(darray), dims, hdf_err)
    call h5dclose_f(dset_id, hdf_err)
    call h5sclose_f(dspace_id, hdf_err)
    if (grp_id .ne. ifile) then
      call h5gclose_f(grp_id, hdf_err)
    endif
  end subroutine hdf5_write_data_2d_complex4

  subroutine hdf5_write_data_3d_complex4(ifile, dsetfull, darray)
    integer(hid_t), intent(in)   :: ifile
    character(len=*), intent(in) :: dsetfull
    complex(4), dimension(:,:,:) :: darray

    integer(hid_t)               :: grp_id, dset_id, dspace_id
    integer, parameter           :: rank = 3
    integer(8), dimension(3)     :: dims

    character(len=150) :: gname, dset

    dims(1) = size(darray,1)
    dims(2) = size(darray,2)
    dims(3) = size(darray,3)

    call hdf5_help_separate_dsetname(dsetfull, gname, dset)
    call hdf5_create_group(ifile, gname, iforce_error=.false.)

    if (trim(adjustl(gname)) == '') then
      grp_id = ifile
    else
      call h5gopen_f(ifile, trim(adjustl(gname)), grp_id, hdf_err)
    endif
    call h5screate_f(h5s_simple_f, dspace_id, hdf_err)
    call h5sset_extent_simple_f(dspace_id, rank, dims, dims, hdf_err)
    call h5dcreate_f(grp_id, trim(adjustl(dset)), complex_id_sp, dspace_id, dset_id, hdf_err)
    call h5dwrite_f(dset_id, complex_id_r_sp, real(darray), dims, hdf_err)
    call h5dwrite_f(dset_id, complex_id_i_sp, aimag(darray), dims, hdf_err)
    call h5dclose_f(dset_id, hdf_err)
    call h5sclose_f(dspace_id, hdf_err)
    if (grp_id .ne. ifile) then
      call h5gclose_f(grp_id, hdf_err)
    endif
  end subroutine hdf5_write_data_3d_complex4

  subroutine hdf5_write_data_4d_complex4(ifile, dsetfull, darray)
    integer(hid_t), intent(in)     :: ifile
    character(len=*), intent(in)   :: dsetfull
    complex(4), dimension(:,:,:,:) :: darray

    integer(hid_t)                 :: grp_id, dset_id, dspace_id
    integer, parameter             :: rank = 4
    integer(8), dimension(4)       :: dims

    character(len=150) :: gname, dset

    dims(1) = size(darray,1)
    dims(2) = size(darray,2)
    dims(3) = size(darray,3)
    dims(4) = size(darray,4)

    call hdf5_help_separate_dsetname(dsetfull, gname, dset)
    call hdf5_create_group(ifile, gname, iforce_error=.false.)

    if (trim(adjustl(gname)) == '') then
      grp_id = ifile
    else
      call h5gopen_f(ifile, trim(adjustl(gname)), grp_id, hdf_err)
    endif
    call h5screate_f(h5s_simple_f, dspace_id, hdf_err)
    call h5sset_extent_simple_f(dspace_id, rank, dims, dims, hdf_err)
    call h5dcreate_f(grp_id, trim(adjustl(dset)), complex_id_sp, dspace_id, dset_id, hdf_err)
    call h5dwrite_f(dset_id, complex_id_r_sp, real(darray), dims, hdf_err)
    call h5dwrite_f(dset_id, complex_id_i_sp, aimag(darray), dims, hdf_err)
    call h5dclose_f(dset_id, hdf_err)
    call h5sclose_f(dspace_id, hdf_err)
    if (grp_id .ne. ifile) then
      call h5gclose_f(grp_id, hdf_err)
    endif
  end subroutine hdf5_write_data_4d_complex4

  subroutine hdf5_write_data_5d_complex4(ifile, dsetfull, darray)
    integer(hid_t), intent(in)       :: ifile
    character(len=*), intent(in)     :: dsetfull
    complex(4), dimension(:,:,:,:,:) :: darray

    integer(hid_t)                   :: grp_id, dset_id, dspace_id
    integer, parameter               :: rank = 5
    integer(8), dimension(5)         :: dims

    character(len=150) :: gname, dset

    dims(1) = size(darray,1)
    dims(2) = size(darray,2)
    dims(3) = size(darray,3)
    dims(4) = size(darray,4)
    dims(5) = size(darray,5)

    call hdf5_help_separate_dsetname(dsetfull, gname, dset)
    call hdf5_create_group(ifile, gname, iforce_error=.false.)

    if (trim(adjustl(gname)) == '') then
      grp_id = ifile
    else
      call h5gopen_f(ifile, trim(adjustl(gname)), grp_id, hdf_err)
    endif
    call h5screate_f(h5s_simple_f, dspace_id, hdf_err)
    call h5sset_extent_simple_f(dspace_id, rank, dims, dims, hdf_err)
    call h5dcreate_f(grp_id, trim(adjustl(dset)), complex_id_sp, dspace_id, dset_id, hdf_err)
    call h5dwrite_f(dset_id, complex_id_r_sp, real(darray), dims, hdf_err)
    call h5dwrite_f(dset_id, complex_id_i_sp, aimag(darray), dims, hdf_err)
    call h5dclose_f(dset_id, hdf_err)
    call h5sclose_f(dspace_id, hdf_err)
    if (grp_id .ne. ifile) then
      call h5gclose_f(grp_id, hdf_err)
    endif
  end subroutine hdf5_write_data_5d_complex4

  subroutine hdf5_write_data_6d_complex4(ifile, dsetfull, darray)
    integer(hid_t), intent(in)         :: ifile
    character(len=*), intent(in)       :: dsetfull
    complex(4), dimension(:,:,:,:,:,:) :: darray

    integer(hid_t)                     :: grp_id, dset_id, dspace_id
    integer, parameter                 :: rank = 6
    integer(8), dimension(6)           :: dims

    character(len=150) :: gname, dset

    dims(1) = size(darray,1)
    dims(2) = size(darray,2)
    dims(3) = size(darray,3)
    dims(4) = size(darray,4)
    dims(5) = size(darray,5)
    dims(6) = size(darray,6)

    call hdf5_help_separate_dsetname(dsetfull, gname, dset)
    call hdf5_create_group(ifile, gname, iforce_error=.false.)

    if (trim(adjustl(gname)) == '') then
      grp_id = ifile
    else
      call h5gopen_f(ifile, trim(adjustl(gname)), grp_id, hdf_err)
    endif
    call h5screate_f(h5s_simple_f, dspace_id, hdf_err)
    call h5sset_extent_simple_f(dspace_id, rank, dims, dims, hdf_err)
    call h5dcreate_f(grp_id, trim(adjustl(dset)), complex_id_sp, dspace_id, dset_id, hdf_err)
    call h5dwrite_f(dset_id, complex_id_r_sp, real(darray), dims, hdf_err)
    call h5dwrite_f(dset_id, complex_id_i_sp, aimag(darray), dims, hdf_err)
    call h5dclose_f(dset_id, hdf_err)
    call h5sclose_f(dspace_id, hdf_err)
    if (grp_id .ne. ifile) then
      call h5gclose_f(grp_id, hdf_err)
    endif
  end subroutine hdf5_write_data_6d_complex4

  subroutine hdf5_write_data_7d_complex4(ifile, dsetfull, darray)
    integer(hid_t), intent(in)           :: ifile
    character(len=*), intent(in)         :: dsetfull
    complex(4), dimension(:,:,:,:,:,:,:) :: darray

    integer(hid_t)                       :: grp_id, dset_id, dspace_id
    integer, parameter                   :: rank = 7
    integer(8), dimension(7)             :: dims

    character(len=150) :: gname, dset

    dims(1) = size(darray,1)
    dims(2) = size(darray,2)
    dims(3) = size(darray,3)
    dims(4) = size(darray,4)
    dims(5) = size(darray,5)
    dims(6) = size(darray,6)
    dims(7) = size(darray,7)

    call hdf5_help_separate_dsetname(dsetfull, gname, dset)
    call hdf5_create_group(ifile, gname, iforce_error=.false.)

    if (trim(adjustl(gname)) == '') then
      grp_id = ifile
    else
      call h5gopen_f(ifile, trim(adjustl(gname)), grp_id, hdf_err)
    endif
    call h5screate_f(h5s_simple_f, dspace_id, hdf_err)
    call h5sset_extent_simple_f(dspace_id, rank, dims, dims, hdf_err)
    call h5dcreate_f(grp_id, trim(adjustl(dset)), complex_id_sp, dspace_id, dset_id, hdf_err)
    call h5dwrite_f(dset_id, complex_id_r_sp, real(darray), dims, hdf_err)
    call h5dwrite_f(dset_id, complex_id_i_sp, aimag(darray), dims, hdf_err)
    call h5dclose_f(dset_id, hdf_err)
    call h5sclose_f(dspace_id, hdf_err)
    if (grp_id .ne. ifile) then
      call h5gclose_f(grp_id, hdf_err)
    endif
  end subroutine hdf5_write_data_7d_complex4

  subroutine hdf5_write_data_0d_complex8(ifile, dsetfull, darray)
    integer(hid_t), intent(in)   :: ifile
    character(len=*), intent(in) :: dsetfull
    complex(8)                   :: darray

    integer(hid_t)               :: grp_id, dset_id, dspace_id
    integer, parameter           :: rank = 0
    integer(8), dimension(0)     :: dims

    character(len=150) :: gname, dset

    call hdf5_help_separate_dsetname(dsetfull, gname, dset)
    call hdf5_create_group(ifile, gname, iforce_error=.false.)

    if (trim(adjustl(gname)) == '') then
      grp_id = ifile
    else
      call h5gopen_f(ifile, trim(adjustl(gname)), grp_id, hdf_err)
    endif
    call h5screate_simple_f(0, dims, dspace_id, hdf_err)
    call h5dcreate_f(grp_id, trim(adjustl(dset)), complex_id_dp, dspace_id, dset_id, hdf_err)
    call h5dwrite_f(dset_id, complex_id_r_dp, real(darray), dims, hdf_err)
    call h5dwrite_f(dset_id, complex_id_i_dp, aimag(darray), dims, hdf_err)
    call h5dclose_f(dset_id, hdf_err)
    call h5sclose_f(dspace_id, hdf_err)
    if (grp_id .ne. ifile) then
      call h5gclose_f(grp_id, hdf_err)
    endif
  end subroutine hdf5_write_data_0d_complex8

  subroutine hdf5_write_data_1d_complex8(ifile, dsetfull, darray)
    integer(hid_t), intent(in)   :: ifile
    character(len=*), intent(in) :: dsetfull
    complex(8), dimension(:)     :: darray

    integer(hid_t)               :: grp_id, dset_id, dspace_id
    integer, parameter           :: rank = 1
    integer(8), dimension(1)     :: dims

    character(len=150) :: gname, dset

    dims(1) = size(darray,1)

    call hdf5_help_separate_dsetname(dsetfull, gname, dset)
    call hdf5_create_group(ifile, gname, iforce_error=.false.)

    if (trim(adjustl(gname)) == '') then
      grp_id = ifile
    else
      call h5gopen_f(ifile, trim(adjustl(gname)), grp_id, hdf_err)
    endif
    call h5screate_f(h5s_simple_f, dspace_id, hdf_err)
    call h5sset_extent_simple_f(dspace_id, rank, dims, dims, hdf_err)
    call h5dcreate_f(grp_id, trim(adjustl(dset)), complex_id_dp, dspace_id, dset_id, hdf_err)
    call h5dwrite_f(dset_id, complex_id_r_dp, real(darray), dims, hdf_err)
    call h5dwrite_f(dset_id, complex_id_i_dp, aimag(darray), dims, hdf_err)
    call h5dclose_f(dset_id, hdf_err)
    call h5sclose_f(dspace_id, hdf_err)
    if (grp_id .ne. ifile) then
      call h5gclose_f(grp_id, hdf_err)
    endif
  end subroutine hdf5_write_data_1d_complex8

  subroutine hdf5_write_data_2d_complex8(ifile, dsetfull, darray)
    integer(hid_t), intent(in)   :: ifile
    character(len=*), intent(in) :: dsetfull
    complex(8), dimension(:,:)   :: darray

    integer(hid_t)               :: grp_id, dset_id, dspace_id
    integer, parameter           :: rank = 2
    integer(8), dimension(2)     :: dims

    character(len=150) :: gname, dset

    dims(1) = size(darray,1)
    dims(2) = size(darray,2)

    call hdf5_help_separate_dsetname(dsetfull, gname, dset)
    call hdf5_create_group(ifile, gname, iforce_error=.false.)

    if (trim(adjustl(gname)) == '') then
      grp_id = ifile
    else
      call h5gopen_f(ifile, trim(adjustl(gname)), grp_id, hdf_err)
    endif
    call h5screate_f(h5s_simple_f, dspace_id, hdf_err)
    call h5sset_extent_simple_f(dspace_id, rank, dims, dims, hdf_err)
    call h5dcreate_f(grp_id, trim(adjustl(dset)), complex_id_dp, dspace_id, dset_id, hdf_err)
    call h5dwrite_f(dset_id, complex_id_r_dp, real(darray), dims, hdf_err)
    call h5dwrite_f(dset_id, complex_id_i_dp, aimag(darray), dims, hdf_err)
    call h5dclose_f(dset_id, hdf_err)
    call h5sclose_f(dspace_id, hdf_err)
    if (grp_id .ne. ifile) then
      call h5gclose_f(grp_id, hdf_err)
    endif
  end subroutine hdf5_write_data_2d_complex8

  subroutine hdf5_write_data_3d_complex8(ifile, dsetfull, darray)
    integer(hid_t), intent(in)   :: ifile
    character(len=*), intent(in) :: dsetfull
    complex(8), dimension(:,:,:) :: darray

    integer(hid_t)               :: grp_id, dset_id, dspace_id
    integer, parameter           :: rank = 3
    integer(8), dimension(3)     :: dims

    character(len=150) :: gname, dset

    dims(1) = size(darray,1)
    dims(2) = size(darray,2)
    dims(3) = size(darray,3)

    call hdf5_help_separate_dsetname(dsetfull, gname, dset)
    call hdf5_create_group(ifile, gname, iforce_error=.false.)

    if (trim(adjustl(gname)) == '') then
      grp_id = ifile
    else
      call h5gopen_f(ifile, trim(adjustl(gname)), grp_id, hdf_err)
    endif
    call h5screate_f(h5s_simple_f, dspace_id, hdf_err)
    call h5sset_extent_simple_f(dspace_id, rank, dims, dims, hdf_err)
    call h5dcreate_f(grp_id, trim(adjustl(dset)), complex_id_dp, dspace_id, dset_id, hdf_err)
    call h5dwrite_f(dset_id, complex_id_r_dp, real(darray), dims, hdf_err)
    call h5dwrite_f(dset_id, complex_id_i_dp, aimag(darray), dims, hdf_err)
    call h5dclose_f(dset_id, hdf_err)
    call h5sclose_f(dspace_id, hdf_err)
    if (grp_id .ne. ifile) then
      call h5gclose_f(grp_id, hdf_err)
    endif
  end subroutine hdf5_write_data_3d_complex8

  subroutine hdf5_write_data_4d_complex8(ifile, dsetfull, darray)
    integer(hid_t), intent(in)     :: ifile
    character(len=*), intent(in)   :: dsetfull
    complex(8), dimension(:,:,:,:) :: darray

    integer(hid_t)                 :: grp_id, dset_id, dspace_id
    integer, parameter             :: rank = 4
    integer(8), dimension(4)       :: dims

    character(len=150) :: gname, dset

    dims(1) = size(darray,1)
    dims(2) = size(darray,2)
    dims(3) = size(darray,3)
    dims(4) = size(darray,4)

    call hdf5_help_separate_dsetname(dsetfull, gname, dset)
    call hdf5_create_group(ifile, gname, iforce_error=.false.)

    if (trim(adjustl(gname)) == '') then
      grp_id = ifile
    else
      call h5gopen_f(ifile, trim(adjustl(gname)), grp_id, hdf_err)
    endif
    call h5screate_f(h5s_simple_f, dspace_id, hdf_err)
    call h5sset_extent_simple_f(dspace_id, rank, dims, dims, hdf_err)
    call h5dcreate_f(grp_id, trim(adjustl(dset)), complex_id_dp, dspace_id, dset_id, hdf_err)
    call h5dwrite_f(dset_id, complex_id_r_dp, real(darray), dims, hdf_err)
    call h5dwrite_f(dset_id, complex_id_i_dp, aimag(darray), dims, hdf_err)
    call h5dclose_f(dset_id, hdf_err)
    call h5sclose_f(dspace_id, hdf_err)
    if (grp_id .ne. ifile) then
      call h5gclose_f(grp_id, hdf_err)
    endif
  end subroutine hdf5_write_data_4d_complex8

  subroutine hdf5_write_data_5d_complex8(ifile, dsetfull, darray)
    integer(hid_t), intent(in)       :: ifile
    character(len=*), intent(in)     :: dsetfull
    complex(8), dimension(:,:,:,:,:) :: darray

    integer(hid_t)                   :: grp_id, dset_id, dspace_id
    integer, parameter               :: rank = 5
    integer(8), dimension(5)         :: dims

    character(len=150) :: gname, dset

    dims(1) = size(darray,1)
    dims(2) = size(darray,2)
    dims(3) = size(darray,3)
    dims(4) = size(darray,4)
    dims(5) = size(darray,5)

    call hdf5_help_separate_dsetname(dsetfull, gname, dset)
    call hdf5_create_group(ifile, gname, iforce_error=.false.)

    if (trim(adjustl(gname)) == '') then
      grp_id = ifile
    else
      call h5gopen_f(ifile, trim(adjustl(gname)), grp_id, hdf_err)
    endif
    call h5screate_f(h5s_simple_f, dspace_id, hdf_err)
    call h5sset_extent_simple_f(dspace_id, rank, dims, dims, hdf_err)
    call h5dcreate_f(grp_id, trim(adjustl(dset)), complex_id_dp, dspace_id, dset_id, hdf_err)
    call h5dwrite_f(dset_id, complex_id_r_dp, real(darray), dims, hdf_err)
    call h5dwrite_f(dset_id, complex_id_i_dp, aimag(darray), dims, hdf_err)
    call h5dclose_f(dset_id, hdf_err)
    call h5sclose_f(dspace_id, hdf_err)
    if (grp_id .ne. ifile) then
      call h5gclose_f(grp_id, hdf_err)
    endif
  end subroutine hdf5_write_data_5d_complex8

  subroutine hdf5_write_data_6d_complex8(ifile, dsetfull, darray)
    integer(hid_t), intent(in)         :: ifile
    character(len=*), intent(in)       :: dsetfull
    complex(8), dimension(:,:,:,:,:,:) :: darray

    integer(hid_t)                     :: grp_id, dset_id, dspace_id
    integer, parameter                 :: rank = 6
    integer(8), dimension(6)           :: dims

    character(len=150) :: gname, dset

    dims(1) = size(darray,1)
    dims(2) = size(darray,2)
    dims(3) = size(darray,3)
    dims(4) = size(darray,4)
    dims(5) = size(darray,5)
    dims(6) = size(darray,6)

    call hdf5_help_separate_dsetname(dsetfull, gname, dset)
    call hdf5_create_group(ifile, gname, iforce_error=.false.)

    if (trim(adjustl(gname)) == '') then
      grp_id = ifile
    else
      call h5gopen_f(ifile, trim(adjustl(gname)), grp_id, hdf_err)
    endif
    call h5screate_f(h5s_simple_f, dspace_id, hdf_err)
    call h5sset_extent_simple_f(dspace_id, rank, dims, dims, hdf_err)
    call h5dcreate_f(grp_id, trim(adjustl(dset)), complex_id_dp, dspace_id, dset_id, hdf_err)
    call h5dwrite_f(dset_id, complex_id_r_dp, real(darray), dims, hdf_err)
    call h5dwrite_f(dset_id, complex_id_i_dp, aimag(darray), dims, hdf_err)
    call h5dclose_f(dset_id, hdf_err)
    call h5sclose_f(dspace_id, hdf_err)
    if (grp_id .ne. ifile) then
      call h5gclose_f(grp_id, hdf_err)
    endif
  end subroutine hdf5_write_data_6d_complex8

  subroutine hdf5_write_data_7d_complex8(ifile, dsetfull, darray)
    integer(hid_t), intent(in)           :: ifile
    character(len=*), intent(in)         :: dsetfull
    complex(8), dimension(:,:,:,:,:,:,:) :: darray

    integer(hid_t)                       :: grp_id, dset_id, dspace_id
    integer, parameter                   :: rank = 7
    integer(8), dimension(7)             :: dims

    character(len=150) :: gname, dset

    dims(1) = size(darray,1)
    dims(2) = size(darray,2)
    dims(3) = size(darray,3)
    dims(4) = size(darray,4)
    dims(5) = size(darray,5)
    dims(6) = size(darray,6)
    dims(7) = size(darray,7)

    call hdf5_help_separate_dsetname(dsetfull, gname, dset)
    call hdf5_create_group(ifile, gname, iforce_error=.false.)

    if (trim(adjustl(gname)) == '') then
      grp_id = ifile
    else
      call h5gopen_f(ifile, trim(adjustl(gname)), grp_id, hdf_err)
    endif
    call h5screate_f(h5s_simple_f, dspace_id, hdf_err)
    call h5sset_extent_simple_f(dspace_id, rank, dims, dims, hdf_err)
    call h5dcreate_f(grp_id, trim(adjustl(dset)), complex_id_dp, dspace_id, dset_id, hdf_err)
    call h5dwrite_f(dset_id, complex_id_r_dp, real(darray), dims, hdf_err)
    call h5dwrite_f(dset_id, complex_id_i_dp, aimag(darray), dims, hdf_err)
    call h5dclose_f(dset_id, hdf_err)
    call h5sclose_f(dspace_id, hdf_err)
    if (grp_id .ne. ifile) then
      call h5gclose_f(grp_id, hdf_err)
    endif
  end subroutine hdf5_write_data_7d_complex8

  subroutine hdf5_write_attribute_logical(ifile, location, attrname, attr)
    integer(hid_t), intent(in)   :: ifile
    character(len=*), intent(in) :: location
    character(len=*), intent(in) :: attrname
    logical                      :: attr

    integer :: attri
    if (attr) then
      attri = 1
    else
      attri = 0
    endif
    call hdf5_write_attribute_int4(ifile, location, attrname, attri)
  end subroutine

  subroutine hdf5_write_attribute_int4(ifile, location, attrname, attr)
    integer(hid_t), intent(in)   :: ifile
    character(len=*), intent(in) :: location
    character(len=*), intent(in) :: attrname
    integer(4)                   :: attr

    integer(8), dimension(0)     :: dims
    integer(hid_t) :: obj_id, dspace_id, attr_id

    call h5oopen_f(ifile, trim(adjustl(location)), obj_id, hdf_err)
    call h5screate_simple_f(0, dims, dspace_id, hdf_err)
    call h5acreate_f(obj_id, trim(adjustl(attrname)), h5t_native_integer, dspace_id, attr_id, hdf_err)
    call h5awrite_f(attr_id, h5t_native_integer, attr, dims, hdf_err)
    call h5aclose_f(attr_id, hdf_err)
    call h5sclose_f(dspace_id, hdf_err)
    call h5oclose_f(obj_id, hdf_err)
  end subroutine

  subroutine hdf5_write_attribute_real4(ifile, location, attrname, attr)
    integer(hid_t), intent(in)   :: ifile
    character(len=*), intent(in) :: location
    character(len=*), intent(in) :: attrname
    real(4)                      :: attr

    integer(8), dimension(0)     :: dims
    integer(hid_t) :: obj_id, dspace_id, attr_id

    call h5oopen_f(ifile, trim(adjustl(location)), obj_id, hdf_err)
    call h5screate_simple_f(0, dims, dspace_id, hdf_err)
    call h5acreate_f(obj_id, trim(adjustl(attrname)), h5t_native_real, dspace_id, attr_id, hdf_err)
    call h5awrite_f(attr_id, h5t_native_real, attr, dims, hdf_err)
    call h5aclose_f(attr_id, hdf_err)
    call h5sclose_f(dspace_id, hdf_err)
    call h5oclose_f(obj_id, hdf_err)
  end subroutine

  subroutine hdf5_write_attribute_real8(ifile, location, attrname, attr)
    integer(hid_t), intent(in)   :: ifile
    character(len=*), intent(in) :: location
    character(len=*), intent(in) :: attrname
    real(8)                      :: attr

    integer(8), dimension(0)     :: dims
    integer(hid_t) :: obj_id, dspace_id, attr_id

    call h5oopen_f(ifile, trim(adjustl(location)), obj_id, hdf_err)
    call h5screate_simple_f(0, dims, dspace_id, hdf_err)
    call h5acreate_f(obj_id, trim(adjustl(attrname)), h5t_native_double, dspace_id, attr_id, hdf_err)
    call h5awrite_f(attr_id, h5t_native_double, attr, dims, hdf_err)
    call h5aclose_f(attr_id, hdf_err)
    call h5sclose_f(dspace_id, hdf_err)
    call h5oclose_f(obj_id, hdf_err)
  end subroutine

  subroutine hdf5_write_attribute_string(ifile, location, attrname, attr)
    integer(hid_t), intent(in)   :: ifile
    character(len=*), intent(in) :: location
    character(len=*), intent(in) :: attrname
    character(len=*)             :: attr

    integer(8), dimension(0) :: dims
    integer(hid_t)           :: type_id, dspace_id, attr_id, obj_id
    character(len=150)       :: strim
    integer(8)               :: string_size

    strim = adjustl(trim(attr))
    string_size = int(len_trim(strim),8)

    call h5tcopy_f(h5t_native_character, type_id, hdf_err)
    call h5tset_size_f(type_id, string_size, hdf_err)
    call h5screate_f(h5s_scalar_f, dspace_id, hdf_err)

    call h5oopen_f(ifile, trim(adjustl(location)), obj_id, hdf_err)
    call h5acreate_f(obj_id, trim(adjustl(attrname)), type_id, dspace_id, attr_id, hdf_err)
    call h5awrite_f(attr_id, type_id, trim(attr), dims, hdf_err)
    call h5aclose_f(attr_id, hdf_err)
    call h5sclose_f(dspace_id, hdf_err)
    call h5oclose_f(obj_id, hdf_err)
    call h5tclose_f(type_id, hdf_err)

  end subroutine

  subroutine hdf5_write_attribute_complex4(ifile, location, attrname, attr)
    integer(hid_t), intent(in)   :: ifile
    character(len=*), intent(in) :: location
    character(len=*), intent(in) :: attrname
    complex(4)                   :: attr

    integer(8), dimension(0)     :: dims
    integer(hid_t) :: obj_id, dspace_id, attr_id

    call h5oopen_f(ifile, trim(adjustl(location)), obj_id, hdf_err)
    call h5screate_simple_f(0, dims, dspace_id, hdf_err)
    call h5acreate_f(obj_id, trim(adjustl(attrname)), complex_id_sp, dspace_id, attr_id, hdf_err)
    ! for whatever reason ... this works ... I really dont know why
    ! and it has to be in this order
    call h5awrite_f(attr_id, complex_id_sp, aimag(attr), dims, hdf_err)
    call h5awrite_f(attr_id, complex_id_sp, real(attr), dims, hdf_err)
    call h5aclose_f(attr_id, hdf_err)
    call h5sclose_f(dspace_id, hdf_err)
    call h5oclose_f(obj_id, hdf_err)
  end subroutine

  subroutine hdf5_write_attribute_complex8(ifile, location, attrname, attr)
    integer(hid_t), intent(in)   :: ifile
    character(len=*), intent(in) :: location
    character(len=*), intent(in) :: attrname
    complex(8)                   :: attr

    integer(8), dimension(0)     :: dims
    integer(hid_t) :: obj_id, dspace_id, attr_id

    call h5oopen_f(ifile, trim(adjustl(location)), obj_id, hdf_err)
    call h5screate_simple_f(0, dims, dspace_id, hdf_err)
    call h5acreate_f(obj_id, trim(adjustl(attrname)), complex_id_dp, dspace_id, attr_id, hdf_err)
    ! same story here
    call h5awrite_f(attr_id, complex_id_dp, aimag(attr), dims, hdf_err)
    call h5awrite_f(attr_id, complex_id_dp, real(attr), dims, hdf_err)
    call h5aclose_f(attr_id, hdf_err)
    call h5sclose_f(dspace_id, hdf_err)
    call h5oclose_f(obj_id, hdf_err)
  end subroutine

  subroutine hdf5_read_attribute_logical(ifile, location, attrname, attr)
    integer(hid_t), intent(in)   :: ifile
    character(len=*), intent(in) :: location
    character(len=*), intent(in) :: attrname
    logical, intent(out)         :: attr

    integer(4) :: attri

    call hdf5_read_attribute_int4(ifile, location, attrname, attri)
    if (attri == 0) then
      attr = .false.
    else
      attr = .true.
    endif

  end subroutine

  subroutine hdf5_read_attribute_int4(ifile, location, attrname, attr)
    integer(hid_t), intent(in)   :: ifile
    character(len=*), intent(in) :: location
    character(len=*), intent(in) :: attrname
    integer(4), intent(out)      :: attr

    integer(8), dimension(0)     :: dims
    integer(hid_t) :: obj_id, dspace_id, attr_id

    call h5oopen_f(ifile, trim(adjustl(location)), obj_id, hdf_err)
    call h5aopen_f(obj_id, trim(adjustl(attrname)), attr_id, hdf_err)
    call h5aread_f(attr_id, h5t_native_integer, attr, dims, hdf_err)
    call h5aclose_f(attr_id, hdf_err)
    call h5oclose_f(obj_id, hdf_err)
  end subroutine

  subroutine hdf5_read_attribute_real4(ifile, location, attrname, attr)
    integer(hid_t), intent(in)   :: ifile
    character(len=*), intent(in) :: location
    character(len=*), intent(in) :: attrname
    real(4), intent(out)      :: attr

    integer(8), dimension(0)     :: dims
    integer(hid_t) :: obj_id, dspace_id, attr_id

    call h5oopen_f(ifile, trim(adjustl(location)), obj_id, hdf_err)
    call h5aopen_f(obj_id, trim(adjustl(attrname)), attr_id, hdf_err)
    call h5aread_f(attr_id, h5t_native_real, attr, dims, hdf_err)
    call h5aclose_f(attr_id, hdf_err)
    call h5oclose_f(obj_id, hdf_err)
  end subroutine

  subroutine hdf5_read_attribute_real8(ifile, location, attrname, attr)
    integer(hid_t), intent(in)   :: ifile
    character(len=*), intent(in) :: location
    character(len=*), intent(in) :: attrname
    real(8), intent(out)      :: attr

    integer(8), dimension(0)     :: dims
    integer(hid_t) :: obj_id, dspace_id, attr_id

    call h5oopen_f(ifile, trim(adjustl(location)), obj_id, hdf_err)
    call h5aopen_f(obj_id, trim(adjustl(attrname)), attr_id, hdf_err)
    call h5aread_f(attr_id, h5t_native_double, attr, dims, hdf_err)
    call h5aclose_f(attr_id, hdf_err)
    call h5oclose_f(obj_id, hdf_err)
  end subroutine

  subroutine hdf5_read_attribute_string(ifile, location, attrname, attr)
    integer(hid_t), intent(in)      :: ifile
    character(len=*), intent(in)    :: location
    character(len=*), intent(in)    :: attrname
    character(len=*), intent(inout) :: attr

    integer(hsize_t) :: strg_size
    integer(8), dimension(0)     :: dims
    integer(hid_t) :: obj_id, dspace_id, attr_id, type_id

    attr = ''

    call h5oopen_f(ifile, trim(adjustl(location)), obj_id, hdf_err)
    call h5aopen_f(obj_id, trim(adjustl(attrname)), attr_id, hdf_err)

    ! create dspace
    call h5aget_storage_size_f(attr_id, strg_size, hdf_err)
    call h5tcopy_f(h5t_native_character, type_id, hdf_err)
    call h5tset_size_f(type_id, strg_size, hdf_err)
    call h5screate_f(h5s_scalar_f, dspace_id, hdf_err)

    call h5aread_f(attr_id, type_id, attr, dims, hdf_err)

    call h5aclose_f(attr_id, hdf_err)
    call h5sclose_f(dspace_id, hdf_err)
    call h5tclose_f(type_id, hdf_err)
    call h5oclose_f(obj_id, hdf_err)
  end subroutine

  subroutine hdf5_read_attribute_complex4(ifile, location, attrname, attr)
    integer(hid_t), intent(in)   :: ifile
    character(len=*), intent(in) :: location
    character(len=*), intent(in) :: attrname
    complex(4), intent(out)      :: attr

    complex(4) :: ci = (0.0, 1.0)
    integer(8), dimension(0)     :: dims
    integer(hid_t) :: obj_id, dspace_id, attr_id
    real(4) :: tmp_r, tmp_i

    call h5oopen_f(ifile, trim(adjustl(location)), obj_id, hdf_err)
    call h5aopen_f(obj_id, trim(adjustl(attrname)), attr_id, hdf_err)
    call h5aread_f(attr_id, complex_id_r_sp, tmp_r, dims, hdf_err)
    call h5aread_f(attr_id, complex_id_i_sp, tmp_i, dims, hdf_err)
    call h5aclose_f(attr_id, hdf_err)
    call h5oclose_f(obj_id, hdf_err)

    attr = tmp_r + ci*tmp_i
  end subroutine

  subroutine hdf5_read_attribute_complex8(ifile, location, attrname, attr)
    integer(hid_t), intent(in)   :: ifile
    character(len=*), intent(in) :: location
    character(len=*), intent(in) :: attrname
    complex(8), intent(out)      :: attr

    complex(8) :: ci = (0.d0, 1.d0)
    integer(8), dimension(0)     :: dims
    integer(hid_t) :: obj_id, dspace_id, attr_id
    real(8) :: tmp_r, tmp_i

    call h5oopen_f(ifile, trim(adjustl(location)), obj_id, hdf_err)
    call h5aopen_f(obj_id, trim(adjustl(attrname)), attr_id, hdf_err)
    call h5aread_f(attr_id, complex_id_r_dp, tmp_r, dims, hdf_err)
    call h5aread_f(attr_id, complex_id_i_dp, tmp_i, dims, hdf_err)
    call h5aclose_f(attr_id, hdf_err)
    call h5oclose_f(obj_id, hdf_err)

    attr = tmp_r + ci*tmp_i
  end subroutine

  integer function hdf5_get_number_groups(ifile, location)
    integer(hid_t), intent(in)                  :: ifile
    character(len=*), intent(in)                :: location

    integer :: nmembers
    integer :: ngroups

    integer :: idx
    character(len=150) :: strim
    integer(hid_t) :: obj_id
    integer :: obj_type

    ngroups = 0
    call h5gn_members_f(ifile, trim(adjustl(location)), nmembers, hdf_err)
    do idx = 0, nmembers-1
      call h5gget_obj_info_idx_f(ifile, trim(adjustl(location)), idx, strim, obj_type, hdf_err)
      if (obj_type .eq. h5g_group_f) then
        ngroups = ngroups + 1
      endif
    enddo

    hdf5_get_number_groups = ngroups
  end function

  integer function hdf5_get_number_datasets(ifile, location)
    integer(hid_t), intent(in)                  :: ifile
    character(len=*), intent(in)                :: location

    integer :: nmembers
    integer :: ndsets

    character(len=150) :: strim
    integer(hid_t) :: obj_id
    integer :: obj_type
    integer :: idx

    ndsets = 0
    call h5gn_members_f(ifile, trim(adjustl(location)), nmembers, hdf_err)
    do idx = 0,nmembers-1
      call h5gget_obj_info_idx_f(ifile, trim(adjustl(location)), idx, strim, obj_type, hdf_err)
      if (obj_type .eq. h5g_dataset_f) then
        ndsets = ndsets + 1
      endif
    enddo

    hdf5_get_number_datasets = ndsets
  end function

  integer function hdf5_get_number_attributes(ifile, location)
    integer(hid_t), intent(in)                  :: ifile
    character(len=*), intent(in)                :: location

    integer(hid_t) :: obj_id
    integer attr_num

    call h5oopen_f(ifile, trim(adjustl(location)), obj_id, hdf_err)
    call h5aget_num_attrs_f(obj_id, attr_num, hdf_err)
    call h5oclose_f(obj_id, hdf_err)

    hdf5_get_number_attributes = attr_num
  end function

  subroutine hdf5_list_groups(ifile, location, groups)
    integer(hid_t), intent(in)                  :: ifile
    character(len=*), intent(in)                :: location
    character(len=*), allocatable, dimension(:) :: groups

    integer :: nmembers, ngroups
    character(len=150) :: strim
    character(len=150), allocatable, dimension(:) :: members
    integer(hid_t) :: obj_id
    integer :: obj_type
    integer :: idx

    ngroups = 0

    call h5gn_members_f(ifile, trim(adjustl(location)), nmembers, hdf_err)
    if (nmembers .eq. 0) return
    call h5oopen_f(ifile, trim(adjustl(location)), obj_id, hdf_err)
    allocate(members(nmembers))
    do idx = 0, nmembers-1
      call h5gget_obj_info_idx_f(ifile, trim(adjustl(location)), idx, strim, obj_type, hdf_err)
      if (obj_type .eq. h5g_group_f) then
        ngroups = ngroups + 1
        members(ngroups) = trim(strim)
      endif
    enddo
    call h5oclose_f(obj_id, hdf_err)
    if (ngroups .eq. 0) return
    allocate(groups(ngroups))
    groups(:) = members(:ngroups)
    deallocate(members)
  end subroutine

  subroutine hdf5_list_datasets(ifile, location, dsets)
    integer(hid_t), intent(in)                  :: ifile
    character(len=*), intent(in)                :: location
    character(len=*), allocatable, dimension(:) :: dsets

    integer :: nmembers, ndsets
    character(len=150) :: strim
    character(len=150), allocatable, dimension(:) :: members
    integer(hid_t) :: obj_id
    integer :: obj_type
    integer :: idx

    ndsets = 0

    call h5gn_members_f(ifile, trim(adjustl(location)), nmembers, hdf_err)
    if (nmembers .eq. 0) return
    call h5oopen_f(ifile, trim(adjustl(location)), obj_id, hdf_err)
    allocate(members(nmembers))
    do idx = 0, nmembers-1
      call h5gget_obj_info_idx_f(ifile, trim(adjustl(location)), idx, strim, obj_type, hdf_err)
      if (obj_type .eq. h5g_dataset_f) then
        ndsets = ndsets + 1
        members(ndsets) = trim(strim)
      endif
    enddo
    call h5oclose_f(obj_id, hdf_err)
    if (ndsets .eq. 0) return
    allocate(dsets(ndsets))
    dsets(:) = members(:ndsets)
    deallocate(members)
  end subroutine

  subroutine hdf5_list_attributes(ifile, location, attrs)
    integer(hid_t), intent(in)                  :: ifile
    character(len=*), intent(in)                :: location
    character(len=*), allocatable, dimension(:) :: attrs

    integer :: nattrs
    integer(hid_t) :: obj_id
    integer :: idx
    character(len=150) :: strim

    call h5oopen_f(ifile, trim(adjustl(location)), obj_id, hdf_err)
    call h5aget_num_attrs_f(obj_id, nattrs, hdf_err)
    call h5oclose_f(obj_id, hdf_err)

    if (nattrs .eq. 0) return
    allocate(attrs(nattrs))
    do idx= 0, nattrs-1
      call h5aget_name_by_idx_f(ifile, trim(adjustl(location)), h5_index_name_f, h5_iter_inc_f, int(idx, hsize_t), strim, hdf_err)
      attrs(idx+1) = strim
    enddo
  end subroutine

  subroutine hdf5_delete_attribute(ifile, location, attr)
    integer(hid_t), intent(in)   :: ifile
    character(len=*), intent(in) :: location
    character(len=*), intent(in) :: attr

    integer(hid_t) :: obj_id
    call h5oopen_f(ifile, trim(adjustl(location)), obj_id, hdf_err)
    call h5adelete_f(obj_id, attr, hdf_err)
    call h5oclose_f(obj_id, hdf_err)
  end subroutine hdf5_delete_attribute

  subroutine hdf5_delete(ifile, location)
    integer(hid_t), intent(in)                  :: ifile
    character(len=*), intent(in)                :: location
    call h5ldelete_f(ifile, location, hdf_err)
  end subroutine hdf5_delete

  logical function hdf5_group_exists(ifile, location)
    integer(hid_t), intent(in)                  :: ifile
    character(len=*), intent(in)                :: location

    integer(hid_t) :: grp_parent_id, grp_id
    character(len=150) :: strim
    integer :: grpcnt, pos, i
    character(len=128), allocatable :: ngroup(:)

    call h5eset_auto_f(0,hdf_err) ! deactivate error printing
    call h5gopen_f(ifile, trim(adjustl(location)), grp_id, hdf_err)
    if (hdf_err .ne. 0) then
      hdf5_group_exists = .false.
    else
      hdf5_group_exists = .true.
      call h5gclose_f(grp_id, hdf_err)
    endif
    call h5eclear_f(hdf_err)      ! clear the error
    call h5eset_auto_f(1,hdf_err) ! acitvate it again
  end function

  logical function hdf5_dataset_exists(ifile, location)
    integer(hid_t), intent(in)                  :: ifile
    character(len=*), intent(in)                :: location

    integer(hid_t) :: dset_id

    call h5eset_auto_f(0,hdf_err) ! deactivate error printing
    call h5dopen_f(ifile, trim(adjustl(location)), dset_id, hdf_err)
    if (hdf_err .ne. 0) then
      hdf5_dataset_exists = .false.
    else
      hdf5_dataset_exists = .true.
      call h5dclose_f(dset_id, hdf_err)
    endif
    call h5eclear_f(hdf_err)      ! clear the error
    call h5eset_auto_f(1,hdf_err) ! acitvate it again
  end function

  logical function hdf5_attribute_exists(ifile, location, attrname)
    integer(hid_t), intent(in)                  :: ifile
    character(len=*), intent(in)                :: location
    character(len=*), intent(in)                :: attrname

    integer(hid_t) :: obj_id, attr_id

    call h5eset_auto_f(0,hdf_err) ! deactivate error printing
    call h5oopen_f(ifile, trim(adjustl(location)), obj_id, hdf_err)
    if (hdf_err .ne. 0) then
      hdf5_attribute_exists = .false.
    else
      call h5aopen_f(obj_id, trim(adjustl(attrname)), attr_id, hdf_err)
      if (hdf_err .ne. 0) then
        hdf5_attribute_exists = .false.
      else
        call h5aclose_f(attr_id, hdf_err)
        hdf5_attribute_exists = .true.
      endif
    endif

    call h5eclear_f(hdf_err)      ! clear the error
    call h5eset_auto_f(1,hdf_err) ! acitvate it again
  end function

  ! help function with separates /group1/group2/dset -> /group1/group2 and dset
  subroutine hdf5_help_separate_dsetname(dsetfull, gname, dset)
    character(len=*), intent(in)    :: dsetfull
    character(len=150), intent(out) :: dset, gname

    character(len=150) :: strim
    integer            :: pos

    ! adjusting string and removing possible leading and trailing '/'
    strim = adjustl(trim(dsetfull))
    if (index(strim, '/', back=.true.) .eq. len_trim(strim)) then
      strim = strim(:(len_trim(strim)-1))
    endif
    if (index(strim, '/') .eq. 1) then
      strim = strim(2:)
    endif

    pos = index(strim, '/', back=.true.)
    if (pos .eq. 0) then
      dset  = strim
      gname = ''
    else
      gname = strim(:(pos-1))
      dset  = strim((pos+1):)
    endif
  end subroutine hdf5_help_separate_dsetname

end module hdf5_wrapper
