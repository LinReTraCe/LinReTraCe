program hdf5test
  use hdf5
  implicit none

  write(*,*) 'begin HDF5 tests'
  call hdf5_interface()
  call hdf5_file()

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

  write(*,*) 'end HDF5 Tests'
  end subroutine

end program
