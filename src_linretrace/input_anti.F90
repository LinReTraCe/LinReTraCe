module Minput_hk
  use Mparams
  use Mtypes
  use Mauxiliary
  use hdf5
  use hdf5_wrapper
  implicit none

contains

! reads in the "Berry Curvature" for the current k step

subroutine read_berry_curv(ifile, edisp, info)
  implicit none
  integer(hid_t)   :: ifile
  type(energydisp) :: edisp
  type(runinfo)    :: info

  real(8), allocatable :: darr3(:,:,:)
  character(len=128)   :: string

  if (edisp%ispin == 1) then
    if (allocated(darr3)) deallocate(darr3)
    write(string,'("kPoint/",I10.10,"/BerryCurv")') info%ik
    call hdf5_read_data(ifile, string, darr3)
    edisp%BerryCurvk(:,:,:,1) = darr3
    deallocate(darr3)
  else
    if (allocated(darr3)) deallocate(darr3)
    write(string,'("up/kPoint/",I10.10,"/BerryCurv")') info%ik
    call hdf5_read_data(ifile, string, darr3)
    edisp%BerryCurvk(:,:,:,1) = darr3
    deallocate(darr3)
    write(string,'("dn/kPoint/",I10.10,"/BerryCurv")') info%ik
    call hdf5_read_data(ifile, string, darr3)
    edisp%BerryCurvk(:,:,:,2) = darr3
    deallocate(darr3)
  endif

end subroutine

end module Minput_hk