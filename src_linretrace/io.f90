module Mio
  use Mparams
  use Mtypes
  use Maux
  use hdf5
  use hdf5_wrapper
  implicit none

contains

subroutine read_preproc_energy_data(algo, kmesh, edisp)
  implicit none
  type(algorithm)      :: algo
  type(kpointmesh)     :: kmesh
  type(energydisp)     :: edisp

  integer(hid_t)       :: ifile
  integer              :: i, is, locderivatives
  integer              :: locgapped
  integer              :: nshape(1)

  integer, allocatable :: irank1arr(:)
  real(8), allocatable :: drank1arr(:)
  real(8), allocatable :: drank2arr(:,:)
  real(8), allocatable :: drank3arr(:,:,:)

  call hdf5_open_file(trim(adjustl(algo%input_energies)), ifile, rdonly=.true.)

  ! mesh information
  call hdf5_read_data(ifile, "/.kmesh/nkp",         kmesh%nkp)
  call hdf5_read_data(ifile, "/.kmesh/weightsum",   kmesh%weightsum)
  call hdf5_read_data(ifile, "/.kmesh/weights",     kmesh%weight)


  ! band information + charge
  call hdf5_read_data(ifile, "/.bands/charge",         edisp%nelect)
  call hdf5_read_data(ifile, "/.bands/energyBandMax",  edisp%nband_max)
  call hdf5_read_data(ifile, "/.bands/opticalBandMin", edisp%nbopt_min)
  call hdf5_read_data(ifile, "/.bands/opticalBandMax", edisp%nbopt_max)
  call hdf5_read_data(ifile, "/.bands/ispin",          edisp%ispin)

  allocate(edisp%gapped(edisp%ispin))
  allocate(edisp%gap(edisp%ispin))

  ! read band gap information
  if (edisp%ispin == 1) then
    call hdf5_read_data(ifile, "/.bands/bandgap/gapped", locgapped)
    if (locgapped == 1) then
      edisp%gapped(1) = .true.
      call hdf5_read_data(ifile, "/.bands/bandgap/gapsize", edisp%gap(1))
      call hdf5_read_data(ifile, "/.bands/bandgap/vband", edisp%valenceBand)
    else
      edisp%gapped(1) = .false.
    endif
  else
    call hdf5_read_data(ifile, "/.bands/bandgap/up/gapped", locgapped)
    if (locgapped == 1) then
      edisp%gapped(1) = .true.
      call hdf5_read_data(ifile, "/.bands/bandgap/gapsize", edisp%gap(1))
      call hdf5_read_data(ifile, "/.bands/bandgap/vband", edisp%valenceBand)
    else
      edisp%gapped(1) = .false.
    endif

    call hdf5_read_data(ifile, "/.bands/bandgap/dn/gapped", locgapped)
    if (locgapped == 1) then
      edisp%gapped(2) = .true.
      call hdf5_read_data(ifile, "/.bands/bandgap/gapsize", edisp%gap(2))
      call hdf5_read_data(ifile, "/.bands/bandgap/vband", edisp%valenceBand)
    else
      edisp%gapped(2) = .false.
    endif
  endif

  if (algo%lScissors) then
    nshape = shape(edisp%scissors)
    if (nshape(1) /= edisp%iSpin) then
      call stop_with_message(stderr, 'Must have as many scissor parameters as spins')
    endif

    do is=1,edisp%ispin
      if (.not. edisp%gapped(is) .and. abs(edisp%scissors(is)) > 0.d0) then
        call log_master(stdout, 'Warning: Cannot apply scissors to gapless energies')
      endif
    enddo
  endif

  if (algo%lScissors .and. .not. algo%lScatteringFile) then
    allocate(edisp%band_shift(edisp%nband_max, kmesh%nkp, edisp%ispin))
    do is=1,edisp%ispin
      edisp%band_shift(:edisp%valenceBand, :, is)   = 0.d0
      edisp%band_shift(edisp%valenceBand+1:, :, is) = edisp%scissors(is)
    enddo
  endif

  ! unit cell information
  call hdf5_read_data(ifile, "/.unitcell/volume", kmesh%vol)

  ! number of saved k-points
  if (edisp%ispin == 2) then
    kmesh%nkp = hdf5_get_number_groups(ifile, "/up/kPoint")
    if (hdf5_group_exists(ifile, "up/derivatives") .and. &
        hdf5_group_exists(ifile, "up/curvatures")) then
      edisp%lDerivatives = .true.
    else
      edisp%lDerivatives = .false.
    endif
  else
    kmesh%nkp = hdf5_get_number_groups(ifile, "/kPoint")
    if (hdf5_group_exists(ifile, "/kPoint/derivatives") .and. &
        hdf5_group_exists(ifile, "/kPoint/curvatures")) then
      edisp%lDerivatives = .true.
    else
      edisp%lDerivatives = .false.
    endif
  endif

  ! now we load the energy data into the according arrays
  ! please be aware here about the implicit Fortran memory transposition
  ! which is happening when loading hdf5 files

  allocate(edisp%band_original(edisp%nband_max, kmesh%nkp, edisp%ispin))
  allocate(edisp%band(edisp%nband_max, kmesh%nkp, edisp%ispin))
  if (edisp%lDerivatives) then
    allocate(edisp%band_dk(3, edisp%nband_max, kmesh%nkp, edisp%ispin))
    allocate(edisp%band_d2k(6, edisp%nband_max, kmesh%nkp, edisp%ispin))
  endif


  if (edisp%ispin == 2) then
    call hdf5_get_shape(ifile, "/up/kPoint/000001/moments", irank1arr)
  else
    call hdf5_get_shape(ifile, "/kPoint/000001/moments", irank1arr)
  endif
  edisp%iOptical = irank1arr(1)
  deallocate(irank1arr)

  ! the optical elements get loaded only for one k-point each time
  allocate(edisp%Mopt(edisp%iOptical,edisp%nbopt_min:edisp%nbopt_max, &
                                     edisp%nbopt_min:edisp%nbopt_max, edisp%ispin))


  if (edisp%ispin == 1) then
    call hdf5_read_data(ifile, "/energies", drank2arr)
    edisp%band_original(:,:,1)     = drank2arr
    deallocate(drank2arr)
    if (edisp%lDerivatives) then
      call hdf5_read_data(ifile, "/derivatives",   drank3arr)
      edisp%band_dk(:,:,:,1) = drank3arr
      deallocate(drank3arr)
      call hdf5_read_data(ifile, "/curvatures",  drank3arr)
      edisp%band_d2k(:,:,:,1) = drank3arr
      deallocate(drank3arr)
    endif
  else if (edisp%ispin == 2) then
    call hdf5_read_data(ifile, "/up/energies", drank2arr)
    edisp%band_original(:,:,1)     = drank2arr
    deallocate(drank2arr)
    call hdf5_read_data(ifile, "/dn/energies", drank2arr)
    edisp%band_original(:,:,2)     = drank2arr
    deallocate(drank2arr)

    if (edisp%lDerivatives) then
      call hdf5_read_data(ifile, "/up/derivatives",  drank3arr)
      edisp%band_dk(:,:,:,1) = drank3arr
      deallocate(drank3arr)
      call hdf5_read_data(ifile, "/up/curvatures",    drank3arr)
      edisp%band_d2k(:,:,:,1) = drank3arr
      deallocate(drank3arr)
      call hdf5_read_data(ifile, "/dn/derivatives", drank3arr)
      edisp%band_dk(:,:,:,2) = drank3arr
      deallocate(drank3arr)
      call hdf5_read_data(ifile, "/dn/curvatures",   drank3arr)
      edisp%band_d2k(:,:,:,2) = drank3arr
      deallocate(drank3arr)
    endif
  endif

  call hdf5_close_file(ifile)

end subroutine

subroutine read_preproc_scattering_data(algo, kmesh, edisp, sct, temp)
  implicit none
  type(algorithm)              :: algo
  type(kpointmesh)             :: kmesh
  type(energydisp)             :: edisp
  type(scattering)             :: sct
  type(temperature)            :: temp

  integer                      :: kpoints
  integer                      :: nbands
  integer                      :: iSpin
  integer(hid_t)               :: ifile
  real(8), allocatable         :: drank1arr(:)
  real(8), allocatable         :: drank2arr(:,:)

  integer :: iT, ik

  call hdf5_init()
  call hdf5_open_file(trim(adjustl(algo%input_scattering)), ifile, rdonly=.true.)

  ! sanity check
  call hdf5_read_data(ifile, "/.quantities/nkp",    kpoints)
  call hdf5_read_data(ifile, "/.quantities/nbands", nbands)
  call hdf5_read_data(ifile, "/.quantities/iSpin",  iSpin)

  if ( kpoints /= kmesh%nkp ) then
     call stop_with_message(stderr, "Number of k-points in preprocessed scattering data &
     do not match")
  endif
  if ( nbands /= edisp%nband_max ) then
     call stop_with_message(stderr, "Number of bands in preprocessed scattering data &
     do not match")
  endif
  if ( iSpin /= edisp%iSpin ) then
     call stop_with_message(stderr, "Number of spins in preprocessed scattering data &
     do not match")
  endif


  ! temperature grid ... already pre-computated
  call hdf5_read_data(ifile, "/.quantities/Tmin", temp%Tmin)
  call hdf5_read_data(ifile, "/.quantities/Tmax", temp%Tmax)
  call hdf5_read_data(ifile, "/.quantities/nT",   temp%nT)

  call hdf5_read_data(ifile, "/.quantities/tempAxis", temp%TT)
  call hdf5_read_data(ifile, "/.quantities/betaAxis", temp%beta)

  ! scattering rates
  ! and quasi particle renormalizations
  allocate(sct%gam(edisp%nband_max, kmesh%nkp, edisp%ispin))
  allocate(sct%zqp(edisp%nband_max, kmesh%nkp, edisp%ispin))

  if (edisp%iSpin == 1) then
    if (hdf5_dataset_exists(ifile, "/tPoint/000001/bandshift")) then
       edisp%lBandShift = .true.
       allocate(edisp%band_shift(edisp%nband_max, kmesh%nkp, edisp%iSpin))
    else
       edisp%lBandShift = .false.
    endif

  else if (edisp%iSpin == 2) then
    if (hdf5_dataset_exists(ifile, "/up/tPoint/000001/bandshift")) then
       edisp%lBandShift = .true.
       allocate(edisp%band_shift(edisp%nband_max, kmesh%nkp, edisp%iSpin))
    else
       edisp%lBandShift = .false.
    endif
  endif

  call hdf5_close_file(ifile)

end subroutine

subroutine output_auxiliary(algo, info, temp, kmesh)
  implicit none
  type(algorithm)   :: algo
  type(runinfo)     :: info
  type(temperature) :: temp
  type(kpointmesh)  :: kmesh

  character(len=128) :: string
  integer(hid_t)     :: ifile

  call hdf5_open_file(algo%output_file, ifile)

  call hdf5_write_data(ifile, '.quantities/tempAxis', temp%TT)
  call hdf5_write_data(ifile, '.quantities/betaAxis', temp%beta)
  call hdf5_write_data(ifile, '.quantities/weights',  kmesh%weight)

  call hdf5_close_file(ifile)

end subroutine

subroutine output_energies(mu, algo, edisp, kmesh, sct, info)
  implicit none
  type(algorithm)  :: algo
  type(energydisp) :: edisp
  type(kpointmesh) :: kmesh
  type(scattering) :: sct
  type(runinfo)    :: info

  real(8), intent(in) :: mu

  integer(hid_t)     :: ifile
  real(8), allocatable :: enrgy(:,:,:)
  character(len=128) :: string

  allocate(enrgy(edisp%nbopt_min:edisp%nbopt_max,kmesh%nkp,edisp%ispin))

  enrgy = sct%zqp(:,:,:) * (edisp%band(edisp%nbopt_min:edisp%nbopt_max,:,:) - mu)

  write(string,'(I6.6,"/energies")') info%iT
  call hdf5_open_file(algo%output_file, ifile)
  call hdf5_write_data(ifile, string, enrgy)

  write(string,'(I6.6)') info%iT
  call hdf5_write_attribute(ifile, string, "temperature", info%temp)
  call hdf5_write_attribute(ifile, string, "invtemperature", info%beta)

  call hdf5_close_file(ifile)

  deallocate(enrgy)

end subroutine

subroutine read_scattering_data(ifile, edisp, sct, info)
  implicit none
  integer(hid_t)   :: ifile
  type(energydisp) :: edisp
  type(scattering) :: sct
  type(runinfo)    :: info

  real(8), allocatable :: darr2(:,:)
  character(len=128)   :: string

  if (edisp%ispin == 1) then
    write(string,'("tPoint/",I6.6,"/scatrate")') info%iT
    call hdf5_read_data(ifile, string, darr2)
    sct%gam(:,:,1) = darr2
    deallocate(darr2)

    write(string,'("tPoint/",I6.6,"/qpweight")') info%iT
    call hdf5_read_data(ifile, string, darr2)
    sct%zqp(:,:,1) = darr2
    deallocate(darr2)

    if (edisp%lBandShift) then
      write(string,'("tPoint/",I6.6,"/bandshift")') info%iT
      call hdf5_read_data(ifile, string, darr2)
      edisp%band_shift(:,:,1) = darr2
      deallocate(darr2)

      edisp%band = edisp%band_original + edisp%band_shift
    endif
  else
    write(string,'("up/tPoint/",I6.6,"/scatrate")') info%iT
    call hdf5_read_data(ifile, string, darr2)
    sct%gam(:,:,1) = darr2
    deallocate(darr2)

    write(string,'("dn/tPoint/",I6.6,"/scatrate")') info%iT
    call hdf5_read_data(ifile, string, darr2)
    sct%gam(:,:,2) = darr2
    deallocate(darr2)

    write(string,'("up/tPoint/",I6.6,"/qpweight")') info%iT
    call hdf5_read_data(ifile, string, darr2)
    sct%zqp(:,:,1) = darr2
    deallocate(darr2)

    write(string,'("dn/tPoint/",I6.6,"/qpweight")') info%iT
    call hdf5_read_data(ifile, string, darr2)
    sct%zqp(:,:,2) = darr2
    deallocate(darr2)


    if (edisp%lBandShift) then
      write(string,'("up/tPoint/",I6.6,"/bandshift")') info%iT
      call hdf5_read_data(ifile, string, darr2)
      edisp%band_shift(:,:,1) = darr2
      deallocate(darr2)

      write(string,'("dn/tPoint/",I6.6,"/bandshift")') info%iT
      call hdf5_read_data(ifile, string, darr2)
      edisp%band_shift(:,:,2) = darr2
      deallocate(darr2)

      edisp%band = edisp%band_original + edisp%band_shift
    endif
  endif
  sct%gam = sct%gam + sct%gamimp ! so we have access to a constant shift right from the config file
  sct%gam = sct%gam * sct%zqp    ! convention

end subroutine

subroutine read_optical_elements(ifile, edisp, sct, info)
  implicit none
  integer(hid_t)   :: ifile
  type(energydisp) :: edisp
  type(scattering) :: sct
  type(runinfo)    :: info

  real(8), allocatable :: darr3(:,:,:)
  character(len=128)   :: string

  if (edisp%ispin == 1) then
    if (allocated(darr3)) deallocate(darr3)
    write(string,'("kPoint/",I6.6,"/moments")') info%ik
    call hdf5_read_data(ifile, string, darr3)
    edisp%Mopt(:,:,:,1) = darr3
    deallocate(darr3)
  else
    if (allocated(darr3)) deallocate(darr3)
    write(string,'("up/kPoint/",I6.6,"/moments")') info%ik
    call hdf5_read_data(ifile, string, darr3)
    edisp%Mopt(:,:,:,1) = darr3
    deallocate(darr3)
    write(string,'("dn/kPoint/",I6.6,"/moments")') info%ik
    call hdf5_read_data(ifile, string, darr3)
    edisp%Mopt(:,:,:,2) = darr3
    deallocate(darr3)
  endif

end subroutine

end module Mio