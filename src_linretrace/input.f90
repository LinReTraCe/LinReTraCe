module Minput
  use Mtypes
  use Maux
  use hdf5
  use hdf5_wrapper
  implicit none

contains

subroutine read_config(algo, edisp)
  implicit none
  type(algorithm)    :: algo
  type(energydisp)   :: edisp
  character(len=256) :: config_file

  real(8) :: mu ! local mu variable

  if (iargc() .ne. 1) then
    er = 1
    erstr = 'The program has to be executed with exactly one argument. (Name of config file)'
    return
  end if
  call getarg(1,config_file)

  open(unit=10,file=trim(config_file),action='read',iostat=stat)
  if (stat .ne. 0) then
    call stop_with_message(0, 'Input file cannot be opened') ! send to stderr
  endif

  open(10,file=config_file,status='old')
  read(10,*) algo%input_bands
  read(10,*) algo%input_scattering
  read(10,*) algo%mumethod, mu
  read(10,*) algo%rootmethod
  read(10,*) algo%lBfield

  if (algo%mumethd == 1) then
    edisp%efer = mu
  endif


  close(10)

end subroutine


subroutine read_preproc_band_data(fname, kmesh, edisp, lat)
  implicit none
  character(len=*), intent(in) :: fname
  type(kpointmesh)             :: kmesh
  type(energydisp)             :: edisp
  type(lattice)                :: lat

  integer(hid_t)               :: ifile
  logical                      :: derivatives_exist
  integer                      :: locortho, i
  character(len=6)             :: nmbstring
  real(8), allocatable         :: drank1arr(:)
  real(8), allocatable         :: drank3arr(:,:,:)
  complex(8), allocatable      :: crank3arr(:,:,:)

  call hdf5_init()
  call hdf5_open_file(trim(adjustl(fname)), ifile, rdonly=.true.)

  ! mesh
  call hdf5_read_data(ifile, "/.kmesh/ktot",         kmesh%ktot)
  call hdf5_read_data(ifile, "/.kmesh/multiplicity", kmesh%multiplicity)
  call hdf5_read_data(ifile, "/.kmesh/weight",       kmesh%weight)

  ! edisp
  call hdf5_read_data(ifile, "/.edisp/nelect",   edisp%band_fill_value)
  call hdf5_read_data(ifile, "/.edisp/nelect",   edisp%nelect)
  call hdf5_read_data(ifile, "/.edisp/efer",     edisp%efer)


  ! crystal
  call hdf5_read_data(ifile, "/.crystal/vol",    lat%vol)
  call hdf5_read_data(ifile, "/.crystal/lortho", locortho) ! there is no boolean type in hdf5
  call hdf5_read_data(ifile, "/.crystal/nalpha", lat%nalpha

  if (locortho == 1) then
     lat%lortho = .true.
     lat%nalpha = 3 ! polarization directions
  else
     lat%lortho = .false.
     lat%nalpha = 6 ! polarization directions
  endif

  ! bands
  call hdf5_read_data(ifile, "/.bands/band_max",         edisp%nband_max)
  call hdf5_read_data(ifile, "/.bands/optical_band_min", edisp%nbopt_min)
  call hdf5_read_data(ifile, "/.bands/optical_band_max", edisp%nbopt_max)

  ! number of saved k-points
  kmesh%ktot = hdf5_get_number_groups(ifile, "/kpoint")
  write(*,*) "MAIN: found ", kmesh%ktot, " kpoints in the preprocessed file."

  ! k-point information
  allocate(edisp%band(edisp%nband_max, kmesh%ktot))
  if (lat%lortho) then
     allocate(edisp%Mopt(3, edisp%nbopt_min:edisp%nbopt_max, edisp%nbopt_min:edisp%nbopt_max))
  else
     allocate(edisp%Mopt(6, edisp%nbopt_min:edisp%nbopt_max, edisp%nbopt_min:edisp%nbopt_max))
  endif

  if (hdf5_get_number_groups(ifile, "/kpoint/000001") > 2) then
     derivatives_exist = .true.
     allocate(edisp%band_dk(3,  edisp%nbopt_min:edisp%nbopt_max, kmesh%ktot))
     allocate(edisp%band_d2k(6, edisp%nbopt_min:edisp%nbopt_max, kmesh%ktot))
  else
     derivatives_exist = .false.
  endif

  do i=1,kmesh%ktot
    write(nmbstring,'(I6.6)') i
    call hdf5_read_data(ifile, "/kpoint/"//nmbstring//"/energies", rank1arr)
    edisp%band(:,i)     = rank1arr
    deallocate(rank1arr)
    call hdf5_read_data(ifile, "/kpoint/"//nmbstring//"/optical",  crank3arr)
    edisp%Mopt(:,:,:,i) = crank3arr
    deallocate(crank3arr)
    if (derivatives_exist) then
      call hdf5_read_data(ifile, "/kpoint/"//nmbstring//"/energies_dk",   drank3arr)
      edisp%band_dk(:,:,:,i) = drank3arr
      deallocate(drank3arr)
      call hdf5_read_data(ifile, "/kpoint/"//nmbstring//"/energies_d2k",  drank3arr)
      edisp%Mopt(:,:,:,i) = drank3arr
      deallocate(drank3arr)
    endif
  enddo

end subroutine

subroutine read_preproc_scat_data(fname, kmesh, edisp, scat)
  implicit none
  character(len=*), intent(in) :: fname
  type(kpointmesh)             :: kmesh
  type(energydisp)             :: edisp
  type(scattering)             :: scat

  real(8)                      :: kpoints
  real(8)                      :: nbands
  integer(hid_t)               :: ifile
  character(len=6)             :: nmbstring
  real(8), allocatable         :: drank2arr(:)

  call hdf5_init()
  call hdf5_open_file(trim(adjustl(fname)), ifile, rdonly=.true.)

  ! mesh
  call hdf5_read_data(ifile, "/.quantities/ktot",   kpoints)
  call hdf5_read_data(ifile, "/.quantities/nbands", nbands)

  if ( kpoints /= kmesh%ktot ) then
     call stop_with_message(0, "Number of k-points in preprocessed scattering data &
     does not match")
  endif
  if ( nbands /= edisp%nband_max ) then
     call stop_with_message(0, "Number of bands in preprocessed scattering data &
     does not match")
  endif


  ! temperature grid
  call hdf5_read_data(ifile, "/.quantities/Tmin", scat%Tmin)
  call hdf5_read_data(ifile, "/.quantities/Tmax", scat%Tmax)
  call hdf5_read_data(ifile, "/.quantities/nT",   scat%nT)

  if ((Tmin .eq. Tmax) .and. (nT .ne. 1)) then
     call stop_with_message(0, "Temperature grid is not properly defined")
  endif

  allocate(scat%TT(nT))
  allocate(scat%mu(nT))
  allocate(scat%d1(nT))
  allocate(scat%d2(nT))
  allocate(scat%d0(nT))


  ! scattering rates
  ! and quasi particle renormalizations
  allocate(scat%gam(edisp%nband_max, kmesh%ktot, scat%nT))
  allocate(scat%zqp(edisp%nband_max, kmesh%ktot, scat%nT))

  if (hdf5_get_number_groups(ifile, "/kpoint/000001") > 2) then
     shift_exist = .true.
     allocate(edisp%band_shift(edisp%nband_max, kmesh%ktot, nT)
  else
     shift_exist = .false.
  endif

  do i=1,kmesh%ktot
     write(nmbstring,'(I6.6)') i
     call hdf5_read_data(ifile, "/Tpoint/"//nmbstring//"/gamma", rank2arr)
     scat%gam(:,:,i)  = rank2arr
     deallocate(rank2arr)
     call hdf5_read_data(ifile, "/Tpoint/"//nmbstring//"/zqp",   rank2arr)
     scat%zqp(:,:,i)   = rank2arr
     deallocate(rank2arr)
     if (shift_exist) then
       call hdf5_read_data(ifile, "/kpoint/"//nmbstring//"/shift",   rank2arr)
       edisp%band_shift(:,:) = rank2arr
       deallocate(drank2arr)
     endif
  enddo

end subroutine

end module Minput
