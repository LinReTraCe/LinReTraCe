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

  integer :: stat ! file status

  if (iargc() .ne. 1) then
    call stop_with_message(5, 'The program has to be executed with exactly one argument. (Name of config file)')
  end if

  call getarg(1,config_file)

  open(unit=10,file=trim(config_file),action='read',iostat=stat)
  if (stat .ne. 0) then
    call stop_with_message(0, 'Input file cannot be opened') ! send to stderr
  endif

  open(10,file=config_file,status='old')
  read(10,*) algo%lDebug
  read(10,*) algo%lBfield
  read(10,*) algo%rootmethod
  read(10,*) algo%mumethod, edisp%mu
  read(10,*)
  read(10,*) algo%input_energies       ! preprocessing of e(k) and such
  read(10,*) algo%input_scattering     ! preprocessing of Gamma(k,n,T)

  close(10)

end subroutine


subroutine read_preproc_energy_data(algo, kmesh, edisp, lat)
  implicit none
  type(algorithm)              :: algo
  type(kpointmesh)             :: kmesh
  type(energydisp)             :: edisp
  type(lattice)                :: lat

  integer(hid_t)               :: ifile
  integer                      :: locortho, i, is, locderivatives
  character(len=6)             :: nmbstring
  real(8), allocatable         :: drank1arr(:)
  real(8), allocatable         :: drank2arr(:,:)

  call hdf5_init()
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


  ! unit cell information
  call hdf5_read_data(ifile, "/.unitcell/volume", lat%vol)
  call hdf5_read_data(ifile, "/.unitcell/ortho",  locortho)

  if (locortho == 1) then
     lat%lOrtho = .true.
     lat%nalpha = 3 ! polarization directions
  else
     lat%lortho = .false.
     lat%nalpha = 6 ! polarization directions
  endif

  ! number of saved k-points
  if (edisp%ispin == 1) then
    kmesh%nkp = hdf5_get_number_groups(ifile, "/up/kpoint")
    locderivatives = hdf5_get_number_groups(ifile, "/up/kpoint/000001")
  else
    kmesh%nkp = hdf5_get_number_groups(ifile, "/kpoint")
    locderivatives = hdf5_get_number_groups(ifile, "/kpoint/000001")
  endif


  ! given by the number of hdf5 groups per kpoint
  ! we can identify whether the derivatives are saved or not
  if (locderivatives == 2) then
    edisp%lDerivatives = .false.
  else
    edisp%lDerivatives = .true.
  endif

  write(*,*) "MAIN: found ", kmesh%nkp, " kpoints in the preprocessed file."
  if (edisp%lDerivatives) write(*,*) "MAIN: found energy derivatives"


  ! now we load the energy data into the according arrays
  ! please be aware here about the implicit Fortran memory transposition
  ! which is happening when loading hdf5 files

  allocate(edisp%band(edisp%nband_max, kmesh%nkp, edisp%ispin))
  if (edisp%lDerivatives) then
    allocate(edisp%band_dk(3, edisp%nband_max, kmesh%nkp, edisp%ispin))
    allocate(edisp%band_d2k(6, edisp%nband_max, kmesh%nkp, edisp%ispin))
  endif


  ! the optical elements get loaded only for one k-point each time
  allocate(edisp%Mopt(6, edisp%nbopt_min:edisp%nbopt_max, edisp%nbopt_min:edisp%nbopt_max, edisp%ispin))


  if (edisp%ispin == 1) then
    do i=1,kmesh%nkp
      write(nmbstring,'(I6.6)') i
      call hdf5_read_data(ifile, "/kpoint/"//nmbstring//"/energies", drank1arr)
      edisp%band(:,i,1)     = drank1arr
      deallocate(drank1arr)
      if (edisp%lDerivatives) then
        call hdf5_read_data(ifile, "/kpoint/"//nmbstring//"/energies_dk",   drank2arr)
        edisp%band_dk(:,:,i,1) = drank2arr
        deallocate(drank2arr)
        call hdf5_read_data(ifile, "/kpoint/"//nmbstring//"/energies_d2k",  drank2arr)
        edisp%Mopt(:,:,i,1) = drank2arr
        deallocate(drank2arr)
      endif
    enddo
  else if (edisp%ispin == 2) then
    do is=1,edisp%ispin
      do i=1,kmesh%nkp
        write(nmbstring,'(I6.6)') i

        if (is==1) then
          call hdf5_read_data(ifile, "/up/kpoint/"//nmbstring//"/energies", drank1arr)
          edisp%band(:,i,1)     = drank1arr
        else
          call hdf5_read_data(ifile, "/dn/kpoint/"//nmbstring//"/energies", drank1arr)
          edisp%band(:,i,2)     = drank1arr
        endif
        deallocate(drank1arr)

        if (edisp%lDerivatives) then
          if (is==1) then
            call hdf5_read_data(ifile, "/up/kpoint/"//nmbstring//"/derivatives",  drank2arr)
            edisp%band_dk(:,:,i,1) = drank2arr
            deallocate(drank2arr)
            call hdf5_read_data(ifile, "/up/kpoint/"//nmbstring//"/curvature",    drank2arr)
            edisp%Mopt(:,:,i,1) = drank2arr
            deallocate(drank2arr)
          else
            call hdf5_read_data(ifile, "/dn//kpoint/"//nmbstring//"/derivatives", drank2arr)
            edisp%band_dk(:,:,i,2) = drank2arr
            deallocate(drank2arr)
            call hdf5_read_data(ifile, "/dn//kpoint/"//nmbstring//"/curvature",   drank2arr)
            edisp%Mopt(:,:,i,2) = drank2arr
            deallocate(drank2arr)
          endif
        endif
      enddo
    enddo
  endif

end subroutine

! subroutine read_preproc_scat_data(algo, kmesh, edisp, scat)
!   implicit none
!   type(algorithm)              :: algo
!   type(kpointmesh)             :: kmesh
!   type(energydisp)             :: edisp
!   type(scattering)             :: scat

!   real(8)                      :: kpoints
!   real(8)                      :: nbands
!   integer(hid_t)               :: ifile
!   character(len=6)             :: nmbstring
!   real(8), allocatable         :: drank2arr(:)

!   call hdf5_init()
!   call hdf5_open_file(trim(adjustl(algo%input_scattering)), ifile, rdonly=.true.)

!   ! mesh
!   call hdf5_read_data(ifile, "/.quantities/ktot",   kpoints)
!   call hdf5_read_data(ifile, "/.quantities/nbands", nbands)

!   if ( kpoints /= kmesh%ktot ) then
!      call stop_with_message(0, "Number of k-points in preprocessed scattering data &
!      does not match")
!   endif
!   if ( nbands /= edisp%nband_max ) then
!      call stop_with_message(0, "Number of bands in preprocessed scattering data &
!      does not match")
!   endif


!   ! temperature grid
!   call hdf5_read_data(ifile, "/.quantities/Tmin", scat%Tmin)
!   call hdf5_read_data(ifile, "/.quantities/Tmax", scat%Tmax)
!   call hdf5_read_data(ifile, "/.quantities/nT",   scat%nT)

!   if ((Tmin .eq. Tmax) .and. (nT .ne. 1)) then
!      call stop_with_message(0, "Temperature grid is not properly defined")
!   endif

!   allocate(scat%TT(nT))
!   allocate(scat%mu(nT))
!   allocate(scat%d1(nT))
!   allocate(scat%d2(nT))
!   allocate(scat%d0(nT))


!   ! scattering rates
!   ! and quasi particle renormalizations
!   allocate(scat%gam(edisp%nband_max, kmesh%ktot, scat%nT))
!   allocate(scat%zqp(edisp%nband_max, kmesh%ktot, scat%nT))

!   if (hdf5_get_number_groups(ifile, "/kpoint/000001") > 2) then
!      shift_exist = .true.
!      allocate(edisp%band_shift(edisp%nband_max, kmesh%ktot, nT)
!   else
!      shift_exist = .false.
!   endif

!   do i=1,kmesh%ktot
!      write(nmbstring,'(I6.6)') i
!      call hdf5_read_data(ifile, "/Tpoint/"//nmbstring//"/gamma", rank2arr)
!      scat%gam(:,:,i)  = rank2arr
!      deallocate(rank2arr)
!      call hdf5_read_data(ifile, "/Tpoint/"//nmbstring//"/zqp",   rank2arr)
!      scat%zqp(:,:,i)   = rank2arr
!      deallocate(rank2arr)
!      if (shift_exist) then
!        call hdf5_read_data(ifile, "/kpoint/"//nmbstring//"/shift",   rank2arr)
!        edisp%band_shift(:,:) = rank2arr
!        deallocate(drank2arr)
!      endif
!   enddo

! end subroutine

end module Minput
