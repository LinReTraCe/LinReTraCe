module Mio
  use Mparams
  use Mtypes
  use Maux
  use hdf5
  use hdf5_wrapper
  implicit none

contains

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
  if (edisp%ispin == 2) then
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
        call hdf5_read_data(ifile, "/kpoint/"//nmbstring//"/derivatives",   drank2arr)
        edisp%band_dk(:,:,i,1) = drank2arr
        deallocate(drank2arr)
        call hdf5_read_data(ifile, "/kpoint/"//nmbstring//"/curvatures",  drank2arr)
        edisp%band_d2k(:,:,i,1) = drank2arr
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
            call hdf5_read_data(ifile, "/up/kpoint/"//nmbstring//"/curvatures",    drank2arr)
            edisp%band_d2k(:,:,i,1) = drank2arr
            deallocate(drank2arr)
          else
            call hdf5_read_data(ifile, "/dn//kpoint/"//nmbstring//"/derivatives", drank2arr)
            edisp%band_dk(:,:,i,2) = drank2arr
            deallocate(drank2arr)
            call hdf5_read_data(ifile, "/dn//kpoint/"//nmbstring//"/curvatures",   drank2arr)
            edisp%band_d2k(:,:,i,2) = drank2arr
            deallocate(drank2arr)
          endif
        endif
      enddo
    enddo
  endif

end subroutine

subroutine read_preproc_scattering_data(algo, kmesh, edisp, sct)
  implicit none
  type(algorithm)              :: algo
  type(kpointmesh)             :: kmesh
  type(energydisp)             :: edisp
  type(scattering)             :: sct

  integer                      :: kpoints
  integer                      :: nbands
  integer                      :: iSpin
  integer(hid_t)               :: ifile
  character(len=6)             :: nmbstring
  real(8), allocatable         :: drank1arr(:)
  real(8), allocatable         :: drank2arr(:,:)

  integer :: iT, ik

  call hdf5_init()
  call hdf5_open_file(trim(adjustl(algo%input_scattering)), ifile, rdonly=.true.)

  ! mesh
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


  ! temperature grid
  call hdf5_read_data(ifile, "/.quantities/Tmin", sct%Tmin)
  call hdf5_read_data(ifile, "/.quantities/Tmax", sct%Tmax)
  call hdf5_read_data(ifile, "/.quantities/nT",   sct%nT)

  if (abs(sct%Tmin - sct%Tmax) < 1e-5 .and. (sct%nT .ne. 1)) then
     call stop_with_message(stderr, "Temperature grid is not properly defined")
  endif

  allocate(sct%TT(sct%nT))
  allocate(sct%mu(sct%nT))
  allocate(sct%d1(sct%nT))
  allocate(sct%d2(sct%nT))
  allocate(sct%d0(sct%nT))

  ! define Temperature grid
  do iT=1,sct%nT
     sct%TT(iT)=real(iT-1,8)*sct%dT+sct%Tmin
  enddo
  sct%TT(sct%nT) = sct%Tmax ! to avoid numerical errors at the last point
  sct%beta = 1.d0/(sct%TT * kB)

  ! scattering rates
  ! and quasi particle renormalizations
  allocate(sct%gam(edisp%nband_max, sct%nT))
  allocate(sct%zqp(edisp%nband_max, sct%nT))

  if (edisp%iSpin == 1) then
    if (hdf5_get_number_groups(ifile, "/kPoint/000001") > 2) then
       edisp%lBandShift = .true.
       allocate(edisp%band_shift(edisp%nband_max, kmesh%nkp, edisp%iSpin))
    else
       edisp%lBandShift = .false.
    endif

    if (edisp%lBandShift) then
      do ik=1,kmesh%nkp
         write(nmbstring,'(I6.6)') ik
         call hdf5_read_data(ifile, "/kPoint/"//nmbstring//"/bandshift", drank1arr)
         edisp%band_shift(:,ik,1)  = drank1arr
         deallocate(drank1arr)
      enddo
    endif

  else if (edisp%iSpin == 2) then
    if (hdf5_get_number_groups(ifile, "/up/kPoint/000001") > 2) then
       edisp%lBandShift = .true.
       allocate(edisp%band_shift(edisp%nband_max, kmesh%nkp, edisp%iSpin))
    else
       edisp%lBandShift = .false.
    endif

    if (edisp%lBandShift) then
      do ik=1,kmesh%nkp
         write(nmbstring,'(I6.6)') ik
         call hdf5_read_data(ifile, "/up/kPoint/"//nmbstring//"/bandshift", drank1arr)
         edisp%band_shift(:,ik,1)  = drank1arr
         deallocate(drank1arr)
         call hdf5_read_data(ifile, "/dn/kPoint/"//nmbstring//"/bandshift", drank1arr)
         edisp%band_shift(:,ik,2)  = drank1arr
         deallocate(drank1arr)
      enddo
    endif
  endif

  call log_master(stdout, 'Detected band shifts in ScatteringFile')

end subroutine

end module Mio
