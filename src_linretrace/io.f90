module Mio
  use Mparams
  use Mtypes
  use Maux
  use hdf5
  use hdf5_wrapper
  implicit none

contains

subroutine read_preproc_energy_data(algo, kmesh, edisp, sct, pot, imp)
  implicit none
  type(algorithm)      :: algo
  type(kpointmesh)     :: kmesh
  type(energydisp)     :: edisp
  type(scattering)     :: sct
  type(potential)      :: pot
  type(impurity)       :: imp

  integer(hid_t)       :: ifile
  integer              :: i, is, locderivatives, iimp, ik
  integer              :: n1shape(1)
  integer              :: n2shape(2)

  integer, allocatable :: irank1arr(:)
  real(8), allocatable :: drank1arr(:)
  real(8), allocatable :: drank2arr(:,:)
  real(8), allocatable :: drank3arr(:,:,:)
  real(8), allocatable :: drank5arr(:,:,:,:,:)

  call hdf5_open_file(trim(adjustl(algo%input_energies)), ifile, rdonly=.true.)

  ! mesh information
  call hdf5_read_data(ifile, "/.kmesh/nkp",         kmesh%nkp)
  call hdf5_read_data(ifile, "/.kmesh/weightsum",   kmesh%weightsum)
  call hdf5_read_data(ifile, "/.kmesh/weights",     kmesh%weight)
  call hdf5_read_data(ifile, "/.kmesh/multiplicity",kmesh%multiplicity)

  ! messy way to achieve proper quad precision in weightQ
  allocate(kmesh%weightQ(kmesh%nkp))
  do ik=1,kmesh%nkp
    kmesh%weightQ(ik) = int(kmesh%multiplicity(ik)) * int(kmesh%weightsum) / real(int(sum(kmesh%multiplicity)),16)
  enddo


  ! band information + charge
  call hdf5_read_data(ifile, "/.bands/mu",             pot%mu_dft_file)
  call hdf5_read_data(ifile, "/.bands/charge",         edisp%nelect_file)
  call hdf5_read_data(ifile, "/.bands/energyBandMax",  edisp%nband_max)
  call hdf5_read_data(ifile, "/.bands/opticalBandMin", edisp%nbopt_min)
  call hdf5_read_data(ifile, "/.bands/opticalBandMax", edisp%nbopt_max)
  call hdf5_read_data(ifile, "/.bands/ispin",          edisp%ispin)

  ! reset electron occupation from config file if it is invalid ... earlierst point to check is here
  if (edisp%nelect_config >= kmesh%weightsum * edisp%nband_max) then ! too many electrons
    edisp%nelect_config = -1.d0
    call stop_with_message(stderr, 'Error: Number of electrons provided is above maximum')
  endif

  if ((edisp%nelect_config > 0.d0)) then
    edisp%nelect = edisp%nelect_config
  else
    edisp%nelect = edisp%nelect_file
  endif

  allocate(edisp%gapped(edisp%ispin))
  allocate(edisp%gap(edisp%ispin))
  allocate(edisp%valenceBand(edisp%ispin))
  allocate(edisp%conductionBand(edisp%ispin))
  allocate(edisp%ene_valenceBand(edisp%ispin))
  allocate(edisp%ene_conductionBand(edisp%ispin))

  edisp%gap_min = 0.d0
  ! read band gap information
  if (edisp%ispin == 1) then
    call hdf5_read_data(ifile, "/.bands/bandgap/gapped", edisp%gapped(1))
    if (edisp%gapped(1)) then
      call hdf5_read_data(ifile, "/.bands/bandgap/gapsize", edisp%gap(1))
      call hdf5_read_data(ifile, "/.bands/bandgap/vband", edisp%valenceBand(1))
      call hdf5_read_data(ifile, "/.bands/bandgap/ene_vband", edisp%ene_valenceBand(1))
      call hdf5_read_data(ifile, "/.bands/bandgap/cband", edisp%conductionBand(1))
      call hdf5_read_data(ifile, "/.bands/bandgap/ene_cband", edisp%ene_conductionBand(1))
    else
      edisp%gap(1) = 0.d0
    endif

  else
    call hdf5_read_data(ifile, "/.bands/bandgap/up/gapped", edisp%gapped(1))
    if (edisp%gapped(1)) then
      call hdf5_read_data(ifile, "/.bands/bandgap/up/gapsize", edisp%gap(1))
      call hdf5_read_data(ifile, "/.bands/bandgap/up/vband", edisp%valenceBand(1))
      call hdf5_read_data(ifile, "/.bands/bandgap/up/ene_vband", edisp%ene_valenceBand(1))
      call hdf5_read_data(ifile, "/.bands/bandgap/up/cband", edisp%conductionBand(1))
      call hdf5_read_data(ifile, "/.bands/bandgap/up/ene_cband", edisp%ene_conductionBand(1))
    else
      edisp%gap(1) = 0.d0
    endif

    call hdf5_read_data(ifile, "/.bands/bandgap/dn/gapped", edisp%gapped(2))
    if (edisp%gapped(2)) then
      call hdf5_read_data(ifile, "/.bands/bandgap/dn/gapsize", edisp%gap(2))
      call hdf5_read_data(ifile, "/.bands/bandgap/dn/vband", edisp%valenceBand(2))
      call hdf5_read_data(ifile, "/.bands/bandgap/dn/ene_vband", edisp%ene_valenceBand(2))
      call hdf5_read_data(ifile, "/.bands/bandgap/dn/cband", edisp%conductionBand(2))
      call hdf5_read_data(ifile, "/.bands/bandgap/dn/ene_cband", edisp%ene_conductionBand(2))
    else
      edisp%gap(2) = 0.d0
    endif
  endif

  ! check if the gap is complete
  ! i.e. check if we are spin-dependent, that a gap exists in both spins
  edisp%gapped_complete = .true.
  do is=1,edisp%ispin
    if (.not. edisp%gapped(is)) edisp%gapped_complete = .false.
  enddo
  edisp%gap_min = minval(edisp%gap)

  if (algo%lImpurities) then
    do iimp = 1, imp%nimp
      ! coming from the config file we might have 0 1 or 2 -> set to 1 for internal access
      if (edisp%ispin == 1) then
        imp%inputspin(iimp) = 1
      endif

      ! we have to explicitly provide the spin
      if (edisp%ispin == 2 .and. .not. (imp%inputspin(iimp) == 1 .or. imp%inputspin(iimp) == 2)) then
        call stop_with_message(stderr, 'Error: Spin type must be 1 or 2')
      endif

      if ((imp%inputtype(iimp) > 1) .and. (edisp%gapped(imp%inputspin(iimp)) .eqv. .false.)) then
        call stop_with_message(stderr, 'Error: Relative impurity position not available in gapless system')
      endif
    enddo
  endif

  if (.not. algo%lScatteringFile) then
    if ( edisp%ispin==2 .and. size(sct%gamcoeff,dim=1)==1 ) then
      ! increase the spin axis by one
      n2shape = shape(sct%gamcoeff)
      allocate(drank2arr(1,n2shape(2)))
      drank2arr(:,:) = sct%gamcoeff
      deallocate(sct%gamcoeff)
      allocate(sct%gamcoeff(2,n2shape(2)))
      sct%gamcoeff(1,:) = drank2arr(1,:)
      sct%gamcoeff(2,:) = drank2arr(1,:)
      deallocate(drank2arr)
    endif
    if ( edisp%ispin==1 .and. size(sct%gamcoeff,dim=1)==2) then
      call stop_with_message(stderr, 'Error: Provided Scattering rates do not match system (one spin type)')
    endif
    if ( edisp%ispin==2 .and. size(sct%zqpcoeff,dim=1)==1 ) then
      ! increase the spin axis by one
      n2shape = shape(sct%zqpcoeff)
      allocate(drank2arr(1,n2shape(2)))
      drank2arr(:,:) = sct%zqpcoeff
      deallocate(sct%zqpcoeff)
      allocate(sct%zqpcoeff(2,n2shape(2)))
      sct%zqpcoeff(1,:) = drank2arr(1,:)
      sct%zqpcoeff(2,:) = drank2arr(1,:)
      deallocate(drank2arr)
    endif
    if ( edisp%ispin==1 .and. size(sct%zqpcoeff,dim=1)==2) then
      call stop_with_message(stderr, 'Error: Provided QuasiParticle weights do not match system (one spin type)')
    endif
  endif

  if (algo%lScissors) then
    ! check for inconsitencies
    n1shape = shape(edisp%scissors)
    if (n1shape(1) /= edisp%iSpin) then
      call stop_with_message(stderr, 'Must have as many scissor parameters as spins')
    endif

    ! save to scissors operator
    do is=1,edisp%ispin
      edisp%scissors(is) = edisp%scissors(is) - edisp%gap(is) ! transform bandgap input to scissors
      if (.not. edisp%gapped(is) .and. abs(edisp%scissors(is)) > 0.d0) then
        call stop_with_message(stdout, 'Error: Cannot apply scissors to gapless band structure')
      endif
    enddo
  endif

  ! unit cell information
  call hdf5_read_data(ifile, "/.unitcell/volume", kmesh%vol)
  call hdf5_read_data(ifile, "/.unitcell/ndim",   kmesh%ndim)
  call hdf5_read_data(ifile, "/.unitcell/dims",   kmesh%dims)

  ! number of saved k-points
  if (edisp%ispin == 2) then
    if (hdf5_dataset_exists(ifile, "up/momentsDiagonalBfield")) then
      edisp%lDerivatives = .true.
    else
      edisp%lDerivatives = .false.
    endif
  else
    if (hdf5_dataset_exists(ifile, "momentsDiagonalBfield")) then
      edisp%lDerivatives = .true.
    else
      edisp%lDerivatives = .false.
    endif
  endif

  ! now we load the energy data into the according arrays
  ! please be aware here about the implicit Fortran memory transposition
  ! which is happening when loading hdf5 files


  ! LOAD THE DIAGONAL ELEMENTS -> these are always here
  if (edisp%ispin == 2) then
    ! get the shape for the number of optical directions
    call hdf5_get_shape(ifile, "up/momentsDiagonal", irank1arr)
  else
    call hdf5_get_shape(ifile, "momentsDiagonal", irank1arr)
  endif
  edisp%iOptical = irank1arr(1)
  deallocate(irank1arr)
  allocate(edisp%MoptDiag(edisp%iOptical, edisp%nbopt_min:edisp%nbopt_max, edisp%ispin, kmesh%nkp))

  if (edisp%ispin == 2) then
    call hdf5_read_data(ifile, "up/momentsDiagonal", drank3arr)
    edisp%MoptDiag(:,:,1,:) = drank3arr
    deallocate(drank3arr)
    call hdf5_read_data(ifile, "dn/momentsDiagonal", drank3arr)
    edisp%MoptDiag(:,:,2,:) = drank3arr
    deallocate(drank3arr)
  else
    call hdf5_read_data(ifile, "momentsDiagonal", drank3arr)
    edisp%MoptDiag(:,:,1,:) = drank3arr
    deallocate(drank3arr)
  endif

  ! FULL ELEMENTS -> if they are here we detect it and are able to calculate inter band contributions
  if (edisp%ispin == 2) then
    if (hdf5_dataset_exists(ifile, "up/kPoint/0000000001/moments") .and. &
        hdf5_dataset_exists(ifile, "dn/kPoint/0000000001/moments")) then
      edisp%lFullMoments = .true.
      call hdf5_get_shape(ifile, "up/kPoint/0000000001/moments", irank1arr)
    else
      edisp%lFullMoments = .false.
    endif
  else
    if (hdf5_dataset_exists(ifile, "kPoint/0000000001/moments")) then
      edisp%lFullMoments = .true.
      call hdf5_get_shape(ifile, "/kPoint/0000000001/moments", irank1arr)
    else
      edisp%lFullMoments = .false.
    endif
  endif

  ! the optical elements get loaded only for one k-point each time
  if (edisp%lFullMoments) then
    allocate(edisp%Moptk(edisp%iOptical,edisp%nbopt_min:edisp%nbopt_max, &
                                       edisp%nbopt_min:edisp%nbopt_max, edisp%ispin))
  endif


  ! ENERGIES & DERIVATIVES
  allocate(edisp%band_file(edisp%nband_max, kmesh%nkp, edisp%ispin))
  allocate(edisp%band(edisp%nband_max, kmesh%nkp, edisp%ispin))
  if (edisp%lDerivatives) then
    allocate(edisp%band_dk(3, edisp%nband_max, kmesh%nkp, edisp%ispin))
    allocate(edisp%band_d2k(6, edisp%nband_max, kmesh%nkp, edisp%ispin))
    allocate(edisp%MBoptdiag(3, 3, 3, edisp%nband_max, kmesh%nkp, edisp%ispin))
  endif

  if (edisp%ispin == 1) then
    call hdf5_read_data(ifile, "/energies", drank2arr)
    edisp%band_file(:,:,1)     = drank2arr
    deallocate(drank2arr)
    if (edisp%lDerivatives) then
      call hdf5_read_data(ifile, "/derivatives",   drank3arr)
      edisp%band_dk(:,:,:,1) = drank3arr
      deallocate(drank3arr)
      call hdf5_read_data(ifile, "/curvatures",  drank3arr)
      edisp%band_d2k(:,:,:,1) = drank3arr
      deallocate(drank3arr)

      call hdf5_read_data(ifile, "/momentsDiagonalBfield",   drank5arr)
      edisp%MBoptDiag(:,:,:,:,:,1) = drank5arr
      deallocate(drank5arr)
    endif

  else if (edisp%ispin == 2) then
    call hdf5_read_data(ifile, "/up/energies", drank2arr)
    edisp%band_file(:,:,1)     = drank2arr
    deallocate(drank2arr)
    call hdf5_read_data(ifile, "/dn/energies", drank2arr)
    edisp%band_file(:,:,2)     = drank2arr
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

      call hdf5_read_data(ifile, "/up/momentsDiagonalBfield",   drank5arr)
      edisp%MBoptDiag(:,:,:,:,:,1) = drank5arr
      deallocate(drank5arr)
      call hdf5_read_data(ifile, "/dn/momentsDiagonalBfield",   drank5arr)
      edisp%MBoptDiag(:,:,:,:,:,2) = drank5arr
      deallocate(drank5arr)
    endif
  endif

  call hdf5_close_file(ifile)

end subroutine

subroutine set_impurities(edisp, imp)
  implicit none
  type(energydisp) :: edisp
  type(impurity)   :: imp

  integer :: iimp

  ! now that we have all the information we can adjust the energies of the impurity levels
  do iimp = 1, imp%nimp
    select case (imp%inputtype(iimp))
      ! case 1 -> already absolute
      case (2) ! relative upwards shift from top of valence band
        imp%Energy(iimp) = imp%Energy(iimp) + edisp%ene_valenceBand(imp%inputspin(iimp))
      case (3) ! relative downwards shift from bottom of conduction band
        imp%Energy(iimp) = -imp%Energy(iimp) + edisp%ene_conductionBand(imp%inputspin(iimp))
      case (4) ! percentage gap shift from top of valence band
        imp%Energy(iimp) = edisp%ene_valenceBand(imp%inputspin(iimp)) &
                         + edisp%gap(imp%inputspin(iimp)) * imp%Energy(iimp)
    end select
  enddo

end subroutine

subroutine read_preproc_scattering_data_hdf5(algo, kmesh, edisp, sct, pot, temp)
  implicit none
  type(algorithm)              :: algo
  type(kpointmesh)             :: kmesh
  type(energydisp)             :: edisp
  type(scattering)             :: sct
  type(potential)              :: pot
  type(temperature)            :: temp

  integer                      :: kpoints
  integer                      :: nbands
  integer                      :: iSpin
  integer(hid_t)               :: ifile
  real(8), allocatable         :: drank1arr(:)
  real(8), allocatable         :: drank2arr(:,:)

  integer :: iT, ik
  logical :: ltemp

  call hdf5_init()
  call hdf5_open_file(trim(adjustl(algo%input_scattering_hdf5)), ifile, rdonly=.true.)

  ! sanity check
  call hdf5_read_data(ifile, "/.quantities/nkp",    kpoints)
  call hdf5_read_data(ifile, "/.quantities/nbands", nbands)
  call hdf5_read_data(ifile, "/.quantities/iSpin",  iSpin)

  call hdf5_read_data(ifile, "/.quantities/tempmode",  ltemp)

  if (algo%lTMODE .neqv. ltemp) then
    call stop_with_message(stderr, "Scattering File mode does not coincide with Config file mode")
  endif

  ! allow identical k-meshes / bands
  ! or kpoints == 1 or nbands == 1 ... this way we are able to save storage
  if ( (kpoints /= kmesh%nkp) .and. .not. (kpoints == 1)) then
     call stop_with_message(stderr, "Number of k-points in preprocessed scattering data &
     do not match")
  endif
  if ( (nbands /= edisp%nband_max) .and. .not. (nbands == 1)) then
     call stop_with_message(stderr, "Number of bands in preprocessed scattering data &
     do not match")
  endif
  if ( iSpin /= edisp%iSpin .and. .not. (iSpin ==1)) then
     call stop_with_message(stderr, "Number of spins in preprocessed scattering data &
     do not match")
  endif

  ! temperature grid ... already pre-computated
  call hdf5_read_data(ifile, "/.quantities/Tmin", temp%Tmin)
  call hdf5_read_data(ifile, "/.quantities/Tmax", temp%Tmax)
  call hdf5_read_data(ifile, "/.quantities/nT",   algo%steps)

  call hdf5_read_data(ifile, "/.quantities/tempAxis", temp%TT)
  call hdf5_read_data(ifile, "/.quantities/betaAxis", temp%BB)

  call hdf5_read_data(ifile, "/.quantities/muAxis", pot%MM)
  allocate(pot%QMM(algo%steps))
  pot%QMM = pot%MM


  ! scattering rates
  ! and quasi particle renormalizations
  allocate(sct%gam(edisp%nband_max, kmesh%nkp, edisp%ispin))
  allocate(sct%zqp(edisp%nband_max, kmesh%nkp, edisp%ispin))

  edisp%lBandShift = .false.
  if (edisp%iSpin == 1) then
    if (hdf5_dataset_exists(ifile, "/step/000001/bandshift")) then
      edisp%lBandShift = .true.
    endif
  else if (edisp%iSpin == 2) then
    if (hdf5_dataset_exists(ifile, "/up/step/000001/bandshift")) then
      edisp%lBandShift = .true.
    endif
  endif

  if (edisp%lBandShift) then
    allocate(edisp%band_shift(edisp%nband_max, kmesh%nkp, edisp%iSpin))
  endif

  pot%mumin = minval(pot%MM)
  pot%mumax = maxval(pot%MM)

  call hdf5_close_file(ifile)

end subroutine

subroutine read_preproc_scattering_data_text(algo, kmesh, edisp, sct, pot, temp)
  implicit none
  type(algorithm)              :: algo
  type(kpointmesh)             :: kmesh
  type(energydisp)             :: edisp
  type(scattering)             :: sct
  type(potential)              :: pot
  type(temperature)            :: temp

  character(len=256) :: str_temp
  integer            :: kpoints
  integer            :: nbands
  integer            :: iT, ik, i,j
  integer            :: lines
  integer            :: stat
  integer            :: empty, pst
  integer            :: cnt
  real(8), allocatable            :: float_array(:)
  character(len=256), allocatable :: file_temp(:), file_save(:)
  character(len=1)                :: cmnt = '#'
  character(len=1)                :: sprt = ' '
  real(8) :: fdum1, fdum2


  open(unit=10,file=trim(algo%input_scattering_text),action='read',iostat=stat)
  if (stat .ne. 0) then
    call stop_with_message(stderr, 'ScatteringText Input file cannot be opened') ! send to stderr
  endif

  ! line counting
  lines=0
  read_count: do
    read(10,'(A)',END=200) str_temp ! read whole line as string, doesnt skip empty lines
    lines=lines+1
  enddo read_count


  200 continue
  rewind 10

  allocate(file_temp(lines))

  ! remove empty lines and comment strings
  empty=0
  read_temp: do i=1,lines
    read(10,'(A)') str_temp
      str_temp = trim(adjustl(str_temp))
      pst=scan(str_temp,cmnt) ! find out possible comment symbol
      if (pst .eq. 1) then ! whole line is commented
        file_temp(i) = ''
      elseif (pst .eq. 0) then ! no comment symbol found
        file_temp(i) = str_temp
      else  ! getting left side of comment
        file_temp(i) = trim(str_temp(1:pst-1))
      endif

      if (len_trim(file_temp(i)) .eq. 0) then ! filter out empty lines
        empty=empty+1
      endif
  enddo read_temp

  ! rewrite everything to a new clean string array
  allocate(file_save(lines-empty))
  j=1
  read_save: do i=1,lines
    if(len_trim(file_temp(i)) .ne. 0) then
      file_save(j)=file_temp(i)
      j=j+1
    endif
  enddo read_save
  deallocate(file_temp)

  lines=lines-empty
  close(unit=10)


  ! check if data has right amount of columns via first row
  cnt = 0
  str_temp = file_save(1)
  do
    pst=scan(str_temp,sprt)
    if (pst == 1) then ! we find the empty space in an empty string at the first position
      exit
    else
      cnt = cnt + 1
      str_temp = trim(adjustl(str_temp(pst+1:)))
    endif
  enddo

  if ((cnt /= 3) .and. (cnt /= 5)) then
    call stop_with_message(stderr, 'Error: ScatteringText file does not have correct number of columns (3 or 5)')
  endif
  if ((edisp%ispin==1) .and. (cnt==5)) then
    call stop_with_message(stderr, 'Error: ScatteringText file does not match energy file (number of spins == 1)')
  endif

  ! now we read this string array into the according data arrays
  algo%steps = lines
  allocate(temp%TT(algo%steps))
  allocate(temp%BB(algo%steps))
  allocate(sct%gamtext(algo%steps, edisp%ispin))
  allocate(sct%zqptext(algo%steps, edisp%ispin))
  allocate(sct%gam(edisp%nband_max, kmesh%nkp, edisp%ispin))
  allocate(sct%zqp(edisp%nband_max, kmesh%nkp, edisp%ispin))

  if (cnt == 3) then
    do i=1,algo%steps
      read(file_save(i),*) temp%TT(i), sct%gamtext(i,1), sct%zqptext(i,1)
    enddo
  else
    do i=1,algo%steps
      read(file_save(i),*) temp%TT(i), sct%gamtext(i,1), sct%gamtext(i,2), sct%zqptext(i,1), sct%zqptext(i,2)
    enddo
  endif
  if ((edisp%ispin==2) .and. (cnt==3)) then
    sct%gamtext(:,2) = sct%gamtext(:,1)
    sct%zqptext(:,2) = sct%zqptext(:,1)
  endif

  temp%tmin = minval(temp%TT)
  temp%tmax = maxval(temp%TT)
  temp%BB = 1.d0/(temp%TT * kB)

  sct%gamtext = sct%gamtext + sct%gamimp

end subroutine

subroutine output_auxiliary(algo, info, pot, temp, kmesh, edisp, sct, imp)
  implicit none
  type(algorithm)   :: algo
  type(runinfo)     :: info
  type(potential)   :: pot
  type(temperature) :: temp
  type(kpointmesh)  :: kmesh
  type(energydisp)  :: edisp
  type(scattering)  :: sct
  type(impurity)    :: imp

  character(len=256) :: string
  integer(hid_t)     :: ifile
  integer            :: iimp

  call hdf5_open_file(algo%output_file, ifile)

  call hdf5_create_group(ifile, '.config')
  call hdf5_write_attribute(ifile, '.config', 'tmode', algo%lTMODE)
  call hdf5_write_attribute(ifile, '.config', 'mumode',algo%lMUMODE)
  call hdf5_write_attribute(ifile, '.config', 'doping',algo%ldoping)
  call hdf5_write_attribute(ifile, '.config', 'quad',algo%lQuad)
  call hdf5_write_attribute(ifile, '.config', 'debug', algo%lDebug)
  call hdf5_write_attribute(ifile, '.config', 'bfield', algo%lBfield)
  call hdf5_write_attribute(ifile, '.config', 'rootmethod', algo%rootMethod)
  call hdf5_write_attribute(ifile, '.config', 'musearch', algo%muSearch)
  call hdf5_write_attribute(ifile, '.config', 'mufermi', algo%muFermi)
  call hdf5_write_attribute(ifile, '.config', 'oldmu', algo%lOldmu)
  call hdf5_write_attribute(ifile, '.config', 'oldmutext', algo%lOldmuText)
  call hdf5_write_attribute(ifile, '.config', 'scatteringfile', algo%lScatteringFile)
  call hdf5_write_attribute(ifile, '.config', 'scatteringtext', algo%lScatteringText)
  call hdf5_write_attribute(ifile, '.config', 'interbandquantities', algo%lInterBandQuantities)
  call hdf5_write_attribute(ifile, '.config', 'intrabandquantities', algo%lIntraBandQuantities)
  call hdf5_write_attribute(ifile, '.config', 'fulloutput', algo%fullOutput)
  call hdf5_write_attribute(ifile, '.config', 'energyOutput', algo%lEnergyOutput)
  call hdf5_write_attribute(ifile, '.config', 'boltzmann', algo%lBoltzmann)
  call hdf5_write_attribute(ifile, '.config', 'scissors', algo%lScissors)
  call hdf5_write_attribute(ifile, '.config', 'steps', algo%steps)
  call hdf5_write_attribute(ifile, '.config', 'steps_dir', algo%step_dir)
  call hdf5_write_attribute(ifile, '.config', 'impurities', algo%lImpurities)
  call hdf5_write_attribute(ifile, '.config', 'input_energies', trim(algo%input_energies))
  if (len(trim(algo%input_scattering_hdf5)) == 0) then
    call hdf5_write_attribute(ifile, '.config', 'input_scattering_hdf5', "-")
  else
    call hdf5_write_attribute(ifile, '.config', 'input_scattering_hdf5', trim(algo%input_scattering_hdf5))
  endif
  if (len(trim(algo%input_scattering_text)) == 0) then
    call hdf5_write_attribute(ifile, '.config', 'input_scattering_text', "-")
  else
    call hdf5_write_attribute(ifile, '.config', 'input_scattering_text', trim(algo%input_scattering_text))
  endif
  if (len(trim(algo%old_output_file)) == 0) then
    call hdf5_write_attribute(ifile, '.config', 'old_output_file', "-")
  else
    call hdf5_write_attribute(ifile, '.config', 'old_output_file', trim(algo%old_output_file))
  endif
  if (len(trim(algo%input_mu_text)) == 0) then
    call hdf5_write_attribute(ifile, '.config', 'input_mu_text', "-")
  else
    call hdf5_write_attribute(ifile, '.config', 'input_mu_text', trim(algo%input_mu_text))
  endif
  if (len(trim(algo%dbgstr)) == 0) then
    call hdf5_write_attribute(ifile, '.config', 'dbgstr', "-")
  else
    call hdf5_write_attribute(ifile, '.config', 'dbgstr', trim(algo%dbgstr))
  endif

  call hdf5_create_group(ifile, '.scattering')
  if (allocated(sct%gamcoeff)) then
    call hdf5_write_data(ifile, '.scattering/gamcoeff', sct%gamcoeff)
  endif
  if (allocated(sct%zqpcoeff)) then
    call hdf5_write_data(ifile, '.scattering/zqpcoeff', sct%zqpcoeff)
  endif
  if (allocated(sct%gamtext)) then
    call hdf5_write_data(ifile, '.scattering/gamtext', sct%gamtext)
  endif
  if (allocated(sct%zqptext)) then
    call hdf5_write_data(ifile, '.scattering/zqptext', sct%zqptext)
  endif
  if (algo%lScatteringFile .or. algo%lScatteringText) then
    call hdf5_write_data(ifile, '.scattering/gamimp', sct%gamimp)
  endif

  call hdf5_write_data(ifile, ".quantities/ispin", edisp%ispin)
  call hdf5_write_data(ifile, ".quantities/charge", edisp%nelect) ! this might have been changed by config
  call hdf5_write_data(ifile, ".quantities/doping", edisp%doping)
  call hdf5_write_data(ifile, '.quantities/mudft', pot%mu_dft)    ! this also might have changed
  call hdf5_write_data(ifile, '.quantities/tempAxis', temp%TT)
  call hdf5_write_data(ifile, '.quantities/betaAxis', temp%BB)
  call hdf5_write_data(ifile, '.quantities/weights',  kmesh%weight)
  if (algo%lTMODE) then
    call hdf5_write_attribute(ifile, '.quantities', 'mode', 'temp')
  else if (algo%lMUMODE) then
    call hdf5_write_attribute(ifile, '.quantities', 'mode', 'mu')
  endif
  call hdf5_write_attribute(ifile, '.quantities', 'identifier', 'LRTC')

  ! output bandgap information to have access to it
  ! in the general output file
  if (edisp%ispin == 1) then
    call hdf5_write_data(ifile, "/.quantities/bandgap/gapped", edisp%gapped(1))
    if (edisp%gapped(1)) then
      call hdf5_write_data(ifile, "/.quantities/bandgap/gapsize", edisp%gap(1))
      call hdf5_write_data(ifile, "/.quantities/bandgap/ene_vband", edisp%ene_valenceBand(1))
      call hdf5_write_data(ifile, "/.quantities/bandgap/ene_cband", edisp%ene_conductionBand(1))
      call hdf5_write_data(ifile, "/.quantities/bandgap/vband", edisp%valenceBand(1))
      call hdf5_write_data(ifile, "/.quantities/bandgap/cband", edisp%conductionBand(1))
    endif
  else
    call hdf5_write_data(ifile, "/.quantities/bandgap/up/gapped", edisp%gapped(1))
    if (edisp%gapped(1)) then
      call hdf5_write_data(ifile, "/.quantities/bandgap/up/gapsize", edisp%gap(1))
      call hdf5_write_data(ifile, "/.quantities/bandgap/up/ene_vband", edisp%ene_valenceBand(1))
      call hdf5_write_data(ifile, "/.quantities/bandgap/up/ene_cband", edisp%ene_conductionBand(1))
      call hdf5_write_data(ifile, "/.quantities/bandgap/up/vband", edisp%valenceBand(1))
      call hdf5_write_data(ifile, "/.quantities/bandgap/up/cband", edisp%conductionBand(1))
    endif

    call hdf5_write_data(ifile, "/.quantities/bandgap/dn/gapped", edisp%gapped(2))
    if (edisp%gapped(2)) then
      call hdf5_write_data(ifile, "/.quantities/bandgap/dn/gapsize", edisp%gap(2))
      call hdf5_write_data(ifile, "/.quantities/bandgap/dn/ene_vband", edisp%ene_valenceBand(2))
      call hdf5_write_data(ifile, "/.quantities/bandgap/dn/ene_cband", edisp%ene_conductionBand(2))
      call hdf5_write_data(ifile, "/.quantities/bandgap/dn/vband", edisp%valenceBand(1))
      call hdf5_write_data(ifile, "/.quantities/bandgap/dn/cband", edisp%conductionBand(1))
    endif
  endif

  call hdf5_write_data(ifile, "/.quantities/impurities/nimp", imp%nimp)
  if (algo%lImpurities) then
    do iimp = 1, imp%nimp
      write(string,'("/.quantities/impurities/imp-",I3.3,"/energy")') iimp
      call hdf5_write_data(ifile, string, imp%Energy(iimp))
      write(string,'("/.quantities/impurities/imp-",I3.3,"/density")') iimp
      call hdf5_write_data(ifile, string, imp%Density(iimp))
      write(string,'("/.quantities/impurities/imp-",I3.3,"/degeneracy")') iimp
      call hdf5_write_data(ifile, string, imp%Degeneracy(iimp))
      write(string,'("/.quantities/impurities/imp-",I3.3,"/dopant")') iimp
      call hdf5_write_data(ifile, string, imp%Dopant(iimp))
      write(string,'("/.quantities/impurities/imp-",I3.3,"/width")') iimp
      call hdf5_write_data(ifile, string, imp%Bandwidth(iimp))
    enddo
  endif

  call hdf5_write_data(ifile, "/.unitcell/dims", kmesh%dims)
  call hdf5_write_data(ifile, "/.unitcell/ndim", kmesh%ndim)
  call hdf5_write_data(ifile, "/.unitcell/vol",  kmesh%vol)

  if (edisp%ispin == 1) then
    call hdf5_write_data(ifile, "/.energies", edisp%band(:,:,1))
  else
    call hdf5_write_data(ifile, "/.energies/up", edisp%band(:,:,1))
    call hdf5_write_data(ifile, "/.energies/dn", edisp%band(:,:,2))
  endif

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

  allocate(enrgy(edisp%nband_max,kmesh%nkp,edisp%ispin))

  enrgy = sct%zqp * (edisp%band - mu)

  write(string,'(I6.6,"/energies")') info%iStep
  call hdf5_open_file(algo%output_file, ifile)
  call hdf5_write_data(ifile, string, enrgy)

  write(string,'(I6.6)') info%iStep
  call hdf5_write_attribute(ifile, string, "temperature", info%temp)
  call hdf5_write_attribute(ifile, string, "invtemperature", info%beta)

  call hdf5_close_file(ifile)

  deallocate(enrgy)

end subroutine

subroutine read_scattering_data_hdf5(ifile, edisp, kmesh, sct, info)
  implicit none
  integer(hid_t)   :: ifile
  type(energydisp) :: edisp
  type(kpointmesh) :: kmesh
  type(scattering) :: sct
  type(runinfo)    :: info

  real(8), allocatable :: darr2_1(:,:), darr2_2(:,:), darr2_3(:,:)
  character(len=128)   :: string

  integer :: iband, ik
  integer :: kpoints, nbands

  logical :: k_dependence
  logical :: band_dependence

  integer, allocatable :: hdf5shape(:)

  ! sanity check
  ! call hdf5_read_data(ifile, "/.quantities/nkp",    kpoints)
  ! call hdf5_read_data(ifile, "/.quantities/nbands", nbands)

  if (edisp%ispin == 1) then
    call hdf5_get_shape(ifile, "step/000001/scatrate", hdf5shape)
  else
    call hdf5_get_shape(ifile, "up/step/000001/scatrate", hdf5shape)
  endif

  ! deduce from the array shape if we have a band or momentum dependence
  if (hdf5shape(1) == 1) then
    band_dependence = .false.
  else
    band_dependence = .true.
  endif
  if (hdf5shape(2) == 1) then
    k_dependence = .false.
  else
    k_dependence = .true.
  endif
  deallocate(hdf5shape)

  if (edisp%ispin == 1) then
    write(string,'("step/",I6.6,"/scatrate")') info%iStep
    call hdf5_read_data(ifile, string, darr2_1)

    write(string,'("step/",I6.6,"/qpweight")') info%iStep
    call hdf5_read_data(ifile, string, darr2_2)

    if (edisp%lBandShift) then
      write(string,'("step/",I6.6,"/bandshift")') info%iStep
      call hdf5_read_data(ifile, string, darr2_3)
    endif

    if (k_dependence .and. band_dependence) then
      sct%gam(:,:,1) = darr2_1
      sct%zqp(:,:,1) = darr2_2
      if (edisp%lBandShift) then
        edisp%band_shift(:,:,1) = darr2_3
      endif
    else if (k_dependence .and. .not. band_dependence) then
      do iband=1,edisp%nband_max
        sct%gam(iband,:,1) = darr2_1(1,:)
        sct%zqp(iband,:,1) = darr2_2(1,:)
        if (edisp%lBandShift) then
          edisp%band_shift(iband,:,1) = darr2_3(1,:)
        endif
      enddo
    else if (.not. k_dependence .and. band_dependence) then
      do ik=1,kmesh%nkp
        sct%gam(:,ik,1) = darr2_1(:,1)
        sct%zqp(:,ik,1) = darr2_2(:,1)
        if (edisp%lBandShift) then
          edisp%band_shift(:,ik,1) = darr2_3(:,1)
        endif
      enddo
    else
      sct%gam(:,:,1) = darr2_1(1,1)
      sct%zqp(:,:,1) = darr2_2(1,1)
      if (edisp%lBandShift) then
        edisp%band_shift(:,:,1) = darr2_3(1,1)
      endif
    endif

    if (allocated(darr2_1)) deallocate(darr2_1)
    if (allocated(darr2_2)) deallocate(darr2_2)
    if (allocated(darr2_3)) deallocate(darr2_3)

  else

    write(string,'("up/step/",I6.6,"/scatrate")') info%iStep
    call hdf5_read_data(ifile, string, darr2_1)

    write(string,'("up/step/",I6.6,"/qpweight")') info%iStep
    call hdf5_read_data(ifile, string, darr2_2)

    if (edisp%lBandShift) then
      write(string,'("up/step/",I6.6,"/bandshift")') info%iStep
      call hdf5_read_data(ifile, string, darr2_3)
    endif

    if (k_dependence .and. band_dependence) then
      sct%gam(:,:,1) = darr2_1
      sct%zqp(:,:,1) = darr2_2
      if (edisp%lBandShift) then
        edisp%band_shift(:,:,1) = darr2_3
      endif
    else if (k_dependence .and. .not. band_dependence) then
      do iband=1,edisp%nband_max
        sct%gam(iband,:,1) = darr2_1(1,:)
        sct%zqp(iband,:,1) = darr2_2(1,:)
        if (edisp%lBandShift) then
          edisp%band_shift(iband,:,1) = darr2_3(1,:)
        endif
      enddo
    else if (.not. k_dependence .and. band_dependence) then
      do ik=1,kmesh%nkp
        sct%gam(:,ik,1) = darr2_1(:,1)
        sct%zqp(:,ik,1) = darr2_2(:,1)
        if (edisp%lBandShift) then
          edisp%band_shift(:,ik,1) = darr2_3(:,1)
        endif
      enddo
    else
      sct%gam(:,:,1) = darr2_1(1,1)
      sct%zqp(:,:,1) = darr2_2(1,1)
      if (edisp%lBandShift) then
        edisp%band_shift(:,:,1) = darr2_3(1,1)
      endif
    endif

    write(string,'("dn/step/",I6.6,"/scatrate")') info%iStep
    call hdf5_read_data(ifile, string, darr2_1)

    write(string,'("dn/step/",I6.6,"/qpweight")') info%iStep
    call hdf5_read_data(ifile, string, darr2_2)

    if (edisp%lBandShift) then
      write(string,'("dn/step/",I6.6,"/bandshift")') info%iStep
      call hdf5_read_data(ifile, string, darr2_3)
    endif

    if (k_dependence .and. band_dependence) then
      sct%gam(:,:,2) = darr2_1
      sct%zqp(:,:,2) = darr2_2
      if (edisp%lBandShift) then
        edisp%band_shift(:,:,2) = darr2_3
      endif
    else if (k_dependence .and. .not. band_dependence) then
      do iband=1,edisp%nband_max
        sct%gam(iband,:,2) = darr2_1(1,:)
        sct%zqp(iband,:,2) = darr2_2(1,:)
        if (edisp%lBandShift) then
          edisp%band_shift(iband,:,2) = darr2_3(1,:)
        endif
      enddo
    else if (.not. k_dependence .and. band_dependence) then
      do ik=1,kmesh%nkp
        sct%gam(:,ik,2) = darr2_1(:,1)
        sct%zqp(:,ik,2) = darr2_2(:,1)
        if (edisp%lBandShift) then
          edisp%band_shift(:,ik,2) = darr2_3(:,1)
        endif
      enddo
    else
      sct%gam(:,:,2) = darr2_1(1,1)
      sct%zqp(:,:,2) = darr2_2(1,1)
      if (edisp%lBandShift) then
        edisp%band_shift(:,:,2) = darr2_3(1,1)
      endif
    endif

    if (allocated(darr2_1)) deallocate(darr2_1)
    if (allocated(darr2_2)) deallocate(darr2_2)
    if (allocated(darr2_3)) deallocate(darr2_3)

  endif

  ! apply the shift
  if (edisp%lBandShift) then
    edisp%band = edisp%band_file + edisp%band_shift
  endif

  ! add possible impurity offsets
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
    write(string,'("kPoint/",I10.10,"/moments")') info%ik
    call hdf5_read_data(ifile, string, darr3)
    edisp%Moptk(:,:,:,1) = darr3
    deallocate(darr3)
  else
    if (allocated(darr3)) deallocate(darr3)
    write(string,'("up/kPoint/",I10.10,"/moments")') info%ik
    call hdf5_read_data(ifile, string, darr3)
    edisp%Moptk(:,:,:,1) = darr3
    deallocate(darr3)
    write(string,'("dn/kPoint/",I10.10,"/moments")') info%ik
    call hdf5_read_data(ifile, string, darr3)
    edisp%Moptk(:,:,:,2) = darr3
    deallocate(darr3)
  endif

end subroutine

subroutine read_muT_hdf5(algo, temp, mu)
  implicit none
  type(algorithm)   :: algo
  type(temperature) :: temp
  real(8)           :: mu(algo%steps)

  integer(hid_t)       :: ifile
  real(8), allocatable :: mutemp(:)
  integer              :: shapemu(1)

  call hdf5_open_file(algo%old_output_file, ifile, rdonly=.true.)
  call hdf5_read_data(ifile, '.quantities/mu', mutemp)

  shapemu = shape(mutemp)

  if (.not. (algo%steps == shapemu(1))) then
    call stop_with_message(stdout, 'chemical potential array from old file does not match user input')
  endif

  mu = mutemp

  call hdf5_close_file(ifile)

end subroutine read_muT_hdf5

subroutine read_muT_text(algo, temp, mu)
  implicit none
  type(algorithm)   :: algo
  type(temperature) :: temp
  real(8)           :: mu(algo%steps)

  character(len=256) :: str_temp
  integer            :: i,j
  integer            :: lines
  integer            :: stat
  integer            :: empty, pst
  character(len=256), allocatable :: file_temp(:), file_save(:)
  character(len=1) :: cmnt = '#'
  real(8) :: fdum1, fdum2


  open(unit=10,file=trim(algo%input_mu_text),action='read',iostat=stat)
  if (stat .ne. 0) then
    call stop_with_message(stderr, 'ScatteringText Input file cannot be opened') ! send to stderr
  endif

  ! line counting
  lines=0
  read_count: do
    read(10,'(A)',END=200) str_temp ! read whole line as string, doesnt skip empty lines
    lines=lines+1
  enddo read_count


  200 continue
  rewind 10

  allocate(file_temp(lines))

  ! remove empty lines and comment strings
  empty=0
  read_temp: do i=1,lines
    read(10,'(A)') str_temp
      str_temp = trim(adjustl(str_temp))
      pst=scan(str_temp,cmnt) ! find out possible comment symbol
      if (pst .eq. 1) then ! whole line is commented
        file_temp(i) = ''
      elseif (pst .eq. 0) then ! no comment symbol found
        file_temp(i) = str_temp
      else  ! getting left side of comment
        file_temp(i) = trim(str_temp(1:pst-1))
      endif

      if (len_trim(file_temp(i)) .eq. 0) then ! filter out empty lines
        empty=empty+1
      endif
  enddo read_temp

  ! rewrite everything to a new clean string array
  allocate(file_save(lines-empty))
  j=1
  read_save: do i=1,lines
    if(len_trim(file_temp(i)) .ne. 0) then
      file_save(j)=file_temp(i)
      j=j+1
    endif
  enddo read_save
  deallocate(file_temp)

  lines=lines-empty
  close(unit=10)

  if (lines /= algo%steps) then
    call stop_with_message(stderr, 'Error: Number of provided mu values does not coincide with temp grid')
  endif
  ! now we read this string array into the according data arrays

  do i=1,algo%steps
    read(file_save(i),*) fdum1, mu(i)
  enddo

end subroutine

end module Mio
