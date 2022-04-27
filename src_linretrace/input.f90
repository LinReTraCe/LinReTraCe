module Minput
  use Mparams
  use Mtypes
  use Mauxiliary
  use hdf5
  use hdf5_wrapper
  implicit none

contains

subroutine read_preproc_energy(algo, kmesh, edisp, sct, pot, imp)
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

  call hdf5_open_file(trim(adjustl(algo%input_energies)), ifile, rdonly=.true.)

  ! mesh information
  call hdf5_read_data(ifile, "/.kmesh/nkp",         kmesh%nkp)
  call hdf5_read_data(ifile, "/.kmesh/nkx",         kmesh%nkx)
  call hdf5_read_data(ifile, "/.kmesh/nky",         kmesh%nky)
  call hdf5_read_data(ifile, "/.kmesh/nkz",         kmesh%nkz)
  call hdf5_read_data(ifile, "/.kmesh/weightsum",   kmesh%weightsum)
  ! call hdf5_read_data(ifile, "/.kmesh/weights",     kmesh%weight)
  call hdf5_read_data(ifile, "/.kmesh/multiplicity",kmesh%multiplicity)
  kmesh%minimal_weight = minval(kmesh%multiplicity) / real(kmesh%nkx*kmesh%nky*kmesh%nkz,8)

  ! band information + charge
  call hdf5_read_data(ifile, "/.bands/mu",             pot%mu_dft_file)
  call hdf5_read_data(ifile, "/.bands/charge",         edisp%nelect_file)
  call hdf5_read_data(ifile, "/.bands/energyBandMax",  edisp%nband_max)
  call hdf5_read_data(ifile, "/.bands/opticalBandMin", edisp%nbopt_min)
  call hdf5_read_data(ifile, "/.bands/opticalBandMax", edisp%nbopt_max)
  call hdf5_read_data(ifile, "/.bands/ispin",          edisp%ispin)

  ! reset electron occupation from config file if it is invalid ... earlierst point to check is here
  if (edisp%nelect_config >= kmesh%weightsum * edisp%nband_max) then ! too many electrons
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
  edisp%gapped_complete = .true.
  do is=1,edisp%ispin
    if (.not. edisp%gapped(is)) edisp%gapped_complete = .false.
  enddo

  ! set common band gap
  if (edisp%gapped_complete) then
    edisp%gap_min = minval(edisp%ene_conductionBand) - maxval(edisp%ene_valenceBand)
  else
    edisp%gap_min = 0.d0
  endif

  ! check impurities for consistency
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

  ! check scattering / zqp / energy coefficient for consistency
  if (.not. algo%lScatteringFile) then
    if (edisp%ispin==2 .and. size(sct%gamcoeff,dim=2)==1) then
      ! increase the spin axis by one
      n2shape = shape(sct%gamcoeff)
      allocate(drank2arr(n2shape(1),n2shape(2)))
      drank2arr(:,:) = sct%gamcoeff
      deallocate(sct%gamcoeff)
      allocate(sct%gamcoeff(n2shape(1),2))
      sct%gamcoeff(:,1) = drank2arr(:,1)
      sct%gamcoeff(:,2) = drank2arr(:,1)
      deallocate(drank2arr)
    endif
    if (edisp%ispin==1 .and. size(sct%gamcoeff,dim=2)==2) then
      call stop_with_message(stderr, 'Error: Provided Scattering Coefficients do not match system (one spin type)')
    endif
    if (edisp%ispin==2 .and. size(sct%zqpcoeff,dim=2)==1) then
      ! increase the spin axis by one
      n2shape = shape(sct%zqpcoeff)
      allocate(drank2arr(n2shape(2),n2shape(2)))
      drank2arr(:,:) = sct%zqpcoeff
      deallocate(sct%zqpcoeff)
      allocate(sct%zqpcoeff(n2shape(1),2))
      sct%zqpcoeff(:,1) = drank2arr(:,1)
      sct%zqpcoeff(:,2) = drank2arr(:,1)
      deallocate(drank2arr)
    endif
    if (edisp%ispin==1 .and. size(sct%zqpcoeff,dim=2)==2) then
      call stop_with_message(stderr, 'Error: Provided QuasiParticle Coefficients do not match system (one spin type)')
    endif

    if (sct%enescaling) then
      if (edisp%ispin==1 .and. size(sct%enecoeff,dim=2)==2) then
        call stop_with_message(stderr, 'Error: Provided Energy Coefficients do not match system (one spin type)')
      endif
      if (edisp%ispin==2 .and. size(sct%enecoeff,dim=2)==1) then
        ! increase the spin axis by one
        n2shape = shape(sct%enecoeff)
        allocate(drank2arr(n2shape(2),n2shape(2)))
        drank2arr(:,:) = sct%enecoeff
        deallocate(sct%enecoeff)
        allocate(sct%enecoeff(n2shape(1),2))
        sct%enecoeff(:,1) = drank2arr(:,1)
        sct%enecoeff(:,2) = drank2arr(:,1)
        deallocate(drank2arr)
      endif
      if (edisp%ispin==1 .and. size(sct%enecoeff,dim=2)==2) then
        call stop_with_message(stderr, 'Error: Provided Energy Coefficients do not match system (one spin type)')
      endif
    endif
  endif

  ! check for Bandgap inconsistency
  if (algo%lScissors) then
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

  ! transform the densities ( doping / impurities ) from cm-3 to electrons
  if (.not. algo%lNominalDoping) then
    if (algo%ldoping) then
      edisp%doping = edisp%doping * 1d-24 * kmesh%vol
    endif
    if (algo%lImpurities) then
      do iimp = 1, imp%nimp
        imp%Density(iimp) = imp%Density(iimp) * 1d-24 * kmesh%vol
      enddo
    endif
  endif

  ! set logical flags whether specific dataset exist in the energy file
  edisp%lIntraMoments  = .false.
  edisp%lFullMoments   = .false.
  edisp%lBIntraMoments = .false.
  edisp%lBFullMoments  = .false.

  if (edisp%ispin == 1) then
    if (.not. hdf5_dataset_exists(ifile, "/energies")) then
      call stop_with_message(stderr, 'Could not find energies in EnergyFile')
    endif
    if (hdf5_dataset_exists(ifile, "/momentsDiagonal")) then
      edisp%lIntraMoments = .true.
      call hdf5_get_shape(ifile, "/momentsDiagonal", irank1arr)
      edisp%iOptical = irank1arr(1) ! number of optical elements 3 6 or 9
      deallocate(irank1arr)
    endif
    if (hdf5_dataset_exists(ifile, "/momentsDiagonalBfield")) then
      edisp%lBIntraMoments = .true.
    endif
    if (hdf5_dataset_exists(ifile, "/kPoint/0000000001/moments")) then
      edisp%lFullMoments = .true.
    endif
    if (hdf5_dataset_exists(ifile, "/kPoint/0000000001/momentsBfield")) then
      edisp%lBFullMoments = .true.
    endif
  else if (edisp%ispin == 2) then
    if (.not. (hdf5_dataset_exists(ifile, "/up/energies") .and. &
               hdf5_dataset_exists(ifile, "/dn/energies"))) then
      call stop_with_message(stderr, 'Could not find energies in EnergyFile')
    endif
    if (hdf5_dataset_exists(ifile, "/up/momentsDiagonal") .and. &
        hdf5_dataset_exists(ifile, "/dn/momentsDiagonal")) then
      edisp%lIntraMoments = .true.
      call hdf5_get_shape(ifile, "/up/momentsDiagonal", irank1arr)
      edisp%iOptical = irank1arr(1) ! number of optical elements 3 6 or 9
      deallocate(irank1arr)
    endif
    if (hdf5_dataset_exists(ifile, "/up/momentsDiagonalBfield") .and. &
        hdf5_dataset_exists(ifile, "/dn/momentsDiagonalBfield")) then
      edisp%lBIntraMoments = .true.
    endif
    if (hdf5_dataset_exists(ifile, "/up/kPoint/0000000001/moments") .and. &
        hdf5_dataset_exists(ifile, "/up/kPoint/0000000001/moments")) then
      edisp%lFullMoments = .true.
    endif
    if (hdf5_dataset_exists(ifile, "/up/kPoint/0000000001/momentsBfield") .and. &
        hdf5_dataset_exists(ifile, "/up/kPoint/0000000001/momentsBfield")) then
      edisp%lBFullMoments = .true.
    endif
  else
    call stop_with_message(stderr, 'Provided number of spins in EnergyFile invalid')
  endif

  call hdf5_close_file(ifile)

end subroutine

subroutine read_energy(algo, edisp)
  implicit none
  type(algorithm)  :: algo
  type(energydisp) :: edisp

  integer(hid_t)   :: ifile
  real(8), allocatable :: drank2arr(:,:)
  real(8), allocatable :: drank3arr(:,:,:)
  complex(8), allocatable :: zrank5arr(:,:,:,:,:)

  allocate(edisp%band_file(edisp%nband_max, ikstr:ikend, edisp%ispin))
  allocate(edisp%band(edisp%nband_max, ikstr:ikend, edisp%ispin))
  if (edisp%lIntraMoments) then
    allocate(edisp%MoptDiag(edisp%iOptical, edisp%nbopt_min:edisp%nbopt_max, edisp%ispin, ikstr:ikend))
  endif
  if (edisp%lBIntraMoments) then
    allocate(edisp%MBoptdiag(3, 3, 3, edisp%nband_max, edisp%ispin, ikstr:ikend))
  endif

  call hdf5_open_file(trim(adjustl(algo%input_energies)), ifile, rdonly=.true.)

  if (edisp%ispin == 1) then
    call hdf5_read_data(ifile, "/energies", drank2arr)
    edisp%band_file(:,:,1)     = drank2arr(:,ikstr:ikend)
    deallocate(drank2arr)
    if (edisp%lBIntraMoments) then
      ! call hdf5_read_data(ifile, "/derivatives",   drank3arr)
      ! edisp%band_dk(:,:,:,1) = drank3arr
      ! deallocate(drank3arr)
      ! call hdf5_read_data(ifile, "/curvatures",  drank3arr)
      ! edisp%band_d2k(:,:,:,1) = drank3arr
      ! deallocate(drank3arr)
      call hdf5_read_data(ifile, "/momentsDiagonalBfield",   zrank5arr)
      edisp%MBoptDiag(:,:,:,:,1,:) = zrank5arr(:,:,:,:,ikstr:ikend)
      deallocate(zrank5arr)
    endif
    if (edisp%lIntraMoments) then
      call hdf5_read_data(ifile, "momentsDiagonal", drank3arr)
      edisp%MoptDiag(:,:,1,:) = drank3arr(:,:,ikstr:ikend)
      deallocate(drank3arr)
    endif
  else if (edisp%ispin == 2) then
    call hdf5_read_data(ifile, "/up/energies", drank2arr)
    edisp%band_file(:,:,1)     = drank2arr(:,ikstr:ikend)
    deallocate(drank2arr)
    call hdf5_read_data(ifile, "/dn/energies", drank2arr)
    edisp%band_file(:,:,2)     = drank2arr(:,ikstr:ikend)
    deallocate(drank2arr)

    if (edisp%lBIntraMoments) then
      ! call hdf5_read_data(ifile, "/up/derivatives",  drank3arr)
      ! edisp%band_dk(:,:,:,1) = drank3arr
      ! deallocate(drank3arr)
      ! call hdf5_read_data(ifile, "/up/curvatures",    drank3arr)
      ! edisp%band_d2k(:,:,:,1) = drank3arr
      ! deallocate(drank3arr)
      ! call hdf5_read_data(ifile, "/dn/derivatives", drank3arr)
      ! edisp%band_dk(:,:,:,2) = drank3arr
      ! deallocate(drank3arr)
      ! call hdf5_read_data(ifile, "/dn/curvatures",   drank3arr)
      ! edisp%band_d2k(:,:,:,2) = drank3arr
      ! deallocate(drank3arr)
      call hdf5_read_data(ifile, "/up/momentsDiagonalBfield",   zrank5arr)
      edisp%MBoptDiag(:,:,:,:,1,:) = zrank5arr(:,:,:,:,ikstr:ikend)
      deallocate(zrank5arr)
      call hdf5_read_data(ifile, "/dn/momentsDiagonalBfield",   zrank5arr)
      edisp%MBoptDiag(:,:,:,:,2,:) = zrank5arr(:,:,:,:,ikstr:ikend)
      deallocate(zrank5arr)
    endif
    if (edisp%lIntraMoments) then
      call hdf5_read_data(ifile, "/up/momentsDiagonal", drank3arr)
      edisp%MoptDiag(:,:,1,:) = drank3arr(:,:,ikstr:ikend)
      deallocate(drank3arr)
      call hdf5_read_data(ifile, "/dn/momentsDiagonal", drank3arr)
      edisp%MoptDiag(:,:,2,:) = drank3arr(:,:,ikstr:ikend)
      deallocate(drank3arr)
    endif
  endif

  call hdf5_close_file(ifile)

end subroutine

subroutine read_preproc_scattering_hdf5(algo, kmesh, edisp, sct, pot, temp)
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
  allocate(sct%gam(edisp%nband_max, ikstr:ikend, edisp%ispin))
  allocate(sct%zqp(edisp%nband_max, ikstr:ikend, edisp%ispin))

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
    allocate(edisp%band_shift(edisp%nband_max, ikstr:ikend, edisp%iSpin))
  endif

  pot%mumin = minval(pot%MM)
  pot%mumax = maxval(pot%MM)

  call hdf5_close_file(ifile)

end subroutine

subroutine read_preproc_scattering_text(algo, kmesh, edisp, sct, pot, temp)
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
  allocate(sct%gam(edisp%nband_max, ikstr:ikend, edisp%ispin))
  allocate(sct%zqp(edisp%nband_max, ikstr:ikend, edisp%ispin))

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

subroutine read_scattering_hdf5(ifile, edisp, kmesh, sct, info)
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
      sct%gam(:,:,1) = darr2_1(:,ikstr:ikend)
      sct%zqp(:,:,1) = darr2_2(:,ikstr:ikend)
      if (edisp%lBandShift) then
        edisp%band_shift(:,:,1) = darr2_3(:,ikstr:ikend)
      endif
    else if (k_dependence .and. .not. band_dependence) then
      do iband=1,edisp%nband_max
        sct%gam(iband,:,1) = darr2_1(1,ikstr:ikend)
        sct%zqp(iband,:,1) = darr2_2(1,ikstr:ikend)
        if (edisp%lBandShift) then
          edisp%band_shift(iband,:,1) = darr2_3(1,ikstr:ikend)
        endif
      enddo
    else if (.not. k_dependence .and. band_dependence) then
      do ik=ikstr,ikend
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
      sct%gam(:,:,1) = darr2_1(:,ikstr:ikend)
      sct%zqp(:,:,1) = darr2_2(:,ikstr:ikend)
      if (edisp%lBandShift) then
        edisp%band_shift(:,:,1) = darr2_3(:,ikstr:ikend)
      endif
    else if (k_dependence .and. .not. band_dependence) then
      do iband=1,edisp%nband_max
        sct%gam(iband,:,1) = darr2_1(1,ikstr:ikend)
        sct%zqp(iband,:,1) = darr2_2(1,ikstr:ikend)
        if (edisp%lBandShift) then
          edisp%band_shift(iband,:,1) = darr2_3(1,ikstr:ikend)
        endif
      enddo
    else if (.not. k_dependence .and. band_dependence) then
      do ik=ikstr,ikend
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
      sct%gam(:,:,2) = darr2_1(:,ikstr:ikend)
      sct%zqp(:,:,2) = darr2_2(:,ikstr:ikend)
      if (edisp%lBandShift) then
        edisp%band_shift(:,:,2) = darr2_3(:,ikstr:ikend)
      endif
    else if (k_dependence .and. .not. band_dependence) then
      do iband=1,edisp%nband_max
        sct%gam(iband,:,2) = darr2_1(1,ikstr:ikend)
        sct%zqp(iband,:,2) = darr2_2(1,ikstr:ikend)
        if (edisp%lBandShift) then
          edisp%band_shift(iband,:,2) = darr2_3(1,ikstr:ikend)
        endif
      enddo
    else if (.not. k_dependence .and. band_dependence) then
      do ik=ikstr,ikend
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

subroutine read_full_optical_elements(ifile, edisp, sct, info)
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

subroutine read_full_magnetic_optical_elements(ifile, edisp, sct, info)
  ! currently unused subroutine which will be used to load in the full magnetic
  ! field optical elements
  implicit none
  integer(hid_t)   :: ifile
  type(energydisp) :: edisp
  type(scattering) :: sct
  type(runinfo)    :: info

  complex(8), allocatable :: zarr5(:,:,:,:,:)
  character(len=128)   :: string

  if (edisp%ispin == 1) then
    if (allocated(zarr5)) deallocate(zarr5)
    write(string,'("kPoint/",I10.10,"/momentsBfield")') info%ik
    call hdf5_read_data(ifile, string, zarr5)
    edisp%MBoptk(:,:,:,:,:,1) = zarr5
    deallocate(zarr5)
  else
    if (allocated(zarr5)) deallocate(zarr5)
    write(string,'("up/kPoint/",I10.10,"/momentsBfield")') info%ik
    call hdf5_read_data(ifile, string, zarr5)
    edisp%MBoptk(:,:,:,:,:,1) = zarr5
    deallocate(zarr5)
    write(string,'("dn/kPoint/",I10.10,"/momentsBfield")') info%ik
    call hdf5_read_data(ifile, string, zarr5)
    edisp%MBoptk(:,:,:,:,:,2) = zarr5
    deallocate(zarr5)
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

  call hdf5_open_file(algo%input_mu_hdf5, ifile, rdonly=.true.)
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

end module Minput
