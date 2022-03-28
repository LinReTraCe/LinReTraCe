program main
  use Mparams
  use Maux
  use Mroot
  use Mtypes
  use Mmpi_org
  use Mconfig
  use Minput
  use Moutput
  use Mresponse
  use hdf5_wrapper
  use hdf5
  implicit none

  type(algorithm)   :: algo
  type(kpointmesh)  :: kmesh   ! contains k-point mesh specifiers and logical switches on how to get the mesh from
  type(energydisp)  :: edisp   ! contains the band dispersion energy and the optical matrix elements
  type(dosgrid)     :: dos     ! DOS, integrated DOS, Fermi level
  type(scattering)  :: sct     ! scattering rates and quasi particle weights
  type(temperature) :: temp    ! temperature quantities
  type(potential)   :: pot     ! chemical potential quantities
  type(runinfo)     :: info    ! runtime information for the calculation routines, temps, betas, etc.
  type(impurity)    :: imp

  type(response_dp) :: resp_intra
  type(response_dp) :: resp_intra_Boltzmann
  type(response_dp) :: resp_inter
  type(response_dp) :: resp_inter_Boltzmann

  type(response_qp) :: qresp_intra
  type(response_qp) :: qresp_intra_Boltzmann
  type(response_qp) :: qresp_inter
  type(response_qp) :: qresp_inter_Boltzmann

  integer(hid_t)    :: ifile_scatter_hdf5
  integer(hid_t)    :: ifile_energy
  integer(hid_t)    :: ifile_output


  integer :: is, ig, iT, ik, iimp, imu
  integer :: niitact
  integer :: iStart, iEnd, iStep
  real(8) :: ndevact
  real(16):: ndevactQ

  character(len=128) :: string
  real(8) :: tstart, tfinish, timings(5)
  real(8) :: time
  real(8), allocatable :: gap_file(:)
  logical :: gapped_file

  ! quantities saved on the Temperature grid
  ! and derived quantities
  real(8), allocatable :: energy(:)    ! total energy
  real(8), allocatable :: carrier(:)   ! carrier concentration
  real(8), allocatable :: electrons(:) ! thermally activated electrons w.r.t chemical potential
  real(8), allocatable :: holes(:)     ! thermally activated holes w.r.t chemical potential
  real(8), allocatable :: imp_contribution(:)

  complex(8), allocatable  :: PolyGamma(:,:,:,:)
  complex(16), allocatable :: PolyGammaQ(:,:,:,:)


  ! work arrays
  real(8), allocatable    :: darr1(:)
  real(8), allocatable    :: darr2(:,:)
  real(8), allocatable    :: darr3(:,:,:)
  real(8), allocatable    :: darr4(:,:,:,:)

  complex(8), allocatable :: zarr1(:)
  complex(8), allocatable :: zarr2(:,:)
  complex(8), allocatable :: zarr3(:,:,:)
  complex(8), allocatable :: zarr4(:,:,:,:)

  call mpi_initialize()
  if (myid.eq.master) call main_greeting(stdout)

#ifdef MPI
  call mpi_barrier(mpi_comm_world, mpierr)
#endif

  call read_config(algo, edisp, sct, temp, pot, imp)
  call check_files(algo)

  timings = 0.d0        ! reset timings
  call cpu_time(tstart) ! start timer

  call hdf5_init()
  ! read the energies, derivatives, diagonal optical elements
  ! + crosscheck config input (mainly spin dependencies)
  call read_preproc_energy_data(algo, kmesh, edisp, sct, pot, imp)

  ! quick checks if run-options are in agreement with provided data
  if (algo%lBfield .and. .not. edisp%lDerivatives) then
    call stop_with_message(stderr, 'Energy derivatives required for Bfield quantities')
  endif

  if (algo%lInterbandQuantities .and. .not. edisp%lFullMoments) then
    call stop_with_message(stderr, 'Full optical elements required for Interband quantities')
  endif

  if (.not. algo%lIntrabandQuantities .and. .not. algo%lInterbandQuantities) then
    ! I will keep this running here
    ! One may just want to look at the chemical potential // total energy
    call log_master(stdout, 'WARNING: Neither Intra nor Interband responses will be calculated')
  endif

  if (algo%ldoping) then
    edisp%doping = edisp%doping * 1.d-24 * kmesh%vol ! cm-3 -> AA-3 -> filling
  endif

  ! distribute k-grid
  call mpi_genkstep(kmesh%nkp)

#ifdef MPI
  if (algo%ldebug .and. (index(algo%dbgstr,"Mpi") .ne. 0)) then
     write(stdout,*) "MPI: myid: ", myid, "ikstr: ", ikstr, "ikend: ", ikend
     call mpi_barrier(mpi_comm_world, mpierr)
  endif
#endif

  allocate(gap_file(edisp%ispin))
  gap_file = edisp%gap
  gapped_file = edisp%gapped_complete

  ! apply scissors here
  ! this has to be saved in the band_file variable
  ! since it is reused if we apply possible band_shifts
  if (algo%lScissors) then
    do is=1,edisp%ispin
      edisp%band_file(edisp%conductionBand(is):, :, is) = &
      edisp%band_file(edisp%conductionBand(is):, :, is) + edisp%scissors(is)
    enddo
  endif
  edisp%band = edisp%band_file

  if (algo%lRedoMudft) then
    ! redo the DFT mu calculation -  if we apply a scissor operator or change the electron occupation
    ! redefine mu and adjust flags / gap sizes etc
    call find_mu_DFT(edisp,kmesh,pot)
  else
    ! use the value from the energy file
    pot%mu_dft = pot%mu_dft_file
  endif

  if (algo%muFermi .and. edisp%gapped_complete) then
    algo%rootMethod = 3
  else
    algo%rootMethod = 2
  endif

  ! debug strings to overwrite the standard behavior
  if (algo%ldebug) then
    if (index(algo%dbgstr,"Bisection") .ne. 0) then
      algo%rootMethod = 3
    else if (index(algo%dbgstr,"Ridders") .ne. 0) then
      algo%rootMethod = 2
    else if (index(algo%dbgstr,"Secant") .ne. 0) then
      algo%rootMethod = 0
    endif
  endif

  if (algo%lImpurities) then
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
  endif

  if (algo%lTMODE) then
    ! construct temperature grid
    if (algo%lScatteringFile) then
      ! at this point we just extract the temperature ranges
      ! and check wether bandshifts exist
      call read_preproc_scattering_data_hdf5(algo, kmesh, edisp, sct, pot, temp)
    else if (algo%lScatteringText) then
      call read_preproc_scattering_data_text(algo, kmesh, edisp, sct, pot, temp)
    else
      if (algo%steps .gt. 1) then
        if (temp%tlogarithmic) then
          temp%dT = (temp%Tmax/temp%Tmin) ** (1.d0 / dble(algo%steps-1)) ! logarithmic step
        else
          temp%dT = (temp%Tmax-temp%Tmin)/dble(algo%steps-1) ! linear step
        endif
      else
        temp%dT = 0 ! 0/0 ...
      endif
      allocate(temp%TT(algo%steps))
      allocate(temp%BB(algo%steps))
      allocate(sct%gam(edisp%nband_max, kmesh%nkp, edisp%ispin))
      allocate(sct%zqp(edisp%nband_max, kmesh%nkp, edisp%ispin))

      ! define Temperature grid
      if (temp%tlogarithmic) then
        do iT=1,algo%steps
           temp%TT(iT)=temp%Tmin * temp%dT**(real(iT-1,8))
        enddo
      else
        do iT=1,algo%steps
           temp%TT(iT)=real(iT-1,8)*temp%dT+temp%Tmin
        enddo
      endif
      temp%BB = 1.d0/(temp%TT * kB)
    endif

    allocate(energy(algo%steps))
    allocate(carrier(algo%steps))
    allocate(electrons(algo%steps))
    allocate(holes(algo%steps))
    allocate(imp_contribution(algo%steps))
    energy = 0.d0
    electrons = 0.d0
    imp_contribution = 0.d0
    holes = 0.d0
    if (.not. allocated(pot%MM))  allocate(pot%MM(algo%steps))
    if (.not. allocated(pot%QMM)) allocate(pot%QMM(algo%steps))
    allocate(pot%occ(algo%steps))
    pot%occ = 0.d0

    if (algo%muSearch) then
      pot%MM  = pot%mu_dft ! initialize MM array with DFT mu (energy file)
      pot%QMM = pot%mu_dft
    else
      pot%MM  = pot%mu_config ! use fixed config mu for all points (config file)
      pot%QMM = pot%mu_config
    endif

    if (algo%lOldmu) then
      call read_muT_hdf5(algo, temp, pot%MM) ! this is saved in double
      pot%QMM = pot%MM
    endif

    if (algo%lOldmuText) then
      call read_muT_text(algo, temp, pot%MM)   ! this usually comes from lprint -> double
      pot%QMM = pot%MM
    endif
  endif

  if (algo%lMUMODE) then
    if (algo%lScatteringFile) then
      call read_preproc_scattering_data_hdf5(algo, kmesh, edisp, sct, pot, temp)
    else
      ! construct mu grid
      allocate(temp%TT(algo%steps))
      allocate(temp%BB(algo%steps))
      allocate(sct%gam(edisp%nband_max, kmesh%nkp, edisp%ispin))
      allocate(sct%zqp(edisp%nband_max, kmesh%nkp, edisp%ispin))
      allocate(pot%MM(algo%steps))
      allocate(pot%QMM(algo%steps))
      temp%TT = temp%temp_config
      temp%BB = 1.d0/(temp%temp_config * kB)

      ! shift
      pot%MuMin = pot%MuMin + pot%mu_dft
      pot%MuMax = pot%MuMax + pot%mu_dft
      ! define Chemical potential grid
      if (algo%steps .gt. 1) then
        if (pot%mlogarithmic) then
          pot%dMu= (pot%Mumax/pot%Mumin) ** (1.q0 / real((algo%steps -1),16))
        else
          pot%dMu= (pot%Mumax-pot%Mumin)/real((algo%steps-1),16)
        endif
      else
        pot%dMu = 0.q0
      endif

      if (pot%mlogarithmic) then
        do imu=1,algo%steps
           pot%QMM(imu)=pot%Mumin * pot%dMu**(real(imu-1,16))
        enddo
      else
        do imu=1,algo%steps
           pot%QMM(imu)=real(imu-1,16)*pot%dMu+pot%MuMin
        enddo
      endif
      pot%MM = real(pot%QMM,8)
    endif


    allocate(energy(algo%steps))
    allocate(electrons(algo%steps))
    allocate(carrier(algo%steps))
    allocate(holes(algo%steps))
    allocate(imp_contribution(algo%steps))
    allocate(pot%occ(algo%steps))
    carrier = 0.d0
    energy = 0.d0
    electrons = 0.d0
    holes = 0.d0
    imp_contribution = 0.d0
    pot%occ = 0.d0
  endif

  if (.not. edisp%lBandShift) then
    deallocate(edisp%band_file)
  endif


  ! for the responses we need psi_1, psi_2 and psi_3
  if (algo%lIntrabandQuantities .or. algo%lInterbandQuantities) then
    if (algo%lQuad) then
      allocate(PolyGammaQ(3, edisp%nbopt_min:edisp%nbopt_max, ikstr:ikend, edisp%ispin))
    else
      allocate(PolyGamma(3, edisp%nbopt_min:edisp%nbopt_max, ikstr:ikend, edisp%ispin))
    endif
  endif

  ! allocate the arrays once outside of the main (temperature) loop
  if (algo%lQuad) then
    if (algo%lIntrabandQuantities) then
      call allocate_response(algo, edisp, temp, qresp_intra)
      if (algo%lBoltzmann) then
        call allocate_response(algo, edisp, temp, qresp_intra_Boltzmann)
      endif
    endif
    if (algo%lInterbandquantities) then
      call allocate_response(algo, edisp, temp, qresp_inter)
      if (algo%lBoltzmann) then
        call allocate_response(algo, edisp, temp, qresp_inter_Boltzmann)
      endif
    endif
  else
    if (algo%lIntrabandQuantities) then
      call allocate_response(algo, edisp, temp, resp_intra)
      if (algo%lBoltzmann) then
        call allocate_response(algo, edisp, temp, resp_intra_Boltzmann)
      endif
    endif
    if (algo%lInterbandquantities) then
      call allocate_response(algo, edisp, temp, resp_inter)
      if (algo%lBoltzmann) then
        call allocate_response(algo, edisp, temp, resp_inter_Boltzmann)
      endif
    endif
  endif

  call hdf5_open_file(algo%input_energies,   ifile_energy,  rdonly=.true.)

  if (algo%lScatteringFile) then
    call hdf5_open_file(algo%input_scattering_hdf5, ifile_scatter_hdf5, rdonly=.true.)
  endif


  ! load the optical elements for each processes k-range
  if (algo%lInterbandQuantities) then
    allocate(edisp%Mopt(edisp%iOptical,edisp%nbopt_min:edisp%nbopt_max, &
                                       edisp%nbopt_min:edisp%nbopt_max, edisp%ispin, ikstr:ikend))
    do ik = ikstr,ikend
      info%ik = ik ! save into the runinformation datatype
      call read_optical_elements(ifile_energy, edisp, sct, info)  ! load them into edisp%Moptk
      edisp%Mopt(:,:,:,:,ik) = edisp%Moptk
    enddo
    deallocate(edisp%Moptk)
  endif

  if (myid .eq. master) then
    write(stdout,*) 'ENERGY INPUT'
    write(stdout,*) '  energy-file:      ', trim(algo%input_energies)
    write(stdout,*) '  mu (file):        ', pot%mu_dft_file
    write(stdout,*) '  electrons (file): ', edisp%nelect_file
    write(stdout,*) '  gapped (file):              ', gapped_file
    if (gapped_file) then
    write(stdout,*) '  gap (file):       ', gap_file
    endif
    write(stdout,*) '  dimensions:       ', kmesh%ndim
    write(stdout,*) '  k-Points:         ', kmesh%nkp
    write(stdout,*) '  bands:            ', edisp%nband_max
    write(stdout,*) '  spins:            ', edisp%ispin
    if (algo%lRedoMudft .or. algo%ldoping) then
      write(stdout,*)
      write(stdout,*) 'ENERGY ADJUSTMENTS'
      if (algo%ldoping) then
      write(stdout,*) '  doping:            ', edisp%doping
      endif
      if (edisp%nelect_config > 0.d0) then
      write(stdout,*) '  new electrons:     ', edisp%nelect
      endif
      if (edisp%gapped_complete .and. (int(edisp%nelect) /= edisp%nelect)) then
      write(stdout,*)
      write(stdout,*) '  WARNING: Detected gapped system without integer filling. Correct input?'
      write(stdout,*)
      endif
      if(edisp%gapped_complete .neqv. gapped_file) then
      write(stdout,*)
      write(stdout,*) '  WARNING - gapped :',  gapped_file, ' -> ', edisp%gapped_complete
      write(stdout,*)
      endif
      if (algo%lscissors) then
      write(stdout,*) '  new gap:           ', edisp%gap
      endif
      if (algo%lRedoMudft) then
      write(stdout,*) '  new mu:            ', pot%mu_dft
      endif
    endif
    if (algo%lTMODE) then
      write(stdout,*)
      write(stdout,*) 'CHEMICAL POTENTIAL'
      if (algo%muSearch) then
        if (algo%muFermi) then
          write(stdout,*) '  via: Fermi function'
        else
          write(stdout,*) '  via: Digamma function'
        endif
      else
        if (algo%lOldmu) then
          write(stdout,*) '  from file: ', trim(adjustl(algo%old_output_file))
        else
          write(stdout,*) '  constant: ', pot%MM(1)
        endif
      endif
    endif
    if (algo%lImpurities) then
      write(stdout,*) '  impurity levels: ', imp%nimp
      write(stdout,*) '    _________________________________________________________________'
      write(stdout,*) '    imp  dopant      density      energy [eV]  degeneracy  width [eV]'
      do iimp = 1,imp%nimp
        if (int(imp%Dopant(iimp)) == 1) then
          string = 'donor'
        else
          string = 'acceptor'
        endif
        write(stdout,'(5X,I1,4X,A8,E14.4,E13.4,I4,7X,E14.4)') &
            iimp, string, imp%Density(iimp), &
            imp%Energy(iimp), int(imp%Degeneracy(iimp)), imp%Bandwidth(iimp)
      enddo
    endif
    write(stdout,*)
    write(stdout,*) 'SCATTERING'
    if (algo%lScatteringFile) then
      write(stdout,*) '  HDF5 scattering-file: ', trim(algo%input_scattering_hdf5)
      write(stdout,*) '  additional impurity offset: ', sct%gamimp
    else if (algo%lScatteringText) then
      write(stdout,*) '  TEXT scattering-file: ', trim(algo%input_scattering_text)
      write(stdout,*) '  additional impurity offset: ', sct%gamimp
    else
      if (edisp%ispin == 1) then
        write(stdout,*) '  scattering coefficients: ', sct%gamcoeff(1,:)
        write(stdout,*) '  quasi particle weight coefficients: ', sct%zqpcoeff(1,:)
      else
        write(stdout,*) '  UP scattering coefficients: ', sct%gamcoeff(1,:)
        write(stdout,*) '  DN scattering coefficients: ', sct%gamcoeff(2,:)
        write(stdout,*) '  UP quasi particle weight coefficients: ', sct%zqpcoeff(1,:)
        write(stdout,*) '  DN quasi particle weight coefficients: ', sct%zqpcoeff(2,:)
      endif
    endif
    write(stdout,*)
    write(stdout,*) 'OUTPUT'
    write(stdout,*) '  output-file: ', trim(algo%output_file)
    select case (algo%fulloutput)
      case (0)
        write(stdout,*) '  full output: false'
      case (1)
        write(stdout,*) '  full output: full dependence'
      case (2)
        write(stdout,*) '  full output: momentum summation'
      case (3)
        write(stdout,*) '  full output: band summation'
    end select
    write(stdout,*)
    if (algo%lTMODE) then
    write(stdout,*) 'TEMPERATURE MODE'
    write(stdout,*) '  Temperature range:'
    write(stdout,*) '  Tmin: ', temp%Tmin
    write(stdout,*) '  Tmax: ', temp%Tmax
    write(stdout,*) '  Temperature points: ', algo%steps
    else if (algo%lMUMODE) then
    write(stdout,*) 'MU MODE'
    write(stdout,*) '  Temperature: ', temp%TT(1)
    write(stdout,*) '  Chemical Potential range:'
    write(stdout,*) '  Mumin: ', pot%Mumin
    write(stdout,*) '  Mumax: ', pot%Mumax
    write(stdout,*) '  Chemical Potential points:   ', algo%steps
    endif
    write(stdout,*)
    write(stdout,*) 'TRANSPORT RESPONSES'
    write(stdout,*) '  Intraband      : ', algo%lIntrabandQuantities
    write(stdout,*) '  Interband      : ', algo%lInterbandQuantities
    write(stdout,*) '  Boltzmann      : ', algo%lBoltzmann
    write(stdout,*) '  Magnetic field : ', algo%lBfield
    write(stdout,*) '  Quad Precision : ', algo%lQuad
    write(stdout,*)
    if (algo%lDebug) then
      write(stdout,*) 'DEBUG MODE'
      write(stdout,*) '  DEBUG STRING: "', trim(adjustl(algo%dbgstr)),'"'
      write(stdout,*)
    endif
    write(stdout,*) 'Starting calculation...'
    write(stdout,*) '____________________________________________________________________________'
    write(stdout,*) 'Temperature[K], invTemperature[1/eV], chemicalPotential[eV], rootFindingSteps'
  endif

  if (myid .eq. master) then
    call hdf5_create_file(algo%output_file)
    call output_auxiliary(algo, info, pot, temp, kmesh, edisp, sct, imp)
  endif


  call cpu_time(tfinish)
  timings(1) = timings(1) + (tfinish - tstart)
  tstart = tfinish

#ifdef MPI
  call mpi_barrier(mpi_comm_world, mpierr)
#endif


  ! MAIN LOOP
  if (algo%lTMODE  .and. .not. (algo%lDebug .and. (index(algo%dbgstr,"LoopReverse") .ne. 0)) .or. &
     (algo%lMUMODE .and. algo%lDebug .and. (index(algo%dbgstr,"LoopReverse") .ne. 0))) then
    iStart = algo%steps
    iEnd   = 1
    algo%step_dir  = -1
  else
    iStart = 1
    iEnd   = algo%steps
    algo%step_dir  = +1
  endif

  do iStep=iStart,iEnd,algo%step_dir
    ! run time information
    info%iStep   = iStep
    info%Temp    = temp%TT(iStep)
    info%beta    = 1.d0/(kB*info%Temp)
    info%beta2p  = info%beta/(2.d0*pi)

    info%TempQ   = real(info%Temp,16)
    info%betaQ   = 1.q0/(kB*info%TempQ)
    info%beta2pQ = info%betaQ/(2.q0*piQ)

    ! define scattering rates and quasi particle weights
    ! for the current temperature
    if (algo%lScatteringFile) then
      ! read in the scattering data for the current temperature
      ! scattering rates; quasi-particle weights and possible band-shifts
      call read_scattering_data_hdf5(ifile_scatter_hdf5, edisp, kmesh, sct, info)
      ! in case we have band shifts we might open or close the DFT gap
      ! for that reason : redo the mudft calculation and set the flags for each step
      ! also change the algorithm if necessary to speed things up / make things more robust
      ! but only when there is no debugstring overwriting it
      if (algo%lTMODE .and. edisp%lBandShift) then
        call find_mu_DFT(edisp,kmesh,pot)
        if (.not. (algo%ldebug .and. ((index(algo%dbgstr,"Bisection") .ne. 0) .or. &
                                      (index(algo%dbgstr,"Ridders") .ne. 0) .or. &
                                      (index(algo%dbgstr,"Secant") .ne. 0)))) then
          if (algo%muFermi .and. edisp%gapped_complete) then
            algo%rootMethod = 3
          else
            algo%rootMethod = 2
          endif
        endif
      endif
    else if (algo%lTMODE .and. algo%lScatteringText) then
      do is=1,edisp%ispin
        sct%gam(:,:,is) = sct%gamtext(iStep,is)
        sct%zqp(:,:,is) = sct%zqptext(iStep,is)
        if (sct%zqp(1,1,is) > 1.d0) then ! since its a constant array
          call log_master(stdout, 'ERROR: Zqp is bigger than 1 ... truncating to 1')
          sct%zqp(:,:,is) = 1.d0
        endif
      enddo
      sct%gam = sct%zqp * sct%gam
    else ! this is entered for both the MuMode and the TempMode
      ! here we are wasting memory however
      ! otherwise we would have to write every single response twice
      sct%gam = 0.d0
      sct%zqp = 0.d0
      do is=1,edisp%ispin
        do ig=1,size(sct%gamcoeff,dim=2)
           sct%gam(:,:,is) = sct%gam(:,:,is) + sct%gamcoeff(is,ig)*(temp%TT(iStep)**(ig-1))
        enddo
        do ig=1,size(sct%zqpcoeff,dim=2)
           sct%zqp(:,:,is) = sct%zqp(:,:,is) + sct%zqpcoeff(is,ig)*(temp%TT(iStep)**(ig-1))
        enddo
        if (sct%zqp(1,1,is) > 1.d0) then ! since its a constant array
          call log_master(stdout, 'WARNING: Zqp is bigger than 1 ... truncating to 1')
          sct%zqp(:,:,is) = 1.d0
        endif
      enddo
      sct%gam = sct%zqp * sct%gam  ! convention we use
    endif


    ! root finding (mu)
    niitact = 0
    if (algo%lTMODE .and. algo%muSearch) then

      ! initialize the new chemical potential
      if (iStep /= iStart) then
        pot%QMM(iStep) = pot%QMM(iStep-algo%step_dir) ! we initialize it to the value from the previous iteration
                                         ! this does nothing if mu is constant
      endif

      call cpu_time(tstart)

      if (edisp%gapped_complete) then
        call find_mu(pot%QMM(iStep),ndevVQ,ndevactQ,niitact, edisp, sct, kmesh, imp, algo, info)
      else
        call find_mu(pot%QMM(iStep),ndevQ,ndevactQ,niitact, edisp, sct, kmesh, imp, algo, info)
      endif

      if (niitact > niitQ) then
        if (iStep == iStart) then
          call stop_with_message(stderr, 'ERROR: Cannot determine mu')
        else if (algo%step_dir == (-1) .and. iStep <= iStart-2) then
          call log_master(stdout, 'Warning: mu determination aborted, using dMu/dT')
          pot%QMM(iStep) = pot%QMM(iStep+1) + (pot%QMM(iStep+2)-pot%QMM(iStep+1))/(temp%TT(iStep+2)-temp%TT(iStep+1)) * &
                                      (temp%TT(iStep)-temp%TT(iStep+1))
        else
          call log_master(stdout, 'Warning: mu determination aborted, mu from previous step')
          pot%QMM(iStep) = pot%QMM(iStep-algo%step_dir)
        endif
      endif

      call cpu_time(tfinish)
      timings(2) = timings(2) + (tfinish - tstart)
      tstart = tfinish

      pot%MM(iStep) = real(pot%QMM(iStep),8)
    endif

    ! after this step the mu stays for the evaluation -> save it
    ! all the response routines then access it via the runinfo data structure
    info%mu  = pot%MM(iStep)
    info%muQ = pot%QMM(iStep)

    ! calculating total energy according to the occuption above
    if (algo%muFermi) then
      call calc_total_energy_fermi(energy(iStep), edisp, sct, kmesh, imp, algo, info)
      call calc_elecholes_fermi(electrons(iStep), holes(iStep), edisp, sct, kmesh, imp, algo, info)
    else
      call calc_total_energy_digamma(energy(iStep), edisp, sct, kmesh, imp, algo, info)
      call calc_elecholes_digamma(electrons(iStep), holes(iStep), edisp, sct, kmesh, imp, algo, info)
    endif

    if (algo%lTMODE .and. algo%lImpurities) then
      call occ_impurity(ndevactQ, pot%QMM(iStep), imp, info)
      imp_contribution(iStep) = real(ndevactQ, 8)
    endif

    ! calculates the difference to the demanded electron number
    call ndeviation_Q(pot%QMM(iStep), ndevactQ, edisp, sct, kmesh, imp, algo, info)
    ! deviation rescaled to volume is the carrier concentration
    carrier(iStep) = -real(ndevactQ / kmesh%vol * 1q24, 8)! 1 / A**3 -> 1 / cm**3
    ! actual occupation is then the difference to the charge neutrality
    pot%occ(iStep) = edisp%nelect - real(ndevactQ,8)

    if (myid.eq.master) then
      write(stdout,'(3X,2F15.7,F18.12,I7)') info%Temp, info%beta, pot%MM(iStep), niitact
    endif

    ! calculate the polygamma function (1...3)
    ! for all optical bands, spins and each core's kpoints
    ! once and use it later for all the different response types
    if (algo%lIntraBandQuantities .or. algo%lInterBandQuantities) then
      if (algo%lQuad) then
        call calc_polygamma(PolyGammaQ, edisp, sct, kmesh, algo, info)
        call initialize_response(algo, qresp_intra)
        if (algo%lInterBandQuantities) then
          call initialize_response(algo, qresp_inter)
        endif
        if (algo%lBoltzmann) then
          call initialize_response(algo, qresp_intra_Boltzmann)
          if (algo%lInterBandQuantities) then
            call initialize_response(algo, qresp_inter_Boltzmann)
          endif
        endif
      else
        call calc_polygamma(PolyGamma, edisp, sct, kmesh, algo, info)
        ! initialize the already allocated arrays to 0
        if (algo%lIntraBandQuantities) then
          call initialize_response(algo, resp_intra)
          if(algo%lBoltzmann) then
            call initialize_response(algo, resp_intra_Boltzmann)
          endif
        endif
        if (algo%lInterBandQuantities) then
          call initialize_response(algo, resp_inter)
          if (algo%lBoltzmann) then
            call initialize_response(algo, resp_inter_Boltzmann)
          endif
        endif
      endif
    endif

    call cpu_time(tfinish)
    timings(3) = timings(3) + (tfinish - tstart)
    tstart = tfinish

    ! do the k-point loop and calculate the response
    do ik = ikstr,ikend
      info%ik = ik ! save into the runinformation datatype
      if (algo%lQuad) then
        ! quad precision routines
        if (algo%lIntrabandQuantities) then
          call response_intra_km_Q(qresp_intra, PolyGammaQ, edisp, sct, kmesh, algo, info)
          if (algo%lBoltzmann) then
            call response_intra_Boltzmann_km_Q(qresp_intra_Boltzmann, edisp, sct, kmesh, algo, info)
          endif
        endif
        if (algo%lInterBandQuantities) then
          call response_inter_km_Q(qresp_inter, PolyGammaQ, edisp, sct, kmesh, algo, info)
          if (algo%lBoltzmann) then
            call response_inter_Boltzmann_km_Q(qresp_inter_Boltzmann, edisp, sct, kmesh, algo, info)
          endif
        endif
      else
        ! double precision routines
        if (algo%lIntrabandQuantities) then
          call response_intra_km(resp_intra,  PolyGamma, edisp, sct, kmesh, algo, info)
          if (algo%lBoltzmann) then
            call response_intra_Boltzmann_km(resp_intra_Boltzmann, edisp, sct, kmesh, algo, info)
          endif
        endif

        if (algo%lInterbandquantities) then
          call response_inter_km(resp_inter, PolyGamma, edisp, sct, kmesh, algo, info)
          if (algo%lBoltzmann) then
            call response_inter_Boltzmann_km(resp_inter_Boltzmann, edisp, sct, kmesh, algo, info)
          endif
        endif
      endif
    enddo

    call cpu_time(tfinish)
    timings(4) = timings(4) + (tfinish - tstart)
    tstart = tfinish

    ! output the response
    ! this subroutines include the summations necessary
    if (algo%lQuad) then
      if (algo%lIntrabandQuantities) then
        call output_response_Q(qresp_intra, "intra", edisp, algo, info, temp, kmesh)
        if (algo%lBoltzmann) then
          call output_response_Q(qresp_intra_Boltzmann, "intraBoltzmann", edisp, algo, info, temp, kmesh)
        endif
      endif
      if (algo%lInterbandQuantities) then
        call output_response_Q(qresp_inter, "inter", edisp, algo, info, temp, kmesh, .false.)
        if (algo%lBoltzmann) then
          call output_response_Q(qresp_inter_Boltzmann, "interBoltzmann", edisp, algo, info, temp, kmesh, .false.)
        endif
      endif
    else
      if (algo%lIntrabandQuantities) then
        call output_response_D(resp_intra, "intra", edisp, algo, info, temp, kmesh)
        if (algo%lBoltzmann) then
          call output_response_D(resp_intra_Boltzmann, "intraBoltzmann", edisp, algo, info, temp, kmesh)
        endif
      endif
      if (algo%lInterbandQuantities) then
        ! here we don't have the Bfield quantities ... no optical elements yet
        call output_response_D(resp_inter, "inter", edisp, algo, info, temp, kmesh, .false.)
        if (algo%lBoltzmann) then
          call output_response_D(resp_inter_Boltzmann, "interBoltzmann", edisp, algo, info, temp, kmesh, .false.)
        endif
      endif
    endif

    ! output the renormalized energies defined by Z*(ek - mu)
    ! for each temperature point
    if (myid.eq.master .and. algo%lEnergyOutput) then
      call output_energies(algo, edisp, kmesh, sct, info)
    endif

    call cpu_time(tfinish)
    timings(5) = timings(5) + (tfinish - tstart)
    tstart = tfinish

  enddo ! end of the outer temperature / chemical potential loop

  if (allocated(edisp%Mopt)) deallocate(edisp%Mopt)

  if (myid.eq.master) then
    call hdf5_open_file(algo%output_file, ifile_output)
    call hdf5_write_data(ifile_output, '.quantities/mu', pot%MM)
    call hdf5_write_data(ifile_output, '.quantities/occupation', pot%occ)
    call hdf5_write_data(ifile_output, '.quantities/carrier', carrier)
    call hdf5_write_data(ifile_output, '.quantities/energy', energy)
    call hdf5_write_data(ifile_output, '.quantities/electrons', electrons)
    call hdf5_write_data(ifile_output, '.quantities/holes', holes)
    if (algo%lImpurities) then
      call hdf5_write_data(ifile_output, '.quantities/imp_contribution', imp_contribution)
    endif
    call hdf5_close_file(ifile_output)
  endif


   ! gather the timings
#ifdef MPI
   time = sum(timings)
   call MPI_allreduce(MPI_IN_PLACE, timings, size(timings), MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpierr) ! total CPU time
   call MPI_allreduce(MPI_IN_PLACE, time, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, mpierr)      ! total real time
   timings = timings / dble(nproc)
#endif
  if (myid .eq. master) then
    write(stdout,*)
    write(stdout,*) '  Timings (average) [s]:'
    write(stdout,'(A21,F19.6)') '    readin + alloc:  ', timings(1)
    write(stdout,'(A21,F19.6)') '    mu-search:       ', timings(2)
    write(stdout,'(A21,F19.6)') '    polygamma-eval:  ', timings(3)
    write(stdout,'(A21,F19.6)') '    response-eval:   ', timings(4)
    write(stdout,'(A21,F19.6)') '    mpi + summation: ', timings(5)
    write(stdout,*)
    write(stdout,'(A23,F17.6)') 'Total real time [s]:', time
    write(stdout,'(A23,F17.6)') 'Total CPU  time [s]:', sum(timings) * dble(nproc)
    write(stdout,*)
  endif

  if (kmesh%ndim /= 3) then
    call log_master(stdout, 'Warning: Reduced dimensions detected')
    call log_master(stdout, 'Warning: Onsager coefficient units may change')
  endif



  call hdf5_close_file(ifile_energy)
  if (algo%lScatteringFile) then
    call hdf5_close_file(ifile_scatter_hdf5)
  endif

  call mpi_close()
  call hdf5_finalize()
end program main
