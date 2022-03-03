program main
  use Mparams
  use Maux
  use Mroot
  use Mtypes
  use Mmpi_org
  use Mconfig
  use Mio
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


  integer :: is, ig, iT, ik, iband, iimp, imu
  integer :: niitact
  integer :: iTstart, iTend, iTstep
  logical :: igap
  real(8) :: maxgap
  real(8) :: ndevact
  real(16):: ndevactQ

  character(len=128) :: string
  real(8) :: tstart, tfinish, timings(5)
  real(8) :: time, maxtime

  ! quantities saved on the Temperature grid
  ! and derived quantities
  real(8), allocatable :: energy(:) ! total energy
  real(8), allocatable :: electrons(:) ! thermally activated electrons w.r.t chemical potential
  real(8), allocatable :: holes(:) ! thermally activated holes w.r.t chemical potential
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
  call read_preproc_energy_data(algo, kmesh, edisp, pot, imp)

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
    if (index(algo%dbgstr,"Verbose") .ne.0) then
      call log_master(stdout, 'Warning: Neither Intra nor Interband responses will be calculated')
    endif
  endif

  if (algo%lTMODE) then
    ! construct temperature grid
    if (algo%lScatteringFile) then
      ! at this point we just extract the temperature ranges
      ! and check wether bandshifts exist
      ! the scattering rates then gets loaded for each temperature-point
      call read_preproc_scattering_data_hdf5(algo, kmesh, edisp, sct, temp)
    else if (algo%lScatteringText) then
      call read_preproc_scattering_data_text(algo, kmesh, edisp, sct, temp)
    else
      if (temp%nT .gt. 1) then
        if (temp%tlogarithmic) then
          temp%dT = (temp%Tmax/temp%Tmin) ** (1.d0 / dble(temp%nT-1)) ! logarithmic step
        else
          temp%dT = (temp%Tmax-temp%Tmin)/dble(temp%nT-1) ! linear step
        endif
      else
        temp%dT = 0 ! 0/0 ...
      endif
      allocate(temp%TT(temp%nT))
      allocate(temp%BB(temp%nT))
      allocate(sct%gam(edisp%nband_max, kmesh%nkp, edisp%ispin))
      allocate(sct%zqp(edisp%nband_max, kmesh%nkp, edisp%ispin))

      ! define Temperature grid
      if (temp%tlogarithmic) then
        do iT=1,temp%nT
           temp%TT(iT)=temp%Tmin * temp%dT**(real(iT-1,8))
        enddo
      else
        do iT=1,temp%nT
           temp%TT(iT)=real(iT-1,8)*temp%dT+temp%Tmin
        enddo
      endif
      temp%BB = 1.d0/(temp%TT * kB)
    endif
    pot%nMu = temp%nT


    allocate(energy(temp%nT))
    allocate(electrons(temp%nT))
    allocate(holes(temp%nT))
    allocate(imp_contribution(temp%nT))
    energy = 0.d0
    electrons = 0.d0
    imp_contribution = 0.d0
    holes = 0.d0
    allocate(pot%MM(pot%nMu))
    allocate(pot%occ(pot%nMu))
    pot%occ = 0.d0
    pot%MM = pot%mu ! here we either have the fixed mu
                    ! or the mu_dft initialized from the preprocessed energy file

    if (algo%lOldmu) then
      call read_muT_hdf5(temp, algo%old_output_file, pot%MM)
    endif

    if (algo%lOldmuText) then
      call read_muT_text(temp, algo%input_mu_text, pot%MM)
    endif
  endif

  if (algo%lMUMODE) then
    ! construct mu grid
    temp%nT = pot%nMu
    allocate(temp%TT(pot%nMu))
    allocate(temp%BB(pot%nMu))
    allocate(sct%gam(edisp%nband_max, kmesh%nkp, edisp%ispin))
    allocate(sct%zqp(edisp%nband_max, kmesh%nkp, edisp%ispin))

    allocate(energy(pot%nMu))
    allocate(electrons(pot%nMu))
    allocate(holes(pot%nMu))
    allocate(imp_contribution(pot%nMu))
    energy = 0.d0
    electrons = 0.d0
    holes = 0.d0
    imp_contribution = 0.d0
    allocate(pot%MM(pot%nMu))
    allocate(pot%occ(pot%nMu))

    temp%TT = temp%temp
    temp%BB = 1.d0/(temp%temp * kB)
    pot%occ = 0.d0

    ! define Chemical potential grid
    pot%dMu= (pot%Mumax-pot%Mumin)/dble(pot%nMu-1)
    do imu=1,pot%nMu
       pot%MM(imu)=real(imu-1,8)*pot%dMu+pot%MuMin
    enddo

  endif

  call mpi_genkstep(kmesh%nkp)

#ifdef MPI
  if (algo%ldebug .and. (index(algo%dbgstr,"Mpi") .ne. 0)) then
     write(stdout,*) "MPI: myid: ", myid, "ikstr: ", ikstr, "ikend: ", ikend
     call mpi_barrier(mpi_comm_world, mpierr)
  endif
#endif


  ! allocate the arrays once outside of the main (temperature) loop
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


  ! for the responses we need psi_1, psi_2 and psi_3
  allocate(PolyGamma(3, edisp%nbopt_min:edisp%nbopt_max, ikstr:ikend, edisp%ispin))

  if (algo%lDebug .and. (index(algo%dbgstr,"Quad") .ne. 0)) then
    call allocate_response(algo, edisp, temp, qresp_intra)
    if (algo%lInterbandquantities) then
      call allocate_response(algo, edisp, temp, qresp_inter)
    endif

    if (algo%lBoltzmann) then
      call allocate_response(algo, edisp, temp, qresp_intra_Boltzmann)
      if (algo%lInterbandquantities) then
        call allocate_response(algo, edisp, temp, qresp_inter_Boltzmann)
      endif
    endif
    allocate(PolyGammaQ(3, edisp%nbopt_min:edisp%nbopt_max, ikstr:ikend, edisp%ispin))
  endif

  call hdf5_open_file(algo%input_energies,   ifile_energy,  rdonly=.true.)

  if (algo%lTMODE .and. algo%lScatteringFile) then
    call hdf5_open_file(algo%input_scattering_hdf5, ifile_scatter_hdf5, rdonly=.true.)
  endif

  ! for this option we can do the shifts once
  ! and deallocate the unused arrays
  ! otherwise we have to recalculate the shifts for every T-point
  if (.not. edisp%lBandShift) then
    if (algo%lScissors) then
      do is=1,edisp%ispin
        edisp%band_shift(:edisp%valenceBand(is), :, is)    = 0.d0
        edisp%band_shift(edisp%conductionBand(is):, :, is) = edisp%scissors(is)
      enddo
      ! apply once and deallocate the unused arrays
      edisp%band = edisp%band_original + edisp%band_shift
      deallocate(edisp%band_original)
      deallocate(edisp%band_shift)
    else
      edisp%band = edisp%band_original ! the energies are constant throughout
      deallocate(edisp%band_shift)
      deallocate(edisp%band_original)
    endif
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

  ! final preparation step for chemical potential mode
  if (algo%lMUMODE) then
    info%iT = 1
    info%Temp=temp%TT(info%iT)
    info%beta=1.d0/(kB*info%Temp)
    info%beta2p=info%beta/(2.d0*pi)

    info%TempQ=real(info%Temp,16)
    info%betaQ=1.q0/(kB*info%TempQ)
    info%beta2pQ=info%betaQ/(2.q0*piQ)

    sct%gam = sct%gamcoeff(1)
    sct%zqp = sct%zqpcoeff(1)
    sct%gam = sct%gam * sct%zqp

    ! find the equilibrium chemical potential
    call find_mu(pot%mu,ndevQ,ndevactQ,niitact, edisp, sct, kmesh, imp, algo, info)

    ! shift the chemical potential ranges by the calculated chemical potential
    ! for the given temperature

    pot%MuMin = pot%MuMin + pot%mu
    pot%MuMax = pot%MuMax + pot%mu
    pot%MM    = pot%MM + pot%mu
  endif


  if (myid .eq. master) then
    write(stdout,*)
    if (algo%lTMODE) then
    write(stdout,*) 'TEMPERATURE MODE'
    write(stdout,*) '  Temperature range:'
    if (.not. algo%lScatteringFile .and. .not. algo%lScatteringText) then
    write(stdout,*) '  Tmin: ', temp%Tmin
    write(stdout,*) '  Tmax: ', temp%Tmax
    endif
    write(stdout,*) '  Temperature points:   ', temp%nT
    else if (algo%lMUMODE) then
    write(stdout,*) 'MU MODE'
    write(stdout,*) '  Temperature: ', temp%TT(1)
    write(stdout,*) '  Chemical Potential for given Temp: ', pot%mu
    write(stdout,*) '  Chemical Potential range:'
    write(stdout,*) '  Mumin: ', pot%Mumin
    write(stdout,*) '  Mumax: ', pot%Mumax
    write(stdout,*) '  Chemical Potential points:   ', pot%nMu
    endif
    write(stdout,*)
    write(stdout,*) 'INPUT'
    write(stdout,*) '  dimensions: ', kmesh%ndim
    write(stdout,*) '  k-Points: ', kmesh%nkp
    write(stdout,*) '  spins: ', edisp%ispin
    write(stdout,*) '  electrons: ', edisp%nelect
    write(stdout,*)
    write(stdout,*) '  energy-file: ', trim(algo%input_energies)
    if (algo%lScissors) then
      write(stdout,*) '  old gap: ', edisp%gap - edisp%scissors
      write(stdout,*) '  new gap: ', edisp%gap
      write(stdout,*) '  scissors: ', edisp%scissors
    else
      write(stdout,*) '  gap: ', edisp%gap
    endif
    if (algo%lTMODE) then
      write(stdout,*)
      if (algo%muSearch) then
        if (algo%muFermi) then
          write(stdout,*) '  chemical potential: determined via Fermi function'
        else
          write(stdout,*) '  chemical potential: determined via Digamma function'
        endif
      else
        if (algo%lOldmu) then
          write(stdout,*) '  chemical potential: from file: ', trim(adjustl(algo%old_output_file))
        else
          write(stdout,*) '  chemical potential: constant: ', pot%MM(1)
        endif
      endif
    endif
    write(stdout,*)
    if (algo%lScatteringFile) then
      write(stdout,*) '  scattering-file: ', trim(algo%input_scattering_hdf5)
      write(stdout,*) '  additional impurity offset: ', sct%gamimp
    else if (algo%lScatteringText) then
      write(stdout,*) '  scattering-file: ', trim(algo%input_scattering_text)
      write(stdout,*) '  additional impurity offset: ', sct%gamimp
    else
      write(stdout,*) '  scattering coefficients: ', sct%gamcoeff
      write(stdout,*) '  quasi particle weight coefficients: ', sct%zqpcoeff
    endif
    write(stdout,*)
    if (algo%lImpurities) then
      write(stdout,*) '  impurity levels: ', imp%nimp
      write(stdout,*) '    ______________________________________________'
      write(stdout,*) '    iimp, dopant, density, energy [eV], degeneracy, width [eV], cutoff [sigma]'
      do iimp = 1,imp%nimp
        write(stdout,'(2X,I5,I5,5F15.10)') iimp, int(imp%Dopant(iimp)), imp%Density(iimp), &
                                           imp%Energy(iimp), imp%Degeneracy(iimp), imp%Bandwidth(iimp), &
                                           imp%Bandcutoff(iimp)
      enddo
    endif
    write(stdout,*)
    write(stdout,*) 'OUTPUT'
    write(stdout,*) '  output-file: ', trim(algo%output_file)
    write(stdout,*) '  full-output: ', algo%lFullOutput
    write(stdout,*)
    write(stdout,*) '  run-options:'
    write(stdout,*) '  interband quantities: ', algo%lIntrabandQuantities
    write(stdout,*) '  interband quantities: ', algo%lInterbandQuantities
    write(stdout,*) '  Boltzmann quantities: ', algo%lBoltzmann
    write(stdout,*) '    Boltzmann with Fermi approximation: ', algo%lBoltzmannFermi
    write(stdout,*) '  B-field   quantities: ', algo%lBfield
    write(stdout,*)
    if (algo%lDebug) then
      write(stdout,*) 'DEBUG MODE'
      write(stdout,*) '  DEBUG STRING: "', trim(adjustl(algo%dbgstr)),'"'
      write(stdout,*)
    endif
    write(stdout,*) 'Starting calculation...'
    write(stdout,*) '____________________________________________________________________________'
    write(stdout,*) 'Temperature[K], invTemperature[1/eV], chemicalPotential[eV]'
  endif

  if (myid .eq. master) then
    call hdf5_create_file(algo%output_file)
    call output_auxiliary(algo, info, temp, kmesh, edisp, sct, imp)
  endif


  call cpu_time(tfinish)
  timings(1) = timings(1) + (tfinish - tstart)
  tstart = tfinish

#ifdef MPI
  call mpi_barrier(mpi_comm_world, mpierr)
#endif



  ! MAIN LOOP
  if (algo%lDebug .and. (index(algo%dbgstr,"TempReverse") .ne. 0)) then
    iTstart = 1
    iTend   = temp%nT
    iTstep  = 1
  else
    iTstart = temp%nT
    iTend   = 1
    iTstep  = -1
  endif

  temp%Tstep = iTstep

  do iT=iTstart,iTend,iTstep
    ! run time information
    info%iT = iT
    info%Temp=temp%TT(iT)
    info%beta=1.d0/(kB*info%Temp)
    info%beta2p=info%beta/(2.d0*pi)

    info%TempQ=real(info%Temp,16)
    info%betaQ=1.q0/(kB*info%TempQ)
    info%beta2pQ=info%betaQ/(2.q0*piQ)

    ! define scattering rates and quasi particle weights
    ! for the current temperature
    if (algo%lTMODE .and. algo%lScatteringFile) then
      ! as mentioned above
      ! apply the scissors:
      if (algo%lScissors .and. edisp%lBandShift) then
        edisp%band_shift = 0.d0
        do is=1,edisp%ispin
          edisp%band_shift(:edisp%valenceBand(is), :, is)    = 0.d0
          edisp%band_shift(edisp%conductionBand(is):, :, is) = edisp%scissors(is)
        enddo
      endif
      ! read in the scattering data for the current temperature
      ! scattering rates; quasi-particle weights and possible band-shifts
      ! this are in ADDITION to the scissors
      call read_scattering_data_hdf5(ifile_scatter_hdf5, edisp, kmesh, sct, info)
    else if (algo%lTMODE .and. algo%lScatteringText) then
      sct%gam = sct%gamtext(iT)
      sct%zqp = sct%zqptext(iT)
      if (sct%zqp(1,1,1) > 1.d0) then ! since its a constant array
        call log_master(stdout, 'ERROR: Zqp is bigger than 1 ... truncating to 1')
        sct%zqp = 1.d0
      endif
      sct%gam = sct%zqp * sct%gam
    else ! this is entered for both the MuMode and the TempMode
      ! here we are wasting memory however
      ! otherwise we would have to write every single response twice
      sct%gam = 0.d0
      sct%zqp = 0.d0
      do ig=1,size(sct%gamcoeff)
         sct%gam = sct%gam + sct%gamcoeff(ig)*(temp%TT(iT)**(ig-1))
      enddo
      do ig=1,size(sct%zqpcoeff)
         sct%zqp = sct%zqp + sct%zqpcoeff(ig)*(temp%TT(iT)**(ig-1))
      enddo
      if (sct%zqp(1,1,1) > 1.d0) then ! since its a constant array
        call log_master(stdout, 'ERROR: Zqp is bigger than 1 ... truncating to 1')
        sct%zqp = 1.d0
      endif
      sct%gam = sct%zqp * sct%gam  ! convention we use
    endif


    ! root finding (mu)
    niitact = 0
    if (algo%lTMODE .and. algo%muSearch) then

      ! initialize the new chemical potential
      if (iT /= iTstart) then
        pot%MM(iT) = pot%MM(iT-iTstep) ! we initialize it to the value from the previous iteration
                               ! this does nothing if mu is constant
      endif

      call cpu_time(tstart)

      if (edisp%gapped_complete) then
        if (algo%muFermi) then
          ! use quadruple precision
          ! internally use refinement method if in the right temperature range
          call find_mu(pot%MM(iT),ndevVQ,ndevactQ,niitact, edisp, sct, kmesh, imp, algo, info)
        else
          ! also use quadruple precision
          ! however the root-finding can be aborted earlier
          call find_mu(pot%MM(iT),ndevQ,ndevactQ,niitact, edisp, sct, kmesh, imp, algo, info)
        endif
      else
        ! if the system is not gapped, simple double precision is enough
        ! call find_mu(mu(iT),ndev,ndevact,niitact, edisp, sct, kmesh, imp, algo, info)
        ! fuck it
        call find_mu(pot%MM(iT),ndevQ,ndevactQ,niitact, edisp, sct, kmesh, imp, algo, info)
      endif
      call cpu_time(tfinish)
      timings(2) = timings(2) + (tfinish - tstart)
      tstart = tfinish
    endif


    if (algo%lTMODE .and. (.not. (index(algo%dbgstr,"TempReverse") .ne. 0)) .and. info%Temp < edisp%gap_min / 1.95q0 &
        .and. algo%muFermi) then
        if (iT+2 <= temp%nT) then
          pot%MM(iT) = pot%MM(iT+1) + (pot%MM(iT+2)-pot%MM(iT+1))/(temp%TT(iT+2)-temp%TT(iT+1)) * &
                   (temp%TT(iT)-temp%TT(iT+1))
          if (index(algo%dbgstr,"Verbose") .ne.0) call log_master(stdout, 'Debug: Applying dmu/dT')
        else
          call stop_with_message(stderr, 'Debug: Cannot apply dmu/dT')
        endif
    endif

    ! calculating total energy according to the occuption above
    if (algo%muFermi) then
      call calc_total_energy_fermi(pot%MM(iT), energy(iT), edisp, sct, kmesh, imp, algo, info)
      call calc_elecholes_fermi(pot%MM(iT), electrons(iT), holes(iT), edisp, sct, kmesh, imp, algo, info)
    else
      call calc_total_energy_digamma(pot%MM(iT), energy(iT), edisp, sct, kmesh, imp, algo, info)
      call calc_elecholes_digamma(pot%MM(iT), electrons(iT), holes(iT), edisp, sct, kmesh, imp, algo, info)
    endif

    if (algo%lTMODE) then
      call occ_impurity(imp_contribution(iT), pot%MM(iT), imp, info)
    endif

    ! calculates the difference to the demanded electron number
    call ndeviation_D(pot%MM(iT), pot%occ(iT), edisp, sct, kmesh, imp, algo, info)
    pot%occ(iT) = edisp%nelect - pot%occ(iT)

    if (myid.eq.master) then
      write(stdout,'(3X,3F15.7,I7)') info%Temp, info%beta, pot%MM(iT), niitact
    endif

    ! calculate the polygamma function (1...3)
    ! for all optical bands, spins and each core's kpoints
    ! once and use it later for all the different response types
    call calc_polygamma(PolyGamma, pot%MM(iT), edisp, sct, kmesh, algo, info)

    if (algo%lDebug .and. (index(algo%dbgstr, "Quad") .ne. 0)) then
      call calc_polygamma(PolyGammaQ, pot%MM(iT), edisp, sct, kmesh, algo, info)
      call initresp_qp(algo, qresp_intra)
      if (algo%lInterBandQuantities) then
        call initresp_qp(algo, qresp_inter)
      endif
      if (algo%lBoltzmann) then
        call initresp_qp(algo, qresp_intra_Boltzmann)
        if (algo%lInterBandQuantities) then
          call initresp_qp(algo, qresp_inter_Boltzmann)
        endif
      endif
    endif

    call cpu_time(tfinish)
    timings(3) = timings(3) + (tfinish - tstart)
    tstart = tfinish



    ! initialize the already allocated arrays to 0
    call initresp(algo, resp_intra)
    if(algo%lBoltzmann) then
      call initresp(algo, resp_intra_Boltzmann)
    endif
    if (algo%lInterBandQuantities) then
      call initresp(algo, resp_inter)
      if (algo%lBoltzmann) then
        call initresp(algo, resp_inter_Boltzmann)
      endif
    endif

    ! do the k-point loop and calculate the response
    do ik = ikstr,ikend
      info%ik = ik ! save into the runinformation datatype

      ! double precision routines
      if (algo%lIntrabandQuantities) then
        call response_intra_km(resp_intra,  PolyGamma, pot%MM(iT), edisp, sct, kmesh, algo, info)
        if (algo%lBoltzmann) then
          call response_intra_Boltzmann_km(resp_intra_Boltzmann, pot%MM(iT), edisp, sct, kmesh, algo, info)
        endif
      endif

      if (algo%lInterbandquantities) then
        call response_inter_km(resp_inter, PolyGamma, pot%MM(iT), edisp, sct, kmesh, algo, info)
        if (algo%lBoltzmann) then
          call response_inter_Boltzmann_km(resp_inter_Boltzmann, pot%MM(iT), edisp, sct, kmesh, algo, info)
        endif
      endif

      ! quad precision routines
      if (algo%lDebug .and. (index(algo%dbgstr, "Quad") .ne. 0)) then
        if (algo%lIntrabandQuantities) then
          call response_intra_km_Q(qresp_intra, PolyGammaQ, pot%MM(iT), edisp, sct, kmesh, algo, info)
          if (algo%lBoltzmann) then
            call response_intra_Boltzmann_km_Q(qresp_intra_Boltzmann, pot%MM(iT), edisp, sct, kmesh, algo, info)
          endif
        endif
        if (algo%lInterBandQuantities) then
          call response_inter_km_Q(qresp_inter, PolyGammaQ, pot%MM(iT), edisp, sct, kmesh, algo, info)
          if (algo%lBoltzmann) then
            call response_inter_Boltzmann_km_Q(qresp_inter_Boltzmann, pot%MM(iT), edisp, sct, kmesh, algo, info)
          endif
        endif
      endif
    enddo

    call cpu_time(tfinish)
    timings(4) = timings(4) + (tfinish - tstart)
    tstart = tfinish

    ! if (myid.eq. master) then
    !   call intldos(mu(iT), dos, edisp, sct, kmesh, algo, info)
    ! endif


    ! output the response
    ! this subroutines include the summations necessary
    if (algo%lIntrabandQuantities) then
      call response_h5_output(resp_intra, "intra", edisp, algo, info, temp, kmesh)
      if (algo%lBoltzmann) then
        call response_h5_output(resp_intra_Boltzmann, "intraBoltzmann", edisp, algo, info, temp, kmesh)
      endif
    endif
    if (algo%lInterbandQuantities) then
      ! here we don't have the Bfield quantities ... no optical elements yet
      call response_h5_output(resp_inter, "inter", edisp, algo, info, temp, kmesh, .false.)
      if (algo%lBoltzmann) then
        call response_h5_output(resp_inter_Boltzmann, "interBoltzmann", edisp, algo, info, temp, kmesh, .false.)
      endif
    endif

    if (algo%lDebug .and. (index(algo%dbgstr, "Quad") .ne. 0)) then
      if (algo%lIntrabandQuantities) then
        call response_h5_output_Q(qresp_intra, "intraQuad", edisp, algo, info, temp, kmesh)
        if (algo%lBoltzmann) then
          call response_h5_output_Q(qresp_intra_Boltzmann, "intraQuadBoltzmann", edisp, algo, info, temp, kmesh)
        endif
      endif

      if (algo%lInterbandQuantities) then
        call response_h5_output_Q(qresp_inter, "interQuad", edisp, algo, info, temp, kmesh, .false.)
        if (algo%lBoltzmann) then
          call response_h5_output_Q(qresp_inter_Boltzmann, "interQuadBoltzmann", edisp, algo, info, temp, kmesh, .false.)
        endif
      endif
    endif

    ! output the renormalized energies defined by Z*(ek - mu)
    ! for each temperature point
    if (myid.eq.master .and. algo%lEnergyOutput) then
      call output_energies(pot%MM(iT), algo, edisp,kmesh,sct,info)
    endif

    call cpu_time(tfinish)
    timings(5) = timings(5) + (tfinish - tstart)
    tstart = tfinish

  enddo ! end of the outer temperature loop

  if (allocated(edisp%Mopt)) deallocate(edisp%Mopt)

  if (myid.eq.master) then
    call hdf5_open_file(algo%output_file, ifile_output)
    ! FIXME ... keep this or not ?
    if (algo%lMuMode) then
      pot%MM = pot%MM - pot%mu ! for output reasons
    endif
    call hdf5_write_data(ifile_output, '.quantities/mu', pot%MM)
    call hdf5_write_data(ifile_output, '.quantities/occupation', pot%occ)
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
