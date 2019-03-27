program main
  use Mparams
  use Maux
  use Mdos
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
  type(runinfo)     :: info    ! runtime information for the calculation routines, temps, betas, etc.

  type(response_dp) :: resp_intra
  type(response_dp) :: resp_intra_Boltzmann
  type(response_dp) :: resp_inter
  type(response_dp) :: resp_inter_Boltzmann

  ! type(response_qp) :: qpresp  ! response quadruple precision

  integer(hid_t)    :: ifile_scatter
  integer(hid_t)    :: ifile_energy

  integer :: is, ig, iT, ik, iband
  integer :: niitact
  real(8) :: ndevact
  real(16):: ndevactQ

  character(len=128) :: string
  real(8) :: tstart, tfinish, timings(4)

  ! quantities saved on the Temperature grid
  ! and derived quantities
  real(8), allocatable :: drhodT(:) ! resisity derivative
  real(8), allocatable :: energy(:) ! total energy
  real(8), allocatable :: cv(:)     ! specific heat
  real(8), allocatable :: mu(:)     ! chemical potential (temperature dependent)

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
  call mpi_barrier(mpi_comm_world, mpierr)

  call read_config(algo, edisp, sct, temp)
  call check_config(algo)

  call hdf5_init()
  ! read the energies ( and derivatives if they exist )
  call read_preproc_energy_data(algo, kmesh, edisp)

  if (algo%lBfield .and. .not. edisp%lDerivatives) then
    call stop_with_message(stderr, 'Energy derivatives required for Bfield quantities')
  endif

  ! calculate a purely DFT density of states with Laurentzian broadening
  ! call gendosel (kmesh, edisp, dos) ! normalization already taken care of
  ! call findef(dos, edisp)   ! finds the (non-interacting) Fermi level

  ! if (myid.eq.master) then
  !   do is = 1,edisp%ispin
  !     if (dos%gap(is) == 0.d0) then
  !       write(stdout,*) 'Detected no band gap in spin ', is, ' / ', edisp%ispin
  !     else
  !       write(stdout,*) 'Detected band gap of size: ', dos%gap(is), 'in spin ', is, ' / ', edisp%ispin
  !     endif
  !   enddo
  ! endif

  ! starting point for the chemical potential
  ! if (myid.eq.master) then
  !    if (algo%muSearch) then
  !       write(stdout,*) 'initialized LDA mu = ', edisp%efer
  !    else
  !       write(stdout,*) 'Running with fixed mu read from file = ',  edisp%mu
  !    endif
  ! endif


  ! so every single process runs with this given mu
  if (.not. algo%muSearch) then
    edisp%efer = edisp%mu ! we overwrite the calculated fermienergy with the provided chem.pot
  endif

  ! construct temperature grid
  if (algo%lScatteringFile) then
    ! at this point we just extract the temperature ranges
    ! and check wether bandshifts exist
    ! the scattering rates then gets loaded for each temperature-point
    call read_preproc_scattering_data(algo, kmesh, edisp, sct, temp)
  else
    temp%dT= (temp%Tmax-temp%Tmin)/dble(temp%nT-1)
    allocate(temp%TT(temp%nT))
    allocate(temp%beta(temp%nT))

    ! define Temperature grid
    do iT=1,temp%nT
       temp%TT(iT)=real(iT-1,8)*temp%dT+temp%Tmin
    enddo
    temp%beta = 1.d0/(temp%TT * kB)
  endif


  allocate(mu(temp%nT))
  allocate(energy(temp%nT))
  allocate(cv(temp%nT))

  ! either we start with the LDA mu
  ! or with the fixed mu from above
  mu     = edisp%efer ! here we either have the fixed mu or the LDA initialized one from above
  energy = 0.d0
  cv     = 0.d0


  call mpi_genkstep(kmesh%nkp)

  if (algo%ldebug) then
     write(stdout,*) "MPI: myid: ", myid, "ikstr: ", ikstr, "ikend: ", ikend
#ifdef MPI
     call mpi_barrier(mpi_comm_world, mpierr)
#endif
  endif


  ! allocate the arrays once outside of the main (temperature) loop
  call allocate_response(algo, edisp, resp_intra)
  if (algo%lBoltzmann) then
    call allocate_response(algo, edisp, resp_intra_Boltzmann)
  endif
  if (algo%lInterbandquantities) then
    call allocate_response(algo, edisp, resp_inter)
    if (algo%lBoltzmann) then
      call allocate_response(algo, edisp, resp_inter_Boltzmann)
    endif
  endif

  ! for the responses we need psi_1, psi_2 and psi_3
  allocate(PolyGamma(3, edisp%nbopt_min:edisp%nbopt_max, ikstr:ikend, edisp%ispin))

  ! if (.not. algo%ldebug) then
  !   call allocate_response(algo, edisp, qpresp)
  !   allocate(PolyGammaQ(3, edisp%nbopt_min:edisp%nbopt_max, ikstr:ikend, edisp%ispin))
  ! endif

  call hdf5_open_file(algo%input_energies,   ifile_energy,  rdonly=.true.)
  if (algo%lScatteringFile) then
    call hdf5_open_file(algo%input_scattering, ifile_scatter, rdonly=.true.)
  endif

  if (.not. edisp%lBandShift) then
    edisp%band = edisp%band_original ! the energies are constant throughout
    deallocate(edisp%band_original)
  endif

  if (edisp%ispin == 1) then
    edisp%nelect = edisp%nelect * 2
  endif

  if (myid .eq. master) then
    write(stdout,*)
    write(stdout,*)
    write(stdout,*) 'Option summary:'
    write(stdout,*)
    write(stdout,*) '  Temperature range:'
    write(stdout,*) '  Tmin: ', temp%Tmin
    write(stdout,*) '  Tmax: ', temp%Tmax
    write(stdout,*) '  Temperature points:   ', temp%nT
    write(stdout,*)
    write(stdout,*) '  k-Points: ', kmesh%nkp
    write(stdout,*) '  spins: ', edisp%ispin
    write(stdout,*) '  electrons: ', edisp%nelect
    write(stdout,*)
    write(stdout,*) '  energy-file: ', trim(algo%input_energies)
    if (.not. algo%muSearch) then
      write(stdout,*) '  constant chemical potential: ', edisp%efer
    endif
      write(stdout,*)
    if (algo%lScatteringFile) then
      write(stdout,*) '  scattering-file: ', trim(algo%input_scattering)
      write(stdout,*) '  additional impurity offset: ', sct%gamimp
    else
      write(stdout,*) '  scattering coefficients: ', sct%gamcoeff
      write(stdout,*) '  quasi particle weight coefficients: ', sct%zqpcoeff
    endif
    write(stdout,*)
    write(stdout,*) '  output-options:'
    write(stdout,*) '  output-file: ', trim(algo%output_file)
    write(stdout,*) '  full-output: ', algo%lFullOutput
    write(stdout,*)
    write(stdout,*) '  run-options:'
    write(stdout,*) '  interband quantities: ', algo%lInterbandQuantities
    write(stdout,*) '  Boltzmann quantities: ', algo%lBoltzmann
    write(stdout,*) '  B-field   quantities: ', algo%lBfield
    write(stdout,*)
    write(stdout,*)
    write(stdout,*) 'Starting calculation...'
    write(stdout,*) '____________________________________________________________________________'
    write(stdout,*) 'Temperature[K], invTemperature[1/eV], chemicalPotential[eV], totalEnergy[eV]'
  endif

  if (myid .eq. master) then
    call hdf5_create_file(algo%output_file)
    call output_auxiliary(algo, info, temp, kmesh)
  endif


  timings = 0.d0        ! reset timings
  call cpu_time(tstart) ! start timer

  ! MAIN LOOP
  do iT=1,temp%nT
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
    if (algo%lScatteringFile) then
      ! read in the scattering data for the current temperature
      ! scattering rates; quasi-particle weights and possible band-shifts
      call read_scattering_data(ifile_scatter, edisp, sct, info)
    else
      sct%gamscalar = 0.d0
      sct%zqpscalar = 0.d0
      do ig=1,size(sct%gamcoeff)
         sct%gamscalar = sct%gamscalar + sct%gamcoeff(ig)*(temp%TT(iT)**(ig-1))
      enddo
      do ig=1,size(sct%zqpcoeff)
         sct%zqpscalar = sct%zqpscalar + sct%zqpcoeff(ig)*(temp%TT(iT)**(ig-1))
      enddo
      if (sct%zqpscalar > 1.d0) then
        call log_master(stdout, 'WARNING: Zqp is bigger than 1 ... truncating to 1')
        sct%zqpscalar = 1.d0
      endif
      sct%gamscalar = sct%zqpscalar * sct%gamscalar  ! convention we use
    endif


    niitact = 0
    if (algo%muSearch) then
      call cpu_time(tstart)
      call find_mu(mu(iT),ndev,ndevact,niitact, edisp, sct, kmesh, algo, info)
      ! call find_mu(mu(iT),ndevQ,ndevactQ,niitact, edisp, sct, kmesh, algo, info)
      ! call find_mu(mu(iT),ndevVQ,ndevactQ,niitact, edisp, sct, kmesh, algo, info)
      call cpu_time(tfinish)
      timings(1) = timings(1) + (tfinish - tstart)
      tstart = tfinish
    endif

    ! call calc_total_energy(mu(iT), energy(iT), edisp, sct, kmesh, algo, info)

    if (myid.eq.master) then
      write(stdout,*)info%Temp, info%beta, mu(iT), energy(iT), niitact
    endif

    ! calculate the polygamma function (1...3)
    ! for all optical bands, spins and each core's kpoints
    ! once and use it later for all the different response types
    call calc_polygamma(PolyGamma, mu(iT), edisp, sct, kmesh, algo, info)
    ! if (.not. algo%lDebug) then
    !   call calc_polygamma(PolyGammaQ, mu(iT), edisp, sct, kmesh, algo, info)
    ! endif
    call cpu_time(tfinish)
    timings(2) = timings(2) + (tfinish - tstart)
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

    ! if (.not. algo%ldebug) then
    !    call initresp_qp (algo, qpresp)
    ! endif


    do ik = ikstr,ikend
      info%ik = ik ! save into the runinformation datatype

      ! load the moments for the current k-point
      call read_optical_elements(ifile_energy, edisp, sct, info)

      ! calculate the response
      call calc_response(PolyGamma, mu(iT), edisp, sct, kmesh, algo, info, &
                         resp_intra, resp_intra_Boltzmann, &
                         resp_inter, resp_inter_Boltzmann)
    enddo

    call cpu_time(tfinish)
    timings(3) = timings(3) + (tfinish - tstart)
    tstart = tfinish

    ! if (myid.eq. master) then
    !   call intldos(mu(iT), dos, edisp, sct, kmesh, algo, info)
    ! endif


    call response_h5_output(resp_intra, "intra", edisp, algo, info, temp, kmesh)
    if (algo%lBoltzmann) then
      call response_h5_output(resp_intra_Boltzmann, "intraBoltzmann", edisp, algo, info, temp, kmesh)
    endif
    if (algo%lInterbandQuantities) then
      ! here we don't have the Bfield quantities
      call response_h5_output(resp_inter, "inter", edisp, algo, info, temp, kmesh, .false.)
      if (algo%lBoltzmann) then
        call response_h5_output(resp_inter_Boltzmann, "interBoltzmann", edisp, algo, info, temp, kmesh, .false.)
      endif
    endif

    if (myid.eq.master .and. algo%lEnergyOutput) then
      call output_energies(mu(iT), algo, edisp,kmesh,sct,info)
    endif

    call cpu_time(tfinish)
    timings(4) = timings(4) + (tfinish - tstart)
    tstart = tfinish

  enddo ! end of the outer temperature loop


   ! gather the timings
#ifdef MPI
   call MPI_allreduce(MPI_IN_PLACE, timings,size(timings), MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpierr)
   timings(:4) = timings(:4) / dble(nproc)
#endif
  if (myid .eq. master) then
    write(stdout,*)
    write(stdout,*) '  Timings (average) [s]:'
    write(stdout,'(A21,F16.6)') '    mu-search:       ', timings(1)
    write(stdout,'(A21,F16.6)') '    polygamma-eval:  ', timings(2)
    write(stdout,'(A21,F16.6)') '    response-eval:   ', timings(3)
    write(stdout,'(A21,F16.6)') '    mpi + summation: ', timings(4)
    write(stdout,*)
  endif


  call hdf5_close_file(ifile_energy)
  if (algo%lScatteringFile) then
    call hdf5_close_file(ifile_scatter)
  endif

  call mpi_close()
  call hdf5_finalize()
end program main
