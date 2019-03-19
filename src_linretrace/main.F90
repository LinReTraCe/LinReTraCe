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
  type(lattice)     :: lat     ! lattice information (volume)
  type(runinfo)     :: info    ! runtime information for the calculation routines, temps, betas, etc.

  type(response_dp) :: dpresp  ! response double precision
  type(response_qp) :: qpresp  ! response quadruple precision
  type(response_dp) :: respBl  ! response Boltzman regime
  type(response_dp) :: dinter  ! response interband
  type(response_dp) :: dderesp ! response functions (intraband conductivity) derivatives in double precision

  integer(hid_t)    :: ifile_scatter
  integer(hid_t)    :: ifile_energy

  integer :: is, ig, iT, ik, iband
  integer :: niitact
  real(8) :: ndevact
  real(16):: ndevactQ

  real(8) :: criterion
  character(len=128) :: string
  real(8) :: tstart, tfinish, timings(5)

  ! quantities saved on the Temperature grid
  ! and derived quantities
  real(8), allocatable :: drhodT(:) ! resisity derivative
  real(8), allocatable :: energy(:) ! total energy
  real(8), allocatable :: cv(:)     ! specific heat
  real(8), allocatable :: mu(:)     ! chemical potential (temperature dependent)
  real(8), allocatable :: d1(:)     ! square of the 1st derivative of sigma
  real(8), allocatable :: d2(:)     ! product of the 2nd derivative times sigma
  real(8), allocatable :: d0(:)     ! linear combination of the two above, whose zero corresponds to T*
  real(8) :: Tstar                  ! temperature for which (d^2 rho)/(d beta^2)=0
                                    ! in practice it is the T for which (d^2 sigma)/(d beta^2) changes sign
  real(8) :: Tflat                  ! T for which (d sigma)/(d beta) changes sign (onset of saturation)

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
  call read_preproc_energy_data(algo, kmesh, edisp, lat)

  if (algo%lBfield .and. .not. edisp%lDerivatives) then
    call stop_with_message(stderr, 'Energy derivatives required for Bfield quantities')
  endif

  if(.not.(algo%lBfield)) then
     call log_master(stdout, 'LINRETRACE will NOT perform calculations with magnetic field')
  else
     call log_master(stdout, 'LINRETRACE WILL perform calculations with magnetic field')
  endif

  ! calculate a purely DFT density of states with Laurentzian broadening
  call gendosel (kmesh, edisp, dos) ! normalization already taken care of
  call findef(dos, edisp)   ! finds the (non-interacting) Fermi level

  if (myid.eq.master) then
    do is = 1,edisp%ispin
      if (dos%gap(is) == 0.d0) then
        write(stdout,*) 'Detected no band gap in spin ', is, ' / ', edisp%ispin
      else
        write(stdout,*) 'Detected band gap of size: ', dos%gap(is), 'in spin ', is, ' / ', edisp%ispin
      endif
    enddo
  endif

  ! starting point for the chemical potential
  if (myid.eq.master) then
     if (algo%muSearch) then
        write(stdout,*) 'initialized LDA mu = ', edisp%efer
     else
        edisp%efer = edisp%mu ! we overwrite the calculated fermienergy with the provided chem.pot
        write(stdout,*) 'Running with fixed mu read from file = ',  edisp%efer
     endif
  endif


  ! construct temperature grid
  if (algo%lScatteringFile) then
    ! we have everything inside here
    ! at this point we just extract the temperature ranges
    ! the scattering rates then gets loaded for each k-point
    ! (similar to the optical elements)
    call read_preproc_scattering_data(algo, kmesh, edisp, sct, temp)
  else
    temp%dT= (temp%Tmax-temp%Tmin)/temp%nT
    allocate(temp%TT(temp%nT))
    allocate(temp%beta(temp%nT))

    ! define Temperature grid
    do iT=1,temp%nT
       temp%TT(iT)=real(iT-1,8)*temp%dT+temp%Tmin
    enddo
    temp%TT(temp%nT) = temp%Tmax ! to avoid numerical errors at the last point
    temp%beta = 1.d0/(temp%TT * kB)
  endif


  allocate(mu(temp%nT))
  allocate(d0(temp%nT))
  allocate(d1(temp%nT))
  allocate(d2(temp%nT))
  allocate(drhodT(temp%nT))
  allocate(energy(temp%nT))
  allocate(cv(temp%nT))

  mu     = edisp%efer ! here we either have the fixed mu or the LDA initialized one from above
  d0     = 0.d0
  d1     = 0.d0
  d2     = 0.d0
  drhodT = 0.d0
  energy = 0.d0
  cv     = 0.d0
  Tstar  = 0.d0
  Tflat  = 0.d0

  call mpi_genkstep(kmesh%nkp) ! thats definitely correct

  if (algo%ldebug) then
     write(stdout,*) "MPI: myid: ", myid, "ikstr: ", ikstr, "ikend: ", ikend
  endif



  ! allocate the arrays once outside of the main (temperature) loop
  call allocate_response(algo, edisp, dpresp)
  call allocate_response(algo, edisp, respBl)
  call allocate_response(algo, edisp, dderesp)
  call allocate_response(algo, edisp, dinter)
  allocate(PolyGamma(3, edisp%nbopt_min:edisp%nbopt_max, ikstr:ikend, edisp%ispin))

  if (.not. algo%ldebug) then
    call allocate_response(algo, edisp, qpresp)
    allocate(PolyGammaQ(3, edisp%nbopt_min:edisp%nbopt_max, ikstr:ikend, edisp%ispin))
  endif
  ! for the responses we need psi_1, psi_2 and psi_3

  call log_master(stdout, 'DEBUG MODE')

  call hdf5_open_file(algo%input_energies,   ifile_energy,  rdonly=.true.)
  if (algo%lScatteringFile) then
    call hdf5_open_file(algo%input_scattering, ifile_scatter, rdonly=.true.)
  endif

  if (.not. edisp%lBandShift) then
    edisp%band = edisp%band_original ! the energies are constant throughout
    deallocate(edisp%band_original)
  endif


  if (myid .eq. master) then
    write(stdout,*)
    write(stdout,*)
    write(stdout,*) 'Calculation options summary:'
    write(stdout,*)
    write(stdout,*) '  Temperature range:'
    write(stdout,*) '  Tmin: ', temp%Tmin
    write(stdout,*) '  Tmax: ', temp%Tmax
    write(stdout,*) '  Temperature points:   ', temp%nT
    write(stdout,*)
    write(stdout,*) '  k-Points: ', kmesh%nkp
    write(stdout,*) '  spins: ', edisp%ispin
    write(stdout,*)
    write(stdout,*) '  energy-file: ', trim(algo%input_energies)
    if (algo%lScatteringFile) then
      write(stdout,*) '  scattering-file: ', trim(algo%input_scattering)
      write(stdout,*) '  additional impurity offset: ', sct%gamimp
    else
      write(stdout,*) '  scattering coefficients: ', sct%gamcoeff
      write(stdout,*) '  quasi particle weight coefficients: ', sct%zqpcoeff
    endif
    write(stdout,*)
  endif

  if (myid .eq. master) then
    call hdf5_create_file(algo%output_file)
  endif
  call mpi_barrier(mpi_comm_world, mpierr)

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

    criterion=info%beta/20.d0

    ! define scattering rates and quasi particle weights
    ! for the current temperature
    ! move this to a subroutine
    if (.not. algo%lScatteringFile) then
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
    else
      if (edisp%ispin == 1) then
        write(string,'("tPoint/",I6.6,"/scatrate")') iT
        call hdf5_read_data(ifile_scatter, string, darr2)
        sct%gam(:,:,1) = darr2
        deallocate(darr2)

        write(string,'("tPoint/",I6.6,"/qpweight")') iT
        call hdf5_read_data(ifile_scatter, string, darr2)
        sct%zqp(:,:,1) = darr2
        deallocate(darr2)

        if (edisp%lBandShift) then
          write(string,'("tPoint/",I6.6,"/bandshift")') iT
          call hdf5_read_data(ifile_scatter, string, darr2)
          edisp%band_shift(:,:,1) = darr2
          deallocate(darr2)

          write(*,*) 'applying band_shift'
          edisp%band = edisp%band_original + edisp%band_shift
        endif

      else
        write(string,'("up/tPoint/",I6.6,"/scatrate")') iT
        call hdf5_read_data(ifile_scatter, string, darr2)
        sct%gam(:,:,1) = darr2
        deallocate(darr2)

        write(string,'("dn/tPoint/",I6.6,"/scatrate")') iT
        call hdf5_read_data(ifile_scatter, string, darr2)
        sct%gam(:,:,2) = darr2
        deallocate(darr2)

        write(string,'("up/tPoint/",I6.6,"/qpweight")') iT
        call hdf5_read_data(ifile_scatter, string, darr2)
        sct%zqp(:,:,1) = darr2
        deallocate(darr2)

        write(string,'("dn/tPoint/",I6.6,"/qpweight")') iT
        call hdf5_read_data(ifile_scatter, string, darr2)
        sct%zqp(:,:,2) = darr2
        deallocate(darr2)


        if (edisp%lBandShift) then
          write(string,'("up/tPoint/",I6.6,"/bandshift")') iT
          call hdf5_read_data(ifile_scatter, string, darr2)
          edisp%band_shift(:,:,1) = darr2
          deallocate(darr2)

          write(string,'("dn/tPoint/",I6.6,"/bandshift")') iT
          call hdf5_read_data(ifile_scatter, string, darr2)
          edisp%band_shift(:,:,2) = darr2
          deallocate(darr2)

          edisp%band = edisp%band_original + edisp%band_shift
        endif
      endif
      sct%gam = sct%gam + sct%gamimp ! so we have access to a constant shift right from the config file
      sct%gam = sct%gam * sct%zqp    ! convention
    endif



    if (algo%muSearch) then
      call cpu_time(tstart)
      if (criterion.lt.20.d0) then !DP
        call find_mu(mu(iT),ndev,ndevact,niitact, edisp, sct, kmesh, algo, info)
        ! if (myid.eq.master) then
        !    ! write(*,'(1F10.5,4E15.7,2I5)')info%Temp, mu, criterion, ndev, abs(ndevact),niit,niitact
        !    write(*,*)info%Temp, mu(iT), criterion, ndev, abs(ndevact),niit,niitact
        ! endif
      elseif (criterion.lt.80.d0) then !QP
        call find_mu(mu(iT),ndevQ,ndevactQ,niitact, edisp, sct, kmesh, algo, info)
        ! if (myid.eq.master) then
        !    write(*,*) info%Temp, mu(iT), criterion, real(ndevQ,8), abs(real(ndevactQ,8)),niitQ,niitact
        ! endif
      else   ! further refinement
        call find_mu(mu(iT),ndevVQ,ndevactQ,niitact, edisp, sct, kmesh, algo, info)
        ! if (myid.eq.master) then
        !    write(*,*)info%Temp ,mu(iT), criterion, real(ndevVQ,8), abs(real(ndevactQ,8)),niitQ,niitact
        ! endif
      endif !criterion
      call cpu_time(tfinish)
      timings(1) = timings(1) + (tfinish - tstart)
      tstart = tfinish
    endif

    ! call calc_total_energy(mu(iT), energy(iT), edisp, sct, kmesh, algo, info)

    if (myid.eq.master) then
      if (iT == 1) then
        write(stdout,*) 'Temperature, invTemperature, chemicalPotential, totalEnergy'
      endif
      write(stdout,*)info%Temp, info%beta, mu(iT), energy(iT), niitact
    endif

    call calc_polygamma(PolyGamma, mu(iT), edisp, sct, kmesh, algo, info)
    if (.not. algo%lDebug) then
      call calc_polygamma(PolyGammaQ, mu(iT), edisp, sct, kmesh, algo, info)
    endif
    call cpu_time(tfinish)
    timings(2) = timings(2) + (tfinish - tstart)
    tstart = tfinish




    ! initialize the already allocated arrays to 0
    call initresp (algo, dpresp)
    call initresp (algo, respBl)
    call initresp (algo, dinter)
    if (.not. algo%ldebug) then
       call initresp_qp (algo, qpresp)
    endif

    do ik = ikstr,ikend
      info%ik = ik
      ! load the moments
      if (edisp%ispin == 1) then
        if (allocated(darr3)) deallocate(darr3)
        write(string,'("kPoint/",I6.6,"/moments")') ik
        call hdf5_read_data(ifile_energy, string, darr3)
        edisp%Mopt(:,:,:,1) = darr3
        deallocate(darr3)
      else
        if (allocated(darr3)) deallocate(darr3)
        write(string,'("up/kPoint/",I6.6,"/moments")') ik
        call hdf5_read_data(ifile_energy, string, darr3)
        edisp%Mopt(:,:,:,1) = darr3
        deallocate(darr3)
        write(string,'("dn/kPoint/",I6.6,"/moments")') ik
        call hdf5_read_data(ifile_energy, string, darr3)
        edisp%Mopt(:,:,:,2) = darr3
        deallocate(darr3)
      endif

      call calc_response(PolyGamma, mu(iT), edisp, sct, kmesh, algo, info, dpresp, respBl, dinter, qpresp)
    enddo

    call cpu_time(tfinish)
    timings(3) = timings(3) + (tfinish - tstart)
    tstart = tfinish

    if (myid.eq. master) then
      call intldos(mu(iT), dos, edisp, sct, kmesh, algo, info)
    endif


  ! SUMMATIONS
    if (myid .eq. master) then
      allocate(dpresp%s_gather(3,3,edisp%nbopt_min:edisp%nbopt_max,edisp%ispin,kmesh%nkp))
    else
      allocate(dpresp%s_gather(1,1,1,1,1))
    endif
#ifdef MPI
    call MPI_gatherv(dpresp%s_full,(ikend-ikstr+1)*9*(edisp%nbopt_max-edisp%nbopt_min+1)*edisp%ispin, &
                     MPI_DOUBLE_COMPLEX, dpresp%s_gather, rcounts*9*(edisp%nbopt_max-edisp%nbopt_min+1)*edisp%ispin, &
                     displs*9*(edisp%nbopt_max-edisp%nbopt_min+1)*edisp%ispin, MPI_DOUBLE_COMPLEX, master, &
                     MPI_COMM_WORLD, mpierr)
#else
    dpresp%s_gather = dpresp%s_full
#endif

    if (myid .eq. master) then
      do ik = 1,kmesh%nkp
        do iband = edisp%nbopt_min,edisp%nbopt_max
          dpresp%s_sum(:,:,:) = dpresp%s_sum(:,:,:) + dpresp%s_gather(:,:,iband,:,ik) * kmesh%weight(ik)
        enddo
      enddo
    endif


    if (myid .eq. master) then
      allocate(dpresp%a_gather(3,3,edisp%nbopt_min:edisp%nbopt_max,edisp%ispin,kmesh%nkp))
    else
      allocate(dpresp%a_gather(1,1,1,1,1))
    endif
#ifdef MPI
      call MPI_gatherv(dpresp%a_full,(ikend-ikstr+1)*9*(edisp%nbopt_max-edisp%nbopt_min+1)*edisp%ispin, &
                       MPI_DOUBLE_COMPLEX, dpresp%a_gather, rcounts*9*(edisp%nbopt_max-edisp%nbopt_min+1)*edisp%ispin, &
                       displs*9*(edisp%nbopt_max-edisp%nbopt_min+1)*edisp%ispin, MPI_DOUBLE_COMPLEX, master, &
                       MPI_COMM_WORLD, mpierr)
#else
    dpresp%a_gather = dpresp%a_full
#endif

    if (myid .eq. master) then
      do ik = 1,kmesh%nkp
        do iband = edisp%nbopt_min,edisp%nbopt_max
          dpresp%a_sum(:,:,:) = dpresp%a_sum(:,:,:) + dpresp%a_gather(:,:,iband,:,ik) * kmesh%weight(ik)
        enddo
      enddo
    endif


    call cpu_time(tfinish)
    timings(4) = timings(4) + (tfinish - tstart)
    tstart = tfinish

    if (myid .eq. master) then
      call output_data(algo, info, temp, dpresp)
    endif

    call cpu_time(tfinish)
    timings(5) = timings(5) + (tfinish - tstart)
    tstart = tfinish


    deallocate(dpresp%s_gather)
    deallocate(dpresp%a_gather)

  enddo ! end of the outer temperature loop

#ifdef MPI
   call MPI_allreduce(MPI_IN_PLACE, timings,size(timings), MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpierr)
   timings(:3) = timings(:3) / dble(nproc)
#endif
  if (myid .eq. master) then
    write(stdout,*)
    write(stdout,*) '  Timings (average) [s]:'
    write(stdout,'(A21,F16.6)') '    mu-search:       ', timings(1)
    write(stdout,'(A21,F16.6)') '    polygamma-eval:  ', timings(2)
    write(stdout,'(A21,F16.6)') '    response-eval:   ', timings(3)
    write(stdout,'(A21,F16.6)') '    mpi + summation: ', timings(4)
    write(stdout,'(A21,F16.6)') '    hdf5-output:     ', timings(5)
    write(stdout,*)
  endif



!!      if (myid.eq.master) write(*,'(A,3E)') 'XXX ', T,beta,1.1d0/(delta/100.d0-0.64d0*gmax)

  !      !if ( (delta.gt.0.d0).and.(criterion.gt.90.d0).and.(gminall.eq.0.d0) ) then
  !    !if ( (dos%gap .gt. 0.d0).and.(criterion .gt. 90.d0).and.(gminall .eq. 0.d0) ) then

  !    !   if (myid.eq.master) then
  !    !      write(*,*) 'Using low T extrapolation'
  !    !      imeth=3

  !    !      if (iflag_dmudt.eq.0) then
  !    !         !dmudT = ( mu-(emax(iband_valence)+delta/2.d0) ) / T
  !    !         dmudT = ( mu-(dos%vbm+dos%gap/2.d0) ) / T  !eM: let's hope that this actually means the same as what above
  !    !         !         write(*,*) ' XXX '
  !    !         !         write(*,*)mu
  !    !         !         write(*,*) (emax(iband_valence)+delta/2.d0)
  !    !         iflag_dmudt=1
  !    !      endif

  !    !      !      write(*,*)emax(iband_valence),emax(iband_valence)+delta/2.d0,dmudt

  !    !      !mu=emax(iband_valence)+ delta/2.d0 + dmudT * T
  !    !      mu=dos%vbm+ dos%gap/2.d0 + dmudT * T !eM: let's hope that this actually means the same as what above

  !    !   endif

  !   ! endif
  !   ! endif !algo%imurestart==0

  !! old_vers
  !!
  !!      if (1.eq.2) then ! XXX beware... this is meant to work with gamma=0 ... but check...
  !!                    ! XXX for gamma=0, we know that it should be linear down to midgap, so we can actually use that...
  !!                    ! XXX so determine slope and extrapolate... need to think harder for gamma>0
  !!                    ! e.g. evaluate impsi at homo/lumu... and look (analytically?)
  !!
  !!!      if ((delta.gt.0.d0).and.(T/delta.lt.100.d0)) then ! then issues for low gamma... but narrow down by checking more
  !!         !if ((delta.gt.0.d0).and.(criterion.gt.90.d0)) then ! then issues for low gamma... but narrow down by checking more
  !!         if ((dos%gap.gt.0.d0).and.(criterion.gt.90.d0)) then ! then issues for low gamma... but narrow down by checking more
  !!            call ndeviationQ(mu,NE,1,test0)
  !!            call ndeviationQ(mu+1.d-3,NE,1,test1)
  !!            imeth=3
  !!            call find_muQ_lowT(mu,iT)
  !!            !if (myid.eq.master) write(*,'(A,1E20.12)')' using find_muQ_lowT ', T
  !!         endif
  !!
  !!      endif

  !   !endif !master


  !   !if (myid.eq.master) then
  !   !   select case (imurestart)
  !   !   case (0)
  !   !      write(20,'(100E20.12)')T,mu,mutmp,real(imeth,kind=8),criterion
  !!! only needed for gfortran... if one wants to monitor progress by continous readout
  !   !      flush(20)
  !!! maybe need to sync on some systems...    iflush=fsync(fnum(20))
  !   !   case(1)
  !!!      else ! read mu
  !   !      read(20,*)dum1,mu
  !   !      if (dum1.ne.T) STOP 'inconsistent mu.dat'
  !!!      endif
  !   !      niitact=0
  !   !      write(*,'(1F10.3,5E15.7,2I5)')T,mu,beta,criterion/80.d0,real(ndevVQ,8),abs(real(ndevactQ,8)),niitQ,niitact
  !   !   case (2)
  !   !      mu = edisp%efer !this is the original assignment to the value read from file
  !   !                      !unless imurestart == 0, in which case it is the value found from the
  !   !                      !dos (if the tetrahedron methos has been selected)
  !   !      niitact=0
  !   !      write(*,'(1F10.3,5E15.7,2I5)')T,mu,beta,criterion/80.d0,real(ndevVQ,8),abs(real(ndevactQ,8)),niitQ,niitact
  !   !      write(20,'(100E20.12)')T,mu,mutmp,real(imeth,kind=8),criterion
  !   !   end select
  !   !endif !master

!! #ifdef MPI
!!      if (imurestart.ne.0) then ! if read mu, need to BCAST it to all procs.
!!         call MPI_BCAST(mu,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpierr)
!!      endif
!! #endif

  !   !copy the given value of mu into the datastructure
  !   sct%mu(iT) = mu

  !   if (.not. algo%ltetra .and. (myid == master) .and. (iT == sct%nT)) then
  !      call intldos(iT, dos, kmesh, edisp, sct)
  !   endif

  !   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !   ! DONE MU. DO TRANSPORT AT THIS POINT.
  !   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !   if (algo%ldebug) then
  !      call calc_response(mu, iT, drhodT, kmesh, edisp, thdr, sct, dpresp, dderesp, dinter, respBl)
  !   else
  !      call calc_response(mu, iT, drhodT, kmesh, edisp, thdr, sct, dpresp, dderesp, dinter, respBl, qpresp)
  !   endif

  !enddo ! iT temperature
  !! the values of the derivative of the resistivity have been accumulated over
  !! during the temparature loop, now find the critical temperature


!! Finished main loop here
!! output
  !if (myid.eq.master) then
  !  do iT=1,sct%nT-1
  !     write(800,*) sct%TT(iT),sct%d1(iT),sct%d2(iT)
  !     write(801,*) sct%TT(iT),sct%d0(iT)
  !  enddo

  !  dum1=drhodT(1)
  !  do iT=2,sct%nT-1
  !     if(drhodT(iT) > dum1) then
  !        dum1=drhodT(iT)
  !        sct%Tstar=sct%TT(iT)
  !     endif
  !  enddo
  !  ! Lines below to be removed
  !  if (sct%Tstar > 0.0d0) then
  !     write(*,*)'found T* (intraband only)', sct%Tstar
  !  else
  !     write(*,*)'T* (intraband only) not found'
  !  endif

  !  if (sct%Tflat > 0.0d0) then
  !     write(*,*)'found T_s (intraband only)', sct%Tflat
  !  else
  !     write(*,*)'T_s (intraband only) not found'
  !  endif
  !  !end lines to be removed

  !  close(20) ! mu.dat
  !  call response_close_files()
  !endif

  call hdf5_close_file(ifile_energy)
  if (algo%lScatteringFile) then
    call hdf5_close_file(ifile_scatter)
  endif

  call mpi_close()
  call hdf5_finalize()

end program main
