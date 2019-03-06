program main
  use Mparams
  use Maux
  use Mdos
  use Mtypes
  use Mmpi_org
  use Minput
  implicit none

  type(algorithm)  :: algo
  type(kpointmesh) :: kmesh   ! contains k-point mesh specifiers and logical switches on how to get the mesh from
  type(energydisp) :: edisp   ! contains the band dispersion energy and the optical matrix elements
  type(dosgrid)    :: dos     ! DOS, integrated DOS, Fermi level
  type(scattering) :: sct     ! temperature grid and scattering rate (only the coefficients, NO BAND DEPENDENCE)
  type(lattice)    :: lat


  type(response_dp) :: dpresp  ! response double precision
  type(response_qp) :: qpresp  ! response quadruple precision
  type(response_dp) :: respBl  ! response Boltzman regime
  type(response_dp) :: dinter  ! response interband
  type(response_dp) :: dderesp ! response functions (intraband conductivity) derivatives in double precision

  real(8)              :: gmax, gminall
  real(8)              :: mu
  integer              :: nT,iT,ierr,imeth,iflag_dmudt,imurestart,niitact
  real(16)             :: test0,test1,ndevactQ
  real(8)              :: criterion,ndevact
  real(8)              :: dmudT ! used only for gamma=0.d0
  real(8), allocatable :: drhodT(:)

  real(8) :: dum1,dum2,dum3,dum4
  integer :: idum, i
  integer :: nk
  integer :: ig, iband, ik !counters for polynomial gamma(T), band and k-point index
  integer :: vb, ib   !valence band (at Gamma point) and auxiliary counter


  call mpi_initialize()
  if (myid.eq.master) call main_greeting(stdout)
  call mpi_barrier(mpi_comm_world, mpierr)

  call read_config(algo, edisp)
  call read_preproc_energy_data(algo, kmesh, edisp, lat)

  if (algo%lBfield .and. .not. edisp%lDerivatives) then
    call stop_with_message(5, 'Energy derivatives required for Bfield quantities')
  endif

  if(myid .eq. master .and. .not.(algo%lBfield)) then
     write(*,*)'LINRETRACE will not perform calculations with magnetic field'
  else
     write(*,*)'LINRETRACE will perform calculations with magnetic field'
  endif


  call gendosel (kmesh, edisp, dos) ! normalization already taken care of
  call findef(dos, edisp)   ! finds the (non-interacting) Fermi level

  if (myid.eq.master) then
    if (dos%gap == 0.d0) then
      write(*,*) 'Detected no band gap'
    else
      write(*,*) 'Detected band gap of size: ', dos%gap
    endif
  endif

  ! starting point for the chemical potential
  if (myid.eq.master) then
     if (algo%muMethod == 0) then
        write(*,*) 'initialized LDA mu = ', edisp%efer
     else
        edisp%efer = edisp%mu ! we overwrite the calculated fermienergy with the provided chem.pot
        write(*,*) 'mu read from file = ',  edisp%efer
     endif
  endif

  ! ! construct temperature grid
  ! sct%nT=int((sct%Tmax-sct%Tmin)/sct%dT)+1
  ! allocate(sct%TT(sct%nT))
  ! allocate(sct%mu(sct%nT))
  ! allocate(sct%gam(sct%nT))
  ! do iT=1,sct%nT
  !    sct%TT(iT)=real(iT-1,8)*sct%dT+sct%Tmin
  ! enddo


  ! sct%gam=0.d0
  ! do iT=1,sct%nT
  !    do ig=0,sct%ng
  !       sct%gam(iT)=sct%gam(iT) + sct%gc(ig)*(sct%TT(iT)**ig)
  !    enddo
  ! enddo

  call mpi_genkstep(kmesh%nkp)

  if (algo%ldebug) then
     write(*,*) "MPI: myid: ", myid, "iqstr: ", iqstr, "iqend: ", iqend
  endif

  !! allocate the arrays once outside of the main (temperature) loop
  !call dpresp_alloc(algo%lBfield, dpresp, edisp%nband_max)
  !call dpresp_alloc(algo%lBfield, respBl, edisp%nband_max)
  !call dpresp_alloc(.false., dderesp, edisp%nband_max)
  !call dpresp_alloc(.false., dinter, edisp%nband_max)
  !if (.not. algo%ldebug) then
  !   call qpresp_alloc(algo%lBfield, qpresp, edisp%nband_max)
  !endif

  !if (myid.eq.master) then
  !   if (algo%ldebug) write(*,*)'MAIN: debug mode'
  !   write(*,*)
  !   write(*,*)
  !   write(*,*)'T, mu,beta,criterion,ndev,ndevact,niit,niitact'
  !   write(*,*)
  !   write(*,*) 'gap', dos%gap
  !endif


  !!initialise the Tstar variable
  !sct%Tstar = 0.0d0
  !sct%Tflat = 0.0d0
  !allocate(drhodT(sct%nT))
  !drhodT = 0.0d0


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! START TEMPERATURE LOOP
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !do iT=sct%nT,1,-1
  !   T=sct%TT(iT)
  !   beta=1.d0/(kB*T)
  !   betaQ=1.q0/(kBQ*T)
  !   beta2p=beta/(2.d0*pi)
  !   beta2pQ=betaQ/(2.q0*piQ)
  !   !determine maximal gamma of bands involved in gap
  !   !gmax=max(gam(it,iband_valence),gam(it,iband_valence+1))
  !   gmax=sct%gam(iT)

  !   ! if (algo%ldmft) gmax=max(sct%ykb(iT,1,vb),sct%ykb(iT,1,vb+1))
  !   ! criterion=1.1d0/beta+0.64d0*gmax ! considerations of FWHM, etc...
  !   !write(*,*) 'criterion/1', criterion
  !   !write(*,*) 'gap', dos%gap

  !   !if (dos%gap.gt.0.d0) then
  !   !   !criterion=dos%gap/criterion
  !   !   criterion=beta/20.d0
  !   !else
  !   !endif

  !   criterion=beta/20.d0

  !   !write(*,*) 'criterion/2', criterion

  !! using FWHM of ImDigamma Funktion... see notes:
  !! require here FWHM >=delta/20 for DP to be sufficient
  !! XXX be careful in the case of metals... then delta makes no sense...

  !   !if (myid.eq.master) then
  !   !also the chemical potential search is not parallelised at the moment
  !   !don't know how difficult it would be to revert to the original implementation
  !   !and/or to implement the tetrahedron method
  !   if (algo%imurestart == 0) then
  !      if (criterion.lt.20.d0) then !DP
  !         call find_mu(mu,iT,ndev,ndevact,niitact, edisp, sct, kmesh, thdr)
  !         if (myid.eq.master) then
  !            write(*,'(1F10.5,5E15.7,2I5)')T,mu,beta,criterion,ndev,abs(ndevact),niit,niitact
  !         endif
  !      elseif (criterion.lt.80.d0) then !QP
  !         call find_mu(mu,iT,ndevQ,ndevactQ,niitact, edisp, sct, kmesh, thdr)! full QUAD on particle number
  !         if (myid.eq.master) then
  !            write(*,'(1F10.5,5E15.7,2I5)')T,mu,beta,criterion,real(ndevQ,8),abs(real(ndevactQ,8)),niitQ,niitact
  !         endif
  !      else   ! further refinement
  !         if (myid.eq.master) write(*,*) 'SUPER QUAD'
  !         call find_mu(mu,iT,ndevVQ,ndevactQ,niitact, edisp, sct, kmesh, thdr)! full QUAD on particle number
  !         if (myid.eq.master) then
  !            write(*,'(1F10.3,5E15.7,2I5)')T,mu,beta,criterion/80.d0,real(ndevVQ,8),abs(real(ndevactQ,8)),niitQ,niitact
  !         endif
  !      endif !criterion
  !   else
  !      write(*,'(1F10.5,5E15.7,2I5)')T,mu,beta
  !   endif

  !   ! mutmp=mu ! possibly unused

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

  call mpi_close()

end program main
