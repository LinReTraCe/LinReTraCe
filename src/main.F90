

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! TRANSPORT CODE for finite static scattering rates.
! Jan M. Tomczak, rewritten from scratch, after N. Berlakovich's project work 2015
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! BUGS
!
! - FIXED. low T mu is different when using several nodes (openMP 12 cores is ok)
! - FIXED. any parallelization fails for sigmaB/peltierB ? values different.
! - FIXED. Seebeck missing a beta, not just constants...
!
! - 1 "ERROR" was caused by NaN problem of the Cluster. Program worked after node in question was rebooted.
! - However: check the DEBUG option in QUAD_REDUCE when called from responseQ. Results reasonable, but error message... does this debug option qp0,qp1 work???
!
!
! TO DO:
! - DONE. NO EXTRA TERM. Nernst coefficient: extra term???
! - automatic detection of precision limits, e.g. precision of mu, differentiation for metals and insulators... in a metal ndev=1.d-15 ok
! - XXXXXXX
! - check debug_? for QUAD precision. numbers from responseQ are at times (metals!!) ENORMOUS and so errors large on absolute scale! Separation into small and big ABSOLUTE value?
! - COMPARE AGAIN mpi vs single proc runs... CASE of METALS
! - TEST initialize mu for 3 bands with partial occupation...
! - in case of T-independent Gamma, compute spectral function
! - yy, zz in sigma Ausgabe !!
!
!
! FUTURE:
! - For comparison: Mott formula description.
! - interface with wien2k, use vol.f90 and read_struct from wientrans/wienkernel for Volumes
! - offdiagonal elements in Akw or vk ? --> ab initio
! - specific heat a la Miyake and Tsuruta JPSJ 84, 094708 (2015) [to be generalized to T-dependent mu]
! - k-mesh !
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program gtmain ! GammaTransport

  use Mparams
  use Mtypes
  use Mmpi_org
  use Mestruct
  use Mresponse
  use Mroot
  implicit none

  !!eM note: it is necessary to declare dderesp with the extended datatype because by doing so in interptra_mu there is
  !! no extra multiplication for some additional factors (required for the conductivity); the additional memory requirement
  !! is negligible
  type(kpointmesh), target :: irrkm   ! contains k-point mesh specifiers and logical switches on how to get the mesh from
  type(kpointmesh), target :: redkm
  type(kpointmesh), target :: fulkm
  type(edisp), target      :: eirrk   ! contains the band dispersion energy and the optical matrix elements
  type(edisp), target      :: eredk
  type(edisp), target      :: efulk
  type(tetramesh)          :: thdr    ! contains the tetrahedra
  type(dosgrid)            :: dos     ! DOS, integrated DOS, Fermi level
  type(scatrate)           :: sct     ! temperature grid and scattering rate (only the coefficients, NO BAND DEPENDENCE)
  type(dp_resp)            :: dpresp  ! response functions in double precision
  type(qp_resp)            :: qpresp  ! response functions in quadruple precision -- same as equations as used in dpresp
  type(dp_resp)            :: respBl  ! response functions in double precision for Boltzmann regime response ( i.e. Gamma -> 0 )
  type(dp_respinter)       :: dinter  ! response functions in double precision for interband transitions
  type(dp_respinter)       :: dderesp ! response functions (intraband conductivity) derivatives in double precision

  ! mP note
  ! after we initialized all the kpointer, optical elements, energy dispersion
  ! we just point to the required variable
  class(kpointmesh), pointer :: kpointer
  class(edisp), pointer      :: epointer

  real(8)              :: gmax, gminall
  real(8)              :: mu,mutmp !chemical potential
  integer              :: nT,iT,ierr,imeth,iflag_dmudt,imurestart,niitact
  real(16)             :: test0,test1,ndevactQ
  real(8)              :: criterion,ndevact
  real(8)              :: dmudT ! used only for gamma=0.d0
  real(8), allocatable :: drhodT(:)

  real(8) :: dum1,dum2,dum3,dum4,muconst !eM 13.03.2018: muconst is likely obsolete
  integer :: idum, i
  integer :: nk
  integer :: ig, iband, ik !counters for polynomial gamma(T), band and k-point index
  integer :: vb, ib   !valence band (at Gamma point) and auxiliary counter

  ! mP: should possibly go into the algo datatype
  ! method 0: secant; method 1: linint; method 2: Riddler; method 3: bisection
  integer, parameter   :: method = 2

  call mpi_initialize()

  if (myid.eq.master) then
     write(*,*)
     write(*,*)
     write(*,*)'#####################################################'
     write(*,*)'#  Lin-ReTraCe -- Linear Response Transport Centre  #'
     write(*,*)'#####################################################'
     write(*,*)'#  E. Maggio, M. Pickem and J.M. Tomczak            #'
     write(*,*)'#####################################################'
     write(*,*)
  endif

  call mpi_barrier(mpi_comm_world, mpierr)

  ! with this flag set to false the quad precision response is computed
  ! currently in developing / debugging mode
  algo%ldebug = .true.
  ! algo%lgenred = .false.
  algo%lgenred = .true.

  ! we read the information into the irreducible datatypes
  call read_config(irrkm, eirrk, sct)
  ! based on the input we setup the other kind of grids
  call setup_algo(irrkm, redkm, fulkm, eirrk, eredk, efulk)
  !read in electronic structure and matrix elements
  call estruct_init(irrkm, redkm, fulkm, eirrk, eredk, efulk, thdr, dos, sct)

  ! now we either have the data in the full BZ or in the reducible one
  if (algo%ltetra) then
     epointer => efulk
     kpointer => fulkm
  else
     epointer => eredk
     kpointer => redkm
  endif

  if(myid .eq. master .and. .not.(algo%lBfield)) then
     write(*,*)'LINRETRACE will not perform calculations with magnetic field'
  endif

  !output of DOS NOS
  if (algo%ldebug .and. myid .eq. master) then
     open(10,file='dos',status='unknown')
     do i=1,size(dos%enrg)
        write (10,*) dos%enrg(i), dos%dos(i), dos%nos(i)
     enddo
     close(10)
  endif

  ! starting point for the chemical potential
  mu = epointer%efer
  if (myid.eq.master) then
     if (algo%imurestart == 0) then
        write(*,*) 'initialized LDA mu = ',mu
     else
        write(*,*) 'mu read from file = ', mu
     endif
  endif

  ! construct temperature grid
  sct%nT=int((sct%Tmax-sct%Tmin)/sct%dT)+1
  allocate(sct%TT(sct%nT))
  allocate(sct%mu(sct%nT))
  do iT=1,sct%nT
     sct%TT(iT)=real(iT-1,8)*sct%dT+sct%Tmin
  enddo

  ! construct scattering rate temperature dependence
  ! (these variables are in the scatrate type)
  ! k-point and band dependence are included in the
  ! intrinsic scattering rate ykb (it inherits them from Im{Sigma})
  ! whereas sct%gam has only T dependence
  allocate(sct%gam(sct%nT)) ! always there, value read from input file
  if (algo%ldmft) then
     if (algo%ltetra) then
        allocate(sct%ykb(sct%nT, thdr%ntet, epointer%nband_max))
     else
        allocate(sct%ykb(sct%nT, kpointer%ktot, epointer%nband_max))
     endif
  endif


  sct%gam=0.d0
  do iT=1,sct%nT
     do ig=0,sct%ng
        sct%gam(iT)=sct%gam(iT) + sct%gc(ig)*(sct%TT(iT)**ig)
     enddo
  enddo

  !intrinsic scattering rate
  !trivial T-dependence at the moment
  !TODO: FIX THIS FOR TETRAHEDRONS
  if (algo%ldmft) then
     sct%ykb=0.0d0
     do iband=1,epointer%nband_max
        do ik=1,kpointer%ktot
           do iT=1,sct%nT
              sct%ykb(iT,ik,iband)=epointer%Im(ik,iband)
           enddo
        enddo
     enddo
  endif

  ! if (algo%ldebug .and. myid.eq.master) then
  !    open(10,file='GamT.dat', status='unknown')
  !    do iT=sct%nT,1,-1
  !       write(10,'(100E20.12)')sct%TT(iT),(sct%gam(iT,iband),iband=1,nband)
  !    enddo
  !    close(10)
  ! endif

  ! distribute k-points
  if (algo%ltetra) then
     call mpi_genkstep(thdr%ntet)
  else
     call mpi_genkstep(kpointer%ktot)
  endif

  if (algo%ldebug) then
     write(*,*) "MPI: myid: ", myid, "iqstr: ", iqstr, "iqend: ", iqend
  endif

  if (myid.eq.master) then
     select case (algo%imurestart)
        case (0)
           open(20,file='mu.dat',status='unknown')
           write(20,'(A,3E15.7,1I10,1E20.12)')'# ', sct%Tmin, sct%Tmax, sct%dT,kpointer%ktot,epointer%nelect
        case (1) ! read and check
           open(20,file='mu.dat',status='old')
           read(20,'(2X,3E15.7,1I10,1E20.12)')dum1,dum2,dum3,idum,dum4
           write(*,'(A,3E15.7,1I10,1E20.12)')' # ',dum1,dum2,dum3,idum,dum4
           if (dum1.ne.sct%Tmin) STOP 'wrong Tmin in mu.dat'
           if (dum2.ne.sct%Tmax) STOP 'wrong Tmax in mu.dat'
           if (dum3.ne.sct%dT) STOP   'wrong dT in mu.dat'
           if (idum.ne.kpointer%ktot) STOP 'wrong nk in mu.dat'
           if (dum4.ne.epointer%nelect) STOP 'wrong NE in mu.dat'
           ! TODO: Here, I should be checking for the Gammas, too...
           write(*,*)'using chemical potential from mu.dat'
        case (2) ! use mu from inp.linretrace
           open(20,file='mu.dat',status='unknown')
           write(20,'(A,3E15.7,1I10,1E20.12)')'# ', sct%Tmin, sct%Tmax, sct%dT, kpointer%ktot, epointer%nelect
           write(20,*)'using constant mu=',mu
           write(*,*)'using constant mu=',mu
     end select
     call response_open_files()
  endif !master


  !find minimal gamma for lowest T. If = 0 then we can do (linear) extrapolation... if gam=O(T,T^2) ... strictly not linear...
  gminall=10.d0
  if (algo%ldmft) then
     do iband=1,epointer%nband_max
        if ((sct%ykb(1,1,iband)<gminall).and.(sct%ykb(1,1,iband)>0.0d0)) &
               gminall=sct%ykb(1,1,iband)
     enddo
  else
  gminall=sct%gam(1)
  endif

  iflag_dmudt=0 !!!what is this for?

  ! allocate the arrays once outside of the main (temperature) loop
  call dpresp_alloc(algo%lBfield, dpresp, epointer%nband_max)
  call dpresp_alloc(algo%lBfield, respBl, epointer%nband_max)
  call dpresp_alloc(.false., dderesp, epointer%nband_max)
  call dpresp_alloc(.false., dinter, epointer%nband_max)
  if (.not. algo%ldebug) then
     call qpresp_alloc(algo%lBfield, qpresp, epointer%nband_max)
  endif

  if (myid.eq.master) then
     if (algo%ldebug) write(*,*)'MAIN: debug mode'
     write(*,*)
     write(*,*)
     write(*,*)'T, mu,beta,criterion,ndev,ndevact,niit,niitact'
     write(*,*)
     write(*,*) 'gap', dos%gap
  endif


  !initialise the Tstar variable
  sct%Tstar = 0.0d0
  sct%Tflat = 0.0d0
  allocate(drhodT(sct%nT))
  drhodT = 0.0d0


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! START TEMPERATURE LOOP
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  do iT=sct%nT,1,-1
     T=sct%TT(iT)
     beta=1.d0/(kB*T)
     betaQ=1.q0/(kBQ*T)
     beta2p=beta/(2.d0*pi)
     beta2pQ=betaQ/(2.q0*piQ)
     !determine maximal gamma of bands involved in gap
     !gmax=max(gam(it,iband_valence),gam(it,iband_valence+1))
     gmax=sct%gam(iT)
     if (algo%ldmft) gmax=max(sct%ykb(iT,1,vb),sct%ykb(iT,1,vb+1))
     criterion=1.1d0/beta+0.64d0*gmax ! considerations of FWHM, etc...
     !write(*,*) 'criterion/1', criterion
     !write(*,*) 'gap', dos%gap

     if (dos%gap.gt.0.d0) then
        !criterion=dos%gap/criterion
        criterion=beta/20.d0
     else
        criterion=beta/20.d0
     endif

     !write(*,*) 'criterion/2', criterion

  ! using FWHM of ImDigamma Funktion... see notes:
  ! require here FWHM >=delta/20 for DP to be sufficient
  ! XXX be careful in the case of metals... then delta makes no sense...

     !if (myid.eq.master) then
     !also the chemical potential search is not parallelised at the moment
     !don't know how difficult it would be to revert to the original implementation
     !and/or to implement the tetrahedron method
     if (algo%imurestart == 0) then
        if (criterion.lt.20.d0) then !DP
           imeth=0
           call find_mu(mu,iT,ndev,ndevact,niitact, epointer, sct, kpointer, thdr, method)
           if (myid.eq.master) then
              write(*,'(1F10.5,5E15.7,2I5)')T,mu,beta,criterion,ndev,abs(ndevact),niit,niitact
           endif
        elseif (criterion.lt.80.d0) then !QP
           imeth=1
           call find_mu(mu,iT,ndevQ,ndevactQ,niitact, epointer, sct, kpointer, thdr, method)! full QUAD on particle number
           if (myid.eq.master) then
              write(*,'(1F10.5,5E15.7,2I5)')T,mu,beta,criterion,real(ndevQ,8),abs(real(ndevactQ,8)),niitQ,niitact
           endif
        else   ! further refinement
           imeth=2
           if (myid.eq.master) write(*,*) 'SUPER QUAD'
           call find_mu(mu,iT,ndevVQ,ndevactQ,niitact, epointer, sct, kpointer, thdr, method)! full QUAD on particle number
           if (myid.eq.master) then
              write(*,'(1F10.3,5E15.7,2I5)')T,mu,beta,criterion/80.d0,real(ndevVQ,8),abs(real(ndevactQ,8)),niitQ,niitact
           endif
        endif !criterion

        mutmp=mu ! possibly unused

  !      if (myid.eq.master) write(*,'(A,3E)') 'XXX ', T,beta,1.1d0/(delta/100.d0-0.64d0*gmax)

           !if ( (delta.gt.0.d0).and.(criterion.gt.90.d0).and.(gminall.eq.0.d0) ) then
        if ( (dos%gap .gt. 0.d0).and.(criterion .gt. 90.d0).and.(gminall .eq. 0.d0) ) then

           if (myid.eq.master) then
              write(*,*) 'Using low T extrapolation'
              imeth=3

              if (iflag_dmudt.eq.0) then
                 !dmudT = ( mu-(emax(iband_valence)+delta/2.d0) ) / T
                 dmudT = ( mu-(dos%vbm+dos%gap/2.d0) ) / T  !eM: let's hope that this actually means the same as what above
                 !         write(*,*) ' XXX '
                 !         write(*,*)mu
                 !         write(*,*) (emax(iband_valence)+delta/2.d0)
                 iflag_dmudt=1
              endif

              !      write(*,*)emax(iband_valence),emax(iband_valence)+delta/2.d0,dmudt

              !mu=emax(iband_valence)+ delta/2.d0 + dmudT * T
              mu=dos%vbm+ dos%gap/2.d0 + dmudT * T !eM: let's hope that this actually means the same as what above

           endif

        endif
     endif !algo%imurestart==0

  ! old_vers
  !
  !      if (1.eq.2) then ! XXX beware... this is meant to work with gamma=0 ... but check...
  !                    ! XXX for gamma=0, we know that it should be linear down to midgap, so we can actually use that...
  !                    ! XXX so determine slope and extrapolate... need to think harder for gamma>0
  !                    ! e.g. evaluate impsi at homo/lumu... and look (analytically?)
  !
  !!      if ((delta.gt.0.d0).and.(T/delta.lt.100.d0)) then ! then issues for low gamma... but narrow down by checking more
  !         !if ((delta.gt.0.d0).and.(criterion.gt.90.d0)) then ! then issues for low gamma... but narrow down by checking more
  !         if ((dos%gap.gt.0.d0).and.(criterion.gt.90.d0)) then ! then issues for low gamma... but narrow down by checking more
  !            call ndeviationQ(mu,NE,1,test0)
  !            call ndeviationQ(mu+1.d-3,NE,1,test1)
  !            imeth=3
  !            call find_muQ_lowT(mu,iT)
  !            !if (myid.eq.master) write(*,'(A,1E20.12)')' using find_muQ_lowT ', T
  !         endif
  !
  !      endif

     !endif !master


     if (myid.eq.master) then
        select case (imurestart)
        case (0)
           write(20,'(100E20.12)')T,mu,mutmp,real(imeth,kind=8),criterion
  ! only needed for gfortran... if one wants to monitor progress by continous readout
           flush(20)
  ! maybe need to sync on some systems...    iflush=fsync(fnum(20))
        case(1)
  !      else ! read mu
           read(20,*)dum1,mu
           if (dum1.ne.T) STOP 'inconsistent mu.dat'
  !      endif
           niitact=0
           write(*,'(1F10.3,5E15.7,2I5)')T,mu,beta,criterion/80.d0,real(ndevVQ,8),abs(real(ndevactQ,8)),niitQ,niitact
        case (2)
           mu = eirrk%efer !this is the original assignment to the value read from file
                           !unless imurestart == 0, in which case it is the value found from the
                           !dos (if the tetrahedron methos has been selected)
           niitact=0
           write(*,'(1F10.3,5E15.7,2I5)')T,mu,beta,criterion/80.d0,real(ndevVQ,8),abs(real(ndevactQ,8)),niitQ,niitact
           write(20,'(100E20.12)')T,mu,mutmp,real(imeth,kind=8),criterion
        end select
     endif !master

#ifdef MPI
     if (imurestart.ne.0) then ! if read mu, need to BCAST it to all procs.
        call MPI_BCAST(mu,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpierr)
     endif
#endif

     !copy the given value of mu into the datastructure
     sct%mu(iT) = mu

     if ((myid == master) .and. (iT == sct%nT)) then
        call intldos(iT, dos, redkm, eredk, sct)
     endif

     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     ! DONE MU. DO TRANSPORT AT THIS POINT.
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     if (algo%ldebug) then
        call calc_response(mu, iT, drhodT, kpointer, epointer, thdr, sct, dpresp, dderesp, dinter, respBl)
     else
        call calc_response(mu, iT, drhodT, kpointer, epointer, thdr, sct, dpresp, dderesp, dinter, respBl, qpresp)
     endif

  enddo ! iT temperature
  ! the values of the derivative of the resistivity have been accumulated over
  ! during the temparature loop, now find the critical temperature


! Finished main loop here
! output
  if (myid.eq.master) then
    do iT=1,sct%nT-1
       write(800,*) sct%TT(iT),sct%d1(iT),sct%d2(iT)
       write(801,*) sct%TT(iT),sct%d0(iT)
    enddo

    dum1=drhodT(1)
    do iT=2,sct%nT-1
       if(drhodT(iT) > dum1) then
          dum1=drhodT(iT)
          sct%Tstar=sct%TT(iT)
       endif
    enddo
    ! Lines below to be removed
    if (sct%Tstar > 0.0d0) then
       write(*,*)'found T* (intraband only)', sct%Tstar
    else
       write(*,*)'T* (intraband only) not found'
    endif

    if (sct%Tflat > 0.0d0) then
       write(*,*)'found T_s (intraband only)', sct%Tflat
    else
       write(*,*)'T_s (intraband only) not found'
    endif
    !end lines to be removed

    close(20) ! mu.dat
    call response_close_files()
  endif

  call mpi_close()

end program gtmain
