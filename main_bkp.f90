

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





! somehow mpif90 on hclm doesn like the more elegant iso way...
!#define qp 16
!#define dp 8

program gtmain ! GammaTransport






use params
use mpi_org
use estruct, only: delta,gam,gmax,iband_valence,gminall,emin,emax!,ek,vka,vkab,zqp ! ek,zqp,vka,vkab needed for testing only!!! move that
use response
implicit none

real(8) mu,mutmp !chemical potential
integer nT,iT,ierr,imeth,iband,iflag_dmudt,imurestart,niitact
real(8), allocatable :: Tvek(:)
real(16) test0,test1,ndevactQ
real(8) criterion,ndevact
real(8) dmudT ! used only for gamma=0.d0

real(8) dum1,dum2,dum3,dum4,muconst
integer idum

!TEST
!real(8)emin0,emax0,emin,emax,ntest
!integer iband


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11


open(10,file='inp.linretrace',status='old')
read(10,*)nband,nk,vol,icubic
read(10,*)NE
read(10,*)Tmin,Tmax,dT
read(10,*)imurestart, muconst
close(10)



idiag=1 ! only do diagonal components of the tensors (in absence of B)
!icubic=0 ! only do 1st element of the diagonal of the tensor (in absence of B, with B, only do xy)


!sets computer environement. need to test 1 version...
! eM: I think that this subroutine sets the values of ikstart and ikend for k-point parallelisation. Possibly this will have to be replaced by parallelisation over tetrahedra
call mpi_env(nk,small,threshold,smallQ,thresholdQ)

if (myid.eq.master) then
   open(700,file='control.dat',status='unknown')



   write(*,*)
   write(*,*)
   write(*,*)'#####################################################'
   write(*,*)'#  Lin-ReTraCe -- Linear Response Transport Centre  #'
   write(*,*)'#####################################################'
   write(*,*)'#  Jan M. Tomczak                                   #'
   write(*,*)'#####################################################'
   write(*,*)
endif




!read in electronic structure and matrix elements
!I assume (although not critical) that ek is sorted in ASCENDING order in iband
call call_estruct_init()


!determines value of gap in case there is one... requires out.dos output, and ek sorted in ascending order with iband
if (myid.eq.master) then
   if (ivkab.eq.0) then
      write(*,*)'ACHTUNG: no vkab present! will not perform B.ne.0 calculations!'
   endif

   call estruct_determine_gap()
endif

   call MPI_BARRIER( MPI_COMM_WORLD, mpierr )
   call MPI_BCAST(delta,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpierr)
   call MPI_BARRIER( MPI_COMM_WORLD, mpierr )
   call MPI_BCAST(iband_valence,1,MPI_INTEGER8,master,MPI_COMM_WORLD,mpierr)

!   write(*,*)'delta ',delta

!mu=0.d0
! or

if (imurestart.eq.0) then
!select case (imurestart)
!   case (0) ! no starting mu
   if (myid.eq.master) then ! only master knows about first k-point
      call initialize_mu(mu)
      write(*,*)'initialized mu = ',mu
   endif

   call MPI_BARRIER( MPI_COMM_WORLD, mpierr )
   call MPI_BCAST(mu,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpierr)

endif
!   case (2) ! reading existing mu.dat
!   mu=muconst
! #ifdef 1
!   call MPI_BARRIER( MPI_COMM_WORLD, mpierr )
!   call MPI_BCAST(mu,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpierr)
! #endif
!   end select

nT=int((Tmax-Tmin)/dT)+1
allocate(Tvek(nT))
do iT=1,nT
   Tvek(iT)=real(iT-1,8)*dT+Tmin
enddo
call call_estruct_read_gamma(nT,Tvek)
call call_estruct_read_Z()


if (myid.eq.master) then

   select case (imurestart)
      case (0)
!   if (imurestart.eq.0) then
         open(20,file='mu.dat',status='unknown')
         write(20,'(A,3E15.7,1I10,1E20.12)')'# ', Tmin,Tmax,dT,nk,ne
!   else
      case (1) ! read and check
         open(20,file='mu.dat',status='old')
         read(20,'(2X,3E15.7,1I10,1E20.12)')dum1,dum2,dum3,idum,dum4
         write(*,'(A,3E15.7,1I10,1E20.12)')' # ',dum1,dum2,dum3,idum,dum4
         if (dum1.ne.Tmin) STOP 'wrong Tmin in mu.dat'
         if (dum2.ne.Tmax) STOP 'wrong Tmax in mu.dat'
         if (dum3.ne.dT) STOP 'wrong dT in mu.dat'
         if (idum.ne.nk) STOP 'wrong nk in mu.dat'
         if (dum4.ne.ne) STOP 'wrong NE in mu.dat'

! Here, I should be checking for the Gammas, too...
         write(*,*)'using chemical potential from mu.dat'
      case (2) ! use mu from inp.linretrace
         open(20,file='mu.dat',status='unknown')
         write(20,'(A,3E15.7,1I10,1E20.12)')'# ', Tmin,Tmax,dT,nk,ne
         write(20,*)'using constant mu=',mu
         write(*,*)'using constant mu=',mu
      end select
      call response_open_files()
   endif

!find minimal gamma for lowest T. If = 0 then we can do (linear) extrapolation... if gam=O(T,T^2) ... strictly not linear...
gminall=10.d0
do it=1,1 ! lowest T
do iband=1,nband
   if (gam(it,iband).lt.gminall) gminall=gam(it,iband)
enddo
enddo

iflag_dmudt=0
call response_allocate_variables()

if (myid.eq.master) write(*,*)
if (myid.eq.master) write(*,*)
if (myid.eq.master) write(*,*)'T, mu,beta,criterion,ndev,ndevact,niit,niitact'
if (myid.eq.master) write(*,*)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! START TEMPERATURE LOOP
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
do iT=nT,1,-1
!do iT=nT,1,-8
   T=Tvek(iT) !real(iT-1,8)*dT+Tmin
   beta=1.d0/(kB*T)
   beta2p=beta/(2.d0*pi)
!determine maximal gamma of bands involved in gap
   gmax=max(gam(it,iband_valence),gam(it,iband_valence+1))
   criterion=1.1d0/beta+0.64d0*gmax ! considerations of FWHM, etc...

   if (delta.gt.0.d0) then
      criterion=delta/criterion
   else
      criterion=beta/20.d0
!    write(*,*) 'XXX AAAAAAAAAAAAAAAAAa'
!    STOP     
   endif


!find the chemical potential through secant root-finding

! using FWHM of ImDigamma Funktion... see notes:
! require here FWHM >=delta/20 for DP to be sufficient
! XXX be careful in the case of metals... then delta makes no sense...


   if (imurestart.eq.0) then
!   case (0)
      !if (imurestart.eq.0) then
      if (criterion.lt.20.d0) then !DP
         !   if (T.gt.150.d0) then
         imeth=0
         call find_mu(mu,iT,nT,ndev,ndevact,niitact)
         if (myid.eq.master) then
            write(*,'(1F10.5,5E15.7,2I5)')T,mu,beta,criterion,ndev,abs(ndevact),niit,niitact
            write(700,'(1F10.5,5E15.7,2I5)')T,mu,beta,criterion,ndev,abs(ndevact),niit,niitact
         endif
      elseif (criterion.lt.80.d0) then !QP
         imeth=1
         if (myid.eq.master) write(*,*) 'QUAD'
         !      call find_muDPQ(mu,iT,nT,ndevDQ) !internal QUAD precision
         call find_muQ(mu,iT,nT,ndevQ,ndevactQ,niitact) ! full QUAD on particle number
         if (myid.eq.master) then
            write(*,'(1F10.5,5E15.7,2I5)')T,mu,beta,criterion,real(ndevQ,8),abs(real(ndevactQ,8)),niitQ,niitact
            write(700,'(1F10.5,5E15.7,2I5)')T,mu,beta,criterion,real(ndevQ,8),abs(real(ndevactQ,8)),niitQ,niitact
         endif
      else
         ! further refinement
         !      if (criterion.gt.80.d0) then
         imeth=2
         if (myid.eq.master) write(*,*) 'SUPER QUAD'
!!!      niitQ=niitQ+niitQ/5
         call find_muQ(mu,it,nt,ndevVQ,ndevactQ,niitact)
         if (myid.eq.master) then
            write(*,'(1F10.3,5E15.7,2I5)')T,mu,beta,criterion/80.d0,real(ndevVQ,8),abs(real(ndevactQ,8)),niitQ,niitact
            write(700,'(1F10.3,5E15.7,2I5)')T,mu,beta,criterion/80.d0,real(ndevVQ,8),abs(real(ndevactQ,8)),niitQ,niitact
         endif
         !      endif
      endif


      mutmp=mu

!   if (myid.eq.master) write(*,'(A,3E)') 'XXX ', T,beta,1.1d0/(delta/100.d0-0.64d0*gmax)

      if ( (delta.gt.0.d0).and.(criterion.gt.90.d0).and.(gminall.eq.0.d0) ) then

         if (myid.eq.master) then
            write(*,*) 'Using low T extrapolation'
            imeth=3
            
            if (iflag_dmudt.eq.0) then
               dmudT = ( mu-(emax(iband_valence)+delta/2.d0) ) / T 
               !         write(*,*) ' XXX '
               !         write(*,*)mu
               !         write(*,*) (emax(iband_valence)+delta/2.d0)
               iflag_dmudt=1
            endif
            
            !      write(*,*)emax(iband_valence),emax(iband_valence)+delta/2.d0,dmudt
            
            mu=emax(iband_valence)+ delta/2.d0 + dmudT * T
            
         endif

         call MPI_BCAST(mu,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpierr)



      endif


      if (1.eq.2) then ! XXX beware... this is meant to work with gamma=0 ... but check... 
                    ! XXX for gamma=0, we know that it should be linear down to midgap, so we can actually use that...
                    ! XXX so determine slope and extrapolate... need to think harder for gamma>0
                    ! e.g. evaluate impsi at homo/lumu... and look (analytically?)

!      if ((delta.gt.0.d0).and.(T/delta.lt.100.d0)) then ! then issues for low gamma... but narrow down by checking more
         if ((delta.gt.0.d0).and.(criterion.gt.90.d0)) then ! then issues for low gamma... but narrow down by checking more
            call ndeviationQ(mu,NE,1,test0)
            call ndeviationQ(mu+1.d-3,NE,1,test1)
            imeth=3
            call find_muQ_lowT(mu,iT)
            if (myid.eq.master) write(*,'(A,1E20.12)')' using find_muQ_lowT ', T
         endif
         
      endif

   endif


!XXX ACHTUNG
!mu=0.d0

   if (myid.eq.master) then
      select case (imurestart)
      case (0)
!      if (imurestart.eq.0) then ! write mu
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
         if (myid.eq.master) then
            write(*,'(1F10.3,5E15.7,2I5)')T,mu,beta,criterion/80.d0,real(ndevVQ,8),abs(real(ndevactQ,8)),niitQ,niitact
            write(700,'(1F10.3,5E15.7,2I5)')T,mu,beta,criterion/80.d0,real(ndevVQ,8),abs(real(ndevactQ,8)),niitQ,niitact
         endif

      case (2)
         mu=muconst
         niitact=0

         if (myid.eq.master) then
            write(*,'(1F10.3,5E15.7,2I5)')T,mu,beta,criterion/80.d0,real(ndevVQ,8),abs(real(ndevactQ,8)),niitQ,niitact
            write(700,'(1F10.3,5E15.7,2I5)')T,mu,beta,criterion/80.d0,real(ndevVQ,8),abs(real(ndevactQ,8)),niitQ,niitact
            write(20,'(100E20.12)')T,mu,mutmp,real(imeth,kind=8),criterion
         endif

      end select

   endif
   if (imurestart.ne.0) then ! if read mu, need to BCAST it to all procs.

      call MPI_BCAST(mu,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpierr)

   endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! DONE MU. GO TRANSPORT AT THIS POINT.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! eM: if the tetrahedron method is used at this point it is necessary to re-evaluate the linear interpolation over the thetrahedra for the quantity you need for the response

! put the below into a function... also prepare a QUAD precision routine...
! low-T sigma is very bad... check with constants and absorption errors.

!   mu=1.6d0

!   if (criterion.lt.20.d0) then !DP
      call calc_response(mu,it) ! writes to files
!   else
      call calc_responseQ(mu,it) ! writes to Q-files
!   endif

      call calc_response_BoltzmannQ(mu,it)  ! writes to B-files, but reuses Q-variables from responseQ...
! Boltzmann-approximation in transport kernels, i.e. digamma --> Fermi and using leading order in 1/Gamma
! I could just set gamma=0 in the digamma functions, but its faster when using Fermi functions...      

      if (myid.eq.master) write(*,*)
enddo ! iT temperature


if (myid.eq.master) then
   close(20)
   close(700)



   call response_close_files()
endif



call mpi_close()


end program gtmain



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

