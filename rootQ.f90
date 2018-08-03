
subroutine find_muQ(mu,iT,nT,dev,target_zero,niitact, ek, sct, mesh, thdr, ltetra)
 ! uses QUAD-PRECISION for particle numbers    
  use params
  use types
  use mpi_org

  implicit none
  ! passed variables
  real(16) target_zero, dev !QUAD    
  real(8) mu
  integer iT, nT, niitact
  type(edisp) :: ek
  type(scatrate) :: sct
  type(kpointmesh) :: mesh
  type(tetramesh)  :: thdr
  logical :: ltetra
  ! local variables
  real(16) target_zero1, target_zero2
  real(8) mu1, mu2, dmu
  integer iit, niit0, itest
  logical lsecant  ! selects the secant root finding algorithm
  ! linear interpolation method
  real(16), allocatable :: Y(:) !array containig the function to minimise 
  real(8), allocatable :: X(:) !array containig the chemical potential
  integer :: nmu  ! number of points that sample the mu interval (mu1,mu2)
  real(8) :: a11, a22, a31, a42
  real(8) :: A(4,4), B(4)
  integer :: i, j
  integer :: ipiv(4)
  integer :: ierr
  logical linint  ! selects the linear interpolation method
  ! Ridders' method
  real(16) :: F(4)
  real(8) :: P(4)
  real(16) :: s 
  real(8) :: psave, ptol  
  integer  :: maxiter ! maximum number of iterations
  logical  :: lridd   ! selects Ridders' method

! deviation from set particle number with initial mu
  call ndeviationQ(mu, iT, ek, sct, mesh, thdr, ltetra, target_zero1)

  !!!!!!!!!!!!!!!!!TEST 
  !write(*,*) 'find_muQ, ndeviation1',myid,target_zero1
  !if (ltetra) then
  !   write(*,*) 'kpoint',mesh%k_coord(:,thdr%idtet(2,1))
  !   write(*,*) 'energy',ek%band(thdr%idtet(2,1), 30)
  !   write(*,*) myid,mu
  !else
  !   write(*,*) 'kpoint',mesh%k_coord(:,2)
  !   write(*,*) 'energy',ek%band(2, 30)
  !   write(*,*) myid,mu
  !endif
  !!!!!!!!!!!!!!!!!TEST END
  target_zero2=target_zero1
!coarse initialization of secant bracket mu1, mu2... Secant doesnt need mu to lie within bracket, but here it does                      
  mu1=mu
  mu2=mu

  dmu=0.005d0
  if (iT.eq.nT) dmu=0.025d0

  do while (target_zero2.gt.0.q0)
     mu2=mu2+dmu
     call ndeviationQ(mu2, iT, ek, sct, mesh, thdr, ltetra, target_zero2)
  enddo
  do while (target_zero1.le.0.q0)
     mu1=mu1-dmu
     call ndeviationQ(mu1, iT, ek, sct, mesh, thdr, ltetra, target_zero1)
  enddo

  niit0=niitQ
!  if (iT.eq.nT) niit0=niitQ*2

  if (dev.eq.ndevVQ) niit0=niit0+niit0/2

  lsecant=.false.
  linint =.false.
  lridd  =.true.
  if (lsecant) then
  !Secant root finding                                                                                                                  
    do iit=1,niit0
       mu=mu1-real(target_zero1,8)*(mu2-mu1)/real(target_zero2-target_zero1,8)
       call ndeviationQ(mu, iT, ek, sct, mesh, thdr, ltetra, target_zero)

       if (abs(target_zero).lt.dev) exit    
       if (target_zero.gt.0.q0) then
          mu1=mu
          target_zero1=target_zero 
          call ndeviationQ(mu2, iT, ek, sct, mesh, thdr, ltetra, target_zero2)
       else
          mu2=mu
          target_zero2=target_zero 
          call ndeviationQ(mu1, iT, ek, sct, mesh, thdr, ltetra, target_zero1)
       endif
    enddo
    niitact=iit


  elseif (linint) then
    ! evaluate the target function on an interval and find the root by linear interpolation 
    ! fix the number of points to sample the (mu1,mu2) 
    nmu=30
    allocate(Y(nmu), X(nmu))
    Y(:)=0.0d0 ; X(:)=0.0d0 
    ! construct linear grid 
    X(1)=mu1; X(nmu)=mu2; dmu=(mu2-mu1)/(nmu-1)
    !write(*,*) 'mu1',mu1,'mu2',mu2,'dmu',dmu
    Y(1)=target_zero1
    Y(nmu)=target_zero2
  
    do i=2,nmu-1
      X(i)=X(i-1)+dmu
    enddo 
    ! evaluate target function in the interval
    do i=2,nmu-1
       call ndeviationQ(X(i), iT, ek, sct, mesh, thdr, ltetra, Y(i))
    enddo 
    do i=1,nmu
      Y(i)=Y(i)+X(i) !this is the correct target function for this method
    enddo 
    
    !!!!!!!!!!!!!!!!!TEST 
    !open(666,file='targeT.dat',status='unknown')
    !write(666,'(A,1I10)')'T ',iT
    !write(666,'(2E15.7)') (X(i),Y(i), i=1,nmu)
    !write(666,'(A)')'   '
    !!!!!!!!!!!!!!!!!TEST END 

    ! find root by linear interpolation
    do i = 1, nmu-1
       do j = 1, nmu-1
          A(:,:) = 0.0d0
          a11 = X(i+1)-X(i)
          a22 = X(j+1)-X(j)
          a31 = Y(i+1)-Y(i)
          a42 = X(j+1)-X(j)
  
          !write(*,*) 'assignment to vectors'
          A(1,1)=a11; A(2,2)=a22; A(3,1)=a31; A(4,2)=a42
          A(1,3)=-1.0d0; A(2,3)=-1.0d0
          A(3,4)=-1.0d0; A(4,4)=-1.0d0
          B(1) = -X(i); B(2) = -X(j)
          B(3) = -Y(i); B(4) = -X(j)
  
          !write(*,*) 'LU factorisation begins'
          call dgetrf(4, 4, A, 4, ipiv, ierr )
          if (ierr /= 0) write(*,*) 'lu fact failed', ierr, i, j, a31
  
          !write(*,*) 'solution lin syst begins'
          call dgetrs( 'N', 4, 1, A, 4, ipiv, B, 4, ierr)
          if (ierr /= 0) write(*,*) 'solution of the system has failed', ierr, i, j
  
          ! check if there is any intersection
          if (B(1) < 1.0d0 .and. B(2) < 1.0d0) then
             if (B(1) >= 0.0d0 .and. B(2) >= 0.0d0) then
                !write(*,*) b(3), b(4)   
                ! save the values of the intersection
                mu = B(3)
                call ndeviationQ(mu, iT, ek, sct, mesh, thdr, ltetra, target_zero)
             endif
          endif
       enddo ! over freq. counter j
    enddo ! over freq. counter i 
    deallocate(Y,X)

  elseif(lridd) then   !Ridders' method for root finding
    !write(*,*) 'Ridders search for chemical potential'

! initialise the varibles
    P(1)=mu1 ; P(2)=mu2
    F(1)= target_zero1; F(2)= target_zero2
    !ptol   =  1.0d-6
    ptol   =  dev
    psave  = -1.1d30
    maxiter= 60    

     do j = 1, maxiter
        P(3) = 0.5d0*(P(1)+P(2))
        call ndeviationQ(P(3), iT, ek, sct, mesh, thdr, ltetra, F(3))
        s = sqrt((F(3)**2)-(F(1)*F(2)))
        if (s==0.0d0) then
           write(*,*) 'Error in Ridders search for chemical potential'
           write(*,*) 'ITER', j, 'x1', P(1),'  x2',P(2),'  x3', P(3)
           write(*,*) 'ITER', j, 'F1', F(1),'  F2',F(2),'  F3', F(3)
           goto 400
        endif
        P(4) = P(3)+(P(3)-P(1))*(sign(1.0q0,F(1)-F(2))*F(3)/s)
        if(abs(P(4)-psave)<=ptol) goto 400
        psave= P(4)
        call ndeviationQ(P(4), iT, ek, sct, mesh, thdr, ltetra, F(4))
        if (f(4) ==0.0d0) goto 400
        if (sign(F(3), F(4)) /= F(3)) then
        !change of sign btw x3 and x4 then reduce search interval
           P(1)  = P(3)
           F(1)  = F(3)
           P(2)  = P(4)
           F(2)  = F(4)
        else if (sign(F(1), F(4)) /= F(1)) then
        !change of sign btw x1 and x4 then reduce search interval
           P(2)  = P(4)
           F(2)  = F(4)
        else if (sign(F(2), F(4)) /= F(2)) then
        !change of sign btw x2 and x4 then reduce search interval
           P(1)  = P(4)
           F(1)  = F(4)
        endif
        !condition for termination
        if (abs(P(2)-P(1)) <= ptol) goto 400
     enddo ! over number of iterations 

 400 if (j == maxiter) write(*,*) 'Ridders seach might not have converged'

     ! save the values of the intersection
     mu = P(4)
     niitact = j
     niit0   = maxiter
     target_zero = F(4)
  
  endif ! root finding algorithm

  if (lsecant .or. lridd) then
    if ((niitact.ge.niit0).and.(myid.eq.master)) then
       write(*,'(A,1E20.12)') "WARNING: diminished root precision. ndevQ_actual =",real(target_zero,8)
       write(*,'(A,1F10.3,A,1I5,A,1E20.12)') "at T=",T, " with  niitQ=",niitQ, " ndevQ =", real(dev,8) !ndevQ       
       write(*,*) "increase niitQ, or allow for bigger ndevQ (see params.F90)"
    endif
  endif

  return
end subroutine find_muQ


subroutine find_muQ_lowT(mu,iT, ek, sct, mesh, thdr, ltetra, dos)
  use params
  use types
  use mpi_org !, only: myid,master
  implicit none

  type(edisp) :: ek
  type(scatrate) :: sct
  type(kpointmesh) :: mesh
  type(tetramesh)  :: thdr
  logical :: ltetra
  type(dosgrid) :: dos
  real(16) test,test1,test2
  real(8) mu,mu1,mu2,muold,sgn,tmp
  integer istep,iT,nstep

  nstep=100

  do istep=1,nstep

     mu1=mu-dos%gap/2.05d0 ! should be related to size of gap... ~Delta/2 or so                                                             
     mu2=mu+dos%gap/2.05d0

     call ndeviationQ(mu, 1, ek, sct, mesh, thdr, ltetra, test)
!     if (myid.eq.master)        write(*,*)test                                                                                         
     !call ndeviationQ(mu1, 1, ek, sct, mesh, test1)
     call ndeviationQ(mu1, 1, ek, sct, mesh, thdr, ltetra, test1)
!     if (myid.eq.master)        write(*,*)test1                                                                                        
     call ndeviationQ(mu2, 1, ek, sct, mesh, thdr, ltetra, test2)
!     if (myid.eq.master)        write(*,*)test2


!     if (myid.eq.master) then 
!        write(*,*)
!        write(*,*)mu
!        write(*,*)test-(test1+test2)/2.q0
!        write(*,*)mu+(mu1-mu2)*(test-(test1+test2)/2.q0)
!     endif                                                                                                                             

     muold=mu
     tmp=real((mu1-mu2)*(test-(test1+test2)/2.q0)*10.d0,8)

     sgn=1.d0
     if (tmp.lt.0.d0) sgn=-1.d0

     mu=mu+sgn*min(abs(tmp),5.d-3)
!     mu=mu+tmp                                                                                                                         
     !call ndeviationQ(mu,NE,1,test2)
     !if (myid.eq.master)        write(*,*)
     !if (myid.eq.master)        write(*,*)test2
     !if (myid.eq.master)        write(*,*)
!     if (myid.eq.master)        write(*,*)' NEW MU ',istep, mu
!     if (myid.eq.master)        write(22,*)istep,mu
!     if (myid.eq.master)        write(*,*)
!     if (myid.eq.master)        write(*,*)                                                                                             

     if (abs(mu-muold).lt.1.d-5) exit

  enddo
  !if ((istep.ge.100).and.(myid.eq.master)) write(*,*) 'increase istepmax in find_muQ_lowT'
  if (istep.ge.100) write(*,*) 'increase istepmax in find_muQ_lowT'
  return
end subroutine find_muQ_lowT


function FERMIQ(eps,beta)
  implicit none
  real(16) :: FERMIQ
  real(8) :: eps,beta

  FERMIQ=1.q0 / (1.q0 + EXP(real(beta*eps,16)))
!  FERMIQ=1.q0 / (1.q0 + QEXP(real(beta*eps,16)))

return
end function FERMIQ



subroutine ndeviationQ(mu, iT, ek, sct, mesh, thdr, ltetra, target_zero)
  use params
  use types
  use mpi_org
  implicit none

  !passed variables
  real(8), intent(in)  :: mu
  integer iT
  type(edisp) :: ek
  type(scatrate) :: sct
  type(kpointmesh) :: mesh
  type(tetramesh)  :: thdr
  logical :: ltetra
  real(16) :: occ_tot
  real(16), intent(out) :: target_zero

  !real(16) ninteger,nbig,nsmall,target_zero,NQ !QP   
  !real(16) ninteger_tot,nbig_tot,nsmall_tot !QP 

  if (ltetra) then
     call varoccQ_tet(mu, iT, ek, sct, thdr, target_zero)
  else
     call varoccQ(mu, iT, ek, sct, mesh, occ_tot)
     target_zero=real(ek%nelect,16)-occ_tot
  endif

  !NQ=real(NE,16)
  ! particles per processor                                                                                                             
  !call n_per_procQ(iT,mu,ninteger,nbig,nsmall)
  ! summing particles from all processors                                                                                               
  !call reduce_nQ(ninteger,nbig,nsmall,ninteger_tot,nbig_tot,nsmall_tot)
  !if (myid.eq.master) write(*,*)'tot ',myid,ninteger_tot,nbig_tot,nsmall_tot                                                           


  !if (myid.eq.master) then

! combine large and small numbers                                                                                                       
     !target_zero= ( (NQ-ninteger_tot) - nbig_tot ) - nsmall_tot
!target_zero=real(((NE-ninteger_tot)-nbig_tot)-nsmall_tot,8) ! here transform to DP                                                    
!write(*,*)'master ', target_zero                                                                                                       

  !endif !master                                                                                                                         
!broadcast target_zero. Not very elegant to do this...
!broadcasting QUAD precision data ...
  !call MPI_BARRIER( MPI_COMM_WORLD, mpierr )
  !call MPI_BCAST_QUAD(target_zero)
  !call MPI_BCAST(target_zero,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpierr)                                                      

return
end subroutine ndeviationQ


subroutine varoccQ(mu, iT, ek, sct, mesh, occ_tot)
  use types
  use params
  use mpi_org

  implicit none

  integer, intent(in) :: iT
  real(8), intent(in) :: mu
  type(edisp) ::ek
  type(scatrate) :: sct
  type(kpointmesh) :: mesh
  real(16) :: occ_loc, occ_tot
!local variables
  real(16),parameter :: thr = 1.0q-30
  complex(16) ::z
  real(16) :: nsmall, nbig, ninteger, tmp, tmp2
  real(8) :: eps
  integer :: iband, ik
  integer :: ktot
  real(16) :: cutQ
  !more sophistication                                                                                                                  
  integer(8) ::IEXP
  !parameters                                                                                                                           
  real(16),parameter :: QCUT=1.Q14 ! not the same as cutQ!!!!!                                                                                        
    ! relevant digits  14?                                                                                         
!external variables
  real(16), external :: FERMIQ
  complex(16), external :: wpsipghp

  occ_tot=0.0q0
  occ_loc=0.0q0

  ninteger=0.q0
  nbig=0.q0
  nsmall=0.q0
  ktot = mesh%kx*mesh%ky*mesh%kz
  cutQ=1.q0/thr

  do ik = iqstr, iqend
     do iband=1,ek%nband_max 
        if (ek%band(ik,iband) .gt. 90.0d0) cycle
        if (iband<ek%nbopt_min) cycle
        if (iband>ek%nbopt_max) cycle
        eps=sct%z*ek%band(ik,iband)-mu

        if (eps.lt.0.d0) then ! occupied state                                                                                                        
           ninteger=ninteger+1.q0
           
           if (sct%gam(iT,iband).eq.0.d0) then
              tmp=1.q0-FERMIQ(eps,beta)
           else
              z=0.5q0 + real(sct%z*sct%gam(iT,iband)*beta2p,16) + ciQ*real(eps*beta2p,16) ! eps --> -eps
              tmp=0.5q0+aimag(wpsipghp(z,0))/piQ ! >0 
           endif
     
           !write(*,*) 'occupied state, occ var', tmp
           if (tmp.gt.thr) then
              IEXP=int(log10(abs(tmp)),8)
              tmp2=( real(int((tmp/(10.q0**iEXP))*QCUT,8),16)*10.q0**iEXP ) / QCUT
              tmp=tmp-tmp2
     
              nbig=nbig-tmp2
              nsmall=nsmall-tmp
           else
              nsmall=nsmall-tmp
           endif
     
        else ! unoccupied state
     
           if (sct%gam(iT,iband).eq.0.d0) then
              tmp=FERMIQ(eps,beta)
           else
              z=0.5q0 + real(sct%z*sct%gam(iT,iband)*beta2p,16) - ciQ*real(eps*beta2p,16) ! eps
              tmp=(0.5q0+aimag(wpsipghp(z,0))/piQ) ! >0 !                               
           endif
     
           !write(*,*) 'unoccupied state, occ var', tmp 
           if (tmp.gt.thr) then
              if (tmp.ne.0.q0) then
                 IEXP=int(log10(abs(tmp)),8)
              else
                 IEXP=0
              endif
     
              tmp2=( real(int((tmp/(10.q0**iEXP))*QCUT,8),16)*10.q0**iEXP ) / QCUT
              tmp=tmp-tmp2
              nbig=nbig+tmp2
              nsmall=nsmall+tmp
           else
              nsmall=nsmall+tmp
           endif
        endif
     enddo ! iband       
  enddo                                                                                                                   
  ninteger=2.q0*ninteger/real(ktot,16) ! 2 for spin                                                                                       
  nbig=2.q0*nbig/real(ktot,16)
  nsmall=2.q0*nsmall/real(ktot,16)
  !ek%occ_tot=real(ninteger+nbig+nsmall,8)

  if (nproc == 1) then
     !target_zero=((real(ek%nelect,16)-ninteger)-nbig)-nsmall
     occ_tot=ninteger+nbig+nsmall
  else
     !target_loc=((real(ek%nelect,16)-ninteger)-nbig)-nsmall
     occ_loc=ninteger+nbig+nsmall
  endif

  if (nproc > 1) then
     call mpi_reduce_quad(occ_loc, occ_tot)
  endif

  return

end subroutine varoccQ

subroutine varoccQ_tet(mu, iT, ek, sct, thdr, target_zero)
  use types
  use params
  use mpi_org
  use estruct

  implicit none

  integer, intent(in) :: iT
  real(8), intent(in) :: mu
  type(edisp) ::ek
  type(scatrate) :: sct
  type(tetramesh) :: thdr
  real(16) :: target_loc
  real(16), intent(out) :: target_zero
!local variables
  real(16),parameter :: thr = 1.0q-30
  complex(16) ::z
  real(16) :: nsmall, nbig, ninteger, tmp, tmp2
  real(16) :: target_tet(4), target_intp
  real(8) :: eps
  integer :: iband, ik
  integer :: itet
  real(16) :: cutQ
  !more sophistication
  integer(8) ::IEXP
  !parameters
  real(16),parameter :: QCUT=1.Q14 ! not the same as cutQ!!!!!
!external variables
  real(16), external :: FERMIQ
  complex(16), external :: wpsipghp

  target_zero=0.0q0
  target_loc=0.0q0
  cutQ=1.q0/thr

  do itet = iqstr, iqend
     target_intp=0.0q0
     do ik = 1, 4
        ninteger=0.q0
        nbig=0.q0
        nsmall=0.q0
        do iband=1,ek%nband_max 
           if (ek%band(thdr%idtet(ik,itet),iband) .gt. 90.0d0) cycle
           if (iband<ek%nbopt_min) cycle
           if (iband>ek%nbopt_max) cycle
           eps=sct%z*ek%band(thdr%idtet(ik,itet),iband)-mu
        
           if (eps.lt.0.d0) then ! occupied state
              ninteger=ninteger+1.q0
              
              if (sct%gam(iT,iband).eq.0.d0) then
                 tmp=1.q0-FERMIQ(eps,beta)
              else
                 z=0.5q0 + real(sct%z*sct%gam(iT,iband)*beta2p,16) + ciQ*real(eps*beta2p,16) ! eps --> -eps
                 tmp=0.5q0+aimag(wpsipghp(z,0))/piQ ! >0 
              endif
        
              !write(*,*) 'occupied state, occ var', tmp
              if (tmp.gt.thr) then
                 IEXP=int(log10(abs(tmp)),8)
                 tmp2=( real(int((tmp/(10.q0**iEXP))*QCUT,8),16)*10.q0**iEXP ) / QCUT
                 tmp=tmp-tmp2
        
                 nbig=nbig-tmp2
                 nsmall=nsmall-tmp
              else
                 nsmall=nsmall-tmp
              endif
        
           else ! unoccupied state
        
              if (sct%gam(iT,iband).eq.0.d0) then
                 tmp=FERMIQ(eps,beta)
              else
                 z=0.5q0 + real(sct%z*sct%gam(iT,iband)*beta2p,16) - ciQ*real(eps*beta2p,16) ! eps
                 tmp=(0.5q0+aimag(wpsipghp(z,0))/piQ) ! >0 
              endif
        
              !write(*,*) 'unoccupied state, occ var', tmp
              if (tmp.gt.thr) then
                 if (tmp.ne.0.q0) then
                    IEXP=int(log10(abs(tmp)),8)
                 else
                    IEXP=0
                 endif
        
                 tmp2=( real(int((tmp/(10.q0**iEXP))*QCUT,8),16)*10.q0**iEXP ) / QCUT
                 tmp=tmp-tmp2
                 nbig=nbig+tmp2
                 nsmall=nsmall+tmp
              else
                 nsmall=nsmall+tmp
              endif
           endif
        enddo ! iband  
        ninteger=2.q0*ninteger ! 2 for spin
        nbig=2.q0*nbig
        nsmall=2.q0*nsmall
        target_tet(ik)=((real(ek%nelect,16)-ninteger)-nbig)-nsmall
     enddo ! corners of the tetrahedron     
     call interptra_muQ (thdr%vltet(itet), target_tet, target_intp)
     if (nproc == 1) then
        ! accumulate globally (nproc = 1)
        target_zero=target_zero+target_intp
     else
        ! or locally (nproc > 1)
        target_loc=target_loc+target_intp
     endif
  enddo ! tetrahedra

  ! MPI DISTRIBUTION
  if (nproc > 1) then
     call mpi_reduce_quad(target_loc, target_zero)
  endif

  return

end subroutine varoccQ_tet




!subroutine nikQ(ik,iT,mu,ninteger,nbig,nsmall)
!  use params
!  use estruct
!  implicit none
!  complex(16) wpsipghp,z
!  external wpsipghp
!  integer iband,ik,iT
!  real(8) eps,mu
!  real(16) nsmall,nbig,ninteger,tmp,tmp2!,threshold                                                                                     
!
!  real(16) cutQ
!!  parameter(cutQ=1Q12)                                                                                                                 

  !more sophistication                                                                                                                  
!  integer(8) ::IEXP
!
!  !parameters                                                                                                                           
!  real(16) :: QCUT   ! not the same as cutQ!!!!!                                                                                        
!  parameter(QCUT=1.Q14)  ! relevant digits  14?                                                                                         
!  !test: FERMI function QUAD, from DP arguments                                                                                         
!  real(16) FERMIQ
!  external FERMIQ
!
!  cutQ=1.q0/thresholdQ
!
!  do iband=1,nband
!     eps=ek(ik,iband)-mu
!
!     if (eps.lt.0.d0) then ! occ                                                                                                        
!        ninteger=ninteger+1.q0
!
!        if (gam(iT,iband).eq.0.d0) then
!           tmp=1.q0-FERMIQ(eps,beta)
!        else
!           z=0.5q0 + real(gam(iT,iband)*beta2p,16) + ciQ*real(eps*beta2p,16) ! eps --> -eps       
!           tmp=0.5q0+aimag(wpsipghp(z,0))/piQ ! >0 !
!        endif

!        write(*,*)tmp                                                                                                                  
        ! XXX TESTING GAMMA = 0                                                                                                         
!        tmp=1.q0-FERMIQ(eps,beta)                                                                                                      

!        write(*,*)tmp                                                                                                                  
!        write(*,*)                                                                                                                     

        !      write(*,*)'test'                                                                                                         
        !      write(*,*)ci*eps*beta2p                                                                                                  
        !      write(*,*)ciQ*real(eps*beta2p,16)                                                                                        
        !      write(*,*)z                                                                                                              
        !      write(*,*)wpsipghp(z,0)                                                                                                  
        !      write(*,*)tmp,thresholdQ                                                                                                 
        !      STOP                                                                                                                     

!        write(33,*)tmp                                                                                                                 

!        if (tmp.gt.thresholdQ) then

!           write(33,*)'I ',tmp                                                                                                         

!try something more sophisticated here...                                                                                               
!           if (qp0.ne.0.q0) then                                                                                                       
!           IEXP=int(log10(abs(tmp)),8)
!           write(*,*)iexp                                                                                                              

!           else                                                                                                                        
!              IEXP=0                                                                                                                   
!           endif                                                                                                                       

!           tmp2=( real(int((tmp/(10.q0**iEXP))*QCUT,8),16)*10.q0**iEXP ) / QCUT
!           tmp2=real(int(tmp*cutQ,8),16)/cutQ                                                                                         
!           tmp=tmp-tmp2

!           write(33,'(A,2E)')'i  ',tmp2,tmp                                                                                            
!           write(33,*)                                                                                                                 


!           nbig=nbig-tmp2 !                                                                                                             
           !         nbig=nbig-tmp !               

!          nsmall=nsmall-tmp

           ! do something more advanced here...                                                                                         
           !    tmp2=nint(tmp/threshold)*threshold                                                                                      
           !    tmp=tmp-tmp2                                                                                                            
           !    nbig=nbig-tmp2                                                                                                          
           !    nsmall=nsmall-tmp                                                                                                       
!        else
!           nsmall=nsmall-tmp
!        endif
!
!     else ! unocc                                                                                                                       
!
!        if (gam(iT,iband).eq.0.d0) then
!           tmp=FERMIQ(eps,beta)
!        else
!           z=0.5q0 + real(gam(iT,iband)*beta2p,16) - ciQ*real(eps*beta2p,16) ! eps                                                      
!           tmp=(0.5q0+aimag(wpsipghp(z,0))/piQ) ! >0 !                                                                                  
!        endif

        !      write(*,*)tmp                                                                                                            
! XXX TESTING GAMMA = 0                                                                                                                 
!        tmp=FERMIQ(eps,beta)                                                                                                           

!        if (tmp.gt.thresholdQ) then

!!           write(33,*)'I ',tmp        
!!try something more sophisticated here...      

!           if (tmp.ne.0.q0) then
!              IEXP=int(log10(abs(tmp)),8)
!           else
!              IEXP=0
!           endif

!!           write(*,*)iexp             
!!           else                       
!!              IEXP=0                  
!!           endif                      

!           tmp2=( real(int((tmp/(10.q0**iEXP))*QCUT,8),16)*10.q0**iEXP ) / QCUT
!!           tmp2=real(int(tmp*cutQ,8),16)/cutQ
                                                                                                                                        
!           tmp=tmp-tmp2

!!           write(33,'(A,2E)')'i  ',tmp2,tmp   
!           !           write(33,*)                        
!                                                                                                                                        
!           nbig=nbig+tmp2 !                   
                                                                                                                                        
           !         nbig=nbig-tmp !          

!           nsmall=nsmall+tmp
!           else
!           nsmall=nsmall+tmp
!        endif
!
!     endif !eps occ or unocc                                                                                                            
!
!  enddo !iband                                                                                                                          
!
!  return
!end subroutine nikQ



!subroutine n_per_procQ(iT,mu,ninteger,nbig,nsmall)
 ! use params
 ! use mpi_org
 ! use estruct
 ! implicit none
 ! integer ik,it
 ! real(8) mu
 ! real(16) ninteger, nbig, nsmall
 !
 ! ninteger=0.q0
 ! nbig=0.q0
 ! nsmall=0.q0
 ! nkthis=(ikend-ikstart+1)
 ! do ik=ikstart,ikend
 !    call nikQ(ik,iT,mu,ninteger,nbig,nsmall) ! adds to the n's        
 ! enddo
 ! ninteger=2.q0*ninteger/real(nk,16) ! 2 for spin                      
 ! nbig=2.q0*nbig/real(nk,16)
 ! nsmall=2.q0*nsmall/real(nk,16)
 !
 ! return
 !end subroutine n_per_procQ
 !

!Subroutine reduce_nQ(ninteger,nbig,nsmall,ninteger_tot,nbig_tot,nsmall_tot)
 ! use mpi_org
 ! implicit none
 ! real(16) ninteger_tot,nsmall_tot,nbig_tot
 ! real(16) ninteger,nsmall,nbig,tmp
 !
 ! ninteger_tot=0.q0
 ! nbig_tot=0.q0
 ! nsmall_tot=0.q0
 !
 !
!! openMPI nor MPICH fully support quad precision... so I have to devise a bloody workaround... 
 ! if (nproc.gt.1) then
 ! call MPI_REDUCE_QUAD(ninteger,ninteger_tot)
 ! call MPI_REDUCE_QUAD(nbig,nbig_tot)
 ! call MPI_REDUCE_QUAD(nsmall,nsmall_tot)
 ! else
 ! ninteger_tot=ninteger
 ! nbig_tot=nbig
 ! nsmall_tot=nsmall
 ! endif
 !
! if MPI_REAL16 is supported on your cluster.... it bloody ain't on HCLM  
  !call MPI_REDUCE(ninteger,ninteger_tot,1,MPI_REAL16,MPI_SUM,master,MPI_COMM_WORLD,mpierr) 
  !call MPI_BARRIER( MPI_COMM_WORLD, mpierr )  
  !call MPI_REDUCE(nbig,nbig_tot,1,MPI_REAL16,MPI_SUM,master,MPI_COMM_WORLD,mpierr)
  !call MPI_BARRIER( MPI_COMM_WORLD, mpierr )  
  !call MPI_REDUCE(nsmall,nsmall_tot,1,MPI_REAL16,MPI_SUM,master,MPI_COMM_WORLD,mpierr) 
  !call MPI_BARRIER( MPI_COMM_WORLD, mpierr )                                                                                           






  !write(*,*)'red ',myid,ninteger_tot,nbig_tot,nsmall_tot                                                                               

!end Subroutine reduce_nQ

