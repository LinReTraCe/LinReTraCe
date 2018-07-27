module Mroot
  use Mparams
  use Mtypes

contains

subroutine find_mu(mu,iT,nT,dev,target_zero,niitact, ek, sct, mesh )
  implicit none
  ! passed variables
  real(8) :: mu, dev, target_zero
  integer  :: iT,nT,niitact
  type(edisp) :: ek
  type(scatrate) :: sct
  type(kpointmesh) :: mesh
  ! local variables
  real(8) target_zero1, target_zero2, mu1, mu2
  integer iit,niit0
  logical lsecant  ! selects the secant root finding algorithm
  ! linear interpolation method
  real(8), allocatable :: Y(:), X(:) !arrays containig the function to minimise and the chemical potential
  integer :: nmu  ! number of points that sample the mu interval (mu1,mu2)
  real(8) :: dmu ! increment
  real(8) :: a11, a22, a31, a42
  real(8) :: A(4,4), B(4)
  integer :: i, j
  integer :: ipiv(4)
  integer :: ierr
  logical linint  ! selects the linear interpolation method
  ! Ridders' method
  real(8) :: F(4), P(4)
  real(8) :: s
  real(8) :: psave, ptol
  integer  :: maxiter ! maximum number of iterations
  logical  :: lridd   ! selects Ridders' method

! deviation from set particle number with initial mu
  call ndeviation(mu, iT, ek, sct, mesh, target_zero1)

  target_zero2=target_zero1
!coarse initialization of secant bracket mu1, mu2... Secant doesnt need mu to lie within bracket, but here it does
  mu1=mu
  mu2=mu

  do while (target_zero2.gt.0.d0)
     mu2=mu2+0.005d0
     call ndeviation(mu2, iT, ek, sct, mesh, target_zero2)

  enddo
  do while (target_zero1.le.0.d0)
     mu1=mu1-0.005d0
     call ndeviation(mu1, iT, ek, sct, mesh, target_zero1)

  enddo

  niit0=niit
!  if (iT.eq.nT) niit0=niit*3
  lsecant=.false.
  linint =.false.
  lridd  =.true.

  if (lsecant) then
  !Secant root finding
    do iit=1,niit0
       mu=mu1-target_zero1*(mu2-mu1)/(target_zero2-target_zero1)
       call ndeviation(mu, iT, ek, sct, mesh, target_zero)

       if (abs(target_zero).lt.dev)  exit
       if (target_zero.gt.0.d0) then
          mu1=mu
          target_zero1=target_zero
          call ndeviation(mu2, iT, ek, sct, mesh, target_zero2)
       else
          mu2=mu
          target_zero2=target_zero
          call ndeviation(mu1, iT, ek, sct, mesh, target_zero1)
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
       call ndeviation(X(i), iT, ek, sct, mesh, Y(i))
    enddo
    do i=1,nmu
      Y(i)=Y(i)+X(i) !this is the correct target function for this method
    enddo

!!!!!!!!!!!!!!!!!test
    !open(666,file='targeT.dat',status='unknown')
    !write(666,'(A,1I10)')'T ',iT
    !write(666,'(2E15.7)') (X(i),Y(i), i=1,nmu)
    !write(666,'(A)')'   '
!!!!!!!!!!!!!!!!!test end

    ! find root by linear interpolation
    do i = 1, nmu-1
       do j = 1, nmu-1
          A(:,:) = 0.0d0
          a11 = X(i+1)-X(i)
          a22 = X(j+1)-X(j)
          a31 = Y(i+1)-Y(i)
          a42 = X(j+1)-X(j)

          A(1,1)=a11; A(2,2)=a22; A(3,1)=a31; A(4,2)=a42
          A(1,3)=-1.0d0; A(2,3)=-1.0d0
          A(3,4)=-1.0d0; A(4,4)=-1.0d0
          B(1) = -X(i); B(2) = -X(j)
          B(3) = -Y(i); B(4) = -X(j)

          !write(*,*) 'LU factorisation begins'
          call dgetrf(4, 4, A, 4, ipiv, ierr )
          if (ierr /= 0) write(*,*) 'LU factorisation failed', ierr, i, j, a31

          !write(*,*) 'solution lin syst begins'
          call dgetrs( 'N', 4, 1, A, 4, ipiv, B, 4, ierr)
          if (ierr /= 0) write(*,*) 'solution of the system of linear equations has failed', ierr, i, j

          ! check if there is any intersection
          if (B(1) < 1.0d0 .and. B(2) < 1.0d0) then
             if (B(1) >= 0.0d0 .and. B(2) >= 0.0d0) then
                !write(*,*) b(3), b(4)
                ! save the values of the intersection
                mu = B(3)
                call ndeviation(mu, iT, ek, sct, mesh, target_zero)
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
    !ptol   =  1.0d-18
    ptol   =  dev
    psave  = -1.1d30
    maxiter= 60

     do j = 1, maxiter
        P(3) = 0.5d0*(P(1)+P(2))
        call ndeviation(P(3), iT, ek, sct, mesh, F(3))
        s = dsqrt((F(3)**2)-(F(1)*F(2)))
        if (s==0.0d0) then
           write(*,*) 'Error in Ridders search for chemical potential'
           write(*,*) 'ITER', j, 'x1', P(1),'  x2',P(2),'  x3', P(3)
           write(*,*) 'ITER', j, 'F1', F(1),'  F2',F(2),'  F3', F(3)
           goto 400
        endif
        P(4) = P(3)+(P(3)-P(1))*(SIGN(1.0d0,F(1)-F(2))*F(3)/s)
        if(abs(P(4)-psave)<=ptol) goto 400
        psave= P(4)
        call ndeviation(P(4), iT, ek, sct, mesh, F(4))
        if (F(4) ==0.0d0) goto 400
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
        !write(*,*)'iter, mu, target_zero',j, P(4), F(4)
     enddo ! over number of iterations

 400 if (j == maxiter) write(*,*) 'Ridders seach might not have converged'

     ! save the values of the intersection
     mu = P(4)
     niitact = j
     niit0   = maxiter
     target_zero = F(4)

  endif ! root finding algorithm

  if (lsecant .or. lridd) then
    !if ((niitact.ge.niit0).and.(myid.eq.master)) then
    if (niitact .ge. niit0) then
       write(*,'(A,1E20.12)') "WARNING: diminished root precision. ndev_actual =",target_zero
       write(*,'(A,1F10.3,A,1I5,A,1E20.12)') "at T=",sct%TT(iT), " with  niit=",niit0, " ndev =", dev !ndevQ
       write(*,*) "increase niit, or allow for bigger ndev (see params.F90)"
       !write(*,*) "myid=",myid
    endif
  endif

  return
end subroutine find_mu


!*******************************************************************************
! The RIDDRS_ROOT routine finds the zeroes of the function Y+X
! using Ridders' method. Given an interval [x1, x2] the function is evaluated at
! the midpoint x3 then the unique exponential function that turns the residual
! function into a straight line is factorised, i.e. the equation:
!       F(x1) - 2F(x3)*exp{Q} + F(x2)*exp{2Q} = 0 is solved.
! The false position method is then applied to the factors in the equation
! above. This leads to the guess for the zero:
!   x4 = x3 + (x3-x1)*[f(x3)*sign(f(x1)-f(x2))]/sqrt(f(x3)^2 - f(x1)*f(x2))
! Properties:
! x4 is always in the interval [x1,x2]
! number of significant digits doubles with each two function evaluations
! the method should be computationally robust.
!*******************************************************************************


function FERMI(eps,beta)
  implicit none
  real(8) :: FERMI
  real(8) :: eps,beta

  FERMI=1.d0 / (1.d0 + EXP(beta*eps))

return
end function FERMI


subroutine ndeviation(mu, iT, ek, sct, mesh, target_zero)
  implicit none

  !passed variables
  real(8), intent(in)  :: mu
  integer iT
  type(edisp) :: ek
  type(scatrate) :: sct
  type(kpointmesh) :: mesh
  real(8), intent(out) :: target_zero
  call varocc(mu, iT, ek, sct, mesh )

  !if (myid.eq.master) then

     target_zero=ek%nelect-ek%occ_tot
     !write(*,*)'mu, target_zero, ek%occ_tot ', mu, target_zero, ek%occ_tot

  !endif !master
  ! broadcast target_zero. Not very elegant to do this...
  !call MPI_BCAST(target_zero,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpierr)
!  call MPI_BARRIER( MPI_COMM_WORLD, mpierr )
  return
end subroutine ndeviation

subroutine varocc(mu, iT, ek, sct, mesh)
  implicit none

  integer, intent(in) :: iT
  real(8), intent(in) :: mu
  type(edisp) :: ek
  type(scatrate) :: sct
  type(kpointmesh) :: mesh
!local variables
  real(8),parameter :: thr = 1.0d-30
  complex(8) :: z
  real(8) :: nsmall, nbig, ninteger, eps, tmp
  integer :: iband, ik, ikx, iky, ikz
  integer :: ktot !total number of k-points
  ! RECALL that mesh%ktot includes extra k-points when the
  ! tetrahedron method is selected
!external variables
  complex(8), external :: wpsipg

  ek%occ_tot=0.0d0
  if (.not.allocated(ek%occ)) then
    allocate(ek%occ(mesh%ktot,ek%nband_max))
  endif
  !write(*,*)'varocc: after allocation',mesh%ktot,mesh%kx,mesh%ky,mesh%kz,ek%nband_max
  !write(*,*)'allocated mesh%k_id',allocated(mesh%k_id),size(mesh%k_id)
  !write(*,*)'allocated ek%band',allocated(ek%band),size(ek%band)

  ninteger=0.d0
  nbig=0.d0
  nsmall=0.d0
  ktot=mesh%kx*mesh%ky*mesh%kz

  do ik = 1, ktot
     do iband=1,ek%nband_max
        if (ek%band(ik,iband) .gt. 90.0d0) cycle
        eps=sct%z*ek%band(ik,iband)-mu

        if (eps .lt. 0.d0) then ! occupied state
           ninteger=ninteger+1.d0

           if (sct%gam(iT,iband).eq.0.d0) then
              tmp=1.d0-FERMI(eps,beta)
           else
              z=0.5d0 + (sct%z*sct%gam(iT,iband) + ci*eps ) * beta2p ! eps --> -eps
              tmp=0.5d0+aimag(wpsipg(z,0))/pi ! >0 !
           endif

           if (tmp.gt.thr) then
              nbig=nbig-tmp
           else
              nsmall=nsmall-tmp
           endif

        else ! unoccupied state

           if (sct%gam(iT,iband).eq.0.d0) then
              tmp=FERMI(eps,beta)
           else
              z=0.5d0 + (sct%z*sct%gam(iT,iband) - ci*eps ) * beta2p ! eps
              tmp=(0.5d0+aimag(wpsipg(z,0))/pi) ! >0 !
           endif

           if (tmp.gt.thr) then
              nbig=nbig+tmp
           else
              nsmall=nsmall+tmp
           endif
        endif
        !TODO: remove from datatype?
        !ek%occ(ik,iband) = 2.0d0*(ninteger_tot+nbig_tot+nsmall_tot)/real(mesh%ktot,8)
        !ek%occ(ik,iband) = 2.0d0*(ninteger+nbig+nsmall)/real(mesh%ktot,8)
        !ek%occ_tot = ek%occ_tot + ek%occ(ik,iband)
     enddo ! iband
  enddo
  ninteger=2.d0*ninteger/real(ktot,8) ! 2 for spin
  nbig=2.d0*nbig/real(ktot,8)
  nsmall=2.d0*nsmall/real(ktot,8)
  ek%occ_tot=ninteger+nbig+nsmall

  return

end subroutine varocc

end module Mroot
