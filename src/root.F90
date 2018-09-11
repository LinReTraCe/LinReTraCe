module Mroot
  use Mmpi_org
  use Mparams
  use Mtypes
  use Mestruct

  interface find_mu
    module procedure find_mu_D, find_mu_Q
  end interface find_mu

  interface ndeviation
    module procedure ndeviation_D, ndeviation_Q
  end interface ndeviation

  interface occ
    module procedure occ_D, occ_Q
  end interface occ

  interface occ_tet
    module procedure occ_tet_D, occ_tet_Q
  end interface occ_tet

  interface fermi
    module procedure fermi_D, fermi_Q
  end interface fermi

  contains

subroutine find_mu_D(mu,iT,dev,target_zero,niitact, ek, sct, mesh, thdr, ltetra, method)
  implicit none

  real(8), intent(inout)        :: mu ! chemical potential which is calculated
  integer, intent(in)           :: iT ! integer of temperature loop
  real(8), intent(in)           :: dev ! allowed deviation
  real(8), intent(out)          :: target_zero ! deviation from root after convergence
  integer, intent(out)          :: niitact ! number of iterations
  logical, intent(in)           :: ltetra  ! use the tetrahedron methods
  integer, intent(in), optional :: method  ! choose root finding method
  type(edisp)                   :: ek
  type(scatrate)                :: sct
  type(kpointmesh)              :: mesh
  type(tetramesh)               :: thdr

  ! local variables
  real(8) target_zero1, target_zero2, mu1, mu2
  integer iit,niit0
  ! Secand method
  logical lsecant
  ! linear interpolation method
  real(8), allocatable :: Y(:), X(:) !arrays containig the function to minimise and the chemical potential
  integer :: nmu  ! number of points that sample the mu interval (mu1,mu2)
  real(8) :: dmu ! increment
  real(8) :: a11, a22, a31, a42
  real(8) :: A(4,4), B(4)
  integer :: i, j
  integer :: ipiv(4)
  integer :: ierr
  ! Linear interpolation
  logical linint
  ! Ridders' method
  real(8) :: F(4), P(4)
  real(8) :: s
  integer  :: maxiter ! maximum number of iterations
  logical  :: lridd   ! selects Ridders' method
  ! Bisection method
  logical :: lbisec

  lsecant = .false.
  linint  = .false.
  lridd   = .false.
  lbisec  = .false.
  ! choose method according to input
  ! if method is not provided: default to Riddler
  if (present(method)) then
    select case (method)
      case (0)
        lsecant = .true.
      case (1)
        linint  = .true.
      case (2)
        lridd   = .true.
      case (3)
        lbisec  = .true.
    end select
  else
    lridd = .true.
  endif

! deviation from set particle number with initial mu
! output: target_zero1 = required - actual electrons
  call ndeviation(mu, iT, ek, sct, mesh, thdr, ltetra, target_zero1)

! coarse initialization of secant bracket mu1, mu2
  target_zero2=target_zero1
  mu1=mu
  mu2=mu

  do while (target_zero2.gt.0.d0) ! too few electrons -> increase mu
     mu2=mu2+0.005d0
     call ndeviation(mu2, iT, ek, sct, mesh, thdr, ltetra, target_zero2)
  enddo
  do while (target_zero1.le.0.d0) ! too many electrons -> decrease mu
     mu1=mu1-0.005d0
     call ndeviation(mu1, iT, ek, sct, mesh, thdr, ltetra, target_zero1)
  enddo

  ! maximal steps for double precision calculations
  niit0=niit

  if (lsecant) then
  !Secant root finding
    do iit=1,niit0
       mu=mu1-target_zero1*(mu2-mu1)/(target_zero2-target_zero1)
       call ndeviation(mu, iT, ek, sct, mesh, thdr, ltetra, target_zero)

       if (abs(target_zero).lt.dev)  exit
       if (target_zero.gt.0.d0) then
          mu1=mu
          target_zero1=target_zero
       else
          mu2=mu
          target_zero2=target_zero
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
       call ndeviation(X(i), iT, ek, sct, mesh, thdr, ltetra, Y(i))
    enddo
    do i=1,nmu
      Y(i)=Y(i)+X(i) !this is the correct target function for this method
    enddo

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
                call ndeviation(mu, iT, ek, sct, mesh, thdr, ltetra, target_zero)
             endif
          endif
       enddo ! over freq. counter j
    enddo ! over freq. counter i
    deallocate(Y,X)

  elseif(lridd) then   !Ridders' method for root finding
    P(1)=mu1 ; P(2)=mu2
    F(1)= target_zero1; F(2)= target_zero2
    maxiter=  niit

     do j = 1, maxiter
        P(3) = 0.5d0*(P(1)+P(2))
        call ndeviation(P(3), iT, ek, sct, mesh, thdr, ltetra, F(3))
        s = dsqrt((F(3)**2)-(F(1)*F(2)))
        if (s==0.0d0) then
           write(*,*) 'Error in Ridders search for chemical potential'
           write(*,*) 'ITER', j, 'x1', P(1),'  x2',P(2),'  x3', P(3)
           write(*,*) 'ITER', j, 'F1', F(1),'  F2',F(2),'  F3', F(3)
           goto 400
        endif
        P(4) = P(3)+(P(3)-P(1))*(SIGN(1.0d0,F(1)-F(2))*F(3)/s)
        call ndeviation(P(4), iT, ek, sct, mesh, thdr, ltetra, F(4))
        if (abs(F(4)) .lt. dev) goto 400
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
     enddo ! over number of iterations

 400 if (j == maxiter) write(*,*) 'Ridders seach might not have converged'

     ! save the values of the intersection
     mu = P(4)
     niitact = j
     niit0   = maxiter
     target_zero = F(4)

  elseif(lbisec) then
    ! Bisection root finding
    do iit=1,niit0
       mu = (mu1+mu2)/2.d0
       call ndeviation(mu, iT, ek, sct, mesh, thdr, ltetra, target_zero)
       if (myid.eq.master .and. iit .ge. 50) write(*,*) mu
       if (abs(target_zero).lt.dev) exit
       if (target_zero.gt.0.q0) then
          mu1=mu
          target_zero1=target_zero
       else
          mu2=mu
          target_zero2=target_zero
       endif
    enddo
    niitact = iit
  endif ! root finding algorithm

  if (myid.eq.master .and. (niitact .ge. niit0 .or. abs(target_zero) .ge. dev)) then
     write(*,'(A,1E20.12)') "WARNING: diminished root precision. ndev_actual =",target_zero
     write(*,'(A,1F10.3,A,1I5,A,1E20.12)') "at T=",sct%TT(iT), " with  niit=",niit0, " ndev =", dev
     write(*,*) "increase niit, or allow for bigger ndev (see params.F90)"
     !write(*,*) "myid=",myid
  endif
end subroutine find_mu_D

subroutine find_mu_Q(mu,iT,dev,target_zero,niitact, ek, sct, mesh, thdr, ltetra, method)
  implicit none
  ! passed variables
  real(8), intent(inout)        :: mu
  integer, intent(in)           :: iT
  real(16), intent(in)          :: dev
  real(16), intent(out)         :: target_zero
  integer, intent(out)          :: niitact
  logical, intent(in)           :: ltetra
  integer, intent(in), optional :: method
  type(edisp)                   :: ek
  type(scatrate)                :: sct
  type(kpointmesh)              :: mesh
  type(tetramesh)               :: thdr

  ! local variables
  real(16) mu_qp
  real(16) target_zero1, target_zero2
  real(16) mu1, mu2, dmu
  integer iit, niit0, itest
  logical lsecant  ! selects the secant root finding algorithm
  ! linear interpolation method
  real(16), allocatable :: Y(:) !array containig the function to minimise
  real(16), allocatable :: X(:) !array containig the chemical potential
  integer :: nmu  ! number of points that sample the mu interval (mu1,mu2)
  real(16) :: a11, a22, a31, a42
  real(16) :: A(4,4), B(4)
  integer :: i, j
  integer :: ipiv(4)
  integer :: ierr
  logical linint  ! selects the linear interpolation method
  ! Ridders' method
  real(16) :: F(4)
  real(16) :: P(4)
  real(16) :: s
  integer  :: maxiter ! maximum number of iterations
  logical  :: lridd   ! selects Ridders' method
  ! Bisection
  logical  :: lbisec  ! selects bisection

  lsecant = .false.
  linint  = .false.
  lridd   = .false.
  lbisec  = .false.
  ! choose method according to input
  ! if method is not provided: default to Riddler
  if (present(method)) then
    select case (method)
      case (0)
        lsecant = .true.
      case (1)
        linint  = .true.
      case (2)
        lridd   = .true.
      case (3)
        lbisec  = .true.
    end select
  else
    lridd = .true.
  endif

  mu_qp = real(mu,16) ! save into a local qp number

! deviation from set particle number with initial mu
  call ndeviation(mu_qp, iT, ek, sct, mesh, thdr, ltetra, target_zero1)

  target_zero2=target_zero1
  mu1=mu_qp
  mu2=mu_qp

  do while (target_zero2.gt.0.q0)
     mu2=mu2+0.005q0
     call ndeviation(mu2, iT, ek, sct, mesh, thdr, ltetra, target_zero2)
  enddo
  do while (target_zero1.le.0.q0)
     mu1=mu1-0.005q0
     call ndeviation(mu1, iT, ek, sct, mesh, thdr, ltetra, target_zero1)
  enddo

  niit0=niitQ

  if (lsecant) then
  !Secant root finding
    do iit=1,niit0
       mu_qp=mu1-target_zero1*mu2/(target_zero2-target_zero1)+target_zero1*mu1/(target_zero2-target_zero1)
       call ndeviation(mu_qp, iT, ek, sct, mesh, thdr, ltetra, target_zero)

       if (abs(target_zero).lt.dev) exit
       if (target_zero.gt.0.q0) then
          mu1=mu_qp
          target_zero1=target_zero
          ! call ndeviation(mu2, iT, ek, sct, mesh, thdr, ltetra, target_zero2)
       else
          mu2=mu_qp
          target_zero2=target_zero
          ! call ndeviation(mu1, iT, ek, sct, mesh, thdr, ltetra, target_zero1)
       endif
    enddo
    niitact=iit


  elseif (linint) then
    if (myid.eq.master) then
      write(*,*) 'Linear interpolation root finding method not implemented for quadruple precision'
      write(*,*) 'Exiting.'
      stop
    endif
  elseif(lridd) then   !Ridders' method for root finding
    !write(*,*) 'Ridders search for chemical potential'

! initialise the varibles
    P(1)=mu1 ; P(2)=mu2
    F(1)= target_zero1; F(2)= target_zero2
    maxiter= niitQ

     do j = 1, maxiter
        P(3) = 0.5q0*(P(1)+P(2))
        call ndeviation(P(3), iT, ek, sct, mesh, thdr, ltetra, F(3))
        s = sqrt((F(3)**2)-(F(1)*F(2)))
        if (s==0.0q0) then
           write(*,*) 'Error in Ridders search for chemical potential'
           write(*,*) 'ITER', j, 'x1', P(1),'  x2',P(2),'  x3', P(3)
           write(*,*) 'ITER', j, 'F1', F(1),'  F2',F(2),'  F3', F(3)
           goto 400
        endif
        P(4) = P(3)+(P(3)-P(1))*(sign(1.0q0,F(1)-F(2))*F(3)/s)
        call ndeviation(P(4), iT, ek, sct, mesh, thdr, ltetra, F(4))
        if (abs(F(4)) .le. dev) goto 400
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
     enddo ! over number of iterations

 400 if (j == maxiter) write(*,*) 'Ridders seach might not have converged'

     ! save the values of the intersection
     mu = P(4)
     niitact = j
     niit0   = maxiter
     target_zero = F(4)

  elseif(lbisec) then
    ! Bisection root finding
    do iit=1,niit0
       mu_qp = (mu1+mu2)/2.q0
       call ndeviation(mu_qp, iT, ek, sct, mesh, thdr, ltetra, target_zero)

       if (abs(target_zero).lt.dev) exit
       if (target_zero.gt.0.q0) then
          mu1=mu_qp
          target_zero1=target_zero
       else
          mu2=mu_qp
          target_zero2=target_zero
       endif
    enddo
    niitact = iit
  endif

  if (myid .eq. master .and. (niitact .ge. niit0 .or. abs(target_zero) .ge. dev)) then
     write(*,'(A,1E20.12)') "WARNING: diminished root precision. ndevQ_actual =",real(target_zero,8)
     write(*,'(A,1F10.3,A,1I5,A,1E20.12)') "at T=",T, " with  niitQ=",niitQ, " ndevQ =", real(dev,8)
     write(*,*) "increase niitQ, or allow for bigger ndevQ (see params.F90)"
  endif

  mu = real(mu_qp, 8) ! transform back to dp
end subroutine find_mu_Q



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



! calculate the occupation and return the difference
! required electrons - current electrons
! if positive we have to increase the chemical potential
! if negative we have to decrease the chemical potential
subroutine ndeviation_D(mu, iT, ek, sct, mesh, thdr, ltetra, target_zero)
  implicit none

  real(8), intent(in)  :: mu
  integer, intent(in)  :: iT
  logical, intent(in)  :: ltetra
  real(8), intent(out) :: target_zero
  type(edisp)          :: ek
  type(scatrate)       :: sct
  type(kpointmesh)     :: mesh
  type(tetramesh)      :: thdr

  real(8) :: occ_tot

  if (ltetra) then
     call occ_tet(mu, iT, ek, sct, thdr, occ_tot)
  else
     call occ(mu, iT, ek, sct, mesh, occ_tot)
  endif
  target_zero=ek%nelect-occ_tot
end subroutine ndeviation_D

subroutine ndeviation_Q(mu, iT, ek, sct, mesh, thdr, ltetra, target_zero)
  implicit none

  !passed variables
  real(16), intent(in)   :: mu
  integer, intent(in)   :: iT
  logical, intent(in)   :: ltetra
  real(16), intent(out) :: target_zero
  type(edisp)           :: ek
  type(scatrate)        :: sct
  type(kpointmesh)      :: mesh
  type(tetramesh)       :: thdr

  real(16) :: occ_tot

  if (ltetra) then
     call occ_tet(mu, iT, ek, sct, thdr, occ_tot)
  else
     call occ(mu, iT, ek, sct, mesh, occ_tot)
  endif
  target_zero=real(ek%nelect,16)-occ_tot
end subroutine ndeviation_Q

! for a given chemical potential mu and Temperature iT
! calculate the occupation and save it in ek%occ_tot
subroutine occ_D(mu, iT, ek, sct, mesh, occ_tot)
  implicit none

  integer, intent(in)  :: iT
  real(8), intent(in)  :: mu
  real(8), intent(out) :: occ_tot
  type(edisp)          :: ek
  type(scatrate)       :: sct
  type(kpointmesh)     :: mesh
!local variables
  real(8),parameter :: thr = 1.0d-10
  complex(8) :: z
  real(8) :: nsmall, nbig, eps, tmp
  real(8) :: occ_loc
  integer :: iband, ik, ikx, iky, ikz
  integer :: ktot !total number of k-points
!external variables
  complex(8), external :: wpsipg

  ek%occ_tot=0.0d0
  occ_loc=0.0d0

  nbig=0.d0
  nsmall=0.d0
  ktot=mesh%kx*mesh%ky*mesh%kz

  do ik = iqstr, iqend
     do iband=1,ek%nband_max
        if (ek%band(ik,iband) .gt. band_fill_value) cycle
        if (iband<ek%nbopt_min) cycle
        if (iband>ek%nbopt_max) cycle

        eps=(ek%z(ik,iband)*ek%band(ik,iband))-mu

        ! the digamma function transitions to the fermi function as Gamma -> 0
        if ((sct%gam(iT).eq.0.d0) .and. (.not. allocated(sct%ykb))) then
           tmp=fermi(eps,beta)
        else if (allocated(sct%ykb)) then
           z=0.5d0 + (ek%z(ik,iband)*(sct%gam(iT)+sct%ykb(iT,ik,iband))+ci*eps)*beta2p ! eps --> -eps
           tmp=0.5d0+aimag(wpsipg(z,0))/pi ! >0
        else
           z=0.5d0 + (ek%z(ik,iband)*sct%gam(iT)+ci*eps)*beta2p ! eps --> -eps
           tmp=0.5d0+aimag(wpsipg(z,0))/pi ! >0
        endif

        ! here we separate the 'big' and 'small' values
        ! in order to not mess with the significant digits
        if (tmp.gt.thr) then
           nbig=nbig+tmp
        else
           nsmall=nsmall+tmp
        endif
     enddo ! iband
  enddo !k-points

  !NORMALISATION AND ACCUMULATION
  nbig=2.d0*nbig/real(ktot,8)
  nsmall=2.d0*nsmall/real(ktot,8)

#ifdef MPI
  occ_loc=nbig+nsmall
  call MPI_ALLREDUCE(occ_loc, occ_tot, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, mpierr)
  ek%occ_tot = occ_tot
#else
  occ_tot=nbig+nsmall
  ek%occ_tot=occ_tot
#endif

end subroutine occ_D


subroutine occ_Q(mu, iT, ek, sct, mesh, occ_tot)
  implicit none

  real(16), intent(in)   :: mu
  integer, intent(in)   :: iT
  real(16), intent(out) :: occ_tot
  type(edisp)           :: ek
  type(scatrate)        :: sct
  type(kpointmesh)      :: mesh

  real(16) :: occ_loc
!local variables
  real(16),parameter :: thr = 1.0q-15
  complex(16) ::z
  real(16) :: nsmall, nbig, tmp, tmp2
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
  complex(16), external :: wpsipghp

  occ_tot=0.0q0
  occ_loc=0.0q0

  nbig=0.q0
  nsmall=0.q0
  ktot = mesh%kx*mesh%ky*mesh%kz
  cutQ=1.q0/thr

  do ik = iqstr, iqend
     do iband=1,ek%nband_max
        if (ek%band(ik,iband) .gt. band_fill_value) cycle
        if (iband<ek%nbopt_min) cycle
        if (iband>ek%nbopt_max) cycle
        eps=(ek%z(ik,iband)*ek%band(ik,iband))-mu

        if ((sct%gam(iT).eq.0.d0) .and. (.not. allocated(sct%ykb))) then
           tmp=fermi(eps,betaQ)
        else if (allocated(sct%ykb)) then
           z=0.5q0 + real(ek%z(ik,iband)*(sct%gam(iT)+sct%ykb(iT,ik,iband))*beta2p,16) + &
                      ciQ*real(eps*beta2p,16) ! eps --> -eps
           tmp=0.5q0+aimag(wpsipghp(z,0))/piQ ! >0
        else
           z=0.5q0 + real(ek%z(ik,iband)*sct%gam(iT)*beta2p,16) + &
                      ciQ*real(eps*beta2p,16) ! eps --> -eps
           tmp=0.5q0+aimag(wpsipghp(z,0))/piQ ! >0
        endif

        if (tmp.gt.thr) then
           IEXP=int(log10(abs(tmp)),8)
           tmp2=( real(int((tmp/(10.q0**iEXP))*QCUT,8),16)*10.q0**iEXP ) / QCUT
           tmp=tmp-tmp2
           nbig=nbig+tmp2
           nsmall=nsmall+tmp
        else
           nsmall=nsmall+tmp
        endif

     enddo ! iband
  enddo
  nbig=2.q0*nbig/real(ktot,16)
  nsmall=2.q0*nsmall/real(ktot,16)
  !ek%occ_tot=real(ninteger+nbig+nsmall,8)

#ifdef MPI
  occ_loc=nbig+nsmall
  call MPI_reduce_quad(occ_loc, occ_tot)
  ek%occ_tot = real(occ_tot,8)
#else
  occ_tot=nbig+nsmall
  ek%occ_tot = real(occ_tot,8)
#endif

end subroutine occ_Q


subroutine occ_tet_D(mu, iT, ek, sct, thdr, occ_tot)
  implicit none

  integer, intent(in)  :: iT
  real(8), intent(in)  :: mu
  real(8), intent(out) :: occ_tot
  type(edisp)          :: ek
  type(scatrate)       :: sct
  type(tetramesh)      :: thdr

!local variables
  real(8),parameter :: thr = 1.0d-30
  complex(8) :: z
  real(8) :: nsmall, nbig, ninteger, eps, tmp
  real(8) :: occ_loc, occ_intp, occ_tet(4)
  integer :: iband, ik, itet, kp
!external variables
  complex(8), external :: wpsipg

  ek%occ_tot=0.0d0
  occ_loc=0.0d0
  do itet = iqstr, iqend
     occ_intp=0.0d0
     do ik = 1, 4
        nbig=0.d0
        nsmall=0.d0
        kp=thdr%idtet(ik,itet)
        do iband=1,ek%nband_max
           if (ek%band(thdr%idtet(ik,itet),iband) .gt. 90.0d0) cycle
           if (iband<ek%nbopt_min) cycle
           if (iband>ek%nbopt_max) cycle

           eps=(ek%z(ik,iband)*ek%band(ik,iband))-mu

           if ((sct%gam(iT).eq.0.d0) .and. (.not. allocated(sct%ykb))) then
              tmp=fermi(eps,beta)
           else if (allocated(sct%ykb)) then
              z=0.5d0 + (ek%z(kp,iband)*(sct%gam(iT)+sct%ykb(iT,kp,iband)) + ci*eps ) * beta2p ! eps --> -eps
              tmp=0.5d0+aimag(wpsipg(z,0))/pi ! >0
           else
              z=0.5d0 + (ek%z(kp,iband)*sct%gam(iT) + ci*eps ) * beta2p ! eps --> -eps
              tmp=0.5d0+aimag(wpsipg(z,0))/pi ! >0
           endif

           if (tmp.gt.thr) then
              nbig=nbig+tmp
           else
              nsmall=nsmall+tmp
           endif
        enddo ! iband

        nbig=2.d0*nbig
        nsmall=2.d0*nsmall
        occ_tet(ik)=nbig+nsmall
     enddo ! corners of the tetrahedron

     ! NOTE: this procedure is different from the one used
     ! for the response function, for the latter the interpolation
     ! is done band by band, whereas here I'm summing over bands
     ! first
     call interptra_mu (thdr%vltet(itet), occ_tet, occ_intp)
     occ_loc = occ_loc + occ_intp
  enddo    ! tetrahedra

#ifdef MPI
  call MPI_ALLREDUCE(occ_loc, occ_tot, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, mpierr)
  ek%occ_tot = occ_tot
#else
  occ_tot = occ_loc
  ek%occ_tot = occ_loc
#endif

end subroutine occ_tet_D

subroutine occ_tet_Q(mu, iT, ek, sct, thdr, occ_tot)
  implicit none

  integer, intent(in)   :: iT
  real(16), intent(in)   :: mu
  real(16), intent(out) :: occ_tot
  type(edisp)           :: ek
  type(scatrate)        :: sct
  type(tetramesh)       :: thdr

!local variables
  real(16),parameter :: thr = 1.0q-30
  complex(16) :: z
  real(16) :: nsmall, nbig, ninteger, tmp, tmp2
  real(16) :: occ_loc, occ_tet(4), occ_intp
  real(8) :: eps
  integer :: iband, ik
  integer :: itet, kp
  real(16) :: cutQ
  !more sophistication
  integer(8) ::IEXP
  !parameters
  real(16),parameter :: QCUT=1.Q14 ! not the same as cutQ!!!!!
!external variables
  complex(16), external :: wpsipghp

  occ_loc=0.0q0
  cutQ=1.q0/thr

  do itet = iqstr, iqend
     occ_intp=0.0q0
     do ik = 1, 4 ! corners of the tetrahedron
        ninteger=0.q0
        nbig=0.q0
        nsmall=0.q0
        kp=thdr%idtet(ik,itet)
        do iband=1,ek%nband_max
           if (ek%band(thdr%idtet(ik,itet),iband) .gt. 90.0d0) cycle
           if (iband<ek%nbopt_min) cycle
           if (iband>ek%nbopt_max) cycle
           eps=(ek%z(kp,iband)*ek%band(kp,iband))-mu

           if ((sct%gam(iT).eq.0.d0) .and. (.not. allocated(sct%ykb))) then
              tmp=fermi(eps,betaQ)
           else if (allocated(sct%ykb)) then
              z=0.5q0+real(ek%z(kp,iband)*(sct%gam(iT)+sct%ykb(iT,kp,iband))*beta2p,16) + &
                       ciQ*real(eps*beta2p,16) ! eps --> -eps
              tmp=0.5q0+aimag(wpsipghp(z,0))/piQ ! >0
           else
              z=0.5q0+real(ek%z(kp,iband)*sct%gam(iT)*beta2p,16) + &
                       ciQ*real(eps*beta2p,16) ! eps --> -eps
              tmp=0.5q0+aimag(wpsipghp(z,0))/piQ ! >0
           endif

           if (tmp.gt.thr) then
              IEXP=int(log10(abs(tmp)),8)
              tmp2=( real(int((tmp/(10.q0**iEXP))*QCUT,8),16)*10.q0**iEXP ) / QCUT
              tmp=tmp-tmp2

              nbig=nbig+tmp2
              nsmall=nsmall+tmp
           else
              nsmall=nsmall+tmp
           endif

        enddo ! iband
        nbig=2.q0*nbig
        nsmall=2.q0*nsmall
        occ_tet(ik)=nbig+nsmall
     enddo ! corners of the tetrahedron

     call interptra_muQ (thdr%vltet(itet), occ_tet, occ_intp)
     occ_loc = occ_loc + occ_intp
  enddo ! tetrahedra

#ifdef MPI
  call mpi_reduce_quad(occ_loc, occ_tot)
  ek%occ_tot = real(occ_tot,8)
#else
  occ_tot = occ_loc
  ek%occ_tot = real(occ_tot,8)
#endif

end subroutine occ_tet_Q

pure elemental function fermi_d(eps,beta) result(f)
  implicit none
  real(8), intent(in) :: eps
  real(8), intent(in) :: beta
  real(8)             :: f
  f = 1.d0 / (1.d0 + EXP(beta*eps))
end function fermi_d

pure elemental function fermi_q(eps,beta) result(f)
  implicit none
  real(8), intent(in)  :: eps
  real(16), intent(in) :: beta
  real(16)             :: f
  f = 1.q0 / (1.q0 + EXP(real(beta*eps,16)))
end function fermi_q

end module Mroot
