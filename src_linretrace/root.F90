module Mroot
  use Mmpi_org
  use Mparams
  use Mtypes
  use Mfermi
  use Maux

  interface find_mu
    module procedure find_mu_D, find_mu_Q
  end interface find_mu

  interface ndeviation
    module procedure ndeviation_D, ndeviation_Q
  end interface ndeviation

  interface occ_digamma
    module procedure occ_digamma_D, occ_digamma_Q
  end interface occ_digamma

  interface occ_fermi
    module procedure occ_fermi_D, occ_fermi_Q
  end interface occ_fermi

  contains

subroutine find_mu_D(mu,dev,target_zero,niitact, edisp, sct, kmesh, imp, algo, info)
  implicit none

  real(8), intent(inout)        :: mu ! chemical potential which is calculated
  real(8), intent(in)           :: dev    ! allowed deviation
  real(8), intent(out)          :: target_zero ! deviation from root after convergence
  integer, intent(out)          :: niitact ! number of iterations

  type(algorithm)  :: algo
  type(energydisp) :: edisp
  type(kpointmesh) :: kmesh
  type(scattering) :: sct
  type(impurity)   :: imp
  type(runinfo)    :: info

  ! local variables
  real(8)  target_zero1, target_zero2, mu1, mu2
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
 select case (algo%rootMethod)
   case (0)
     lsecant = .true.
   case (1)
     linint  = .true.
   case (2)
     lridd   = .true.
   case (3)
     lbisec  = .true.
   case default
     call stop_with_message(stderr, "Root-finding method not properly defined")
 end select

! deviation from set particle number with initial mu
! output: target_zero1 = required - actual electrons
  call ndeviation(mu, target_zero1, edisp, sct, kmesh, imp, algo, info)

! coarse initialization of secant bracket mu1, mu2
  target_zero2=target_zero1
  mu1=mu
  mu2=mu

  do while (target_zero2.gt.0.d0) ! too few electrons -> increase mu
     mu2 = mu2 + 0.5d0
     call ndeviation(mu2, target_zero2, edisp, sct, kmesh, imp, algo, info)
  enddo
  do while (target_zero1.le.0.d0) ! too many electrons -> decrease mu
     mu1 = mu1 - 0.5d0
     call ndeviation(mu1, target_zero1, edisp, sct, kmesh, imp, algo, info)
  enddo

  ! maximal steps for double precision calculations
  niit0=niit

  if (lsecant) then
  !Secant root finding
    do iit=1,niit0
       mu=mu1-target_zero1*(mu2-mu1)/(target_zero2-target_zero1)
       call ndeviation(mu, target_zero, edisp, sct, kmesh, imp, algo, info)

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
       call ndeviation(X(i), Y(i), edisp, sct, kmesh, imp, algo, info)
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
                call ndeviation(mu, target_zero, edisp, sct, kmesh, imp, algo, info)
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
        call ndeviation(P(3), F(3), edisp, sct, kmesh, imp, algo, info)
        s = dsqrt((F(3)**2)-(F(1)*F(2)))
        if (s==0.0d0) then
           write(*,*) 'Error in Ridders search for chemical potential'
           write(*,*) 'ITER', j, 'x1', P(1),'  x2',P(2),'  x3', P(3)
           write(*,*) 'ITER', j, 'F1', F(1),'  F2',F(2),'  F3', F(3)
           goto 400
        endif
        P(4) = P(3)+(P(3)-P(1))*(SIGN(1.0d0,F(1)-F(2))*F(3)/s)
        call ndeviation(P(4), F(4), edisp, sct, kmesh, imp, algo, info)
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
       call ndeviation(mu, target_zero, edisp, sct, kmesh, imp, algo, info)
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



  ! if (myid.eq.master .and. (niitact .ge. niit0 .or. abs(target_zero) .ge. dev)) then
  !   write(*,'(A,1E20.12)') "WARNING: diminished root precision. ndev_actual =",target_zero
  !   write(*,'(A,1F10.3,A,1I5,A,1E20.12)') "at T=",sct%TT(iT), " with  niit=",niit0, " ndev =", dev
  !   write(*,*) "increase niit, or allow for bigger ndev (see params.F90)"
  ! endif
end subroutine find_mu_D

subroutine find_mu_Q(mu,dev,target_zero,niitact, edisp, sct, kmesh, imp, algo, info)
  implicit none
  ! passed variables
  real(8), intent(inout)        :: mu
  real(16), intent(in)          :: dev         ! allowed deviation
  real(16), intent(out)         :: target_zero ! actual deviation
  integer, intent(out)          :: niitact

  type(energydisp) :: edisp
  type(scattering) :: sct
  type(kpointmesh) :: kmesh
  type(impurity)   :: imp
  type(algorithm)  :: algo
  type(runinfo)    :: info


  ! local variables
  real(16) mu_qp
  real(16) target_zero1, target_zero2
  real(16) mu1, mu2, dmu
  real(16) :: target_zero_old
  integer iit, niit0, itest
  logical :: switchsign
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
  select case (algo%rootmethod)
     case (0)
        lsecant = .true.
     case (1)
        linint  = .true.
     case (2)
        lridd   = .true.
     case (3)
        lbisec  = .true.
     case default
       call stop_with_message(stderr, "Root-finding method not properly defined")
  end select

  mu_qp = real(mu,16) ! save into a local qp number

! deviation from set particle number with initial mu
  call ndeviation(mu_qp, target_zero1, edisp, sct, kmesh, imp, algo, info)

  target_zero2=target_zero1
  mu1=mu_qp
  mu2=mu_qp

  do while (target_zero2.gt.0.q0)
     mu2=mu2+0.5q0
     call ndeviation(mu2, target_zero2, edisp, sct, kmesh, imp, algo, info)
  enddo
  do while (target_zero1.le.0.q0)
     mu1=mu1-0.5q0
     call ndeviation(mu1, target_zero1, edisp, sct, kmesh, imp, algo, info)
  enddo

  niit0=niitQ

  if (lsecant) then
  !Secant root finding
    do iit=1,niit0
       mu_qp=mu1-target_zero1*mu2/(target_zero2-target_zero1)+target_zero1*mu1/(target_zero2-target_zero1)
       call ndeviation(mu_qp, target_zero, edisp, sct, kmesh, imp, algo, info)

       if (abs(target_zero).lt.dev) exit
       if (target_zero.gt.0.q0) then
          mu1=mu_qp
          target_zero1=target_zero
       else
          mu2=mu_qp
          target_zero2=target_zero
       endif
    enddo
    niitact=iit


  elseif (linint) then
    call stop_with_message(stderr, 'Linear interpolation root finding method not implemented for quadruple precision')

  elseif(lridd) then   !Ridders' method for root finding
    !write(*,*) 'Ridders search for chemical potential'

! initialise the varibles
    P(1)=mu1 ; P(2)=mu2
    F(1)= target_zero1; F(2)= target_zero2
    maxiter= niitQ

     do j = 1, maxiter
        P(3) = 0.5q0*(P(1)+P(2))
        call ndeviation(P(3), F(3), edisp, sct, kmesh, imp, algo, info)
        s = sqrt((F(3)**2)-(F(1)*F(2)))
        if (s==0.0q0) then
           write(*,*) 'Error in Ridders search for chemical potential'
           write(*,*) 'ITER', j, 'x1', P(1),'  x2',P(2),'  x3', P(3)
           write(*,*) 'ITER', j, 'F1', F(1),'  F2',F(2),'  F3', F(3)
           goto 400
        endif
        P(4) = P(3)+(P(3)-P(1))*(sign(1.0q0,F(1)-F(2))*F(3)/s)
        call ndeviation(P(4), F(4), edisp, sct, kmesh, imp, algo, info)
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
     mu_qp = P(4)
     niitact = j
     niit0   = maxiter
     target_zero = F(4)

  elseif(lbisec) then
    ! Bisection root finding
    do iit=1,niit0
       mu_qp = (mu1+mu2)/2.q0
       call ndeviation(mu_qp, target_zero, edisp, sct, kmesh, imp, algo, info)

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



  ! now lets try something fun
  dmu = 0.05q0
  mu1 = mu_qp
  mu2 = mu_qp
  call occ_fermi_Q_refine(mu_qp, target_zero1, edisp, sct, kmesh, algo, info)
  ! if (abs(target_zero).lt. ndevVVQ) return
  if (target_zero1 < 0.q0) then
    dmu = dmu
    target_zero2 = target_zero1
  else
    dmu = -dmu
    target_zero2 = target_zero1
  endif
  do while ((target_zero1 <= 0.q0 .and. target_zero2 <= 0.q0) .or. &
            (target_zero1 >= 0.q0  .and. target_zero2 >= 0.q0))
    mu2 = mu2 + dmu
    call occ_fermi_Q_refine(mu2, target_zero2, edisp, sct, kmesh, algo, info)
  enddo

  ! write(*,*) mu1, mu2
  ! write(*,*) target_zero1, target_zero2

  ! write(*,*)
  ! write(*,*) mu1, mu2
  ! write(*,*) target_zero1, target_zero2

! #if TRUE
  target_zero_old = 0.q0
  do
     mu_qp = (mu1+mu2)/2.q0
     call occ_fermi_Q_refine(mu_qp, target_zero, edisp, sct, kmesh, algo, info)

     if (target_zero == target_zero_old) exit
     target_zero_old = target_zero
     if ((target_zero .gt. 0.q0 .and. target_zero2.gt. 0.q0) &
          .or. (target_zero .lt. 0.q0 .and. target_zero2 .lt. 0.q0)) then
        mu2=mu_qp
        target_zero2=target_zero
     else
        mu1=mu_qp
        target_zero1=target_zero
     endif
     ! write(*,*) mu1, mu2
     ! write(*,*) target_zero1, target_zero2
  enddo
  ! write(*,*) mu_qp, target_zero
  niitact = iit
! #endif



  ! write(*,*) mu1, mu2
  ! write(*,*) target_zero1, target_zero2

  ! do iit=1,10000
  !   mu_qp = mu_qp + dmu
  !   call occ_fermi_Q_refine(mu_qp, target_zero2, edisp, sct, kmesh, algo, info)
  !   if (abs(target_zero2).lt. ndevVVQ) exit
  !   if ((target_zero1 < 0.q0 .and. target_zero2 < 0.q0) &
  !       .or. (target_zero1 > 0.q0 .and. target_zero2 > 0.q0))

  !   if target_zero < 0.q0 then
  !   if (target_zero.gt.0.q0) then
  !      mu1=mu_qp
  !      target_zero1=target_zero
  !   else
  !      mu2=mu_qp
  !      target_zero2=target_zero
  !   endif
  ! enddo

  ! if (myid .eq. master .and. (niitact .ge. niit0 .or. abs(target_zero) .ge. dev)) then
  !    write(*,'(A,1E20.12)') "WARNING: diminished root precision. ndevQ_actual =",real(target_zero,8)
  !    write(*,'(A,1F10.3,A,1I5,A,1E20.12)') "at T=",T, " with  niitQ=",niitQ, " ndevQ =", real(dev,8)
  !    write(*,*) "increase niitQ, or allow for bigger ndevQ (see params.F90)"
  ! endif

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
subroutine ndeviation_D(mu, target_zero, edisp, sct, kmesh, imp, algo, info)
  implicit none

  real(8), intent(in)  :: mu
  real(8), intent(out) :: target_zero

  type(energydisp) :: edisp
  type(scattering) :: sct
  type(kpointmesh) :: kmesh
  type(impurity)   :: imp
  type(algorithm)  :: algo
  type(runinfo)    :: info

  integer :: iimp
  real(8) :: occ_tot

  if (algo%muFermi) then
    call occ_fermi(mu, occ_tot, edisp, sct, kmesh, algo, info)
    ! call occ_fermi_comp_D(mu, occ_tot, edisp, kmesh, info)
  else
    call occ_digamma(mu, occ_tot, edisp, sct, kmesh, algo, info)
    ! call occ_digamma_comp_D(mu, occ_tot, edisp, sct, kmesh, algo, info)
  endif

  ! nvalence = nsearch - N_D^+ + N_A^-
  ! N_D^+ = N_D/(1 + g * exp(beta * (mu - E_D)))
  ! N_A^+ = N_D/(1 + g * exp(-beta * (mu - E_A)))
  if (algo%lImpurities) then
    do iimp = 1,imp%nimp
      occ_tot = occ_tot - imp%Dopant(iimp)*imp%Density(iimp) &
        / (1.d0 + imp%Degeneracy(iimp) * exp(info%beta*imp%Dopant(iimp)*(mu-imp%Energy(iimp))))
    enddo
  endif

  target_zero = edisp%nelect - occ_tot
end subroutine ndeviation_D

subroutine ndeviation_Q(mu, target_zero, edisp, sct, kmesh, imp, algo, info)
  implicit none

  !passed variables
  real(16), intent(in)  :: mu
  real(16), intent(out) :: target_zero

  type(energydisp) :: edisp
  type(scattering) :: sct
  type(kpointmesh) :: kmesh
  type(impurity)   :: imp
  type(algorithm)  :: algo
  type(runinfo)    :: info

  integer  :: iimp
  real(16) :: occ_tot

  if (algo%muFermi) then
    call occ_fermi(mu, occ_tot, edisp, sct, kmesh, algo, info)
  else
    call occ_digamma(mu, occ_tot, edisp, sct, kmesh, algo, info)
  endif

  ! nvalence = nsearch - N_D^+ + N_A^-
  ! N_D^+ = N_D/(1 + g * exp(beta * (mu - E_D)))
  ! N_A^+ = N_D/(1 + g * exp(-beta * (mu - E_A)))
  if (algo%lImpurities) then
    do iimp = 1,imp%nimp
      occ_tot = occ_tot - imp%Dopant(iimp)*imp%Density(iimp) &
        / (1.d0 + imp%Degeneracy(iimp) * exp(info%beta*imp%Dopant(iimp)*(mu-imp%Energy(iimp))))
    enddo
  endif

  target_zero=real(edisp%nelect,16) - occ_tot
end subroutine ndeviation_Q

subroutine occ_digamma_D(mu, occ_tot, edisp, sct, kmesh, algo, info)
  implicit none

  real(8), intent(in)  :: mu
  real(8), intent(out) :: occ_tot

  type(energydisp) :: edisp
  type(scattering) :: sct
  type(kpointmesh) :: kmesh
  type(algorithm)  :: algo
  type(runinfo)    :: info
  !local variables

  real(8) :: occ_loc
  integer :: is, ik, iband
  complex(8), allocatable :: to_evaluate(:,:,:)
  real(8), allocatable    :: occupation(:,:,:)
  !external variables
  complex(8), external :: wpsipg

  allocate(to_evaluate(edisp%nband_max, ikstr:ikend, edisp%ispin))
  allocate(occupation(edisp%nband_max, ikstr:ikend, edisp%ispin))

  to_evaluate = 0.5d0 + info%beta2p * &
                (sct%gam(:,ikstr:ikend,:) - ci*sct%zqp(:,ikstr:ikend,:)*(edisp%band(:,ikstr:ikend,:) - mu))

  ! evaluate the function
  do is = 1,edisp%ispin
    do ik = ikstr, ikend
      do iband=1,edisp%nband_max
        occupation(iband,ik,is) = 0.5d0 + aimag(wpsipg(to_evaluate(iband,ik,is),0))/pi
        occupation(iband,ik,is) = occupation(iband,ik,is) * kmesh%weight(ik)
      enddo
    enddo
  enddo

  deallocate(to_evaluate)
  occ_loc = sum(occupation)
  deallocate(occupation)

#ifdef MPI
  call MPI_ALLREDUCE(occ_loc, occ_tot, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, mpierr)
#else
  occ_tot = occ_loc
#endif

end subroutine occ_digamma_D


subroutine occ_digamma_Q(mu, occ_tot, edisp, sct, kmesh, algo, info)
  implicit none

  real(16), intent(in)  :: mu
  real(16), intent(out) :: occ_tot

  type(energydisp) :: edisp
  type(scattering) :: sct
  type(kpointmesh) :: kmesh
  type(algorithm)  :: algo
  type(runinfo)    :: info

  !local variables
  real(16) :: occ_loc
  integer  :: iband, is, ik

  complex(16), allocatable :: to_evaluate(:,:,:)
  real(16), allocatable    :: occupation(:,:,:)
  !external variables
  complex(16), external :: wpsipghp


  allocate(to_evaluate(edisp%nband_max, ikstr:ikend, edisp%ispin))
  allocate(occupation(edisp%nband_max, ikstr:ikend, edisp%ispin))

  to_evaluate = 0.5q0 + info%beta2pQ * &
                (sct%gam(:,ikstr:ikend,:) - ciQ*sct%zqp(:,ikstr:ikend,:)*(edisp%band(:,ikstr:ikend,:) - mu))

  ! evaluate the function
  do is = 1,edisp%ispin
    do ik = ikstr, ikend
      do iband=1,edisp%nband_max
        occupation(iband,ik,is) = 0.5q0 + aimag(wpsipghp(to_evaluate(iband,ik,is),0))/piQ
        occupation(iband,ik,is) = occupation(iband,ik,is) * kmesh%weight(ik)
      enddo
    enddo
  enddo

  deallocate(to_evaluate)
  occ_loc = sum(occupation)
  deallocate(occupation)

#ifdef MPI
  call MPI_reduce_quad(occ_loc, occ_tot)
#else
  occ_tot = occ_loc
#endif

end subroutine occ_digamma_Q

subroutine occ_fermi_D(mu, occ_tot, edisp, sct, kmesh, algo, info)
  implicit none

  real(8), intent(in)  :: mu
  real(8), intent(out) :: occ_tot

  type(energydisp) :: edisp
  type(scattering) :: sct
  type(kpointmesh) :: kmesh
  type(algorithm)  :: algo
  type(runinfo)    :: info
  !local variables

  real(8) :: occ_loc
  integer :: is, ik, iband

  real(8), allocatable :: occupation(:,:,:)

  allocate(occupation(edisp%nband_max, ikstr:ikend, edisp%ispin))

  ! evaluate the function
  do is = 1,edisp%ispin
    do ik = ikstr, ikend
      do iband=1,edisp%nband_max
        occupation(iband,ik,is) = fermi(sct%zqp(iband,ik,is)*(edisp%band(iband,ik,is)-mu), info%beta)
        occupation(iband,ik,is) = occupation(iband,ik,is) * kmesh%weight(ik)
      enddo
    enddo
  enddo

  occ_loc = sum(occupation)
  occ_loc = occ_loc + edisp%ispin * (ikend-ikstr+1)*edisp%nband_max * 0.5d0
  deallocate(occupation)

#ifdef MPI
  call MPI_ALLREDUCE(occ_loc, occ_tot, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, mpierr)
#else
  occ_tot = occ_loc
#endif

end subroutine occ_fermi_D

! Neumaier algorithm to increase summation accuracy
subroutine occ_fermi_comp_D(mu, occ_tot, edisp, kmesh, algo, info)
  implicit none

  real(8), intent(in)  :: mu
  real(8), intent(out) :: occ_tot

  type(energydisp) :: edisp
  type(kpointmesh) :: kmesh
  type(runinfo)    :: info
  type(algorithm)  :: algo
  !local variables

  real(8) :: occ_sum
  integer :: is, ik, iband

  real(8), allocatable :: occupation(:,:,:)
  real(8) :: t,c

  allocate(occupation(edisp%nband_max, ikstr:ikend, edisp%ispin))

  ! evaluate the function
  occ_sum = 0.d0
  t = 0.d0
  c = 0.d0
  do is = 1,edisp%ispin
    do ik = ikstr, ikend
      do iband=1,edisp%nband_max
        occupation(iband,ik,is) = fermi((edisp%band(iband,ik,is)-mu), info%beta)
        occupation(iband,ik,is) = occupation(iband,ik,is) * kmesh%weight(ik)

        t = occ_sum + occupation(iband,ik,is)
        if (abs(occ_sum) >= abs(occupation(iband,ik,is))) then
          c = c + (occ_sum - t) + occupation(iband,ik,is)
        else
          c = c + (occupation(iband,ik,is) - t) + occ_sum
        endif
        occ_sum = t

      enddo
    enddo
  enddo

  occ_sum = occ_sum + c
  deallocate(occupation)

#ifdef MPI
  call MPI_ALLREDUCE(occ_sum, occ_tot, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, mpierr)
#else
  occ_tot = occ_sum
#endif

end subroutine occ_fermi_comp_D

subroutine occ_fermi_Q(mu, occ_tot, edisp, sct, kmesh, algo, info)
  implicit none

  real(16), intent(in)  :: mu
  real(16), intent(out) :: occ_tot

  type(energydisp) :: edisp
  type(scattering) :: sct
  type(kpointmesh) :: kmesh
  type(runinfo)    :: info
  type(algorithm)  :: algo

  logical:: ingap
  !local variables

  real(16) :: occ_loc
  integer :: is, ik, iband
  real(16) :: locecc, lochole

  real(16), allocatable :: electrons(:,:,:)
  real(16), allocatable :: holes(:,:,:)
  real(16), allocatable :: occupation(:,:,:)

  allocate(occupation(edisp%nband_max, ikstr:ikend, edisp%ispin))
  allocate(electrons(edisp%nband_max, ikstr:ikend, edisp%ispin))
  allocate(holes(edisp%nband_max, ikstr:ikend, edisp%ispin))

  ingap = .false.

  if (.not. ingap) then
    do is = 1,edisp%ispin
      do ik = ikstr, ikend
        do iband=1,edisp%nband_max
          ! directly call the specific fermi function in order to avoid unnecessary many
          ! vtable look-ups
          occupation(iband,ik,is) = fermi(sct%zqp(iband,ik,is)*(edisp%band(iband,ik,is)-mu), info%betaQ)
          occupation(iband,ik,is) = occupation(iband,ik,is) * kmesh%weightQ(ik)

          electrons(iband,ik,is) = fermi(sct%zqp(iband,ik,is)*(edisp%band(iband,ik,is)-mu), info%betaQ)
          electrons(iband,ik,is) = electrons(iband,ik,is) * kmesh%weightQ(ik)

          holes(iband,ik,is) = omfermi(sct%zqp(iband,ik,is)*(edisp%band(iband,ik,is)-mu), info%betaQ)
          holes(iband,ik,is) = holes(iband,ik,is) * kmesh%weightQ(ik)
        enddo
      enddo
    enddo
    occ_loc = sum(occupation)

    do iband=1,edisp%nband_max
      if (abs(2.q0 - sum(electrons(iband,:,:))) /= 2.q0) then! significant digits
        locecc = 0.q0
      else
        locecc = sum(electrons(iband,:,:))
      endif
      if (abs(2.q0 - sum(holes(iband,:,:))) /= 2.q0) then! significant digits
        lochole = 0.q0
      else
        lochole = sum(holes(iband,:,:))
      endif
      ! write(*,*) iband, sum(electrons(iband,:,:)), sum(holes(iband,:,:))
    enddo
  else
    allocate(electrons(edisp%nband_max, ikstr:ikend, edisp%ispin))
    allocate(holes(edisp%nband_max, ikstr:ikend, edisp%ispin))
    electrons = 0.q0
    holes = 0.q0
    do is = 1,edisp%ispin
      do ik = ikstr, ikend
        do iband=1,edisp%nband_max
          electrons(iband,ik,is) = fermi(sct%zqp(iband,ik,is)*(edisp%band(iband,ik,is)-mu), info%betaQ)
          electrons(iband,ik,is) = electrons(iband,ik,is) * kmesh%weightQ(ik)

          holes(iband,ik,is) = omfermi(sct%zqp(iband,ik,is)*(edisp%band(iband,ik,is)-mu), info%betaQ)
          holes(iband,ik,is) = holes(iband,ik,is) * kmesh%weightQ(ik)
        enddo
      enddo
    enddo

    ! do iband=1,edisp%nband_max
    !   write(*,*) iband, sum(occupation(iband,:,:))
    ! enddo
    occ_loc = sum(electrons) - sum(holes)
    ! write(*,*) mu, sum(electrons), sum(holes)
    occ_loc = occ_loc + int(sum(occupation))
    ! write(*,*) mu, occ_loc
    deallocate(electrons)
    deallocate(holes)
  endif

  deallocate(occupation)

#ifdef MPI
  call mpi_reduce_quad(occ_loc, occ_tot) ! custom quad reduction
#else
  occ_tot = occ_loc
#endif

end subroutine occ_fermi_Q

subroutine occ_fermi_Q_refine(mu, deviation, edisp, sct, kmesh, algo, info)
  implicit none

  real(16), intent(in)  :: mu
  real(16), intent(out) :: deviation

  type(energydisp) :: edisp
  type(scattering) :: sct
  type(kpointmesh) :: kmesh
  type(runinfo)    :: info
  type(algorithm)  :: algo

  logical:: ingap
  !local variables

  real(16) :: deviation_loc
  integer :: is, ik, iband
  real(16) :: locelec, lochole
  real(16) :: sumelec, sumhole

  real(16), allocatable :: electrons(:,:,:)
  real(16), allocatable :: holes(:,:,:)

  allocate(electrons(edisp%nband_max, ikstr:ikend, edisp%ispin))
  allocate(holes(edisp%nband_max, ikstr:ikend, edisp%ispin))

  do is = 1,edisp%ispin
    do ik = ikstr, ikend
      do iband=1,edisp%nband_max
        ! directly call the specific fermi function in order to avoid unnecessary many
        ! vtable look-ups
        electrons(iband,ik,is) = fermi(sct%zqp(iband,ik,is)*(edisp%band(iband,ik,is)-mu), info%betaQ)
        electrons(iband,ik,is) = electrons(iband,ik,is) * kmesh%weightQ(ik)

        holes(iband,ik,is) = omfermi(sct%zqp(iband,ik,is)*(edisp%band(iband,ik,is)-mu), info%betaQ)
        holes(iband,ik,is) = holes(iband,ik,is) * kmesh%weightQ(ik)
      enddo
    enddo
  enddo

  sumelec = 0.q0
  sumhole = 0.q0
  do iband=1,edisp%nband_max
    ! large deviation
    ! if (abs(sum(electrons(iband,:,:)) - 2.q0) > 1.q-3) then
    !   lochole = 2.q0 - sum(electrons(iband,:,:)) - sum(holes(iband,:,:))
    !   locelec = 0.q0
    !   goto 400
    ! endif

    locelec = 0.q0
    lochole = 0.q0

    if (abs(2.q0 - sum(electrons(iband,:,:))) /= 2.q0) then! significant digits
      locelec = 0.q0
    else
      locelec = sum(electrons(iband,:,:))
      goto 400
    endif

    if (abs(2.q0 - sum(holes(iband,:,:))) /= 2.q0) then! significant digits
      lochole = 0.q0
    else
      lochole = sum(holes(iband,:,:))
      goto 400
    endif

    ! locelec = sum(electrons(iband,:,:))
    lochole = 0.q0
    locelec = 2.q0 - sum(holes(iband,:,:)) - sum(electrons(iband,:,:))



! 400 write(*,*) iband, locelec, lochole
    ! sumelec = sumelec + locelec
400 sumelec = sumelec + locelec
    sumhole = sumhole + lochole
  enddo

  deallocate(electrons)
  deallocate(holes)

  deviation_loc = sumelec - sumhole


#ifdef MPI
  call mpi_reduce_quad(deviation_loc, deviation) ! custom quad reduction
#else
  deviation = deviation_loc
#endif

end subroutine occ_fermi_Q_refine

subroutine occ_digamma_comp_D(mu, occ_tot, edisp, sct, kmesh, algo, info)
  implicit none

  real(8), intent(in)  :: mu
  real(8), intent(out) :: occ_tot

  type(energydisp) :: edisp
  type(scattering) :: sct
  type(kpointmesh) :: kmesh
  type(algorithm)  :: algo
  type(runinfo)    :: info
  !local variables

  real(8) :: occ_sum
  integer :: is, ik, iband
  complex(8), allocatable :: to_evaluate(:,:,:)
  real(8), allocatable    :: occupation(:,:,:)
  !external variables
  complex(8), external :: wpsipg

  real(8) :: t,c

  allocate(to_evaluate(edisp%nband_max, ikstr:ikend, edisp%ispin))
  allocate(occupation(edisp%nband_max, ikstr:ikend, edisp%ispin))

  to_evaluate = 0.5d0 + info%beta2p * &
                (sct%gam(:,ikstr:ikend,:) - ci*sct%zqp(:,ikstr:ikend,:)*(edisp%band(:,ikstr:ikend,:) - mu))

  ! evaluate the function
  occ_sum = 0.d0
  t = 0.d0
  c = 0.d0
  do is = 1,edisp%ispin
    do ik = ikstr, ikend
      do iband=1,edisp%nband_max
        occupation(iband,ik,is) = 0.5d0 + aimag(wpsipg(to_evaluate(iband,ik,is),0))/pi
        occupation(iband,ik,is) = occupation(iband,ik,is) * kmesh%weight(ik)

        ! what the fuck
        t = occ_sum + occupation(iband,ik,is)
        if (abs(occ_sum) >= abs(occupation(iband,ik,is))) then
          c = c + (occ_sum - t) + occupation(iband,ik,is)
        else
          c = c + (occupation(iband,ik,is) - t) + occ_sum
        endif
        occ_sum = t

      enddo
    enddo
  enddo

  occ_sum = occ_sum + c
  deallocate(to_evaluate)
  deallocate(occupation)

#ifdef MPI
  call MPI_ALLREDUCE(occ_sum, occ_tot, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpierr)
#else
  occ_tot = occ_sum
#endif


end subroutine occ_digamma_comp_D

end module Mroot
