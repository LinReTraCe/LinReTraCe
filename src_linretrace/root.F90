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
  logical :: skipped
  real(16) mu_qp
  real(16) target_zero1, target_zero2
  real(16) test_up, test_dn
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

  ! quadruple precision for digamma function occupation
  ! is enough (even when the system is gapped)
  if (.not. algo%muFermi) then
    mu = real(mu_qp, 8) ! transform back to dp
    return
  endif


  ! mu refinement is numerically unstable below a certain Temperate/Gap ratio
  ! i.e. the fermi function with quadruple precision is not accurate neough
  if (.not. algo%lImpurities .and. info%Temp < edisp%gap_min / 1.95q0) then
    call log_master(stdout, 'Warning: mu-refinement does not work at this temperature')
    mu = real(mu_qp, 8) ! transform back to dp
    return
  endif

  ! perform the mu_refinement if we have a gap
  call ndeviation_Q(mu_qp, target_zero2, edisp, sct, kmesh, imp, algo, info)
  call occ_fermi_Q_refine(mu_qp, target_zero1, edisp, sct, kmesh, imp, algo, info)

  ! write(*,*) 'after qp root finding: ',mu_qp
  ! write(*,*) 'deviation QP: ',target_zero2
  ! write(*,*) 'deviation refine: ',target_zero1
  ! write(*,*)

  ! if (myid.eq.master) write(*,*) target_zero1, target_zero2

  if ( (info%Temp < edisp%gap_min*200) .and. &  ! hard temperature cutoff
       (abs(target_zero1) < 1d-18)) then        ! we have a reasonal value from thisfunction
                                                ! if this value is too high
                                                ! we might run off to another root ..
    ! TODO: more testing here
    ! essentially whats happening is that the deviation is this large
    ! when we are really close to an impurity crossing...
    ! this close around the impurity normal QP is enough
    ! the same thing that happens when the temperature gets too large
    ! -> we get too many contributions and run-off to the wrong root (artificial root)
    ! so this if conditions cuts off the two causes for the same result
    ! -> avoidance of artifical roots
    ! -> results in a smooth chemical potential

    dmu = edisp%gap_min/100.q0
    mu1 = mu_qp
    mu2 = mu_qp
    iit = 1

    if (target_zero1 < 0.q0) then
      dmu = +dmu
    else
      dmu = -dmu
    endif

    ! decrease the step size until we dont cross an impurity level essentially
    ! thats wrong
    do
      call occ_fermi_Q_refine(mu_qp+dmu, test_up, edisp, sct, kmesh, imp, algo, info)
      call occ_fermi_Q_refine(mu_qp-dmu, test_dn, edisp, sct, kmesh, imp, algo, info)
      test_up = test_up - target_zero1
      test_dn = test_dn - target_zero1
      ! debug
      ! check if both point in the same direction
      ! write(*,*) 'test step up: ', mu_qp+dmu, test_up
      ! write(*,*) 'test step dn: ', mu_qp+dmu, test_dn
      ! if they do, decrease step size and try again
      if ( (test_up > 0 .and. test_dn > 0) .or. (test_up < 0 .and. test_dn < 0) ) then
        dmu = dmu/5.q0
      else
        exit
      endif
    enddo


    ! abort if we have a sudden change ( by crossing an impurity level )
    ! if(myid.eq.master) write(*,*) target_test
    if (abs(target_zero) > 1d-15) then
      mu = real(mu_qp, 8)
      return
    endif

    target_zero2 = target_zero1 ! from the top

    ! get the working interval
    do while (((target_zero1 <= 0.q0 .and. target_zero2 <= 0.q0) .or. &
              (target_zero1 >= 0.q0  .and. target_zero2 >= 0.q0)) .and. iit < niitQ)
      mu2 = mu2 + dmu
      call occ_fermi_Q_refine(mu2, target_zero2, edisp, sct, kmesh, imp, algo, info)
      iit = iit + 1
    enddo

    niitact = niitact + iit

    if (iit >= niitQ) then
      call log_master(stdout, 'Warning: mu-refinement did not converge!')
      mu = real(mu_qp, 8) ! transform back to dp
      return
    endif

    mu1 = mu2 - dmu ! one step before
                    ! this ensures the smallest possible working range
                    ! with opposing signs in the root - finding

    ! perform a bisection in the previous working interval
    do iit = 1,niitQ
       mu_qp = (mu1+mu2)/2.q0
       call occ_fermi_Q_refine(mu_qp, target_zero, edisp, sct, kmesh, imp, algo, info)
       niitact = niitact + 1

       if ( abs(mu1-mu2) < 1q-12) exit ! we go all out here

       if ((target_zero .gt. 0.q0 .and. target_zero2.gt. 0.q0) &
            .or. (target_zero .lt. 0.q0 .and. target_zero2 .lt. 0.q0)) then
          mu2=mu_qp
          target_zero2=target_zero ! here we are on the bisection side of mu2
                                   ! we set the middle point as new mu2
       else
          mu1=mu_qp                ! the other way aorund here
          target_zero1=target_zero
       endif
    enddo

    if (iit >= niitQ) then
      call log_master(stdout, 'Warning: mu-refinement did not converge!')
      return
    endif
    niitact = niitact + iit
  endif

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

  integer :: ii
  real(8) :: densii, eneii
  real(8) :: dist
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
      if (.not. imp%Band(iimp)) then
        occ_tot = occ_tot - imp%Dopant(iimp)*imp%Density(iimp) &
          / (1.d0 + imp%Degeneracy(iimp) * exp(info%beta*imp%Dopant(iimp)*(mu-imp%Energy(iimp))))
      else
        if (imp%Bandtype(iimp) == 0) then ! box
          densii = imp%Density(iimp) / 1001.d0
          do ii=-500,500
            eneii  = imp%Energy(iimp) + ii/1000.d0 * imp%Bandwidth(iimp)
            occ_tot = occ_tot - imp%Dopant(iimp)*densii &
              / (1.d0 + imp%Degeneracy(iimp) * exp(info%betaQ*imp%Dopant(iimp)*(mu-eneii)))
          enddo
        else if (imp%Bandtype(iimp) == 1) then ! Lorentzian
                                               ! 1/pi 0.5 Gamma / ((x-xo)**2 + (0.5*Gamma)**2 )
          densii = imp%Density(iimp) / 5001.d0
          do ii=-2500,2500 ! we go to +- 2.5 * Gamma
            eneii  = imp%Energy(iimp) + ii/1000.d0 * imp%Bandwidth(iimp)
            dist   = 1.d0/pi * 0.5d0 * imp%Bandwidth(iimp) / ((eneii - imp%Energy(iimp))**2 + (0.5d0 * imp%Bandwidth(iimp))**2 )

            occ_tot = occ_tot - imp%Dopant(iimp)*densii*dist &
              / (1.d0 + imp%Degeneracy(iimp) * exp(info%betaQ*imp%Dopant(iimp)*(mu-eneii)))
          enddo

        else if (imp%Bandtype(iimp) == 2) then ! Gaussian

          densii = imp%Density(iimp) / 10001.d0
          do ii=-5000,5000 ! we go to +- 2.5 * sigma
            eneii  = imp%Energy(iimp) + ii/1000.d0 * imp%Bandwidth(iimp)
            dist   = 1.d0/(imp%Bandwidth(iimp) * sqrt(2.d0 * pi)) * exp(-0.5d0 * ((eneii - imp%Energy(iimp))/imp%Bandwidth(iimp))**2)
            occ_tot = occ_tot - imp%Dopant(iimp)*densii*dist &
              / (1.d0 + imp%Degeneracy(iimp) * exp(info%betaQ*imp%Dopant(iimp)*(mu-eneii)))
          enddo

        endif
      endif
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

  integer  :: ii
  real(16)  :: eneii, densii
  integer  :: iimp
  real(16) :: occ_tot
  real(16) :: dist

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
      if (.not. imp%Band(iimp)) then
        occ_tot = occ_tot - imp%Dopant(iimp)*imp%Density(iimp) &
          / (1.q0 + imp%Degeneracy(iimp) * exp(info%beta*imp%Dopant(iimp)*(mu-imp%Energy(iimp))))
      else
        if (imp%Bandtype(iimp) == 0) then ! box
          densii = imp%Density(iimp) / 1001.q0
          do ii=-500,500
            eneii  = imp%Energy(iimp) + ii/1000.q0 * imp%Bandwidth(iimp)
            occ_tot = occ_tot - imp%Dopant(iimp)*densii &
              / (1.q0 + imp%Degeneracy(iimp) * exp(info%betaQ*imp%Dopant(iimp)*(mu-eneii)))
          enddo
        else if (imp%Bandtype(iimp) == 1) then ! Lorentzian
                                               ! 1/pi 0.5 Gamma / ((x-xo)**2 + (0.5*Gamma)**2 )
          densii = imp%Density(iimp) / 5001.q0
          do ii=-2500,2500 ! we go to +- 2.5 * Gamma
            eneii  = imp%Energy(iimp) + ii/1000.q0 * imp%Bandwidth(iimp)
            dist   = 1.q0/piQ * 0.5q0 * imp%Bandwidth(iimp) / ((eneii - imp%Energy(iimp))**2 + (0.5q0 * imp%Bandwidth(iimp))**2 )

            occ_tot = occ_tot - imp%Dopant(iimp)*densii*dist &
              / (1.q0 + imp%Degeneracy(iimp) * exp(info%betaQ*imp%Dopant(iimp)*(mu-eneii)))
          enddo

        else if (imp%Bandtype(iimp) == 2) then ! Gaussian

          densii = imp%Density(iimp) / 10001.q0
          do ii=-5000,5000 ! we go to +- 2.5 * sigma
            eneii  = imp%Energy(iimp) + ii/1000.q0 * imp%Bandwidth(iimp)
            dist   = 1.q0/(imp%Bandwidth(iimp) * sqrt(2.q0 * piQ)) * exp(-0.5q0 * ((eneii - imp%Energy(iimp))/imp%Bandwidth(iimp))**2)
            occ_tot = occ_tot - imp%Dopant(iimp)*densii*dist &
              / (1.q0 + imp%Degeneracy(iimp) * exp(info%betaQ*imp%Dopant(iimp)*(mu-eneii)))
          enddo

        endif
      endif
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

subroutine occ_fermi_Q_refine(mu, deviation, edisp, sct, kmesh, imp, algo, info)
  implicit none

  real(16), intent(in)  :: mu
  real(16), intent(out) :: deviation

  type(energydisp) :: edisp
  type(scattering) :: sct
  type(kpointmesh) :: kmesh
  type(impurity)   :: imp
  type(runinfo)    :: info
  type(algorithm)  :: algo

  logical:: ingap
  logical :: found
  !local variables

  integer  :: iimp, ii
  real(16) :: diff = 1q-38
  integer  :: is, ik, iband

  real(16) :: elecmpi, holempi
  real(16) :: sumelec, sumhole
  real(16) :: impelec, imphole

  real(16) :: elec
  real(16) :: hole

  real(16) :: densii, eneii, dist

  sumelec = 0.q0
  sumhole = 0.q0

  do is = 1,edisp%ispin
    do ik = ikstr, ikend
      do iband=1,edisp%nband_max
        ! directly call the specific fermi function in order to avoid unnecessary many
        ! vtable look-ups
        elec = fermi_qp(sct%zqp(iband,ik,is)*(edisp%band(iband,ik,is)-mu), info%betaQ)
        hole = omfermi_qp(sct%zqp(iband,ik,is)*(edisp%band(iband,ik,is)-mu), info%betaQ)

        ! here we take the smaller of the two quantities
        ! and weigh it with the quadruple precision weight
        if (hole > elec) then
          sumelec = sumelec + elec * kmesh%weightQ(ik)
        else
          sumhole = sumhole + hole * kmesh%weightQ(ik)
        endif

      enddo
    enddo
  enddo

#ifdef MPI
  call mpi_reduce_quad(sumelec, elecmpi) ! custom quad reduction
#else
  elecmpi = sumelec
#endif

#ifdef MPI
  call mpi_reduce_quad(sumhole, holempi) ! custom quad reduction
#else
  holempi = sumhole
#endif

  ! deviation purely from the band structure
  deviation =  elecmpi - holempi


  ! now we add the impurity differences
  ! this has to be done after the MPI communication ....
  sumelec = 0.q0
  sumhole = 0.q0

  ! nvalence = nsearch - N_D^+ + N_A^-
  ! N_D^+ = N_D/(1 + g * exp(beta * (mu - E_D)))
  ! N_A^+ = N_D/(1 + g * exp(-beta * (mu - E_A)))
  if (algo%lImpurities) then
    do iimp = 1,imp%nimp
      ! so we are strictly between 0 and 1
      ! for numerical reasons
      if (.not. imp%Band(iimp)) then
        ! for impurity level
        impelec = 1.q0 &
          / (1.q0 + imp%Degeneracy(iimp) * exp(info%betaQ*imp%Dopant(iimp)*(mu-imp%Energy(iimp))))

        imphole = 1.q0 &
          / (1.q0 + imp%Degeneracy(iimp)**(-1.q0) * exp(info%betaQ*imp%Dopant(iimp)*(imp%Energy(iimp)-mu)))

        ! here we apply the signs and the weights (density)
        if (imphole > impelec) then
          sumelec = sumelec - impelec*imp%Dopant(iimp)*imp%Density(iimp)
        else
          sumhole = sumhole - imphole*imp%Dopant(iimp)*imp%Density(iimp)
        endif
      else
        ! for impurity bands
        if (imp%Bandtype(iimp) == 0) then ! box
          densii = imp%Density(iimp) / 1001.d0
          do ii=-500,500
            eneii  = imp%Energy(iimp) + ii/1000.d0 * imp%Bandwidth(iimp)
            impelec = 1.q0 &
              / (1.q0 + imp%Degeneracy(iimp) * exp(info%betaQ*imp%Dopant(iimp)*(mu-eneii)))

            imphole = 1.q0 &
              / (1.q0 + imp%Degeneracy(iimp)**(-1.q0) * exp(info%betaQ*imp%Dopant(iimp)*(eneii-mu)))

            ! here we apply the signs and the weights (density)
            if (imphole > impelec) then
              sumelec = sumelec - impelec*imp%Dopant(iimp)*densii
            else
              sumhole = sumhole - imphole*imp%Dopant(iimp)*densii
            endif
          enddo
        else if (imp%Bandtype(iimp) == 1) then ! Lorentzian
          densii = imp%Density(iimp) / 5001.d0
          do ii=-2500,2500 ! we go to +- 2.5 * Gamma
            eneii  = imp%Energy(iimp) + ii/1000.d0 * imp%Bandwidth(iimp)
            dist   = 1.d0/pi * 0.5d0 * imp%Bandwidth(iimp) / ((eneii - imp%Energy(iimp))**2 + (0.5d0 * imp%Bandwidth(iimp))**2 )

            impelec = 1.q0 &
              / (1.q0 + imp%Degeneracy(iimp) * exp(info%betaQ*imp%Dopant(iimp)*(mu-eneii)))

            imphole = 1.q0 &
              / (1.q0 + imp%Degeneracy(iimp)**(-1.q0) * exp(info%betaQ*imp%Dopant(iimp)*(eneii-mu)))

            ! here we apply the signs and the weights (density)
            if (imphole > impelec) then
              sumelec = sumelec - impelec*imp%Dopant(iimp)*densii*dist
            else
              sumhole = sumhole - imphole*imp%Dopant(iimp)*densii*dist
            endif
          enddo
        else if (imp%Bandtype(iimp) == 2) then ! Gaussian
          densii = imp%Density(iimp) / 10001.d0
          do ii=-5000,5000 ! we go to +- 2.5 * sigma
            eneii  = imp%Energy(iimp) + ii/1000.d0 * imp%Bandwidth(iimp)
            dist   = 1.d0/(imp%Bandwidth(iimp) * sqrt(2.d0 * pi)) * exp(-0.5d0 * ((eneii - imp%Energy(iimp))/imp%Bandwidth(iimp))**2)

            impelec = 1.q0 &
              / (1.q0 + imp%Degeneracy(iimp) * exp(info%betaQ*imp%Dopant(iimp)*(mu-eneii)))

            imphole = 1.q0 &
              / (1.q0 + imp%Degeneracy(iimp)**(-1.q0) * exp(info%betaQ*imp%Dopant(iimp)*(eneii-mu)))

            if (imphole > impelec) then
              sumelec = sumelec - impelec*imp%Dopant(iimp)*densii*dist
            else
              sumhole = sumhole - imphole*imp%Dopant(iimp)*densii*dist
            endif
          enddo
        endif
      endif

    enddo
  endif

  deviation = deviation + sumelec - sumhole
  return

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
