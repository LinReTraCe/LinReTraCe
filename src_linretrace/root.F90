module Mroot
  use Mmpi_org
  use Mparams
  use Mtypes
  use Mfermi
  use Mauxiliary
  use psi_fast, only: psi0_imag_dp, psi0_imag_qp

  interface ndeviation
    module procedure ndeviation_D, ndeviation_Q
  end interface ndeviation

  interface occ_digamma
    module procedure occ_digamma_D, occ_digamma_Q
  end interface occ_digamma

  interface occ_fermi
    module procedure occ_fermi_D, occ_fermi_Q
  end interface occ_fermi

  interface occ_impurity
    module procedure occ_impurity_D, occ_impurity_Q
  end interface occ_impurity

  ! all occupation kernels of a given precision share one signature;
  ! we exploit this to resolve the (muFermi x T=0) choice ONCE per temperature
  ! step -- outside all momentum/band loops -- via procedure pointers.
  ! inside the leaf kernels the distribution functions remain direct calls to
  ! pure elemental module procedures, i.e. they stay inlinable/vectorizable.
  abstract interface
    subroutine occ_kernel_D(mu, occ_tot, edisp, sct, kmesh, algo, info)
      import :: energydisp, scattering, kpointmesh, algorithm, runinfo
      real(8), intent(in)  :: mu
      real(8), intent(out) :: occ_tot
      type(energydisp) :: edisp
      type(scattering) :: sct
      type(kpointmesh) :: kmesh
      type(algorithm)  :: algo
      type(runinfo)    :: info
    end subroutine occ_kernel_D

    subroutine occ_kernel_Q(mu, occ_tot, edisp, sct, kmesh, algo, info)
      import :: energydisp, scattering, kpointmesh, algorithm, runinfo
      real(16), intent(in)  :: mu
      real(16), intent(out) :: occ_tot
      type(energydisp) :: edisp
      type(scattering) :: sct
      type(kpointmesh) :: kmesh
      type(algorithm)  :: algo
      type(runinfo)    :: info
    end subroutine occ_kernel_Q
  end interface

  ! occupation kernels used by ndeviation_D/Q
  ! set via set_occ_kernels once per temperature step (main.F90)
  procedure(occ_kernel_D), pointer :: occ_ptr_D => null()
  procedure(occ_kernel_Q), pointer :: occ_ptr_Q => null()

  contains

! resolve which occupation kernels the chemical potential search uses:
! (fermi | digamma) x (finite T | T=0)
! to be called once per temperature step, after info%lT0 and the scattering
! quantities for the step have been set
subroutine set_occ_kernels(algo, info)
  implicit none
  type(algorithm) :: algo
  type(runinfo)   :: info

  if (info%lT0) then
    if (algo%muFermi) then
      occ_ptr_D => occ_fermi_T0_D
      occ_ptr_Q => occ_fermi_T0_Q
    else
      occ_ptr_D => occ_digamma_T0_D
      occ_ptr_Q => occ_digamma_T0_Q
    endif
  else
    if (algo%muFermi) then
      occ_ptr_D => occ_fermi_D
      occ_ptr_Q => occ_fermi_Q
    else
      occ_ptr_D => occ_digamma_D
      occ_ptr_Q => occ_digamma_Q
    endif
  endif
end subroutine set_occ_kernels

! if band shifts are applied onto the input DFT energy range
! we have to re-determine the chemical potential of the system
subroutine find_mu_DFT(edisp,kmesh,pot)
  implicit none
  type(energydisp) :: edisp
  type(kpointmesh) :: kmesh
  type(potential)  :: pot

  real(8) :: mu1,mu2,mu
  real(8) :: occ1,occ2,occ
  real(8) :: target_zero1, target_zero2, target_zero

  real(8) :: emin, emax
  integer :: is, iband

  ! root finding purely for a DFT band structure
  ! i.e. T = 0
  ! bisection
  mu1 = minval(edisp%band)
  mu2 = maxval(edisp%band)
#ifdef MPI
  call MPI_allreduce(MPI_IN_PLACE, mu1, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, mpierr)
  call MPI_allreduce(MPI_IN_PLACE, mu2, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, mpierr)
#endif
  mu1 = mu1-0.01
  mu2 = mu2+0.01

  call occ_DFT(mu1,occ1,edisp,kmesh)
  target_zero1 = edisp%nelect - occ1

  call occ_DFT(mu2,occ2,edisp,kmesh)
  target_zero2 = edisp%nelect - occ2

  ! Bisection root finding
  do
     mu = (mu1+mu2)/2.d0
     call occ_DFT(mu,occ,edisp,kmesh)
     target_zero = edisp%nelect - occ

     if (abs(target_zero) < kmesh%minimal_weight) exit ! if we hit the band gap
     if ((mu2-mu1)*10.d0 < kmesh%minimal_weight) exit ! if converged

     if (target_zero.gt.0.q0) then
        mu1=mu
        target_zero1=target_zero
     else
        mu2=mu
        target_zero2=target_zero
     endif
  enddo

  ! Detect if we are gapped
  edisp%gapped = .true.
  do is = 1,edisp%ispin
    do iband = 1,edisp%nband_max
      emin = minval(edisp%band(iband,:,is))
      emax = maxval(edisp%band(iband,:,is))
#ifdef MPI
      call MPI_allreduce(MPI_IN_PLACE, emin, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, mpierr)
      call MPI_allreduce(MPI_IN_PLACE, emax, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, mpierr)
#endif
      if ((mu > emin) .and. (mu < emax)) then ! we cut through a band -> not gapped
        edisp%gapped(is) = .false.
        exit
      endif
    enddo
  enddo

  ! Redefine gap sizes
  do is = 1,edisp%ispin
    if (edisp%gapped(is)) then
      do iband = 1, edisp%nband_max
        emin = minval(edisp%band(iband,:,is))
        emax = maxval(edisp%band(iband,:,is))
#ifdef MPI
        call MPI_allreduce(MPI_IN_PLACE, emin, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, mpierr)
        call MPI_allreduce(MPI_IN_PLACE, emax, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, mpierr)
#endif
        if (mu>emax) then
          edisp%ene_valenceBand(is) = emax
          edisp%valenceBand(is) = iband
        endif
        if (mu<emin) then
          edisp%ene_conductionBand(is) = emin
          edisp%conductionBand(is) = iband
          exit
        endif
      enddo

      if ((edisp%conductionBand(is) - edisp%valenceBand(is)) /= 1) then
        call stop_with_message(stderr, 'ERROR: something went wrong in find_mu_DFT')
      endif

      edisp%gap(is) = edisp%ene_conductionBand(is) - edisp%ene_valenceBand(is)
    else
      edisp%gap(is) = 0.d0
    endif
  enddo

  ! Adjust the flags
  edisp%gap_min = 0.d0
  edisp%gapped_complete = .false.
  if (edisp%ispin == 1) then
    if (edisp%gapped(1) .eqv. .true.) then
      edisp%gapped_complete = .true.
      edisp%gap_min = edisp%ene_conductionBand(1) - edisp%ene_valenceBand(1)
      mu = (edisp%ene_conductionBand(1) + edisp%ene_valenceBand(1)) / 2.d0
    endif
  else if (edisp%ispin == 2) then
    if ((edisp%gapped(1) .eqv. .true.) .and. (edisp%gapped(2) .eqv. .true.)) then
      edisp%gapped_complete = .true.
      edisp%gap_min = minval(edisp%ene_conductionBand) - maxval(edisp%ene_valenceBand)
      mu = (minval(edisp%ene_conductionBand) + maxval(edisp%ene_valenceBand)) / 2.d0
    endif
  endif

  pot%mu_dft = mu ! finally, save it into the data structure

end subroutine find_mu_DFT

! old double precision routine
! this is depecrecated since we always ensure that the chemical potential is found in quad precision
! the result is then always transformed back into double precision
! commented out: march 21st 2022 -- due to lapack dependencies

!subroutine find_mu_D(mu,dev,target_zero,niitact, edisp, sct, kmesh, imp, algo, info)
!  implicit none

!  real(8), intent(inout)        :: mu ! chemical potential which is calculated
!  real(8), intent(in)           :: dev    ! allowed deviation
!  real(8), intent(out)          :: target_zero ! deviation from root after convergence
!  integer, intent(out)          :: niitact ! number of iterations

!  type(algorithm)  :: algo
!  type(energydisp) :: edisp
!  type(kpointmesh) :: kmesh
!  type(scattering) :: sct
!  type(impurity)   :: imp
!  type(runinfo)    :: info

!  ! local variables
!  real(8)  target_zero1, target_zero2, mu1, mu2
!  integer iit,niit0
!  ! Secand method
!  logical lsecant
!  ! linear interpolation method
!  real(8), allocatable :: Y(:), X(:) !arrays containig the function to minimise and the chemical potential
!  integer :: nmu  ! number of points that sample the mu interval (mu1,mu2)
!  real(8) :: dmu ! increment
!  real(8) :: a11, a22, a31, a42
!  real(8) :: A(4,4), B(4)
!  integer :: i, j
!  integer :: ipiv(4)
!  integer :: ierr
!  ! Linear interpolation
!  logical linint
!  ! Ridders' method
!  real(8) :: F(4), P(4)
!  real(8) :: s
!  integer  :: maxiter ! maximum number of iterations
!  logical  :: lridd   ! selects Ridders' method
!  ! Bisection method
!  logical :: lbisec

!  lsecant = .false.
!  linint  = .false.
!  lridd   = .false.
!  lbisec  = .false.
!  ! choose method according to input
!  ! if method is not provided: default to Riddler
! select case (algo%rootMethod)
!   case (0)
!     lsecant = .true.
!   case (1)
!     linint  = .true.
!   case (2)
!     lridd   = .true.
!   case (3)
!     lbisec  = .true.
!   case default
!     call stop_with_message(stderr, "Root-finding method not properly defined")
! end select

!! deviation from set particle number with initial mu
!! output: target_zero1 = required - actual electrons
!  call ndeviation(mu, target_zero1, edisp, sct, kmesh, imp, algo, info)

!! coarse initialization of secant bracket mu1, mu2
!  target_zero2=target_zero1
!  mu1=mu
!  mu2=mu

!  do while (target_zero2.gt.0.d0) ! too few electrons -> increase mu
!     mu2 = mu2 + 0.5d0
!     call ndeviation(mu2, target_zero2, edisp, sct, kmesh, imp, algo, info)
!  enddo
!  do while (target_zero1.le.0.d0) ! too many electrons -> decrease mu
!     mu1 = mu1 - 0.5d0
!     call ndeviation(mu1, target_zero1, edisp, sct, kmesh, imp, algo, info)
!  enddo

!  ! maximal steps for double precision calculations
!  niit0=niit

!  if (lsecant) then
!  !Secant root finding
!    do iit=1,niit0
!       mu=mu1-target_zero1*(mu2-mu1)/(target_zero2-target_zero1)
!       call ndeviation(mu, target_zero, edisp, sct, kmesh, imp, algo, info)

!       if (abs(target_zero).lt.dev)  exit
!       if (target_zero.gt.0.d0) then
!          mu1=mu
!          target_zero1=target_zero
!       else
!          mu2=mu
!          target_zero2=target_zero
!       endif
!    enddo
!    niitact=iit

!  elseif (linint) then
!    ! evaluate the target function on an interval and find the root by linear interpolation
!    ! fix the number of points to sample the (mu1,mu2)
!    nmu=30
!    allocate(Y(nmu), X(nmu))
!    Y(:)=0.0d0 ; X(:)=0.0d0
!    ! construct linear grid
!    X(1)=mu1; X(nmu)=mu2; dmu=(mu2-mu1)/(nmu-1)
!    !write(*,*) 'mu1',mu1,'mu2',mu2,'dmu',dmu
!    Y(1)=target_zero1
!    Y(nmu)=target_zero2

!    do i=2,nmu-1
!      X(i)=X(i-1)+dmu
!    enddo
!    ! evaluate target function in the interval
!    do i=2,nmu-1
!       call ndeviation(X(i), Y(i), edisp, sct, kmesh, imp, algo, info)
!    enddo
!    do i=1,nmu
!      Y(i)=Y(i)+X(i) !this is the correct target function for this method
!    enddo

!    ! find root by linear interpolation
!    do i = 1, nmu-1
!       do j = 1, nmu-1
!          A(:,:) = 0.0d0
!          a11 = X(i+1)-X(i)
!          a22 = X(j+1)-X(j)
!          a31 = Y(i+1)-Y(i)
!          a42 = X(j+1)-X(j)

!          A(1,1)=a11; A(2,2)=a22; A(3,1)=a31; A(4,2)=a42
!          A(1,3)=-1.0d0; A(2,3)=-1.0d0
!          A(3,4)=-1.0d0; A(4,4)=-1.0d0
!          B(1) = -X(i); B(2) = -X(j)
!          B(3) = -Y(i); B(4) = -X(j)

!          !write(*,*) 'LU factorisation begins'
!          call dgetrf(4, 4, A, 4, ipiv, ierr )
!          if (ierr /= 0) write(*,*) 'LU factorisation failed', ierr, i, j, a31

!          !write(*,*) 'solution lin syst begins'
!          call dgetrs( 'N', 4, 1, A, 4, ipiv, B, 4, ierr)
!          if (ierr /= 0) write(*,*) 'solution of the system of linear equations has failed', ierr, i, j

!          ! check if there is any intersection
!          if (B(1) < 1.0d0 .and. B(2) < 1.0d0) then
!             if (B(1) >= 0.0d0 .and. B(2) >= 0.0d0) then
!                !write(*,*) b(3), b(4)
!                ! save the values of the intersection
!                mu = B(3)
!                call ndeviation(mu, target_zero, edisp, sct, kmesh, imp, algo, info)
!             endif
!          endif
!       enddo ! over freq. counter j
!    enddo ! over freq. counter i
!    deallocate(Y,X)

!  elseif(lridd) then   !Ridders' method for root finding
!    P(1)=mu1 ; P(2)=mu2
!    F(1)= target_zero1; F(2)= target_zero2
!    maxiter=  niit

!     do j = 1, maxiter
!        P(3) = 0.5d0*(P(1)+P(2))
!        call ndeviation(P(3), F(3), edisp, sct, kmesh, imp, algo, info)
!        s = dsqrt((F(3)**2)-(F(1)*F(2)))
!        if (s==0.0d0) then
!           write(*,*) 'Error in Ridders search for chemical potential'
!           write(*,*) 'ITER', j, 'x1', P(1),'  x2',P(2),'  x3', P(3)
!           write(*,*) 'ITER', j, 'F1', F(1),'  F2',F(2),'  F3', F(3)
!           goto 400
!        endif
!        P(4) = P(3)+(P(3)-P(1))*(SIGN(1.0d0,F(1)-F(2))*F(3)/s)
!        call ndeviation(P(4), F(4), edisp, sct, kmesh, imp, algo, info)
!        if (abs(F(4)) .lt. dev) goto 400
!        if (sign(F(3), F(4)) /= F(3)) then
!        !change of sign btw x3 and x4 then reduce search interval
!           P(1)  = P(3)
!           F(1)  = F(3)
!           P(2)  = P(4)
!           F(2)  = F(4)
!        else if (sign(F(1), F(4)) /= F(1)) then
!        !change of sign btw x1 and x4 then reduce search interval
!           P(2)  = P(4)
!           F(2)  = F(4)
!        else if (sign(F(2), F(4)) /= F(2)) then
!        !change of sign btw x2 and x4 then reduce search interval
!           P(1)  = P(4)
!           F(1)  = F(4)
!        endif
!     enddo ! over number of iterations

! 400 if (j == maxiter) write(*,*) 'Ridders seach might not have converged'

!     ! save the values of the intersection
!     mu = P(4)
!     niitact = j
!     niit0   = maxiter
!     target_zero = F(4)

!  elseif(lbisec) then
!    ! Bisection root finding
!    do iit=1,niit0
!       mu = (mu1+mu2)/2.d0
!       call ndeviation(mu, target_zero, edisp, sct, kmesh, imp, algo, info)
!       if (myid.eq.master .and. iit .ge. 50) write(*,*) mu
!       if (abs(target_zero).lt.dev) exit
!       if (target_zero.gt.0.q0) then
!          mu1=mu
!          target_zero1=target_zero
!       else
!          mu2=mu
!          target_zero2=target_zero
!       endif
!    enddo
!    niitact = iit
!  endif ! root finding algorithm



!  ! if (myid.eq.master .and. (niitact .ge. niit0 .or. abs(target_zero) .ge. dev)) then
!  !   write(*,'(A,1E20.12)') "WARNING: diminished root precision. ndev_actual =",target_zero
!  !   write(*,'(A,1F10.3,A,1I5,A,1E20.12)') "at T=",sct%TT(iT), " with  niit=",niit0, " ndev =", dev
!  !   write(*,*) "increase niit, or allow for bigger ndev (see params.F90)"
!  ! endif
!end subroutine find_mu_D



! determine the chemical potentail with the preset algorithm
! if the system is gapped and it is determined that the standard accuracy is not good enough
! apply the refinement algorithm after the fact.
subroutine find_mu(mu,dev,target_zero,niitact, edisp, sct, kmesh, imp, algo, info)
  implicit none
  ! passed variables
  real(16), intent(inout)       :: mu
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
  logical  :: refine_abort
  real(16) :: mu_refine
  real(16) target_zero1, target_zero2
  real(16) test_up, test_dn
  real(16) mu1, mu2, dmu
  integer iit, niit0
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
  real(8)  :: sctmin
  real(8)  :: gapdistance, distance_cb, distance_vb
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

  ! at T=0 with Heaviside statistics the deviation is a staircase in mu:
  ! secant/Ridders operate on plateaus (equal function values -> ill-conditioned
  ! updates), while bisection converges the bracket robustly onto the step.
  ! enforce bisection for this combination.
  if (info%lT0 .and. algo%muFermi) then
    lsecant = .false.
    linint  = .false.
    lridd   = .false.
    lbisec  = .true.
  endif

! deviation from set particle number with initial mu
  call ndeviation(mu, target_zero1, edisp, sct, kmesh, imp, algo, info)

  target_zero2=target_zero1
  mu1=mu
  mu2=mu

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
       mu=mu1-target_zero1*mu2/(target_zero2-target_zero1)+target_zero1*mu1/(target_zero2-target_zero1)
       call ndeviation(mu, target_zero, edisp, sct, kmesh, imp, algo, info)

       if (abs(target_zero).lt.dev) exit
       if (target_zero.gt.0.q0) then
          mu1=mu
          target_zero1=target_zero
       else
          mu2=mu
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
     mu = P(4)
     niitact = j
     niit0   = maxiter
     target_zero = F(4)

  elseif(lbisec) then
    ! Bisection root finding
    do iit=1,niit0
       mu = (mu1+mu2)/2.q0
       call ndeviation(mu, target_zero, edisp, sct, kmesh, imp, algo, info)

       if (abs(target_zero).lt.dev) exit
       ! T=0 Heaviside statistics: on a discrete k-grid the demanded electron
       ! number may fall between two occupation plateaus (metals), so
       ! |target_zero| can never drop below dev. the bracket nevertheless
       ! collapses geometrically onto the step, i.e. onto the Fermi level ->
       ! accept once the interval is converged. (finite-T behavior untouched)
       if (info%lT0 .and. algo%muFermi .and. abs(mu2-mu1) < 1q-14) exit
       if (target_zero.gt.0.q0) then
          mu1=mu
          target_zero1=target_zero
       else
          mu2=mu
          target_zero2=target_zero
       endif
    enddo
    niitact = iit
  endif

  ! T=0 with pure Fermi statistics (Heaviside occupation):
  ! the carrier-balance refinement below is degenerate here -- inside a gap the
  ! activated electron and hole counts are BOTH identically zero for ANY mu in
  ! the gap, so occ_refine carries no information and find_mu_refine's
  ! bracketing would not terminate meaningfully. never enter it.
  !
  ! instead, in the pure fully gapped case we place mu at the center of the gap,
  ! which is the exact T->0 limit of Fermi statistics (the root finder above
  ! merely returns SOME point of the zero-deviation plateau, dependent on
  ! bracket history). the band edges used here are maintained by find_mu_DFT,
  ! which main.F90 re-evaluates every step on the CURRENT dispersion whenever
  ! band shifts are present (edisp%band = band_file + band_shift) -- static
  ! self-energy renormalizations from a scattering file are therefore honored.
  ! the quasi-particle weight Z > 0 cannot move the T=0 Fermi mu at all, since
  ! it drops out of sign(Z*(eps-mu)) identically (band- and k-dependent or not).
  !
  ! with impurities or doping the root-finder result above stands: mu is then
  ! pinned by the (staircase) charge balance, and mid-gap centering would be
  ! wrong. note that at strict T=0 the position of mu WITHIN a zero-deviation
  ! plateau of the impurity balance is physically undetermined; if the precise
  ! T->0 pinning is required, evaluate at a small finite temperature instead.
  if (info%lT0 .and. algo%muFermi) then
    if (edisp%gapped_complete .and. &
        .not. (algo%lTMODE .and. (algo%lImpurities .or. algo%lDoping))) then
      if (edisp%ispin == 1) then
        mu = (real(edisp%ene_conductionBand(1),16) + real(edisp%ene_valenceBand(1),16)) / 2.q0
      else
        mu = (real(minval(edisp%ene_conductionBand),16) + real(maxval(edisp%ene_valenceBand),16)) / 2.q0
      endif
      ! re-evaluate the deviation at the centered chemical potential
      ! (identically zero in the gap; keeps the reported target_zero consistent)
      call ndeviation(mu, target_zero, edisp, sct, kmesh, imp, algo, info)
    endif
    return
  endif

  if ( (algo%ldebug .and. (index(algo%dbgstr,"NoRefine") .ne. 0)) ) then
    return
  endif

  if ( .not. (algo%ldebug .and. (index(algo%dbgstr,"AlwaysRefine") .ne. 0)) ) then
    sctmin = minval(sct%gam)
#ifdef MPI
    call MPI_allreduce(MPI_IN_PLACE, sctmin, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, mpierr)
#endif
    gapdistance = 0.d0 ! smalles distance to either valence or conduction band
    if (edisp%gapped_complete) then
      distance_cb = minval(edisp%ene_conductionBand) - mu
      if (distance_cb < 0.d0) then
        return ! we are outside gap
      endif
      gapdistance = distance_cb
      distance_vb = mu - maxval(edisp%ene_valenceBand)
      if (distance_vb < 0.d0) then
        return ! we are outside gap
      endif
      if (distance_cb < gapdistance) then
        gapdistance = distance_cb
      endif
    else
      return ! not gapped
    endif
  endif

  ! hard temperature cutoff relative to band gap size for fermi function
  ! scattering rate and temperature cutoff for digamma function
  ! these values here have been determined empirically
  if ( (algo%ldebug .and. (index(algo%dbgstr,"AlwaysRefine") .ne. 0)) .or. &
       (algo%muFermi .and. (info%Temp < edisp%gap_min*300)) .or. &
       (.not. algo%muFermi .and. (sctmin < 3d-10 * edisp%gap_min) .and. (info%Temp < gapdistance*800)) ) then

    call occ_refine(mu, target_zero1, edisp, sct, kmesh, imp, algo, info)
    ! do not enter refinement algorithm if this value is too large
    ! this is for security puropeses, one can somehow end up here by Debug flags
    ! if one enters the refinement routine with larger values, a run-off to another root might happen -> avoid at all costs
    ! so 'AlwaysRefine' is more like always refine wherever safely possible
    if (algo%muFermi .and. abs(target_zero1) > 1d-18) then
      return
    endif
    if (.not. algo%muFermi .and. abs(target_zero1) > 1d-16) then
      return
    endif
    ! refinement
    call find_mu_refine(mu, mu_refine, refine_abort, niitact, edisp, sct, kmesh, imp, algo, info)
    if (.not. refine_abort) then
      ! if the algorithm is not aborted due to numerical problems at the refinmenet stage
      mu = mu_refine
      return
    endif
    niitact = niitQ+1 ! notify main ! this number triggers the dmu/dT calculation in main.F90
  endif

end subroutine find_mu

! refinement algorithm of the chemical potential
! essentially balances the thermally activated electrons / holes in the system
! if the chemical potential is inside a fully gapped system
subroutine find_mu_refine(mu_in, mu_out, refine_abort, niitact, edisp, sct, kmesh, imp, algo, info)
  implicit none
  ! passed variables
  real(16), intent(in)   :: mu_in
  real(16), intent(out)  :: mu_out
  logical, intent(out)   :: refine_abort
  integer, intent(inout) :: niitact

  type(energydisp) :: edisp
  type(scattering) :: sct
  type(kpointmesh) :: kmesh
  type(impurity)   :: imp
  type(algorithm)  :: algo
  type(runinfo)    :: info

  ! local variables
  real(16) :: target_zero1, target_zero2, target_zero
  real(16) :: mu1, mu2, dmu
  real(16) :: test_up, test_dn
  integer  :: iit
  integer  :: is
  integer  :: cnt


  refine_abort = .false.
  ! flag that is only triggered if we ran into numerical problems at the very end
  ! i.e. we are limited by quad precision evaluation or ran or ran into some false root
  ! this flag is only used there because it we need to set the chemical potential to some values
  ! i.e. via dmu/dT

  ! calculate numeric deviations from 0 for the two methods
  call ndeviation_Q(mu_in, target_zero2, edisp, sct, kmesh, imp, algo, info)
  call occ_refine(mu_in, target_zero1, edisp, sct, kmesh, imp, algo, info)

  ! create a small stepping size in relation to the gap size
  dmu = edisp%gap_min/100.q0
  mu1 = mu_in
  mu2 = mu_in
  iit = 1

  ! determine direction in which we have to move
  if (target_zero1 < 0.q0) then
    dmu = +dmu
  else
    dmu = -dmu
  endif

  ! decrease the step size if necessary
  do cnt = 1,10
    call occ_refine(mu_in+dmu, test_up, edisp, sct, kmesh, imp, algo, info)
    call occ_refine(mu_in-dmu, test_dn, edisp, sct, kmesh, imp, algo, info)
    test_up = test_up - target_zero1
    test_dn = test_dn - target_zero1
    ! check if both point in the same direction
    if ( (test_up > 0 .and. test_dn > 0) .or. (test_up < 0 .and. test_dn < 0) ) then
      dmu = dmu/2.q0
    else
      exit
    endif
  enddo

  ! abort if this for some reason fails
  if (cnt >= 10) then
    mu_out = mu_in
    return
  endif

  ! get the working interval for the bisection
  target_zero2 = target_zero1 ! from the top
  do while (((target_zero1 <= 0.q0 .and. target_zero2 <= 0.q0) .or. &
            (target_zero1 >= 0.q0  .and. target_zero2 >= 0.q0)) .and. iit < niitQ)
    mu2 = mu2 + dmu
    call occ_refine(mu2, target_zero2, edisp, sct, kmesh, imp, algo, info)
    iit = iit + 1
  enddo
  ! avoid adding this to the count
  ! this some times becomes large if the initial chemical potential is far away from the solution


  ! minimize the working interval
  ! i.e. interval that _just_ contains the sign change
  mu1 = mu2 - dmu

  ! finally perform a bisection in the previous working interval
  do iit = 1,niitQ
    mu_out = (mu1+mu2)/2.q0
    call occ_refine(mu_out, target_zero, edisp, sct, kmesh, imp, algo, info)
    niitact = niitact + 1

    if ( abs(mu1-mu2) < 1q-15) then
      exit ! we go all out here
    endif

    if ((target_zero .gt. 0.q0 .and. target_zero2.gt. 0.q0) &
         .or. (target_zero .lt. 0.q0 .and. target_zero2 .lt. 0.q0)) then
      mu2=mu_out
      target_zero2=target_zero ! here we are on the bisection side of mu2
                                ! we set the middle point as new mu2
    else
      mu1=mu_out                ! the other way aorund here
      target_zero1=target_zero
    endif
  enddo

  ! too many interations
  if (iit >= niitQ) then
    mu_out = mu_in
    return
  endif

  ! numeric boundary .. undefined ~ quad precision
  ! smallest quad number = 2^{-16494} ~ 10^{-4965}
  if ((abs(target_zero) < 1q-4900) .or. &
      (abs(target_zero1) < 1q-4900) .or. &
      (abs(target_zero2) < 1q-4900)) then
    refine_abort = .true.
    return
  endif

  ! we ran into some direction we dont want
  ! note: we enter this routine with smaller values than this
  if ((abs(target_zero) > 1q-15) .or. &
      (abs(target_zero1) > 1q-15) .or. &
      (abs(target_zero2) > 1q-15)) then
    refine_abort = .true.
    return
  endif

  ! check if the mu is inside a band -- do not use this value
  do is=1,edisp%ispin
    if (mu_out < edisp%ene_valenceBand(is) .or. mu_out > edisp%ene_conductionBand(is)) then
      mu_out = mu_in
    endif
  enddo

  ! if(myid.eq.master) then
  !   write(*,*) target_zero1, target_zero, target_zero2
  ! endif

  ! the final destination (:
  niitact = niitact + iit

end subroutine


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



! determine the deviation of the calculated total occupation from the demanded occupation
! returns the difference between the two numbers
! -- double precision
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
  real(8) :: occimp

  ! kernel choice (fermi | digamma) x (finite T | T=0) resolved once per
  ! temperature step in set_occ_kernels -- no branching in the hot loops
  if (.not. associated(occ_ptr_D)) then
    call stop_with_message(stderr, 'ndeviation_D: occupation kernel not set (set_occ_kernels)')
  endif
  call occ_ptr_D(mu, occ_tot, edisp, sct, kmesh, algo, info)

  if (algo%lTMODE .and. algo%lImpurities) then
    call occ_impurity(occimp, mu, imp, info)
    occ_tot = occ_tot - occimp
  endif

  if (algo%lTMODE .and. algo%ldoping) then
    occ_tot = occ_tot - edisp%doping
  endif

  target_zero = edisp%nelect - occ_tot
end subroutine ndeviation_D

! determine the deviation of the calculated total occupation from the demanded occupation
! returns the difference between the two numbers
! -- quad precision
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
  real(16) :: occimp

  ! kernel choice (fermi | digamma) x (finite T | T=0) resolved once per
  ! temperature step in set_occ_kernels -- no branching in the hot loops
  if (.not. associated(occ_ptr_Q)) then
    call stop_with_message(stderr, 'ndeviation_Q: occupation kernel not set (set_occ_kernels)')
  endif
  call occ_ptr_Q(mu, occ_tot, edisp, sct, kmesh, algo, info)

  if (algo%lTMODE .and. algo%lImpurities) then
    call occ_impurity(occimp, mu, imp, info)
    occ_tot = occ_tot - occimp
  endif

  if (algo%lTMODE .and. algo%ldoping) then
    occ_tot = occ_tot - edisp%doping
  endif

  target_zero=real(edisp%nelect,16) - occ_tot
end subroutine ndeviation_Q

! determine the occupation of the full system via the digamma function
! -- double precision
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
  ! --- deprecated CERNLIB path -------------------------------------------
  ! complex(8), external :: wpsipg

  occ_loc = 0.d0
  ! evaluate the function
  do is = 1,edisp%ispin
    do ik = ikstr, ikend
      do iband=1,edisp%nband_max
        ! psi0_imag_dp (module psi_fast)
        !
        !   x = psi0_imag_dp(z [, ierr])
        !     z    complex(8), argument
        !     x    real(8), Im psi^(0)(z)
        !     ierr optional integer, PSI_ERR_OK / _ORDER / _POLE
        !
        !   This is the most frequently executed floating-point kernel in the
        !   program: nband_max x nkp evaluations per root-finder iteration.
        !   Only psi_0 is needed here -- the higher orders are computed once,
        !   after mu has converged, by calc_polygamma.
        occ_loc = occ_loc + (0.5d0 + &
        psi0_imag_dp(0.5d0 + info%beta2p * &
             (sct%gam(iband,ik,is) - ci*sct%zqp(iband,ik,is)*(edisp%band(iband,ik,is) - mu)))/pi) &
        * kmesh%weight(ik)
        ! --- deprecated CERNLIB path -----------------------------------
        ! occ_loc = occ_loc + (0.5d0 + &
        ! aimag(wpsipg(0.5d0 + info%beta2p * &
        !      (sct%gam(iband,ik,is) - ci*sct%zqp(iband,ik,is)*(edisp%band(iband,ik,is) - mu)),0))/pi) &
        ! * kmesh%weight(ik)
      enddo
    enddo
  enddo

#ifdef MPI
  call MPI_ALLREDUCE(occ_loc, occ_tot, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, mpierr)
#else
  occ_tot = occ_loc
#endif

end subroutine occ_digamma_D


! determine the occupation of the full system via the digamma function
! -- quad precision
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
  ! --- deprecated CERNLIB path -------------------------------------------
  ! complex(16), external :: wpsipghp

  occ_loc = 0.q0
  ! evaluate the function
  do is = 1,edisp%ispin
    do ik = ikstr, ikend
      do iband=1,edisp%nband_max
        ! psi0_imag_qp(z [, ierr]) -> real(16) Im psi^(0)(z); see occ_digamma_D
        occ_loc = occ_loc + (0.5q0 + &
        psi0_imag_qp(0.5q0 + info%beta2pQ * &
             (sct%gam(iband,ik,is) - ciQ*sct%zqp(iband,ik,is)*(edisp%band(iband,ik,is) - mu)))/piQ) &
        * kmesh%weightQ(ik)
        ! --- deprecated CERNLIB path -----------------------------------
        ! occ_loc = occ_loc + (0.5q0 + &
        ! aimag(wpsipghp(0.5q0 + info%beta2pQ * &
        !      (sct%gam(iband,ik,is) - ciQ*sct%zqp(iband,ik,is)*(edisp%band(iband,ik,is) - mu)),0))/piQ) &
        ! * kmesh%weightQ(ik)
      enddo
    enddo
  enddo

#ifdef MPI
  call MPI_reduce_quad(occ_loc, occ_tot)
#else
  occ_tot = occ_loc
#endif

end subroutine occ_digamma_Q

! determine the occupation of the full system via the fermi function
! -- double precision
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

  occ_loc = 0.d0
  ! evaluate the function
  do is = 1,edisp%ispin
    do ik = ikstr, ikend
      do iband=1,edisp%nband_max
        occ_loc = occ_loc + &
        fermi(sct%zqp(iband,ik,is)*(edisp%band(iband,ik,is)-mu), info%beta) * kmesh%weight(ik)
      enddo
    enddo
  enddo

#ifdef MPI
  call MPI_ALLREDUCE(occ_loc, occ_tot, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, mpierr)
#else
  occ_tot = occ_loc
#endif

end subroutine occ_fermi_D


! determine the occupation of the full system via the fermi function
! -- quad precision
subroutine occ_fermi_Q(mu, occ_tot, edisp, sct, kmesh, algo, info)
  implicit none

  real(16), intent(in)  :: mu
  real(16), intent(out) :: occ_tot

  type(energydisp) :: edisp
  type(scattering) :: sct
  type(kpointmesh) :: kmesh
  type(runinfo)    :: info
  type(algorithm)  :: algo

  real(16) :: occ_loc
  integer :: is, ik, iband

  occ_loc = 0.q0
  do is = 1,edisp%ispin
    do ik = ikstr, ikend
      do iband=1,edisp%nband_max
        occ_loc = occ_loc + &
        fermi_qp(sct%zqp(iband,ik,is)*(edisp%band(iband,ik,is)-mu), info%betaQ) * kmesh%weightQ(ik)
      enddo
    enddo
  enddo

#ifdef MPI
  call mpi_reduce_quad(occ_loc, occ_tot) ! custom quad reduction
#else
  occ_tot = occ_loc
#endif

end subroutine occ_fermi_Q

!________________________________________________________
! T=0 occupation kernels
! these mirror the finite-temperature loop nests above one-to-one;
! only the distribution function is replaced by its beta -> infinity limit
! (see Mfermi for the derivations). they contain no beta and no exp/wpsipg.
!
! note on renormalizations:
! - the kernels act on edisp%band, which contains the band shifts of the
!   current step (edisp%band = band_file + band_shift, input.f90), so static
!   real-part renormalizations from a scattering file are fully honored
! - the quasi-particle weight Z > 0 drops out of the Heaviside identically
!   (sign(Z*(eps-mu)) = sign(eps-mu)) but remains in the digamma/arctan kernel

! T=0 occupation via the fermi function (Heaviside)
! -- double precision
subroutine occ_fermi_T0_D(mu, occ_tot, edisp, sct, kmesh, algo, info)
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

  occ_loc = 0.d0
  ! evaluate the function
  do is = 1,edisp%ispin
    do ik = ikstr, ikend
      do iband=1,edisp%nband_max
        occ_loc = occ_loc + &
        fermi_T0(sct%zqp(iband,ik,is)*(edisp%band(iband,ik,is)-mu)) * kmesh%weight(ik)
      enddo
    enddo
  enddo

#ifdef MPI
  call MPI_ALLREDUCE(occ_loc, occ_tot, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, mpierr)
#else
  occ_tot = occ_loc
#endif

end subroutine occ_fermi_T0_D

! T=0 occupation via the fermi function (Heaviside)
! -- quad precision
subroutine occ_fermi_T0_Q(mu, occ_tot, edisp, sct, kmesh, algo, info)
  implicit none

  real(16), intent(in)  :: mu
  real(16), intent(out) :: occ_tot

  type(energydisp) :: edisp
  type(scattering) :: sct
  type(kpointmesh) :: kmesh
  type(algorithm)  :: algo
  type(runinfo)    :: info

  real(16) :: occ_loc
  integer :: is, ik, iband

  occ_loc = 0.q0
  do is = 1,edisp%ispin
    do ik = ikstr, ikend
      do iband=1,edisp%nband_max
        occ_loc = occ_loc + &
        fermi_T0_qp(sct%zqp(iband,ik,is)*(edisp%band(iband,ik,is)-mu)) * kmesh%weightQ(ik)
      enddo
    enddo
  enddo

#ifdef MPI
  call mpi_reduce_quad(occ_loc, occ_tot) ! custom quad reduction
#else
  occ_tot = occ_loc
#endif

end subroutine occ_fermi_T0_Q

! T=0 occupation via the digamma function limit
! n = 0.5 - 1/pi * atan2(Z*(eps-mu), Gamma)
! i.e. the Lorentzian-broadened occupation; recovers the Heaviside for Gamma -> 0
! -- double precision
subroutine occ_digamma_T0_D(mu, occ_tot, edisp, sct, kmesh, algo, info)
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

  occ_loc = 0.d0
  ! evaluate the function
  do is = 1,edisp%ispin
    do ik = ikstr, ikend
      do iband=1,edisp%nband_max
        occ_loc = occ_loc + (0.5d0 - &
        atankern_T0(sct%gam(iband,ik,is), sct%zqp(iband,ik,is)*(edisp%band(iband,ik,is) - mu))) &
        * kmesh%weight(ik)
      enddo
    enddo
  enddo

#ifdef MPI
  call MPI_ALLREDUCE(occ_loc, occ_tot, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, mpierr)
#else
  occ_tot = occ_loc
#endif

end subroutine occ_digamma_T0_D

! T=0 occupation via the digamma function limit
! -- quad precision
subroutine occ_digamma_T0_Q(mu, occ_tot, edisp, sct, kmesh, algo, info)
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

  occ_loc = 0.q0
  ! evaluate the function
  do is = 1,edisp%ispin
    do ik = ikstr, ikend
      do iband=1,edisp%nband_max
        occ_loc = occ_loc + (0.5q0 - &
        atankern_T0(sct%gam(iband,ik,is), sct%zqp(iband,ik,is)*(edisp%band(iband,ik,is) - mu))) &
        * kmesh%weightQ(ik)
      enddo
    enddo
  enddo

#ifdef MPI
  call MPI_reduce_quad(occ_loc, occ_tot)
#else
  occ_tot = occ_loc
#endif

end subroutine occ_digamma_T0_Q

! determine the difference between the thermally activated electrons and holes in the system
! -- quad precision
subroutine occ_refine(mu, deviation, edisp, sct, kmesh, imp, algo, info)
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
  real(16) :: psikern

  real(16) :: densii, eneii, dist
  real(16) :: occimp

  ! --- deprecated CERNLIB path -------------------------------------------
  ! complex(16), external :: wpsipghp

  sumelec = 0.q0
  sumhole = 0.q0

  ! do this if statement outside
  ! code speed > code duplication
  ! T=0 note: the fermi branch is degenerate at T=0 (elec/hole identically 0
  ! inside a gap); find_mu therefore never enters the refinement for
  ! (muFermi .and. lT0). the T=0 Heaviside branch below exists only for
  ! completeness / debug paths (AlwaysRefine).
  if (info%lT0 .and. algo%muFermi) then
    do is = 1,edisp%ispin
      do ik = ikstr, ikend
        do iband=1,edisp%nband_max
          elec = fermi_T0_qp(sct%zqp(iband,ik,is)*(edisp%band(iband,ik,is)-mu))
          hole = omfermi_T0_qp(sct%zqp(iband,ik,is)*(edisp%band(iband,ik,is)-mu))

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
  else if (info%lT0) then
    do is = 1,edisp%ispin
      do ik = ikstr, ikend
        do iband=1,edisp%nband_max
          ! T=0 limit of the psi_0 kernel: balances the Gamma-activated carriers
          psikern = atankern_T0(sct%gam(iband,ik,is), &
                                sct%zqp(iband,ik,is)*(edisp%band(iband,ik,is)-mu))
          elec = 0.5q0 - psikern
          hole = 0.5q0 + psikern
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
  else if (algo%muFermi) then
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
  else
    do is = 1,edisp%ispin
      do ik = ikstr, ikend
        do iband=1,edisp%nband_max
          ! only calculate the '0th' polygamma function
          ! psi0_imag_qp(z [, ierr]) -> real(16) Im psi^(0)(z); see occ_digamma_D
          psikern = 1.q0/piQ * psi0_imag_qp(0.5q0 + info%betaQ/2.q0/piQ &
                             * (sct%gam(iband,ik,is) + ciQ*sct%zqp(iband,ik,is)*(edisp%band(iband,ik,is)-mu)))
          ! --- deprecated CERNLIB path -------------------------------
          ! psikern = 1.q0/piQ * aimag(wpsipghp(0.5q0 + info%betaQ/2.q0/piQ &
          !                    * (sct%gam(iband,ik,is) + ciQ*sct%zqp(iband,ik,is)*(edisp%band(iband,ik,is)-mu)),0))
          elec = 0.5q0 - psikern
          hole = 0.5q0 + psikern
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
  endif

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

  ! nvalence = nsearch - N_D^+ + N_A^-
  ! N_D^+ = N_D/(1 + g * exp(beta * (mu - E_D)))
  ! N_A^+ = N_D/(1 + g * exp(-beta * (mu - E_A)))
  if (algo%lTMODE .and. algo%lImpurities) then
    call occ_impurity_Q(occimp, mu, imp, info)
    deviation = deviation - occimp
  endif

  if (algo%lTMODE .and. algo%ldoping) then
    deviation = deviation - edisp%doping
  endif

  return
end subroutine occ_refine

! determine the correction of the occupation, caused by the impurity states
! -- double precision
subroutine occ_impurity_D(occimp, mu, imp, info)
  real(8), intent(in)    :: mu
  real(8), intent(inout) :: occimp
  type(impurity)         :: imp
  type(runinfo)          :: info

  integer :: ii, iimp
  real(8) :: densii, eneii

  occimp = 0.d0

  ! impurity occupation
  ! = N_D^+ - N_A^-
  ! nvalence = nsearch - impurity occupation

  ! T=0 note: the exponentials below are replaced by their step-function limits;
  ! previously T=0 was handled implicitly via exp overflow of the huge-beta
  ! surrogate -- fragile under floating point exception trapping
  do iimp = 1,imp%nimp
    if (.not. imp%Band(iimp)) then
      if (info%lT0) then
        occimp = occimp + imp%Dopant(iimp)*imp%Density(iimp) &
          * gfermi_T0(imp%Dopant(iimp)*(mu-imp%Energy(iimp)), imp%Degeneracy(iimp))
      else
        occimp = occimp + imp%Dopant(iimp)*imp%Density(iimp) &
          / (1.d0 + imp%Degeneracy(iimp) * exp(info%beta*imp%Dopant(iimp)*(mu-imp%Energy(iimp))))
      endif
    else
      do ii=-500,500 ! this stays hard coded

        ! we segment the interval into 1001 levels
        eneii  = imp%Energy(iimp) + ii/1000.d0 * imp%Bandwidth(iimp)

        ! we weigh the point according to the given shape
        select case (imp%Bandtype(iimp))
          case (1) ! box
            densii = imp%Density(iimp) / 1001.d0
          case (2) ! triangle
            densii = imp%Density(iimp) * (500-abs(ii)) / 250.d0 &
                                         / 1000.d0
          case (3) ! half circle
            densii = imp%Density(iimp) * sqrt(1.d0 - (ii/500.d0)**2.d0) &
                                         / 785.371869250535197031882735354341d0
          case (4) ! cosine
            densii = imp%Density(iimp) * cos(ii/500.d0 * pi / 2.d0) &
                                         / 636.619248768719616210088408079917d0
          case (5) ! cosine^2
            densii = imp%Density(iimp) * cos(ii/500.d0 * pi / 2.d0)**2.d0 &
                                         / 500.d0
          case (6) ! cosine^3
            densii = imp%Density(iimp) * cos(ii/500.d0 * pi / 2.d0)**3.d0 &
                                         / 424.413181578904334542728254502383d0
          case (7) ! cosine^4
            densii = imp%Density(iimp) * cos(ii/500.d0 * pi / 2.d0)**4.d0 &
                                         / 375.d0
        end select
        ! the additional factors at the end were set such that the sum sum_ii densii = density

        ! and finally simply add the contribution
        if (info%lT0) then
          occimp = occimp + imp%Dopant(iimp)*densii &
            * gfermi_T0(imp%Dopant(iimp)*(mu-eneii), imp%Degeneracy(iimp))
        else
          occimp = occimp + imp%Dopant(iimp)*densii &
            / (1.d0 + imp%Degeneracy(iimp) * exp(info%beta*imp%Dopant(iimp)*(mu-eneii)))
        endif
      enddo
    endif
  enddo

end subroutine occ_impurity_D

! determine the correction of the occupation, caused by the impurity states
! -- quad precision
subroutine occ_impurity_Q(occimp, mu, imp, info)
  real(16), intent(in)    :: mu
  real(16), intent(inout) :: occimp
  type(impurity)          :: imp
  type(runinfo)           :: info

  integer :: ii, iimp
  real(16) :: densii, eneii

  occimp = 0.q0

  ! impurity occupation
  ! = N_D^+ - N_A^-
  ! nvalence = nsearch - impurity occupation

  do iimp = 1,imp%nimp
    if (.not. imp%Band(iimp)) then
      if (info%lT0) then
        occimp = occimp + imp%Dopant(iimp)*imp%Density(iimp) &
          * gfermi_T0(real(imp%Dopant(iimp),16)*(mu-imp%Energy(iimp)), imp%Degeneracy(iimp))
      else
        occimp = occimp + imp%Dopant(iimp)*imp%Density(iimp) &
          / (1.q0 + imp%Degeneracy(iimp) * exp(info%betaQ*imp%Dopant(iimp)*(mu-imp%Energy(iimp))))
      endif
    else
      do ii=-500,500 ! this stays hard coded

        ! we segment the interval into 1001 levels
        eneii  = imp%Energy(iimp) + ii/1000.q0 * imp%Bandwidth(iimp)

        ! we weigh the point according to the given shape
        select case (imp%Bandtype(iimp))
          case (1) ! box
            densii = imp%Density(iimp) / 1001.q0
          case (2) ! triangle
            densii = imp%Density(iimp) * (500-abs(ii)) / 250.q0 &
                                         / 1000.q0
          case (3) ! half circle
            densii = imp%Density(iimp) * sqrt(1.q0 - (ii/500.q0)**2.d0) &
                                         / 785.371869250535197031882735354341q0
          case (4) ! cosine .. we call it sine in the config
            densii = imp%Density(iimp) * cos(ii/500.q0 * pi / 2.q0) &
                                         / 636.619248768719616210088408079917q0
          case (5) ! cosine^2
            densii = imp%Density(iimp) * cos(ii/500.q0 * pi / 2.q0)**2.q0 &
                                         / 500.q0
          case (6) ! cosine^3
            densii = imp%Density(iimp) * cos(ii/500.q0 * pi / 2.q0)**3.q0 &
                                         / 424.413181578904334542728254502383q0
          case (7) ! cosine^4
            densii = imp%Density(iimp) * cos(ii/500.q0 * pi / 2.q0)**4.q0 &
                                         / 375.q0
        end select
        ! the additional factors at the end were set such that the sum sum_ii densii = density

        ! and finally simply add the contribution
        if (info%lT0) then
          occimp = occimp + imp%Dopant(iimp)*densii &
            * gfermi_T0(real(imp%Dopant(iimp),16)*(mu-eneii), imp%Degeneracy(iimp))
        else
          occimp = occimp + imp%Dopant(iimp)*densii &
            / (1.q0 + imp%Degeneracy(iimp) * exp(info%betaQ*imp%Dopant(iimp)*(mu-eneii)))
        endif
      enddo
    endif
  enddo

end subroutine occ_impurity_Q

! determine the occupation of the system for DFT (T=0)
! i.e. n(eps < mu) = 1; n(eps > mu) 0
subroutine occ_DFT(mu, occ_tot, edisp, kmesh)
  implicit none

  real(8), intent(in)  :: mu
  real(8), intent(out) :: occ_tot

  type(energydisp) :: edisp
  type(kpointmesh) :: kmesh

  real(8) :: occ_loc
  integer :: is, ik, iband

  occ_loc = 0.d0
  ! evaluate the function
  do is = 1,edisp%ispin
    do ik = ikstr, ikend
      do iband=1,edisp%nband_max
        if (edisp%band(iband,ik,is) < mu) then
          occ_loc = occ_loc + kmesh%weight(ik)
        endif
      enddo
    enddo
  enddo

#ifdef MPI
  call MPI_ALLREDUCE(occ_loc, occ_tot, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, mpierr)
#else
  occ_tot = occ_loc
#endif

end subroutine occ_DFT

end module Mroot
