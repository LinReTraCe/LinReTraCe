!------------------------------------------------------------------------------
! psi_fast.F90
!
! Fast polygamma functions psi^(K)(z), K = 0..4, for complex arguments.
!
! Replacement for the CERNLIB C317 routines WPSIPG (double) and WPSIPGHP
! (quad) as used by LinReTraCe.  The old routines remain in the archive and
! are not touched; call sites migrate one at a time.  Same algorithm --
! recurrence to push the argument up, then a Bernoulli asymptotic series --
! but with four changes:
!
!   1. The recurrence is entered on |z|, not on |Re z|.
!   2. The recurrence walks until |V| reaches a target chosen so that the
!      series actually delivers the precision of the working kind.
!   3. Several orders K are produced from a single pass.
!   4. The number of series terms is chosen at run time from |V|.
!
! See the notes on each below.  Nothing here changes the mathematics; the
! results agree with the old routines to the accuracy the old routines
! achieve, and exceed it near the real axis.
!
! Jan Tomczak's quad-precision conversion of the CERNLIB routine, and its
! extension to 16 Bernoulli terms, are the basis for this module.
!------------------------------------------------------------------------------
module psi_fast

  use, intrinsic :: iso_fortran_env, only: dp => real64, qp => real128, &
                                           int64, error_unit
  implicit none
  private

  public :: psi_range_qp, psi0_qp, psi0_imag_qp, psi123_qp
  public :: psi_range_dp, psi0_dp, psi0_imag_dp, psi123_dp
  public :: PSI_ERR_OK, PSI_ERR_ORDER, PSI_ERR_POLE

  ! -- error codes ------------------------------------------------------------
  integer, parameter :: PSI_ERR_OK    = 0
  integer, parameter :: PSI_ERR_ORDER = 1   ! K outside 0..4
  integer, parameter :: PSI_ERR_POLE  = 2   ! argument is a non-positive integer

  ! -- Bernoulli numbers B_2i, i = 1..16, as exact fractions -------------------
  !
  ! Given as integer numerator/denominator pairs rather than as decimal
  ! literals on purpose.  The CERNLIB table hard-codes all 80 values C(i,K)
  ! as 30-digit literals, and one of them is wrong: C(13,0) is written
  ! 5.4827583333...Q+3 where B_26/26 = 5.4827583333...Q+4, a factor of ten.
  ! With kmax >= 13 -- i.e. in every quad build -- that degrades psi_0, and
  ! hence the chemical-potential search, from ~2e-31 to ~3e-27 relative.
  ! Deriving the table from 16 exact fractions removes the whole class of
  ! transcription error, and every value below can be checked against any
  ! standard Bernoulli table by eye.
  integer(int64), parameter :: BNUM(16) = [ &
        1_int64,            -1_int64,             1_int64,            -1_int64, &
        5_int64,          -691_int64,             7_int64,         -3617_int64, &
    43867_int64,       -174611_int64,        854513_int64,    -236364091_int64, &
  8553103_int64, -23749461029_int64, 8615841276005_int64, -7709321041217_int64 ]

  integer(int64), parameter :: BDEN(16) = [ &
        6_int64,   30_int64,   42_int64,   30_int64, &
       66_int64, 2730_int64,    6_int64,  510_int64, &
      798_int64,  330_int64,  138_int64, 2730_int64, &
        6_int64,  870_int64, 14322_int64, 510_int64 ]

  real(qp), parameter :: BQ(16) = real(BNUM, qp) / real(BDEN, qp)

  ! -- series coefficients C(i,K) = B_2i * (2i+K-1)! / (2i)! ------------------
  !
  ! The factorial ratio is an exact small integer (or its reciprocal):
  !   K=0 : 1/(2i)          K=1 : 1            K=2 : (2i+1)
  !   K=3 : (2i+1)(2i+2)    K=4 : (2i+1)(2i+2)(2i+3)
  ! so the K dependence is built here instead of being typed out.  Largest
  ! integer factor is 33*34*35 = 39270 at i = 16.
  integer, parameter :: KMAX_TAB = 16

  ! 2i for i = 1..16, written out so that the table below is built from
  ! whole-array elemental operations only.  An array constructor with an
  ! implied-do would need its index variable declared under IMPLICIT NONE,
  ! which would mean a module-scope integer shadowed by every local i.
  integer, parameter :: TWOI(KMAX_TAB) = &
      [ 2,  4,  6,  8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32]

  real(qp), parameter :: CQ(KMAX_TAB,0:4) = reshape( [        &
      BQ / real(TWOI, qp),                                    &
      BQ,                                                     &
      BQ * real(TWOI + 1, qp),                                &
      BQ * real((TWOI + 1)*(TWOI + 2), qp),                   &
      BQ * real((TWOI + 1)*(TWOI + 2)*(TWOI + 3), qp) ],      &
      [KMAX_TAB, 5] )

  real(dp), parameter :: CD(KMAX_TAB,0:4) = real(CQ, dp)

  ! -- algorithm constants ----------------------------------------------------
  real(qp), parameter :: PIQ = 3.14159265358979323846264338327950288_qp
  real(dp), parameter :: PID = 3.14159265358979323846_dp

  integer,  parameter :: FCT(-1:4) = [0, 1, 1, 2, 6, 24]
  integer,  parameter :: SGN(0:4)  = [-1, 1, -1, 1, -1]

  ! Modulus the recurrence walks the argument up to.  Chosen so that the
  ! 16-term (quad) / 9-term (double) series is accurate to the working
  ! epsilon at the hardest point, Im z ~ 0.  Measured worst case over
  ! K = 0..4 at Re z = 1/2 (quad eps = 1.9e-34, double eps = 2.2e-16):
  !
  !   quad,   walk to 22: 6.4e-37     walk to 28: 2.1e-40
  !   double, walk to 10: 4.4e-21     walk to 12: 1.5e-21
  !
  ! and, for arguments that need no walk at all, kmax = KMAX_TAB gives
  !   |V| = 22: 8.2e-31   |V| = 26: 2.7e-33   |V| = 28: 2.2e-34  <= 4*eps
  ! so 28 is the smallest target at which the gate and the series agree on
  ! full quad precision.  The old routine walked to Re V ~ 15 and delivered
  ! ~2e-31, i.e. three digits short of quad, precisely at E = mu.
  !
  ! Lowering PSI_VTARGET_QP to 22 reproduces the old accuracy at lower cost;
  ! it is a single constant, deliberately not a preprocessor symbol.
  real(qp), parameter :: PSI_VTARGET_QP = 28.0_qp
  real(dp), parameter :: PSI_VTARGET_DP = 12.0_dp

  ! Hard cap on recurrence steps.  Unreachable for the arguments LinReTraCe
  ! produces (Re z = 1/2 + beta*Gamma/2pi > 0), present only so that a
  ! pathological call cannot spin.
  integer, parameter :: NS_CAP = 64

  real(qp), parameter :: DELTA_QP = 1.0e-20_qp
  real(dp), parameter :: DELTA_DP = 1.0e-12_dp

contains

!------------------------------------------------------------------------------
! Number of series terms needed at a given |V|.
!
! The old routines fixed this at compile time via "-D kmax=", which also
! controlled which DATA statements existed, through a nest of #if blocks.
! Here the full table always exists and only the Horner start index varies,
! so the choice costs a handful of comparisons and no recompilation.
!
! Thresholds are the smallest kmax reaching 4*epsilon with no recurrence,
! worst case over K = 0..4 at Re z = 1/2.
!------------------------------------------------------------------------------
  pure integer function kterms_qp(absv) result(k)
    real(qp), intent(in) :: absv
    if      (absv >= 2000.0_qp) then ; k =  5
    else if (absv >=  500.0_qp) then ; k =  6
    else if (absv >=  200.0_qp) then ; k =  8
    else if (absv >=  100.0_qp) then ; k =  9
    else if (absv >=   60.0_qp) then ; k = 11
    else if (absv >=   40.0_qp) then ; k = 13
    else if (absv >=   30.0_qp) then ; k = 15
    else                             ; k = KMAX_TAB
    end if
  end function kterms_qp

  pure integer function kterms_dp(absv) result(k)
    real(dp), intent(in) :: absv
    if      (absv >= 60.0_dp) then ; k = 4
    else if (absv >= 35.0_dp) then ; k = 5
    else if (absv >= 25.0_dp) then ; k = 6
    else if (absv >= 16.0_dp) then ; k = 7
    else if (absv >= 14.0_dp) then ; k = 8
    else                           ; k = 9
    end if
  end function kterms_dp

!------------------------------------------------------------------------------
! Recurrence step count.
!
! The recurrence adds 1 to Re V per step, so |V| after n steps is
! sqrt((X+n)^2 + Y^2).  Requiring that to reach T gives
!
!     n = ceil( sqrt(T^2 - Y^2) - X )   for |Y| < T,      n = 0 otherwise.
!
! The old routines instead tested |Re z| < 15 and always took 14 - int|Re z|
! steps.  For LinReTraCe arguments Re z = 1/2 + beta*Gamma/2pi is never near
! 15, so the test always fired -- and because the step only moves the real
! part, a state with |Im z| = 185 (a 1 eV band energy at 10 K) had its |V|
! moved from 185.0007 to 185.6 by fourteen quad-precision divisions that
! bought nothing.  Below the recurrence is skipped outright for such states,
! which is where nearly all of the speedup comes from at low temperature and
! for the wide E - mu windows of the topological kernels.
!------------------------------------------------------------------------------
  pure integer function nsteps_qp(x, y, tgt) result(n)
    real(qp), intent(in) :: x, y, tgt
    real(qp) :: rad
    rad = tgt*tgt - y*y
    if (rad <= 0.0_qp) then
      n = 0
    else
      n = ceiling(sqrt(rad) - x)
      n = max(n, 0)
      n = min(n, NS_CAP)
    end if
  end function nsteps_qp

  pure integer function nsteps_dp(x, y, tgt) result(n)
    real(dp), intent(in) :: x, y, tgt
    real(dp) :: rad
    rad = tgt*tgt - y*y
    if (rad <= 0.0_dp) then
      n = 0
    else
      n = ceiling(sqrt(rad) - x)
      n = max(n, 0)
      n = min(n, NS_CAP)
    end if
  end function nsteps_dp

!==============================================================================
! psi_range_qp : psi^(K)(z) for K = klo .. khi, from one recurrence pass.
!
! The recurrence sum is  H_K = sum_j 1/(V_0 + j)^(K+1).  All orders share the
! same V sequence, so one reciprocal u = 1/V per step yields every order by
! multiplication:  u^2, u^3, u^4 for K = 1,2,3.  The old code called the
! routine once per K and recomputed 1/V**(K+1) from scratch each time -- 45
! complex divisions where 15 divisions and a few multiplies suffice.  In quad
! that matters disproportionately: a complex division is ~625 ns against
! ~184 ns for a multiply, both software-emulated.
!
! psi = psi(klo:khi), i.e. psi(K) holds psi^(K)(z).
!==============================================================================
  subroutine psi_range_qp(z, klo, khi, psi, ierr)
    complex(qp), intent(in)            :: z
    integer,     intent(in)            :: klo, khi
    complex(qp), intent(out)           :: psi(klo:khi)
    integer,     intent(out), optional :: ierr

    complex(qp) :: u, v, h(0:4), uu(0:4), p, r, refl, rsq
    real(qp)    :: x, y, absz, absv
    integer     :: k, i, ns, nk, err
    logical     :: reflect

    err = PSI_ERR_OK
    psi = (0.0_qp, 0.0_qp)

    if (klo < 0 .or. khi > 4 .or. klo > khi) then
      err = PSI_ERR_ORDER
      call report(err, 'psi_range_qp: order out of range 0..4')
      if (present(ierr)) ierr = err
      return
    end if

    u = z
    x = real(u, qp)

    ! pole check: non-positive integer argument
    if (abs(aimag(u)) < DELTA_QP .and. abs(x + nint(abs(x))) < DELTA_QP) then
      err = PSI_ERR_POLE
      call report(err, 'psi_range_qp: argument is a non-positive integer')
      if (present(ierr)) ierr = err
      return
    end if

    ! reflection to the right half plane; never taken by LinReTraCe, where
    ! Re z = 1/2 + beta*Gamma/2pi > 0, but kept for generality
    reflect = (x < 0.0_qp)
    if (reflect) u = -u

    v    = u
    h    = (0.0_qp, 0.0_qp)
    x    = real(v, qp)
    y    = abs(aimag(v))
    absz = abs(v)

    if (absz < PSI_VTARGET_QP) then
      ns = nsteps_qp(x, y, PSI_VTARGET_QP)
      do i = 1, ns
        uu(0) = 1.0_qp / v                     ! the only division per step
        if (khi >= 1) uu(1) = uu(0)*uu(0)
        if (khi >= 2) uu(2) = uu(1)*uu(0)
        if (khi >= 3) uu(3) = uu(1)*uu(1)
        if (khi >= 4) uu(4) = uu(3)*uu(0)
        do k = klo, khi
          h(k) = h(k) + uu(k)
        end do
        v = v + 1.0_qp
      end do
    end if

    absv = abs(v)
    nk   = kterms_qp(absv)
    r    = 1.0_qp / (v*v)

    do k = klo, khi
      p = r * CQ(nk, k)
      do i = nk-1, 1, -1
        p = r * (CQ(i, k) + p)
      end do
      psi(k) = SGN(k) * ( FCT(k)*h(k)                                     &
                        + (v*(FCT(k-1) + p) + 0.5_qp*FCT(k)) / v**(k+1) )
      if (k == 0) psi(k) = psi(k) + log(v)
    end do

    if (reflect) then
      refl = cot_pi_qp(u)
      rsq  = refl*refl
      do k = klo, khi
        select case (k)
        case (0)
          psi(k) =  psi(k) + 1.0_qp/u          + PIQ*refl
        case (1)
          psi(k) = -psi(k) + 1.0_qp/u**2       + PIQ**2*(rsq + 1.0_qp)
        case (2)
          psi(k) =  psi(k) + 2.0_qp/u**3       + 2.0_qp*PIQ**3*refl*(rsq + 1.0_qp)
        case (3)
          psi(k) = -psi(k) + 6.0_qp/u**4       + 2.0_qp*PIQ**4*((3.0_qp*rsq + 4.0_qp)*rsq + 1.0_qp)
        case (4)
          psi(k) =  psi(k) + 24.0_qp/u**5      + 8.0_qp*PIQ**5*refl*((3.0_qp*rsq + 5.0_qp)*rsq + 2.0_qp)
        end select
      end do
    end if

    if (present(ierr)) ierr = err
  end subroutine psi_range_qp

!------------------------------------------------------------------------------
! cot(pi*u), written the way the CERNLIB routine does it, i.e. split into
! real trigonometric and hyperbolic parts so that large |Im u| does not
! overflow the complex cotangent.
!------------------------------------------------------------------------------
  pure complex(qp) function cot_pi_qp(u) result(p)
    complex(qp), intent(in) :: u
    complex(qp) :: w
    real(qp)    :: xx, yy, sa, cb, th
    w  = PIQ*u
    xx = real(w, qp)
    yy = aimag(w)
    sa = sin(xx)
    cb = cos(xx)
    th = tanh(yy)
    p  = cmplx(cb, -sa*th, qp) / cmplx(sa, cb*th, qp)
  end function cot_pi_qp

!==============================================================================
! psi0_qp / psi123_qp : named convenience wrappers.
!
! psi123_qp is what calc_polygamma wants; psi0_qp is what a kernel needing
! the full complex digamma wants (the topological transport terms carry a
! Fermi function rather than a df/dw factor, so psi_0 enters the response and
! not only the particle number).
!==============================================================================
  complex(qp) function psi0_qp(z, ierr) result(p)
    complex(qp), intent(in)            :: z
    integer,     intent(out), optional :: ierr
    complex(qp) :: tmp(0:0)
    call psi_range_qp(z, 0, 0, tmp, ierr)
    p = tmp(0)
  end function psi0_qp

  subroutine psi123_qp(z, psi, ierr)
    complex(qp), intent(in)            :: z
    complex(qp), intent(out)           :: psi(1:3)
    integer,     intent(out), optional :: ierr
    call psi_range_qp(z, 1, 3, psi, ierr)
  end subroutine psi123_qp

!==============================================================================
! psi0_imag_qp : Im psi^(0)(z) only.
!
! occ_digamma uses nothing but aimag(psi_0), and psi_0 is the one order that
! needs a complex logarithm at the end.  In quad that log costs ~2046 ns,
! against ~625 ns for a complex reciprocal and ~184 ns for a multiply -- so
! once the recurrence is skipped it is the single most expensive operation in
! the call.  Only its imaginary part is wanted, and
!
!     Im log(V) = atan2(Im V, Re V)     (~663 ns)
!
! so the log is replaced outright.  This routine is called once per band, per
! k-point, per root-finder iteration, which makes it the most frequently
! executed floating-point kernel in the program.
!
! The reflected branch (Re z < 0) falls back to the full complex evaluation;
! it cannot occur for LinReTraCe arguments.
!==============================================================================
  real(qp) function psi0_imag_qp(z, ierr) result(pim)
    complex(qp), intent(in)            :: z
    integer,     intent(out), optional :: ierr

    complex(qp) :: u, v, h, p, r, tmp(0:0)
    real(qp)    :: x, y, absz, absv
    integer     :: i, ns, nk, err

    err = PSI_ERR_OK
    pim = 0.0_qp

    u = z
    x = real(u, qp)

    if (abs(aimag(u)) < DELTA_QP .and. abs(x + nint(abs(x))) < DELTA_QP) then
      err = PSI_ERR_POLE
      call report(err, 'psi0_imag_qp: argument is a non-positive integer')
      if (present(ierr)) ierr = err
      return
    end if

    if (x < 0.0_qp) then                       ! rare path, not taken here
      call psi_range_qp(z, 0, 0, tmp, err)
      pim = aimag(tmp(0))
      if (present(ierr)) ierr = err
      return
    end if

    v    = u
    h    = (0.0_qp, 0.0_qp)
    y    = abs(aimag(v))
    absz = abs(v)

    if (absz < PSI_VTARGET_QP) then
      ns = nsteps_qp(x, y, PSI_VTARGET_QP)
      do i = 1, ns
        h = h + 1.0_qp/v                       ! K = 0 needs only u^1
        v = v + 1.0_qp
      end do
    end if

    absv = abs(v)
    nk   = kterms_qp(absv)
    r    = 1.0_qp / (v*v)
    p    = r * CQ(nk, 0)
    do i = nk-1, 1, -1
      p = r * (CQ(i, 0) + p)
    end do

    ! SGN(0) = -1, FCT(0) = 1, FCT(-1) = 0
    pim = aimag( -( h + (v*p + 0.5_qp)/v ) ) + atan2(aimag(v), real(v, qp))

    if (present(ierr)) ierr = err
  end function psi0_imag_qp

!==============================================================================
! Double precision counterparts.  Identical structure; the coefficient table
! is the quad table rounded, so the two kinds cannot drift apart.
!==============================================================================
  subroutine psi_range_dp(z, klo, khi, psi, ierr)
    complex(dp), intent(in)            :: z
    integer,     intent(in)            :: klo, khi
    complex(dp), intent(out)           :: psi(klo:khi)
    integer,     intent(out), optional :: ierr

    complex(dp) :: u, v, h(0:4), uu(0:4), p, r, refl, rsq
    real(dp)    :: x, y, absz, absv
    integer     :: k, i, ns, nk, err
    logical     :: reflect

    err = PSI_ERR_OK
    psi = (0.0_dp, 0.0_dp)

    if (klo < 0 .or. khi > 4 .or. klo > khi) then
      err = PSI_ERR_ORDER
      call report(err, 'psi_range_dp: order out of range 0..4')
      if (present(ierr)) ierr = err
      return
    end if

    u = z
    x = real(u, dp)

    if (abs(aimag(u)) < DELTA_DP .and. abs(x + nint(abs(x))) < DELTA_DP) then
      err = PSI_ERR_POLE
      call report(err, 'psi_range_dp: argument is a non-positive integer')
      if (present(ierr)) ierr = err
      return
    end if

    reflect = (x < 0.0_dp)
    if (reflect) u = -u

    v    = u
    h    = (0.0_dp, 0.0_dp)
    x    = real(v, dp)
    y    = abs(aimag(v))
    absz = abs(v)

    if (absz < PSI_VTARGET_DP) then
      ns = nsteps_dp(x, y, PSI_VTARGET_DP)
      do i = 1, ns
        uu(0) = 1.0_dp / v                     ! the only division per step
        if (khi >= 1) uu(1) = uu(0)*uu(0)
        if (khi >= 2) uu(2) = uu(1)*uu(0)
        if (khi >= 3) uu(3) = uu(1)*uu(1)
        if (khi >= 4) uu(4) = uu(3)*uu(0)
        do k = klo, khi
          h(k) = h(k) + uu(k)
        end do
        v = v + 1.0_dp
      end do
    end if

    absv = abs(v)
    nk   = kterms_dp(absv)
    r    = 1.0_dp / (v*v)

    do k = klo, khi
      p = r * CD(nk, k)
      do i = nk-1, 1, -1
        p = r * (CD(i, k) + p)
      end do
      psi(k) = SGN(k) * ( FCT(k)*h(k)                                     &
                        + (v*(FCT(k-1) + p) + 0.5_dp*FCT(k)) / v**(k+1) )
      if (k == 0) psi(k) = psi(k) + log(v)
    end do

    if (reflect) then
      refl = cot_pi_dp(u)
      rsq  = refl*refl
      do k = klo, khi
        select case (k)
        case (0)
          psi(k) =  psi(k) + 1.0_dp/u     + PID*refl
        case (1)
          psi(k) = -psi(k) + 1.0_dp/u**2  + PID**2*(rsq + 1.0_dp)
        case (2)
          psi(k) =  psi(k) + 2.0_dp/u**3  + 2.0_dp*PID**3*refl*(rsq + 1.0_dp)
        case (3)
          psi(k) = -psi(k) + 6.0_dp/u**4  + 2.0_dp*PID**4*((3.0_dp*rsq + 4.0_dp)*rsq + 1.0_dp)
        case (4)
          psi(k) =  psi(k) + 24.0_dp/u**5 + 8.0_dp*PID**5*refl*((3.0_dp*rsq + 5.0_dp)*rsq + 2.0_dp)
        end select
      end do
    end if

    if (present(ierr)) ierr = err
  end subroutine psi_range_dp

  pure complex(dp) function cot_pi_dp(u) result(p)
    complex(dp), intent(in) :: u
    complex(dp) :: w
    real(dp)    :: xx, yy, sa, cb, th
    w  = PID*u
    xx = real(w, dp)
    yy = aimag(w)
    sa = sin(xx)
    cb = cos(xx)
    th = tanh(yy)
    p  = cmplx(cb, -sa*th, dp) / cmplx(sa, cb*th, dp)
  end function cot_pi_dp

  complex(dp) function psi0_dp(z, ierr) result(p)
    complex(dp), intent(in)            :: z
    integer,     intent(out), optional :: ierr
    complex(dp) :: tmp(0:0)
    call psi_range_dp(z, 0, 0, tmp, ierr)
    p = tmp(0)
  end function psi0_dp

  subroutine psi123_dp(z, psi, ierr)
    complex(dp), intent(in)            :: z
    complex(dp), intent(out)           :: psi(1:3)
    integer,     intent(out), optional :: ierr
    call psi_range_dp(z, 1, 3, psi, ierr)
  end subroutine psi123_dp

  real(dp) function psi0_imag_dp(z, ierr) result(pim)
    complex(dp), intent(in)            :: z
    integer,     intent(out), optional :: ierr

    complex(dp) :: u, v, h, p, r, tmp(0:0)
    real(dp)    :: x, y, absz, absv
    integer     :: i, ns, nk, err

    err = PSI_ERR_OK
    pim = 0.0_dp

    u = z
    x = real(u, dp)

    if (abs(aimag(u)) < DELTA_DP .and. abs(x + nint(abs(x))) < DELTA_DP) then
      err = PSI_ERR_POLE
      call report(err, 'psi0_imag_dp: argument is a non-positive integer')
      if (present(ierr)) ierr = err
      return
    end if

    if (x < 0.0_dp) then
      call psi_range_dp(z, 0, 0, tmp, err)
      pim = aimag(tmp(0))
      if (present(ierr)) ierr = err
      return
    end if

    v    = u
    h    = (0.0_dp, 0.0_dp)
    y    = abs(aimag(v))
    absz = abs(v)

    if (absz < PSI_VTARGET_DP) then
      ns = nsteps_dp(x, y, PSI_VTARGET_DP)
      do i = 1, ns
        h = h + 1.0_dp/v
        v = v + 1.0_dp
      end do
    end if

    absv = abs(v)
    nk   = kterms_dp(absv)
    r    = 1.0_dp / (v*v)
    p    = r * CD(nk, 0)
    do i = nk-1, 1, -1
      p = r * (CD(i, 0) + p)
    end do

    pim = aimag( -( h + (v*p + 0.5_dp)/v ) ) + atan2(aimag(v), real(v, dp))

    if (present(ierr)) ierr = err
  end function psi0_imag_dp

!------------------------------------------------------------------------------
! Error reporting.
!
! The old routines called MTLPRT, which pulls in the CERNLIB error machinery
! and, for C317.2, aborts.  This module has no dependencies beyond the
! intrinsic modules, so it writes to error_unit and lets the caller decide via
! the optional ierr.  Callers that want the old behaviour check ierr and call
! their own stop_with_message.
!------------------------------------------------------------------------------
  subroutine report(code, msg)
    integer,      intent(in) :: code
    character(*), intent(in) :: msg
    write(error_unit,'(a,i0,a,a)') 'psi_fast: error ', code, ': ', msg
  end subroutine report

end module psi_fast
