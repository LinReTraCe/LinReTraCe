module Mfermi
  use Mparams
  use psi_fast, only: psi_range_dp, psi_range_qp, psi0_dp, psi0_qp, &
                      psi0_imag_dp, psi0_imag_qp, psi123_dp, psi123_qp
  implicit none

  ! fermi function
  interface fermi
    module procedure fermi_dp, fermi_qp, fermi_dpqp, fermi_qpdp
  end interface fermi

  ! one minus fermi function
  interface omfermi
    module procedure omfermi_dp, omfermi_qp, omfermi_dpqp, omfermi_qpdp
  end interface omfermi

  ! fermi function - 0.5
  interface fermimhalf
    module procedure fermi2_dp, fermi2_qp, fermi2_dpqp, fermi2_qpdp
  end interface fermimhalf

  ! substitution of Re psi_1 in the Gamma -> 0 limit
  ! -> Re[psi_1[0.5 + beta/2pi * (i*a)]] ==  pi**2 / 2 / cosh(beta * eps / 2)**2
  interface polygamma2fermi
    module procedure polygamma2fermi_dp, polygamma2fermi_qp
  end interface polygamma2fermi

  ! -> Re[psi_1[0.5 + beta/2pi * (gamma + i*a)]]
  interface polygamma2psi1
    module procedure polygamma2psi1_dp, polygamma2psi1_qp
  end interface polygamma2psi1

  ! T=0 limit of the fermi function: Heaviside theta(-eps) with value 1/2 at eps=0
  interface fermi_T0
    module procedure fermi_T0_dp, fermi_T0_qp
  end interface fermi_T0

  ! T=0 limit of one minus the fermi function: Heaviside theta(eps) with value 1/2 at eps=0
  interface omfermi_T0
    module procedure omfermi_T0_dp, omfermi_T0_qp
  end interface omfermi_T0

  ! T=0 (beta -> infinity) limit of the psi_0 (digamma) occupation kernel
  ! 1/pi * Im[psi_0(0.5 + beta/2pi * (gamma + i*a))]  --beta->inf-->  1/pi * atan2(a, gamma)
  ! since psi_0(z) -> ln(z) for |z| -> infinity and Im[ln(gamma + i*a)] = arg(gamma + i*a);
  ! the subleading term -1/(2z) contributes O(1/(beta*gamma)) and vanishes identically at T=0.
  ! atan2 covers the corner cases gracefully:
  !   gamma  > 0          : Lorentzian-broadened occupation kernel
  !   gamma -> 0, a <> 0  : +-1/2, i.e. the Heaviside limit is recovered
  !   gamma  = 0, a  = 0  : 0, i.e. occupation 1/2 -- no 0/0
  ! occupations built from this kernel:
  !   occ_digamma convention (gamma - i*a): n = 1/2 - atankern(gamma, a)
  !   occ_refine  convention (gamma + i*a): elec = 1/2 - atankern ; hole = 1/2 + atankern
  ! with a = Z*(eps - mu)
  interface atankern_T0
    module procedure atankern_T0_dp, atankern_T0_qp
  end interface atankern_T0

  ! T=0 limit of the degeneracy-weighted fermi factor 1/(1 + g*exp(beta*x)):
  ! theta(-x) with value 1/(1+g) at x=0
  ! (used for the impurity occupation; g=1 recovers fermi_T0)
  interface gfermi_T0
    module procedure gfermi_T0_dp, gfermi_T0_qp
  end interface gfermi_T0

  contains

!________________________________________________________
! fermi function in different precisions

  ! fermi function: energy, beta in double precision
  pure elemental function fermi_dp(eps,beta) result(f)
    implicit none
    real(8), intent(in) :: eps
    real(8), intent(in) :: beta
    real(8)             :: f
    f = 1.d0 / (1.d0 + EXP(beta*eps))
  end function fermi_dp

  ! fermi function: energy, beta in quad precision
  pure elemental function fermi_qp(eps,beta) result(f)
    implicit none
    real(16), intent(in) :: eps
    real(16), intent(in) :: beta
    real(16)             :: f
    f = 1.q0 / (1.q0 + EXP(eps*beta))
  end function fermi_qp

  ! fermi function: energy in double precision, beta in quad precision
  pure elemental function fermi_dpqp(eps,beta) result(f)
    implicit none
    real(8), intent(in)  :: eps
    real(16), intent(in) :: beta
    real(16)             :: f
    f = 1.q0 / (1.q0 + EXP(beta*eps))
  end function fermi_dpqp

  ! fermi function: energy in quad precision, beta in double precision
  pure elemental function fermi_qpdp(eps,beta) result(f)
    implicit none
    real(16), intent(in) :: eps
    real(8), intent(in)  :: beta
    real(16)             :: f
    f = 1.q0 / (1.q0 + EXP(eps*beta))
  end function fermi_qpdp

!________________________________________________________
! one minus fermi function in different precision

  pure elemental function omfermi_dp(eps,beta) result(f)
    implicit none
    real(8), intent(in) :: eps
    real(8), intent(in) :: beta
    real(8)             :: f
    f = 1.d0 / (1.d0 + EXP(-beta*eps))
  end function omfermi_dp

  pure elemental function omfermi_qp(eps,beta) result(f)
    implicit none
    real(16), intent(in) :: eps
    real(16), intent(in) :: beta
    real(16)             :: f
    f = 1.q0 / (1.q0 + EXP(-eps*beta))
  end function omfermi_qp

  pure elemental function omfermi_dpqp(eps,beta) result(f)
    implicit none
    real(8), intent(in)  :: eps
    real(16), intent(in) :: beta
    real(16)             :: f
    f = 1.q0 / (1.q0 + EXP(-beta*eps))
  end function omfermi_dpqp

  pure elemental function omfermi_qpdp(eps,beta) result(f)
    implicit none
    real(16), intent(in) :: eps
    real(8), intent(in)  :: beta
    real(16)             :: f
    f = 1.q0 / (1.q0 + EXP(-eps*beta))
  end function omfermi_qpdp

!________________________________________________________
! fermi function - 1/2 in different precisions

  pure elemental function fermi2_dp(eps,beta) result(f)
    implicit none
    real(8), intent(in) :: eps
    real(8), intent(in) :: beta
    real(8)             :: f
    f = (1.d0 - EXP(beta*eps)) / (2.d0* (1.d0 + EXP(beta*eps)))
  end function fermi2_dp

  pure elemental function fermi2_qp(eps,beta) result(f)
    implicit none
    real(16), intent(in) :: eps
    real(16), intent(in) :: beta
    real(16)             :: f
    f = (1.q0 - EXP(beta*eps)) / (2.q0* (1.q0 + EXP(beta*eps)))
  end function fermi2_qp

  pure elemental function fermi2_dpqp(eps,beta) result(f)
    implicit none
    real(8), intent(in)  :: eps
    real(16), intent(in) :: beta
    real(16)             :: f
    f = (1.q0 - EXP(beta*eps)) / (2.d0* (1.q0 + EXP(beta*eps)))
  end function fermi2_dpqp

  pure elemental function fermi2_qpdp(eps,beta) result(f)
    implicit none
    real(16), intent(in) :: eps
    real(8), intent(in)  :: beta
    real(16)             :: f
    f = (1.q0 - EXP(beta*eps)) / (2.q0* (1.q0 + EXP(beta*eps)))
  end function fermi2_qpdp

!________________________________________________________

  ! these functions represent the leading term of the limit Gamma-> 0
  ! used in the Boltzmann regime equations
  !
  ! Re[psi_1[0.5 + beta/2pi * (i*a)]] == -d/deps fermi * 2 pi**2 / beta
  ! -d/deps = beta / 4 / cosh(beta * eps / 2)**2
  ! -> Re[psi_1[0.5 + beta/2pi * (i*a)]] ==  pi**2 / 2 / cosh(beta * eps / 2)**2
  pure elemental function polygamma2fermi_dp(eps,beta)
    implicit none
    real(8), intent(in) :: eps,beta
    real(8) :: polygamma2fermi_dp
    polygamma2fermi_dp = pi**2 / (2.d0 * cosh(beta*eps/2.d0)**2)
  end function polygamma2fermi_dp

  pure elemental function polygamma2fermi_qp(eps,beta)
    implicit none
    real(16), intent(in) :: eps,beta
    real(16) :: polygamma2fermi_qp
    polygamma2fermi_qp = piQ**2 / (2.q0 * cosh(beta*eps/2.q0)**2)
  end function polygamma2fermi_qp

  ! these functions represent the psi_1 approximation
  ! i.e. all higher order psi_i i>1 are thrown out
  ! psi_range_dp / psi_range_qp (module psi_fast)
  !
  !   call psi_range_xp(z, klo, khi, psi [, ierr])
  !     z    complex(8|16), argument; Re z > 0 for every LinReTraCe argument
  !     klo  integer, lowest derivative order wanted   (0 <= klo <= khi <= 4)
  !     khi  integer, highest derivative order wanted
  !     psi  complex(8|16) psi(klo:khi), psi(K) = psi^(K)(z) on return
  !     ierr optional integer, PSI_ERR_OK / PSI_ERR_ORDER / PSI_ERR_POLE
  !
  !   All requested orders come out of a single recurrence pass, so asking
  !   for a narrow range costs correspondingly less: with klo = khi = 1 only
  !   1/V**2 is accumulated.
  function polygamma2psi1_dp(gamma,eps,beta)
    implicit none
    real(8), intent(in) :: gamma,eps,beta
    real(8) :: polygamma2psi1_dp
    complex(8) :: psi(1:1)
    call psi_range_dp(0.5d0 + beta/2.d0/pi * (gamma + ci*eps), 1, 1, psi)
    polygamma2psi1_dp = real(psi(1))
    ! --- deprecated CERNLIB path (kept until wpsipg is retired) -------------
    ! complex(8), external :: wpsipg
    ! polygamma2psi1_dp = real(wpsipg(0.5d0 + beta/2.d0/pi * (gamma + ci*eps),1))
  end function polygamma2psi1_dp

  function polygamma2psi1_qp(gamma,eps,beta)
    implicit none
    real(16), intent(in) :: eps,beta
    real(8), intent(in)  :: gamma
    real(16) :: polygamma2psi1_qp
    complex(16) :: psi(1:1)
    call psi_range_qp(0.5q0 + beta/2.q0/piQ * (gamma + ciQ*eps), 1, 1, psi)
    polygamma2psi1_qp = real(psi(1))
    ! --- deprecated CERNLIB path (kept until wpsipghp is retired) ----------
    ! complex(16), external :: wpsipghp
    ! polygamma2psi1_qp = real(wpsipghp(0.5q0 + beta/2.q0/piQ * (gamma + ciQ*eps),1))
  end function polygamma2psi1_qp

!________________________________________________________
! T=0 limits of the distribution functions used in the chemical potential search
! note: these are step / arctan functions of the *sign* of the argument;
! for arguments a = Z*(eps-mu) with Z > 0 the quasi-particle weight hence
! drops out of the Heaviside functions identically

  ! T=0 fermi function: theta(-eps), 1/2 at eps=0
  ! the eps=0 case is not measure-zero on a discrete k-grid (bisection can land
  ! exactly on a band energy), so it is handled explicitly
  pure elemental function fermi_T0_dp(eps) result(f)
    implicit none
    real(8), intent(in) :: eps
    real(8)             :: f
    if (eps < 0.d0) then
      f = 1.d0
    else if (eps > 0.d0) then
      f = 0.d0
    else
      f = 0.5d0
    endif
  end function fermi_T0_dp

  pure elemental function fermi_T0_qp(eps) result(f)
    implicit none
    real(16), intent(in) :: eps
    real(16)             :: f
    if (eps < 0.q0) then
      f = 1.q0
    else if (eps > 0.q0) then
      f = 0.q0
    else
      f = 0.5q0
    endif
  end function fermi_T0_qp

  ! T=0 one minus fermi function: theta(eps), 1/2 at eps=0
  pure elemental function omfermi_T0_dp(eps) result(f)
    implicit none
    real(8), intent(in) :: eps
    real(8)             :: f
    if (eps > 0.d0) then
      f = 1.d0
    else if (eps < 0.d0) then
      f = 0.d0
    else
      f = 0.5d0
    endif
  end function omfermi_T0_dp

  pure elemental function omfermi_T0_qp(eps) result(f)
    implicit none
    real(16), intent(in) :: eps
    real(16)             :: f
    if (eps > 0.q0) then
      f = 1.q0
    else if (eps < 0.q0) then
      f = 0.q0
    else
      f = 0.5q0
    endif
  end function omfermi_T0_qp

  ! T=0 digamma occupation kernel: 1/pi * atan2(a, gamma)
  ! (see interface block above for derivation and sign conventions)
  ! gamma is kept double precision in accordance with sct%gam
  pure elemental function atankern_T0_dp(gamma,a) result(f)
    implicit none
    real(8), intent(in) :: gamma
    real(8), intent(in) :: a
    real(8)             :: f
    f = atan2(a, gamma) / pi
  end function atankern_T0_dp

  pure elemental function atankern_T0_qp(gamma,a) result(f)
    implicit none
    real(8), intent(in)  :: gamma
    real(16), intent(in) :: a
    real(16)             :: f
    f = atan2(a, real(gamma,16)) / piQ
  end function atankern_T0_qp

  ! T=0 degeneracy-weighted fermi factor: theta(-x), 1/(1+g) at x=0
  pure elemental function gfermi_T0_dp(x,g) result(f)
    implicit none
    real(8), intent(in) :: x
    real(8), intent(in) :: g
    real(8)             :: f
    if (x < 0.d0) then
      f = 1.d0
    else if (x > 0.d0) then
      f = 0.d0
    else
      f = 1.d0 / (1.d0 + g)
    endif
  end function gfermi_T0_dp

  pure elemental function gfermi_T0_qp(x,g) result(f)
    implicit none
    real(16), intent(in) :: x
    real(8), intent(in)  :: g
    real(16)             :: f
    if (x < 0.q0) then
      f = 1.q0
    else if (x > 0.q0) then
      f = 0.q0
    else
      f = 1.q0 / (1.q0 + real(g,16))
    endif
  end function gfermi_T0_qp

end module Mfermi
