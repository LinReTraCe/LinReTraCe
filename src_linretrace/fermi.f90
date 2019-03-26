module Mfermi
  use Mparams
  implicit none

  interface fermi
    module procedure fermi_dp, fermi_qp, fermi_dpqp, fermi_qpdp
  end interface fermi

  interface polygamma2fermi
    module procedure polygamma2fermi_dp, polygamma2fermi_qp
  end interface polygamma2fermi

  contains

  pure elemental function fermi_dp(eps,beta) result(f)
    implicit none
    real(8), intent(in) :: eps
    real(8), intent(in) :: beta
    real(8)             :: f
    f = 1.d0 / (1.d0 + EXP(beta*eps))
  end function fermi_dp

  pure elemental function fermi_qp(eps,beta) result(f)
    implicit none
    real(16), intent(in) :: eps
    real(16), intent(in) :: beta
    real(16)             :: f
    f = 1.q0 / (1.q0 + EXP(eps*beta))
  end function fermi_qp

  pure elemental function fermi_dpqp(eps,beta) result(f)
    implicit none
    real(8), intent(in)  :: eps
    real(16), intent(in) :: beta
    real(16)             :: f
    f = 1.q0 / (1.q0 + EXP(beta*eps))
  end function fermi_dpqp

  pure elemental function fermi_qpdp(eps,beta) result(f)
    implicit none
    real(16), intent(in) :: eps
    real(8), intent(in)  :: beta
    real(16)             :: f
    f = 1.q0 / (1.q0 + EXP(eps*beta))
  end function fermi_qpdp


  ! these functions represent the analytical limit of
  ! lim Gamma-> 0     Re[psi_1[0.5 + beta/2pi * (Gamma + i*a)]]
  ! used in the Boltzmann regime equations
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
    polygamma2fermi_qp = piQ**2 / (2.d0 * cosh(beta*eps/2.d0)**2)
  end function polygamma2fermi_qp

end module Mfermi
