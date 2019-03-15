module Mfermi
  implicit none

  interface fermi
    module procedure fermi_dp, fermi_qp, fermi_dpqp, fermi_qpdp
  end interface fermi

  interface dfermi
    module procedure dfermi_dp, dfermi_qp
  end interface dfermi

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

  pure elemental function dfermi_dp(eps,beta)
    implicit none
    real(8), intent(in) :: eps,beta
    real(8) :: dfermi_dp
    dfermi_dp = (-1.d0) / ( exp(-beta*eps/2.d0) + exp(beta*eps/2.d0) )**2
  end function dfermi_dp

  pure elemental function dfermi_qp(eps,beta)
    implicit none
    real(16), intent(in) :: eps,beta
    real(16) :: dfermi_qp
    dfermi_qp = (-1.d0) / ( exp(-beta*eps/2.d0) + exp(beta*eps/2.d0) )**2
  end function dfermi_qp

end module Mfermi
