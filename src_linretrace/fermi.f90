module Mfermi
  use Mparams
  implicit none

  ! fermi function
  interface fermi
    module procedure fermi_dp, fermi_qp, fermi_dpqp, fermi_qpdp
  end interface fermi

  ! one minus fermi function
  interface omfermi
    module procedure omfermi_dp, omfermi_qp, omfermi_dpqp, omfermi_qpdp
  end interface omfermi

  ! fermi - 1/2
  interface fermi2
    module procedure fermi2_dp, fermi2_qp, fermi2_dpqp, fermi2_qpdp
  end interface fermi2

  interface polygamma2fermi
    module procedure polygamma2fermi_dp, polygamma2fermi_qp
  end interface polygamma2fermi

  contains

!________________________________________________________
! fermi function in different precisions

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
    polygamma2fermi_qp = piQ**2 / (2.q0 * cosh(beta*eps/2.q0)**2)
  end function polygamma2fermi_qp

end module Mfermi
