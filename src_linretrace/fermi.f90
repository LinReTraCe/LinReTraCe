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
  interface fermimhalf
    module procedure fermi2_dp, fermi2_qp, fermi2_dpqp, fermi2_qpdp
  end interface fermimhalf

  interface polygamma2fermi
    module procedure polygamma2fermi_dp, polygamma2fermi_qp
  end interface polygamma2fermi

  interface polygamma2psi1
    module procedure polygamma2psi1_dp, polygamma2psi1_qp
  end interface polygamma2psi1

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
  function polygamma2psi1_dp(gamma,eps,beta)
    implicit none
    real(8), intent(in) :: gamma,eps,beta
    complex(8), external :: wpsipg
    real(8) :: polygamma2psi1_dp
    polygamma2psi1_dp = real(wpsipg(0.5d0 + beta/2.d0/pi * (gamma + ci*eps),1))
  end function polygamma2psi1_dp

  function polygamma2psi1_qp(gamma,eps,beta)
    implicit none
    real(16), intent(in) :: eps,beta
    real(8), intent(in)  :: gamma
    complex(16), external :: wpsipghp
    real(16) :: polygamma2psi1_qp
    polygamma2psi1_qp = real(wpsipghp(0.5q0 + beta/2.q0/piQ * (gamma + ciQ*eps),1))
  end function polygamma2psi1_qp

end module Mfermi
