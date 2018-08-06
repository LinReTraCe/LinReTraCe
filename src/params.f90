module Mparams
  implicit none

  real(8)  :: NE ! #electrons
  real(8)  :: Tmin,Tmax,dT,T ! Temperature min, max, interval
  real(8)  :: beta, beta2p, threshold
  real(16) :: betaQ,beta2pQ,thresholdQ

  ! mathematical constants
  complex(8), parameter  :: ci = (0.d0,1.d0)
  real(8), parameter     :: pi=3.1415926535897932384626433832795028841971693993751D0

  complex(16), parameter :: ciQ = (0.q0,1.q0)
  real(16), parameter    :: piQ=3.1415926535897932384626433832795028841971693993751Q0

  ! physical constants
  real(8), parameter     :: kB=8.6173324D-5 ! eV/K
  real(8), parameter     :: hbarevs=(4.135667516D-15)/(2.d0*3.1415926535897932385d0) !evs
  real(8), parameter     :: hbarjs=(6.62607D-34)/(2.d0*3.1415926535897932385d0) !Js
  real(8), parameter     :: echarge=1.6021766208d-19 ! C

  real(16), parameter    :: kBQ=8.6173324Q-5 ! eV/K

  ! numerical parameters
  real(8), parameter     :: ndev=5D-14 ! ~1D-12 allowed deviation from set particle number in mu-root finding.
  real(8), parameter     :: ndevhp=5D-15 ! when using find_muDPQ
  real(16), parameter    :: ndevQ=5Q-18  ! when using full QUAD in find_muQ
  real(16), parameter    :: ndevVQ=1Q-21 !
  integer, parameter     :: niit=100 ! maximal number of secant steps
  integer, parameter     :: niitQ=150 ! maximal number of secant steps
  real(8), parameter     :: small=1D-11 ! ~1D-11 for real(8): use 10 significant digits... use 20 for QUAD
  real(16), parameter    :: smallQ=1Q-21 ! ~1Q-18 doesnt seem to matter for QUAD...
end module Mparams
