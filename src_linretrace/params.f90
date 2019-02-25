module Mparams
  implicit none

  ! mathematical constants
  complex(8),  parameter :: ci      = (0.d0,1.d0)
  real(8),     parameter :: pi      = 3.1415926535897932384626433832795028841971693993751D0

  complex(16), parameter :: ciQ     = (0.q0,1.q0)
  real(16),    parameter :: piQ     = 3.1415926535897932384626433832795028841971693993751Q0

  ! physical constants
  real(8),     parameter :: kB      = 8.6173324D-5 ! eV/K
  real(8),     parameter :: hbarevs = (4.135667516D-15)/(2.d0*pi) !evs
  real(8),     parameter :: hbarjs  = (6.62607D-34)/(2.d0*pi) !Js
  real(8),     parameter :: echarge = 1.6021766208d-19 ! C

  ! numerical parameters
  real(8),     parameter :: ndev    = 1D-13
  real(16),    parameter :: ndevQ   = 5Q-14
  real(16),    parameter :: ndevVQ  = 5Q-15
  integer,     parameter :: niit    = 100
  integer,     parameter :: niitQ   = 200
  real(8),     parameter :: small   = 1D-11
  real(16),    parameter :: smallQ  = 1Q-18

  ! fortran specifier
  integer,     parameter :: stdout = 6
  integer,     parameter :: stdin  = 5
  integer,     parameter :: stderr = 0

end module Mparams
