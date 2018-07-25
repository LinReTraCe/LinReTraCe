

!module gfortran_fix ! gfortran doesn't understand iso_fortran_env properly and screws up quad precision...
!integer dp,qp
!parameter(dp=8)
!parameter(qp=16)
!end module gfortran_fix

module params



integer dp,qp
!parameter(dp=8)
!parameter(qp=16)
!parameter(dp=selected_real_kind(8)) !for HCLM
!parameter(qp=selected_real_kind(16))



  !integer nband, nk ! #bands
  real(8) NE, Tmax,Tmin,dT ! #electrons, Temperature max, min, T-interval
  real(8) kB,T,beta,beta2p,pi,ndev,small,threshold,ndevhp,hbarevs,hbarjs,echarge
  real(16) thresholdQ,smallQ,betaQ,beta2pQ,piQ,ndevQ,ndevVQ,kBQ
  integer niit,niitQ
  complex(8) ci
  complex(16) ciQ
  parameter(kB=8.6173324D-5) ! eV/K
  parameter(kBQ=8.6173324Q-5) ! eV/K
  parameter(pi=3.1415926535897932384626433832795028841971693993751D0)
  parameter(piQ=3.1415926535897932384626433832795028841971693993751Q0)
  parameter(hbarevs=(4.135667516D-15)/(2.d0*3.1415926535897932385d0)) !evs
  parameter(hbarjs=(6.62607D-34)/(2.d0*3.1415926535897932385d0)) !Js
  parameter(echarge=1.6021766208d-19) ! C
!complex unit
  parameter(ci=(0.d0,1.d0))
  parameter(ciQ=(0.q0,1.q0))
! numerical parameters 
  parameter(ndev=5D-14) ! ~1D-12 allowed deviation from set particle number in mu-root finding.
  parameter(ndevhp=5D-15) ! when using find_muDPQ
  parameter(ndevQ=5Q-18)  ! when using full QUAD in find_muQ
  parameter(ndevVQ=1Q-21) ! 
  parameter(niit=100) ! maximal number of secant steps
  parameter(niitQ=150) ! maximal number of secant steps
  parameter(small=1D-11) ! ~1D-11 for real(8): use 10 significant digits... use 20 for QUAD
  parameter(smallQ=1Q-21) ! ~1Q-18 doesnt seem to matter for QUAD...
end module params
