! derived types for the main program
! this should be completely agnostic to the type of data we provide

module Mtypes
  implicit none

  contains

  type algorithm
    logical :: lDebug         ! evaluate quad precision data?
    logical :: lBfield        ! calculations in the presence of a magnetic field
                              ! this requires the existance of de(k)/dk d2e(k)/dk2
    logical :: lDerivatives   ! we have the derivatives
    integer :: rootmethod     ! numerical method to find the chemical potential
    integer :: mumethod       ! 0: find_mu
                              ! 1: constant mu
                              ! 2: fixed mu for each temperature
                              ! __used to compare pure Boltzmann with Boltzmann and Kubo mu
    character(len=256) :: input_bands
    character(len=256) :: input_scattering
  end type

  type lattice
    real(8) :: vol            ! volume of the real space unit cell -- for prefactors
    logical :: lortho         ! do we have an orthogonal bravais lattice (cubic,tetragonal,orthorhombic)
    integer :: nalpha         ! number of polarization directions ... 3 or 6
  end type

  ! ( 1 4 6 )
  ! ( - 2 5 )
  ! ( - - 3 )

  type kpointmesh
    integer, allocatable :: multiplicity(:)
    real(8), allocatable :: weight(:)
    integer :: ktot
  end type

  type energydisp
    integer :: nband_max
    integer :: nbopt_min                     ! number of bands (interval) included in the optical matrix elements
    integer :: nbopt_max                     ! number of bands (interval) included in the optical matrix elements
    real(8) :: efer                          ! Fermi energy
    real(8) :: nelect                        ! number of electrons (provided from inputfile at first iteration)
    real(8) :: band_fill_value               ! if we get incomplete bands (i.e. from wien2k)
    real(8), allocatable    :: band(:,:)        ! energy(nband,ik)
    real(8), allocatable    :: band_dk(:,:,:)   ! d/dk_i band(nband,ik)
    real(8), allocatable    :: band_d2k(:,:,:)  ! d2/(dk_i dk_j) band(nband,ik)
    complex(8), allocatable :: Mopt(:,:,:,:)    ! M(xy,n,n',k)= <n,k|p.e_x|n',k> * <n',k|p.e_y|n,k> *
    real(8), allocatable    :: band_shift(:,:)
  end type

  type dosgrid
    integer :: nnrg                ! number of points in the energy window
    real(8) :: emin                ! bottom of the energy window
    real(8) :: emax                ! top of the energy window
    real(8) :: vbm                 ! valence band maximum
    real(8) :: cbm                 ! conduction band minimum
    real(8) :: gap                 ! band gap
    real(8), allocatable :: enrg(:)! energy grid
    real(8), allocatable :: dos(:) ! density of states (as computed in PRB,49,16223, appx C )
    real(8), allocatable :: nos(:) ! number  of states (as computed in PRB,49,16223, appx A )
  end type

  type scattering
    ! temperature grid
    integer :: nT                  ! number of points in the temperature window
    real(8) :: Tmin                ! bottom of the temperature window
    real(8) :: Tmax                ! top of the temperature window
    real(8) :: dT                  ! temperature interval
    real(8), allocatable :: TT(:)  ! temperature grid
    real(8), allocatable :: beta(:)

    ! temperature tempendent quantities
    real(8), allocatable :: mu(:)  ! chemical potential (temperature dependent)
    real(8), allocatable :: d1(:)  ! square of the 1st derivative of sigma
    real(8), allocatable :: d2(:)  ! product of the 2nd derivative times sigma
    real(8), allocatable :: d0(:)  ! linear combination of the two above, whose zero corresponds to T*
    real(8) :: Tstar               ! temperature for which (d^2 rho)/(d beta^2)=0
                                   ! in practice it is the T for which (d^2 sigma)/(d beta^2) changes sign
    real(8) :: Tflat               ! T for which (d sigma)/(d beta) changes sign (onset of saturation)
    real(8), allocatable :: gam(:,:) ! n, k
    real(8), allocatable :: zqp(:,:) ! n, k
  end type

  type response_dp
    !kernels
    real(8) ::  s_ker  ! for conductivity (s)igma
    real(8) ::  sB_ker ! for conductivity in B-field
    real(8) ::  a_ker  ! for Peltier (a)lpha
    real(8) ::  aB_ker ! for Peltier in B-field

    ! band-resolved response functions
    real(8), allocatable :: s_full(:,:,:)   ! [3 or 6], nband, nk -- for conductivity
    real(8), allocatable :: sB_full(:,:,:)  ! for conductivity in B-field
    real(8), allocatable :: a_full(:,:,:)   ! for Peltier
    real(8), allocatable :: aB_full(:,:,:)  ! for Peltier in B-field

    ! global k-summation
    real(8), allocatable :: s_gather(:,:)   ! for conductivity
    real(8), allocatable :: sB_gather(:,:)  ! for conductivity in B-field
    real(8), allocatable :: a_gather(:,:)   ! for Peltier
    real(8), allocatable :: aB_gather(:,:)  ! for Peltier in B-field

    ! local k-summation
    real(8), allocatable :: s_local(:,:)
    real(8), allocatable :: sB_local(:,:)
    real(8), allocatable :: a_local(:,:)
    real(8), allocatable :: aB_local(:,:)

    ! total band and k-summation
    real(8) :: s_tot(:)
    real(8) :: sB_tot(:)
    real(8) :: a_tot(:)
    real(8) :: aB_tot(:)

    ! derived quantities
    real(8) :: Seebeck(:)
    real(8) :: Nernst(:)
    real(8) :: RH(:)
  end type

  type response_qp
    !kernels
    real(16) ::  s_ker  ! for conductivity
    real(16) ::  sB_ker ! for conductivity in B-field
    real(16) ::  a_ker  ! for Peltier
    real(16) ::  aB_ker ! for Peltier in B-field

    ! band-resolved response functions
    real(16), allocatable :: s_full(:,:,:)   ! nk,nband,3,3 for conductivity
    real(16), allocatable :: sB_full(:,:,:)  ! nk,nband,3,3 for conductivity in B-field
    real(16), allocatable :: a_full(:,:,:)   ! nk,nband,3,3 for Peltier
    real(16), allocatable :: aB_full(:,:,:)  ! nk,nband,3,3 for Peltier in B-field

    ! global k-summation
    real(16), allocatable :: s_gather(:,:)   ! nband,3,3 for conductivity
    real(16), allocatable :: sB_gather(:,:)  ! nband,3,3 for conductivity in B-field
    real(16), allocatable :: a_gather(:,:)   ! nband,3,3 for Peltier
    real(16), allocatable :: aB_gather(:,:)  ! nband,3,3 for Peltier in B-field

    ! local k-summation
    real(16), allocatable :: s_local(:,:)
    real(16), allocatable :: sB_local(:,:)
    real(16), allocatable :: a_local(:,:)
    real(16), allocatable :: aB_local(:,:)

    ! total band and k-summation
    real(16) :: s_tot(:)
    real(16) :: sB_tot(:)
    real(16) :: a_tot(:)
    real(16) :: aB_tot(:)

    ! derived quantities
    real(16) :: Seebeck(:)
    real(16) :: Nernst(:)
    real(16) :: RH(:)
  end type

end module Mtypes
