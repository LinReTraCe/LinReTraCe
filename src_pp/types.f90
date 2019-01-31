module Mtypes

  type algorithm
    logical :: ltbind         ! tight binding lattice?
    logical :: lw2k           ! use Wien2k input
    logical :: lvasp          ! use Vasp input
    character(256) :: mysyst  ! label that is used to open the Wien2k files
  end type

  type lattice
    real(8) :: a(3)  ! direct lattice vectors
    real(8) :: angle(3)
    real(8) :: vol   ! volume of the lattice
    logical :: lortho ! extra switch -> if the crystal is orthogonal our matrices reduce in size!
    integer :: nalpha ! number of polarization directions (1 for cubic; 3 for everything else)
    integer :: ibravais
    integer :: spacegroup
    ! order for ibravais:
    ! 1: primitive cubic
    ! 2: body-centered cubic
    ! 3: face-centered cubic
    ! 4: primitive tetragonal
    ! 5: body-contered tetragonal
    ! 6: primitive orthorhombic
    ! 7: base-centered orthorhombic
    ! 8: body-centered orthorhombic
    ! 9: face-centered orthorhombic
    ! 10: primitive hexagonal
    ! 11: primitive rhombohedral
    ! 12: primitive monoclinic
    ! 13: base-centered monoclinic
    ! 14: primitive triclinic
  end type

  type kpointmesh
    integer, allocatable :: multiplicity(:)
    real(8), allocatable :: weight(:)
    integer :: kx                        ! number of reducible k-points in each cartesian direction
    integer :: ky                        ! corresponds to a k value of [0, 1)
    integer :: kz
    integer :: ktot                      ! total number of k-points from read-in (either reducible or irreducible)
    integer :: kred                      ! total number of REDUCIBLE k-points
  end type

  type energydisp
    integer :: nband_max
    integer :: nbopt_min                     ! number of bands (interval) included in the optical matrix elements
    integer :: nbopt_max                     ! number of bands (interval) included in the optical matrix elements
    real(8) :: efer                          ! Fermi energy
    real(8) :: nelect                        ! number of electrons (provided from inputfile at first iteration)
    real(8), allocatable :: band(:,:)        ! band(ik,nband)
    real(8), allocatable :: band_dk(:,:)
    real(8), allocatable :: band_d2k(:,:)
    real(8), allocatable :: Mopt(:,:,:,:)    ! M(x,k,n',n)= | <n',k|p.e_x|n,k> |^2
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

end module Mtypes
