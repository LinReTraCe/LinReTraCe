module Mtypes

  type algorithm
    logical :: ldebug       ! evaluate quad precision data?
    logical :: ltbind       ! tight binding lattice?
    logical :: ltetra       ! use tetrahedron method?
    logical :: lw2k         ! use Wien2k input
    logical :: loptic       ! use also optic matrix elements
    logical :: lBfield      ! calculations in the presence of a magnetic field?
    logical :: ldmft        ! k-point and band dependent renormalisation factor and scattering rate provided
    logical :: lsymm        ! reads in the full BZ or generates it using symmetry operations?
    logical :: lpreproc     ! use preprocessed data
    integer :: imurestart   ! restart with a privided value of the chemical potential?
    character(100) :: mysyst ! label that is used to open the Wien2k files
  end type

  type lattice
    real(8) :: alat  ! lattice constant
    real(8) :: a(3)  ! direct lattice vectors in units of ALAT. needed for derivatives
    real(8) :: vol   ! volume of the lattice
    integer :: nalpha ! number of polarization directions (1 for cubic; 3 for everything else)
    logical :: lcubic
    logical :: ltetragonal
    logical :: lorthorhombic
    ! not used at the moment
    logical :: lhexagonal
    logical :: lrhombohedral
    logical :: lmonoclinic
    logical :: ltriclinic
  end type

  type kpointmesh
    integer, allocatable :: k_id(:,:,:)  ! counter idenfifying a specific k-point through values of ikx, iky, ikz
    real(8), allocatable :: k_coord(:,:) ! coordinates in reciprocal space associated to a k-point (exclunding the BZ endpoints)
    integer :: kx                        ! number of k-points in each cartesian direction
    integer :: ky
    integer :: kz
    integer :: ktot                      ! total number of k-points (which varies depending on the sampling method for the BZ)
  end type

  type edisp
    integer :: nband_max
    integer :: nbopt_min                     ! number of bands (interval) included in the optical matrix elements
    integer :: nbopt_max                     ! number of bands (interval) included in the optical matrix elements
    real(8) :: efer                          ! Fermi energy
    real(8) :: nelect                        ! number of electrons (provided from inputfile at first iteration)
    real(8) :: occ_tot                       ! number of electrons computed for the updated value of chemical potential (non-eq value)
    real(8) :: ztmp                          ! (constant) renormalisation factor -only used if algo%ldmft=.false.-
    real(8), allocatable :: band(:,:)        ! band(k_id,nband) contains nband_max dispersion curves, if for some k-points
                                             ! fewer bands are computed, then the missing elements are set to 1000 (arbitrarily large value)
    real(8), allocatable :: Mopt(:,:,:,:)    ! M(x,k,n',n)= | <n',k|p.e_x|n,k> |^2
    real(8), allocatable :: Z(:,:)           ! renormalisation factor 1/z= 1 - (d/dw)S_{k,n}(w) @ w=0
    real(8), allocatable :: Im(:,:)          ! =-imaginary part of the self-energy at zero frequency
    real(8), allocatable :: Mopt_tetra(:,:,:,:)    ! M(x,t,n',n)= sum_k(j) {w(j)*|<n',k(j)|p.e_x|n,k(j)>|^2 } with w(j) the tetrahedron weights

    real(8)              :: tmax             ! TIGHT BINDING PARAMETER: number of hopping parameter -- currently not implemented properly
    real(8), allocatable :: a(:,:)           ! TIGHT BINDING PARAMETER: lattice spacing to the higher order hoppings (3, tmax)
    real(8), allocatable :: E0(:)            ! TIGHT BINDING PARAMETER: band energy at Gamma point (nbands)
    real(8), allocatable :: t(:,:)           ! TIGHT BINDING PARAMETER: hopping parameter (~band width) (nbands,tmax)
    real(8), allocatable :: M2(:,:,:,:)      ! TIGHT BINDING PARAMETER: 2nd derivative of the energy dispersion (idir,idir2,nk,nbands)
  end type

  type symop
    integer, allocatable :: symop_id(:,:)   ! symop_id(2,redkp) this counter tells me if for a given reducible k-point the corresponding
                                            ! irreducible kpoint (1) and the required symmetry operation (2)
                                            ! produces a new element or a redundant one (1 or 0 in the 1st entry)
    real(8), allocatable :: Msym(:,:,:)     ! Msym(3,3,nsym) matrix containing the 3x3 symmetry transformations
    real(8), allocatable :: Tras(:,:)       ! Tras(3,nsym) matrix containing additional lattice traslations for non-symmorphic groups
    integer :: nsym                         ! number of symmetry operations
    logical :: lnsymmr                      ! .true. for nonsymmorphic space groups
    character(3) :: cntr
  end type

  type tetramesh
    integer :: ntet                            ! number of inequivalent tetrahedra found starting from an equally spaced k-grid
    integer, allocatable :: idtet(:,:)         ! dimension (0:4,ntet); idtet(0,ntet) is the tetrahedron's multiplicity, 1:4 are the identifiers of the vertices
    real(8) :: vltot                  ! vltot is the total tetrahedron's volume (in reciprocal space)
    real(8), allocatable :: vltet(:)  ! dimension (ntet); vltet(itet) is the tetrahedron's volume (in reciprocal space)
    real(8), allocatable :: nocc(:)   ! dimension (nbands); number of occupied tetrahedra at a given temperature
  end type

  type dosgrid
    integer :: nnrg                         ! number of points in the energy window
    real(8) :: emin                ! bottom of the energy window
    real(8) :: emax                ! top of the energy window
    real(8) :: vbm                 ! valence band maximum
    real(8) :: cbm                 ! conduction band minimum
    real(8) :: gap                 ! band gap
    real(8), allocatable :: enrg(:)! energy grid
    real(8), allocatable :: dos(:) ! density of states (as computed in PRB,49,16223, appx C )
    real(8), allocatable :: nos(:) ! number  of states (as computed in PRB,49,16223, appx A )
  end type

  type scatrate
    integer :: nT                           ! number of points in the temperature window
    real(8) :: Tmin                ! bottom of the temperature window
    real(8) :: Tmax                ! top of the temperature window
    real(8) :: dT                  ! temperature interval
    real(8), allocatable :: TT(:)  ! temperature grid
    real(8), allocatable :: mu(:)  ! chemical potential (temperature dependent)
    real(8), allocatable :: d1(:)  ! square of the 1st derivative of sigma
    real(8), allocatable :: d2(:)  ! product of the 2nd derivative times sigma
    real(8), allocatable :: d0(:)  ! linear combination of the two above, whose zero corresponds to T*
    real(8) :: Tstar               ! temperature for which (d^2 rho)/(d beta^2)=0
                                            ! in practice it is the T for which (d^2 sigma)/(d beta^2) changes sign
    real(8) :: Tflat               ! T for which (d sigma)/(d beta) changes sign (onset of saturation)
    !gamma
    integer :: ng                           ! degree of the polynomial temperature dependence of gamma
    real(8), allocatable :: gc(:)  ! coefficient of the polynomial temperature dependence of gamma
    real(8), allocatable :: gam(:) ! temparature dependent scattering rate (extrinsic scattering events), defects
    real(8), allocatable :: ykb(:,:,:) ! intrinsic scattering rate ykb(T,k,band)

  end type

  type dp_resp
    !kernels
    real(8) ::  s_ker  ! for conductivity
    real(8) ::  sB_ker ! for conductivity in B-field
    real(8) ::  a_ker  ! for Peltier
    real(8) ::  aB_ker ! for Peltier in B-field
    !response functions...
    real(8), allocatable :: s_tmp(:,:,:,:)   ! nk,nband,3,3 for conductivity
    real(8), allocatable :: sB_tmp(:,:,:,:)  ! nk,nband,3,3 for conductivity in B-field
    real(8), allocatable :: a_tmp(:,:,:,:)   ! nk,nband,3,3 for Peltier
    real(8), allocatable :: aB_tmp(:,:,:,:)  ! nk,nband,3,3 for Peltier in B-field
    real(8), allocatable :: s(:,:,:)         ! nband,3,3 for conductivity
    real(8), allocatable :: sB(:,:,:)        ! nband,3,3 for conductivity in B-field
    real(8), allocatable :: a(:,:,:)         ! nband,3,3 for Peltier
    real(8), allocatable :: aB(:,:,:)        ! nband,3,3 for Peltier in B-field
    real(8), allocatable :: s_local(:,:,:)
    real(8), allocatable :: sB_local(:,:,:)
    real(8), allocatable :: a_local(:,:,:)
    real(8), allocatable :: aB_local(:,:,:)

    real(8) :: s_tot(3,3),  s_tet(3,3)
    real(8) :: sB_tot(3,3), sB_tet(3,3)
    real(8) :: a_tot(3,3),  a_tet(3,3)
    real(8) :: aB_tot(3,3), aB_tet(3,3)
    real(8) :: Seebeck(3),Nernst(3),RH(3)

    real(8) :: RePolyGamma(0:4),ImPolyGamma(0:4),gamma,aqp,z,tmp
    complex(8) :: ctmp,zarg
  end type

  !interband transitions functions
  type, extends(dp_resp) :: dp_respinter
    !variables for the second band
    real(8) :: RePolyGamma1(0:4),ImPolyGamma1(0:4)
    real(8) :: RePolyGamma2(0:4),ImPolyGamma2(0:4)
    real(8) :: gamma1,gamma2, aqp1,aqp2, z1,z2
  end type

  type qp_resp
    !kernels
    real(16) ::  s_ker  ! for conductivity
    real(16) ::  sB_ker ! nband,3,3 for conductivity in B-field
    real(16) ::  a_ker  ! nband,3,3 for Peltier
    real(16) ::  aB_ker ! nband,3,3 for Peltier in B-field
    !response functions...
    real(16), allocatable :: s_tmp(:,:,:,:)   ! nk,nband,3,3 for conductivity
    real(16), allocatable :: sB_tmp(:,:,:,:)  ! nk,nband,3,3 for conductivity in B-field
    real(16), allocatable :: a_tmp(:,:,:,:)   ! nk,nband,3,3 for Peltier
    real(16), allocatable :: aB_tmp(:,:,:,:)  ! nk,nband,3,3 for Peltier in B-field
    real(16), allocatable :: s(:,:,:)         ! nband,3,3 for conductivity
    real(16), allocatable :: sB(:,:,:)        ! nband,3,3 for conductivity in B-field
    real(16), allocatable :: a(:,:,:)         ! nband,3,3 for Peltier
    real(16), allocatable :: aB(:,:,:)        ! nband,3,3 for Peltier in B-field
    real(16), allocatable :: s_local(:,:,:)
    real(16), allocatable :: sB_local(:,:,:)
    real(16), allocatable :: a_local(:,:,:)
    real(16), allocatable :: aB_local(:,:,:)

    real(16) :: s_tot(3,3),  s_tet(3,3)
    real(16) :: sB_tot(3,3), sB_tet(3,3)
    real(16) :: a_tot(3,3),  a_tet(3,3)
    real(16) :: aB_tot(3,3), aB_tet(3,3)
    real(16) :: Seebeck(3),Nernst(3),RH(3),Nernstpart(2) ! nernstpart is axy sxx bzw axx sxy / sxx^2

    real(16) :: RePolyGamma(0:4),ImPolyGamma(0:4),gamma,aqp,z,tmp
    complex(16) :: ctmp,zarg
  end type

  type, extends(qp_resp) :: qp_respinter
    real(16) :: RePolyGamma1(0:4),ImPolyGamma1(0:4)
    real(16) :: RePolyGamma2(0:4),ImPolyGamma2(0:4)
    real(16) :: gamma1,gamma2, aqp1,aqp2, z1,z2
  end type

end module Mtypes
