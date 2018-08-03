module Mtypes

  type algorithm
    logical :: ldebug       ! evaluate quad precision data?
    logical :: ltbind       ! tight binding lattice?
    logical :: ltetra       ! use tetrahedron method?
    logical :: lw2k         ! use Wien2k input
    logical :: loptic       ! use also optic matrix elements
    logical :: lBfield      ! calculations in the presence of a magnetic field?
    logical :: lsymm        ! reads in the full BZ or generates it using symmetry operations?
    integer :: imurestart   ! restart with a privided value of the chemical potential?
    character(10) :: mysyst ! label that is used to open the Wien2k files
  end type

  type kpointmesh
    integer, allocatable :: k_id(:,:,:)           ! counter idenfifying a specific k-point through values of ikx, iky, ikz
    double precision, allocatable :: k_coord(:,:) ! coordinates in reciprocal space associated to a k-point (exclunding the BZ endpoints)
    double precision, allocatable :: k_weight(:)  ! k-point weights
    double precision, allocatable :: mult(:)      ! k-point multiplicity
    integer :: kx                                 ! number of k-points in each cartesian direction
    integer :: ky
    integer :: kz
    integer :: ktot                               ! total number of k-points (which varies depending on the sampling method for the BZ)
    double precision :: alat                      ! lattice constant
    double precision :: a(3)                      ! direct lattice vectors in units of ALAT. needed for derivatives --> units of vk's... ???
    integer :: nkeq                               ! number of direction with the same number of k-point (nkeq=3 uniform k-mesh)
  end type

  type edisp
    double precision, allocatable :: E0(:)            ! TIGHT BINDING PARAMETER: band energy at Gamma point (nbands)
    double precision, allocatable :: t(:,:)           ! TIGHT BINDING PARAMETER: band width (nbands,tmax), tmax/=1 is used to set anisotropic dispersion
    double precision, allocatable :: M2(:,:,:,:)      ! TIGHT BINDING PARAMETER: 2nd derivative of the energy dispersion (idir,idir2,nk,nbands), not used if ltbind=F
    integer :: nband_max
    integer :: nbopt_min                              ! number of bands (interval) included in the optical matrix elements
    integer :: nbopt_max                              ! number of bands (interval) included in the optical matrix elements
    double precision :: efer                          ! Fermi energy
    double precision :: nelect                        ! number of electrons (provided from inputfile at first iteration)
    double precision :: occ_tot                       ! number of electrons computed for the updated value of chemical potential (non-eq value)
    double precision, allocatable :: band(:,:)        ! band(k_id,nband) contains nband_max dispersion curves, if for some k-points
                                                      ! fewer bands are computed, then the missing elements are set to 1000 (arbitrarily large value)
    double precision, allocatable :: occ(:,:)         ! occ(k_id,nband) contains the occupation numbers
    double precision, allocatable :: Mopt(:,:,:,:)    ! M(x,k,n',n)= | <n',k|p.e_x|n,k> |^2
    double precision, allocatable :: Mopt_tetra(:,:,:,:)    ! M(x,t,n',n)= sum_k(j) {w(j)*|<n',k(j)|p.e_x|n,k(j)>|^2 } with w(j) the tetrahedron weights
  end type

  type symop
    integer, allocatable :: symop_id(:,:)        ! symop_id(3,irrkp*nsym) this counter tells me if for a given k-point (2nd entry) a given symmetry operation (3rd entry)
                                                 ! produces a new element or a redundant one (1 or 0 in the 1st entry)
    integer :: nsym                              ! number of symmetry operations
    double precision, allocatable :: Msym(:,:,:) ! Msym(3,3,nsym) matrix containing the 3x3 symmetry transformations
    double precision, allocatable :: Tras(:,:)   ! Tras(3,nsym) matrix containing additional lattice traslations for non-symmorphic groups
    logical :: lcubic
    logical :: ltetrag
    logical :: lorthor
    logical :: lnsymmr                           ! .true. for nonsymmorphic space groups
    logical :: lBZpnt
  end type

  type tetramesh
    integer :: ntet                            ! number of inequivalent tetrahedra found starting from an equally spaced k-grid
    integer, allocatable :: idtet(:,:)         ! dimension (0:4,ntet); idtet(0,ntet) is the tetrahedron's multiplicity, 1:4 are the identifiers of the vertices
    double precision :: vltot                  ! vltot is the total tetrahedron's volume (in reciprocal space)
    double precision, allocatable :: vltet(:)  ! dimension (ntet); vltet(itet) is the tetrahedron's volume (in reciprocal space)
    double precision, allocatable :: nocc(:)   ! dimension (nbands); number of occupied tetrahedra at a given temperature
  end type

  type dosgrid
    integer :: nnrg                         ! number of points in the energy window
    double precision :: emin                ! bottom of the energy window
    double precision :: emax                ! top of the energy window
    double precision :: vbm                 ! valence band maximum
    double precision :: cbm                 ! conduction band minimum
    double precision :: gap                 ! band gap
    double precision, allocatable :: enrg(:)! energy grid
    double precision, allocatable :: dos(:) ! density of states (as computed in PRB,49,16223, appx C )
    double precision, allocatable :: nos(:) ! number  of states (as computed in PRB,49,16223, appx A )
  end type

  type scatrate
    !temperature grid
    integer :: nT                           ! number of points in the temperature window
    double precision :: Tmin                ! bottom of the temperature window
    double precision :: Tmax                ! top of the temperature window
    double precision :: dT                  ! temperature interval
    double precision, allocatable :: TT(:)  ! temperature grid
    double precision, allocatable :: mu(:)  ! chemical potential (temperature dependent)
    !double precision, allocatable :: drhodT(:) ! 1st derivative of rho w.r.t. T
    double precision, allocatable :: d1(:)  ! square of the 1st derivative of sigma
    double precision, allocatable :: d2(:)  ! product of the 2nd derivative times sigma
    double precision, allocatable :: d0(:)  ! linear combination of the two above, whose zero corresponds to T*
    double precision :: Tstar               ! temperature for which (d^2 rho)/(d beta^2)=0
                                            ! in practice it is the T for which (d^2 sigma)/(d beta^2) changes sign
    double precision :: Tflat               ! T for which (d sigma)/(d beta) changes sign (onset of saturation)
    double precision :: z                   ! renormalisation factor (read from file)
    !gamma
    integer :: ng                           ! degree of the polynomial temperature dependence of gamma
    double precision, allocatable :: gc(:)  ! coefficient of the polynomial temperature dependence of gamma
    double precision, allocatable :: gam(:,:)  ! temparature and band dependent scattering rate
  end type

  ! DOUBLE PRECISION
  type dp_resp
    !integer,parameter :: iq=8
    !integer, kind :: iq=8
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
    real(8), allocatable :: s_local(:,:,:)   ! used only in the mpi version
    real(8), allocatable :: sB_local(:,:,:)  !
    real(8), allocatable :: a_local(:,:,:)   !
    real(8), allocatable :: aB_local(:,:,:)  !

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

  ! QUAD PRECISION
  type qp_resp
    !integer,parameter :: iq=16
    !integer, kind :: iq=16
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
    real(16), allocatable :: s_local(:,:,:)   ! used only in the mpi version
    real(16), allocatable :: sB_local(:,:,:)  !
    real(16), allocatable :: a_local(:,:,:)   !
    real(16), allocatable :: aB_local(:,:,:)  !

    real(16) :: s_tot(3,3),  s_tet(3,3)
    real(16) :: sB_tot(3,3), sB_tet(3,3)
    real(16) :: a_tot(3,3),  a_tet(3,3)
    real(16) :: aB_tot(3,3), aB_tet(3,3)
    real(16) :: Seebeck(3),Nernst(3),RH(3),Nernstpart(2) ! nernstpart is axy sxx bzw axx sxy / sxx^2

    real(16) :: RePolyGamma(0:4),ImPolyGamma(0:4),gamma,aqp,z,tmp
    complex(16) :: ctmp,zarg
  end type

end module Mtypes
