! derived types for the main program
! this should be completely agnostic to the type of data we provide

module Mtypes
  implicit none

  ! contains the methods employed for the root-finding and
  ! the way we treat the chemical potential
  ! also contains the file names of the input files
  type algorithm
    logical :: lDebug         ! debug mode --> evaluate quad precision data?
    logical :: lBfield        ! calculations in the presence of a magnetic field
                              ! this requires the existance of the band derivatives
    integer :: rootMethod     ! numerical method to find the chemical potential
    logical :: muSearch
    logical :: muFermi
    logical :: lScatteringFile
    character(len=256) :: input_energies
    character(len=256) :: input_scattering
    character(len=256) :: output_file
  end type

  ! lattice information which is necessary for us
  ! volume (prefactors) and the orthogonality of the crystal
  ! note: an orthogonal crystal reduces the effective number of polarization directions
  type lattice
    real(8) :: vol            ! volume of the real space unit cell -- for prefactors
    logical :: lOrtho         ! do we have an orthogonal bravais lattice (cubic,tetragonal,orthorhombic)
    integer :: nalpha         ! number of polarization directions ... 3 or 6
  end type

  ! ( 1 4 6 )
  ! ( - 2 5 )
  ! ( - - 3 )

  ! information about the k-points which is necessary for us
  ! that is: number of k-points and its weight
  type kpointmesh
    real(8)              :: weightsum
    real(8), allocatable :: weight(:)
    integer              :: nkp
  end type

  ! energy dispersion and derived quantities
  ! direct derivatives and transition elements
  type energydisp
    integer :: nband_max
    integer :: nbopt_min                     ! number of bands (interval) included in the optical matrix elements
    integer :: nbopt_max                     ! number of bands (interval) included in the optical matrix elements
    integer :: iSpin                         ! number of spins
    logical :: lDerivatives                  ! do we have the derivatives (band_dk, band_d2k)
    logical :: lBandShift
    integer :: iOptical
    real(8) :: efer                          ! Fermi energy
    real(8) :: mu
    real(8) :: nelect

    real(8), allocatable    :: band_original(:,:,:)
    real(8), allocatable    :: band_shift(:,:,:)  ! same as band
    real(8), allocatable    :: band(:,:,:)        ! energy(nband,ik,ispin)
    real(8), allocatable    :: band_dk(:,:,:,:)   ! d/dk_i band(nband,ik,ispin)
    real(8), allocatable    :: band_d2k(:,:,:,:)  ! d2/(dk_i dk_j) band(nband,ik,ispin)

    ! optical elements (because of the double band dependencies)
    ! are loaded for each k-point and each spin
    real(8), allocatable    :: Mopt(:,:,:,:)    ! M(xy,n,n')= <n,k|p.e_x|n',k> * <n',k|p.e_y|n,k> *
                                                ! 3..9, nband,nband
  end type

  type dosgrid
    integer :: nnrg                ! number of points in the energy window
    real(8) :: emin                ! bottom of the energy window
    real(8) :: emax                ! top of the energy window
    real(8), allocatable :: vbm(:)                 ! valence band maximum
    real(8), allocatable :: cbm(:)                 ! conduction band minimum
    real(8), allocatable :: gap(:)                 ! band gap
    real(8), allocatable :: enrg(:)  ! energy grid
    real(8), allocatable :: dos(:,:) ! density of states (as computed in PRB,49,16223, appx C )
    real(8), allocatable :: nos(:,:) ! number  of states (as computed in PRB,49,16223, appx A )
  end type

  type temperature
    integer :: nT                  ! number of points in the temperature window
    real(8) :: Tmin                ! bottom of the temperature window
    real(8) :: Tmax                ! top of the temperature window
    real(8) :: dT                  ! temperature spacing
    real(8), allocatable :: TT(:)  ! temperature grid [K]
    real(8), allocatable :: beta(:)! inverse temperature grid [eV]
  end type

  type runinfo
    integer  :: iT
    real(8)  :: Temp
    real(8)  :: beta
    real(8)  :: beta2p

    real(16) :: TempQ
    real(16) :: betaQ
    real(16) :: beta2pQ

    integer  :: ik
  end type

  type scattering
    ! scattering rates and quasiparticle weights
    real(8), allocatable :: gamcoeff(:)
    real(8), allocatable :: zqpcoeff(:)
    real(8), allocatable :: gam(:,:,:) ! same as band
    real(8), allocatable :: zqp(:,:,:) ! same as band
    real(8)              :: gamscalar
    real(8)              :: zqpscalar
    real(8)              :: gamimp   ! additional additivie impurity term
  end type

  type response_dp
    ! band-resolved response functions
    complex(8), allocatable :: s_full(:,:,:,:,:)
    complex(8), allocatable :: sB_full(:,:,:,:,:)
    complex(8), allocatable :: a_full(:,:,:,:,:)
    complex(8), allocatable :: aB_full(:,:,:,:,:)

    ! gather arrays for MPI
    complex(8), allocatable :: s_gather(:,:,:,:,:)
    complex(8), allocatable :: sB_gather(:,:,:,:,:)
    complex(8), allocatable :: a_gather(:,:,:,:,:)
    complex(8), allocatable :: aB_gather(:,:,:,:,:)

    ! total band and k-summation
    complex(8), allocatable :: s_sum(:,:,:)
    complex(8), allocatable :: sB_sum(:,:,:)
    complex(8), allocatable :: a_sum(:,:,:)
    complex(8), allocatable :: aB_sum(:,:,:)

    ! derived quantities
    complex(8), allocatable :: Seebeck(:)
    complex(8), allocatable :: Nernst(:)
    complex(8), allocatable :: RH(:)
  end type

  type response_qp
    ! band-resolved response functions
    complex(16), allocatable :: s_full(:,:,:,:,:)
    complex(16), allocatable :: sB_full(:,:,:,:,:)
    complex(16), allocatable :: a_full(:,:,:,:,:)
    complex(16), allocatable :: aB_full(:,:,:,:,:)

    ! gather arrays for MPI
    complex(16), allocatable :: s_gather(:,:,:,:,:)
    complex(16), allocatable :: sB_gather(:,:,:,:,:)
    complex(16), allocatable :: a_gather(:,:,:,:,:)
    complex(16), allocatable :: aB_gather(:,:,:,:,:)

    ! band and k-summation
    complex(16), allocatable :: s_sum(:,:,:)
    complex(16), allocatable :: sB_sum(:,:,:)
    complex(16), allocatable :: a_sum(:,:,:)
    complex(16), allocatable :: aB_sum(:,:,:)

    ! derived quantities
    complex(16), allocatable :: Seebeck(:)
    complex(16), allocatable :: Nernst(:)
    complex(16), allocatable :: RH(:)
  end type

end module Mtypes
