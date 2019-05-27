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
    logical :: muSearch              ! mu fixed or find mu?
    logical :: muFermi               ! calculate the occupation with fermi functions instead of digamma functions
    logical :: lScatteringFile       ! do we get the scattering information from another file
    logical :: lInterbandQuantities  ! calc inter band response
    logical :: lFullOutput    ! output full response dependency
    logical :: lEnergyOutput  ! output renormalized energies
    logical :: lBoltzmann     ! calc boltzmann response
    logical :: lScissors      ! apply gap widening
    logical :: lImpurities    ! include impurity levels
    character(len=256) :: input_energies
    character(len=256) :: input_scattering
    character(len=256) :: output_file
  end type

  ! information about the k-points which is necessary for us
  ! that is: number of k-points and its weight
  type kpointmesh
    real(8)               :: weightsum
    real(8), allocatable  :: weight(:)
    real(16), allocatable :: weightQ(:)
    real(8), allocatable  :: multiplicity(:)
    integer               :: nkp
    real(8)               :: vol
  end type

  ! energy dispersion and derived quantities
  ! direct derivatives and transition elements
  type energydisp
    integer :: nband_max
    integer :: nbopt_min                     ! number of bands (interval) included in the optical matrix elements
    integer :: nbopt_max                     ! number of bands (interval) included in the optical matrix elements
    integer :: iSpin                         ! number of spins
    logical :: lDerivatives                  ! do we have the derivatives (band_dk, band_d2k)
    logical :: lBandShift   ! do we get band_shifts from the scattering file?
    logical :: lFullMoments ! do we have the full optical elements (n n' dependence)
    integer :: iOptical     ! number of optical elements 3 6 or 9
    real(8) :: efer         ! Fermi energy       -- only for temporary use in config file
    real(8) :: mu           ! chemical potential -- only for temporary use in config file
    real(8) :: nelect

    ! gap information
    logical              :: gapped_complete       ! is the system completely gapped (false if spin-dependnet gap)
    real(8)              :: gap_min               ! smallest gap -> important for mu-refinement
    logical, allocatable :: gapped(:)
    real(8), allocatable :: gap(:)
    integer, allocatable :: valenceBand(:)        ! band number
    integer, allocatable :: conductionBand(:)     ! band number
    real(8), allocatable :: ene_valenceBand(:)    ! highest energy of valence band
    real(8), allocatable :: ene_conductionBand(:) ! lowest energy of conduction band

    ! gap widening
    real(8), allocatable :: scissors(:)

    real(8), allocatable    :: band_original(:,:,:)
    real(8), allocatable    :: band_shift(:,:,:)  ! same as band
    real(8), allocatable    :: band(:,:,:)        ! energy(nband,ik,ispin)
    real(8), allocatable    :: band_dk(:,:,:,:)   ! d/dk_i band(nband,ik,ispin)
    real(8), allocatable    :: band_d2k(:,:,:,:)  ! d2/(dk_i dk_j) band(nband,ik,ispin)

    ! optical elements (because of the double band dependencies)
    ! are loaded for each k-point and each spin
    real(8), allocatable    :: Mopt(:,:,:,:)     ! M(xy,n,n')= <n,k|p.e_x|n',k> * <n',k|p.e_y|n,k> *
                                                 ! 3..9, nband,nband, spin
    ! the diagonal optical elements
    ! are loaded in one go
    real(8), allocatable    :: MoptDiag(:,:,:,:) ! 3..9, nband, spin, k-points
  end type

  type impurity
    integer :: nimp                          ! number of impurities
    integer, allocatable :: inputtype(:)     ! how the energy is provided
    integer, allocatable :: inputspin(:)
      ! 0: absolute
      ! 1: relative from top of the valence band (+ -> higher)
      ! 2: relative from bottom of conduction band (+ -> lower)
      ! 3: percentage of gap -> added to the top of the valence band
    real(8), allocatable :: Dopant(:)
    real(8), allocatable :: Density(:)
    real(8), allocatable :: Energy(:)
    real(8), allocatable :: Degeneracy(:)
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
