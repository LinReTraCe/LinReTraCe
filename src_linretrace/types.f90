! derived types for the main program
! this should be completely agnostic to the type of data we provide

module Mtypes
  implicit none

  ! contains the methods employed for the root-finding and
  ! the way we treat the chemical potential
  ! also contains the file names of the input files
  type algorithm
    logical :: lTMODE         ! temperature mode
    logical :: lMUMODE        ! mu mode

    logical :: lDebug         ! debug mode --> evaluate quad precision data?
    logical :: lBfield        ! calculations in the presence of a magnetic field
                              ! this requires the existance of the band derivatives

    integer :: rootMethod     ! numerical method to find the chemical potential
    integer :: fullOutput     ! output full response -- 0:none - 1:full - 2:ksum - 3:bsum
    logical :: muSearch       ! mu fixed or find mu?
    logical :: lOldmuHdf5     ! mus from old run
    logical :: lOldmuText     ! mus from text file (lprint or the same format)
    logical :: muFermi               ! calculate the occupation with fermi functions instead of digamma functions
    logical :: lScatteringFile       ! do we get the scattering information from another file (hdf5)
    logical :: lScatteringText       ! do we get the scattering information from a text file
    logical :: lInterBandQuantities  ! calc inter band response
    logical :: lIntraBandQuantities  ! calc intra band response
    logical :: lEnergyOutput  ! output renormalized energies
    logical :: lBoltzmann     ! calc boltzmann response
    logical :: lBoltzmannFermi! calc boltzmann response with Fermi function or PolyGamma
    logical :: lScissors      ! apply gap widening
    logical :: lImpurities    ! include impurity levels
    logical :: lQuad          ! quad precision response
    logical :: lRedoMudft     ! flag to recalculate the provided dft chemical potential
    logical :: lDoping        ! include doping
    logical :: lNominalDoping  ! interpret the provided densities of the config file as nominal electron values

    integer :: steps          ! number of steps
    integer :: step_dir       ! step direction
                              ! +1 [1...steps] -- -1 [steps...1]

    character(len=256) :: input_energies
    character(len=256) :: input_scattering_hdf5
    character(len=256) :: input_scattering_text
    character(len=256) :: output_file
    character(len=256) :: input_mu_hdf5
    character(len=256) :: input_mu_text
    character(len=256) :: dbgstr
  end type

  ! information about the k-points which is necessary for us
  ! that is: number of k-points and their weights
  type kpointmesh
    real(8)               :: weightsum
    real(8), allocatable  :: weight(:)
    real(16), allocatable :: weightQ(:)
    real(8), allocatable  :: multiplicity(:)
    integer               :: nkp
    real(8)               :: vol
    integer               :: ndim
    logical, allocatable  :: dims(:)
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
    real(8) :: nelect_file   ! electrons given by energy file
    real(8) :: nelect_config ! number of electrons given by config file
    real(8) :: nelect

    real(8) :: doping       ! additional eletrons (>0) or holes (<0) w.r.t. nominal filling (nelect)

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

    real(8), allocatable    :: band_file(:,:,:)   ! directly from the hdf5 file
    real(8), allocatable    :: band_shift(:,:,:)  ! real part shifts

    real(8), allocatable    :: band(:,:,:)        ! energy(nband,ik,ispin)
    real(8), allocatable    :: band_dk(:,:,:,:)   ! d/dk_i band(nband,ik,ispin)
    real(8), allocatable    :: band_d2k(:,:,:,:)  ! d2/(dk_i dk_j) band(nband,ik,ispin)

    ! optical elements (because of the double band dependencies)
    ! for one k-point
    real(8), allocatable    :: Moptk(:,:,:,:)    ! M(xy,n,n')= <n,k|p.e_x|n',k> * <n',k|p.e_y|n,k> *
                                                 ! 3..9, nband,nband, spin
    ! for the MPI k-range
    real(8), allocatable    :: Mopt(:,:,:,:,:)   ! M(xy,n,n')= <n,k|p.e_x|n',k> * <n',k|p.e_y|n,k> *
                                                 ! 3..9, nband,nband, spin, krange
    ! the diagonal optical elements
    ! are loaded in one go
    real(8), allocatable    :: MoptDiag(:,:,:,:) ! 3..9, nband, spin, k-points

    ! diagonal magnetic optical elements
    real(8), allocatable    :: MBoptDiag(:,:,:,:,:,:) ! 3, 3, 3, nband, spin, k-points
  end type

  type impurity
    integer :: nimp                          ! number of impurities
    integer, allocatable :: inputtype(:)     ! how the energy is provided
      ! 1: absolute
      ! 2: relative from top of the valence band (+ -> higher)
      ! 3: relative from bottom of conduction band (+ -> lower)
      ! 4: percentage of gap -> added to the top of the valence band
    integer, allocatable :: inputspin(:)
    logical, allocatable :: Band(:) ! True if band - False if Level
    real(8), allocatable :: Bandwidth(:)
    integer, allocatable :: Bandtype(:)
      ! 1: Box , 2: Triangle, 3: Halfcircle, 4: Sine, 5: Sine^2 6: Sine^3, 7: Sine^4

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
    logical :: tlogarithmic        ! logarithmic steps
    real(8) :: Tmin                ! bottom of the temperature window
    real(8) :: Tmax                ! top of the temperature window
    real(8) :: dT                  ! temperature spacing
    real(8), allocatable :: TT(:)  ! temperature grid [K]
    real(8), allocatable :: BB(:)  ! inverse temperature grid [eV]

    real(8) :: temp_config         ! temperature provided by config file
  end type

  type potential
    logical  :: mabsolute            ! absolute values from input
    logical  :: mlogarithmic         ! logarithmic steps
    real(8)  :: MuMin
    real(8)  :: MuMax
    real(16) :: dMu
    real(8), allocatable  :: MM(:)  ! mu array
    real(16), allocatable :: QMM(:) ! same mu array in quad precision
    real(8), allocatable  :: occ(:) ! corresponding occupation

    real(8)  :: mu_config             ! chemical potential provided by config file
    real(8)  :: mu_dft_file           ! chemical potential provided by energy file
    real(8)  :: mu_dft                ! new DFT chemical potential (changed electrons / bandgap)
  end type

  type runinfo
    ! information about the current status of the run

    integer  :: iStep ! current step number

    real(8)  :: Temp  ! current temperature
    real(8)  :: beta
    real(8)  :: beta2p

    real(16) :: TempQ ! current temperature quad
    real(16) :: betaQ
    real(16) :: beta2pQ

    real(8)  :: mu    ! current chemical potential
    real(16) :: muQ

    integer  :: ik
  end type

  type scattering
    ! scattering rates and quasiparticle weights
    real(8), allocatable :: gamcoeff(:,:) ! coeff spin -- temperature coeffs
    real(8), allocatable :: zqpcoeff(:,:) ! coeff spin -- temperature coeffs
    real(8), allocatable :: enecoeff(:,:) ! coeff spin -- gamma energy coeffs
    real(8), allocatable :: gamtext(:,:)  ! step, spin
    real(8), allocatable :: zqptext(:,:)  ! step, spin
    real(8), allocatable :: gam(:,:,:) ! (nband,ik,ispin)
    real(8), allocatable :: zqp(:,:,:) ! (nband,ik,ispin)
    real(8)              :: gamimp   ! additional additivie impurity term

    logical              :: enescaling
  end type

  type response_dp
    ! band-resolved response functions
    complex(8), allocatable :: s_full(:,:,:,:,:)
    complex(8), allocatable :: sB_full(:,:,:,:,:,:)
    complex(8), allocatable :: a_full(:,:,:,:,:)
    complex(8), allocatable :: aB_full(:,:,:,:,:,:)
    complex(8), allocatable :: x_full(:,:,:,:,:)
    complex(8), allocatable :: xB_full(:,:,:,:,:,:)

    ! total band and k-summation
    complex(8), allocatable :: s_sum(:,:,:)
    complex(8), allocatable :: sB_sum(:,:,:,:)
    complex(8), allocatable :: a_sum(:,:,:)
    complex(8), allocatable :: aB_sum(:,:,:,:)
    complex(8), allocatable :: x_sum(:,:,:)
    complex(8), allocatable :: xB_sum(:,:,:,:)

    ! gather arrays for all T-points of the band and k-summed quantities
    complex(8), allocatable :: s_sum_range(:,:,:,:)
    complex(8), allocatable :: sB_sum_range(:,:,:,:,:)
    complex(8), allocatable :: a_sum_range(:,:,:,:)
    complex(8), allocatable :: aB_sum_range(:,:,:,:,:)
    complex(8), allocatable :: x_sum_range(:,:,:,:)
    complex(8), allocatable :: xB_sum_range(:,:,:,:,:)
  end type

  type response_qp
    ! band-resolved response functions
    complex(16), allocatable :: s_full(:,:,:,:,:)
    complex(16), allocatable :: sB_full(:,:,:,:,:,:)
    complex(16), allocatable :: a_full(:,:,:,:,:)
    complex(16), allocatable :: aB_full(:,:,:,:,:,:)
    complex(16), allocatable :: x_full(:,:,:,:,:)
    complex(16), allocatable :: xB_full(:,:,:,:,:,:)

    ! band and k-summation
    complex(16), allocatable :: s_sum(:,:,:)
    complex(16), allocatable :: sB_sum(:,:,:,:)
    complex(16), allocatable :: a_sum(:,:,:)
    complex(16), allocatable :: aB_sum(:,:,:,:)
    complex(16), allocatable :: x_sum(:,:,:)
    complex(16), allocatable :: xB_sum(:,:,:,:)

    ! gather arrays for all T-points of the band and k-summed quantities
    ! these have to be double precision not quadruple ( they are for output only)
    complex(8), allocatable :: s_sum_range(:,:,:,:)
    complex(8), allocatable :: sB_sum_range(:,:,:,:,:)
    complex(8), allocatable :: a_sum_range(:,:,:,:)
    complex(8), allocatable :: aB_sum_range(:,:,:,:,:)
    complex(8), allocatable :: x_sum_range(:,:,:,:)
    complex(8), allocatable :: xB_sum_range(:,:,:,:,:)
  end type

end module Mtypes
