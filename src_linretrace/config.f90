module Mconfig
  use Mlookup
  use Mtypes
  use Mparams
  use Maux
  implicit none

contains

subroutine read_config(algo, edisp, sct, temp, imp)
  implicit none
  type(algorithm)   :: algo
  type(kpointmesh)  :: kmesh
  type(energydisp)  :: edisp
  type(scattering)  :: sct
  type(temperature) :: temp
  type(impurity)    :: imp

  character(len=256) :: config_file, output
  character(len=256) :: str_temp, str_imp
  integer :: i,j,k,l,stat,iimp

  integer :: search_start, search_end
  integer :: subsearch_start, subsearch_end
  integer :: pst, empty

  real(8), allocatable :: impurityinfo(:)
  integer :: nshape(1)

  logical :: found

  if (iargc() .ne. 1) then
    call stop_with_message(stderr, 'The program has to be executed with exactly one argument. (Name of config file)')
  end if

  call getarg(1,config_file)

  open(unit=10,file=trim(config_file),action='read',iostat=stat)
  if (stat .ne. 0) then
    call stop_with_message(stderr, 'Input file cannot be opened') ! send to stderr
  endif

  ! line counting
  lines=0
  read_count: do
    read(10,'(A)',END=200) str_temp ! read whole line as string, doesnt skip empty lines
    lines=lines+1
  enddo read_count

  200 continue
  rewind 10

  allocate(file_temp(lines))

  ! remove empty lines and comment strings
  empty=0
  read_temp: do i=1,lines
    read(10,'(A)') str_temp
      str_temp = trim(adjustl(str_temp))
      pst=scan(str_temp,cmnt) ! find out possible comment symbol
      if (pst .eq. 1) then ! whole line is commented
        file_temp(i) = ''
      elseif (pst .eq. 0) then ! no comment symbol found
        file_temp(i) = str_temp
      else  ! getting left side of comment
        file_temp(i) = trim(str_temp(1:pst-1))
      endif

      if (len_trim(file_temp(i)) .eq. 0) then ! filter out empty lines
        empty=empty+1
      endif
  enddo read_temp

  ! rewrite everything to a new clean string array
  allocate(file_save(lines-empty))
  j=1
  read_save: do i=1,lines
    if(len_trim(file_temp(i)) .ne. 0) then
      file_save(j)=file_temp(i)
      j=j+1
    endif
  enddo read_save
  deallocate(file_temp)

  lines=lines-empty
  close(unit=10)


  ! setting up defaults
  algo%output_file    = 'linretrace-output.hdf5'
  algo%input_energies = ''
  algo%lBField        = .false.
  algo%rootMethod     = 2     ! 0 -> secant; 1 -> linint; 2 -> riddler; 3 -> bisection
  algo%muFermi        = .false. ! we evaluate the occupation with the digamma function
  algo%lInterbandQuantities = .true.
  algo%lEnergyOutput  = .true.
  algo%lBoltzmann     = .true.
  algo%lFullOutput    = .false.
  sct%gamimp          = 0.d0
  imp%nimp            = 0

  !--------------------------------------------------------------------------------
  !--------------------------------------------------------------------------------
  !--------------------------------------------------------------------------------

  ! search for General stuff + Allocation of values
  call group_find('[General]', search_start, search_end)
  if (search_start .eq. 0) then
    call stop_with_message(stderr, 'General Group not found')
  else if (search_start .eq. -1) then
    call stop_with_message(stderr, 'General Group empty')
  endif
  !--------------------------------------------------------------------------------
  call string_find('EnergyFile', algo%input_energies, search_start, search_end, found)
  call string_find('OutputFile', algo%output_file, search_start, search_end, found)

  call bool_find('BFieldMode', algo%lBfield, search_start, search_end, found)
  call bool_find('DebugMode', algo%lDebug, search_start, search_end, found)

  call float_find('ChemicalPotential', edisp%mu, search_start, search_end, found)
  if (found) then
    algo%muSearch = .false.
  else
    algo%muSearch = .true.
  endif

  call bool_find('FermiOccupation', algo%muFermi, search_start, search_end, found)
  call int_find('RootMethod', algo%rootMethod, search_start, search_end, found)

  call bool_find('FullOutput', algo%lFullOutput, search_start, search_end, found)
  call bool_find('EnergyOutput', algo%lEnergyOutput, search_start, search_end, found)
  call bool_find('Boltzmann', algo%lBoltzmann, search_start, search_end, found)
  call bool_find('Interband', algo%lInterbandQuantities, search_start, search_end, found)

  call floatn_find('Scissors', edisp%scissors, search_start, search_end, found)
  if (found) then
    algo%lScissors = .true.
  else
    algo%lScissors = .false.
  endif

  call int_find('NImp', imp%nimp, search_start, search_end, found)
  if (found) then
    if (imp%nimp > 0) then
      algo%lImpurities = .true.
      allocate(imp%inputspin(imp%nimp))
      allocate(imp%inputtype(imp%nimp))
      allocate(imp%Dopant(imp%nimp))
      allocate(imp%Density(imp%nimp))
      allocate(imp%Energy(imp%nimp))
      allocate(imp%Degeneracy(imp%nimp))
    else if (imp%nimp == 0) then
      algo%lImpurities = .false.
    else
      call stop_with_message(stderr, 'Error: Negative number of impurities')
    endif
  else
    algo%lImpurities = .false.
  endif
  !--------------------------------------------------------------------------------
  !--------------------------------------------------------------------------------
  !--------------------------------------------------------------------------------

  call group_find('[Scattering]', search_start, search_end)
  if (search_start .eq. 0) then
    call stop_with_message(stderr, 'Scattering Group not found')
  else if (search_start .eq. -1) then
    call stop_with_message(stderr, 'Scattering Group empty')
  endif
  !--------------------------------------------------------------------------------
  call string_find('ScatteringFile', algo%input_scattering, search_start, search_end, found)
  if (found) then
    algo%lScatteringFile = .true.
  else
    algo%lScatteringFile = .false.
  endif
  call float_find('ScatteringImpurity', sct%gamimp, search_start, search_end, found)

  if (.not. algo%lScatteringFile) then
    call float_find('TMinimum', temp%Tmin, search_start, search_end, found)
    if (.not. found) call stop_with_message(stderr, 'TMinimum in Scattering group not found')
    call float_find('TMaximum', temp%Tmax, search_start, search_end, found)
    if (.not. found) call stop_with_message(stderr, 'TMaximum in Scattering group not found')
    call int_find('TPoints', temp%nT, search_start, search_end, found)
    if (.not. found) call stop_with_message(stderr, 'TPoints in Scattering group not found')
    call floatn_find('ScatteringCoefficients', sct%gamcoeff, search_start, search_end, found)
    if (.not. found) call stop_with_message(stderr, 'ScatteringCoefficients in Scattering group not found')
    call floatn_find('QuasiParticleCoefficients', sct%zqpcoeff, search_start, search_end, found)
    if (.not. found) call stop_with_message(stderr, 'QuasiParticleCoefficients in Scattering group not found')

    edisp%lBandShift = .false. ! only with scattering File
  endif

  call group_find('[Impurities]', search_start, search_end)
  if (algo%lImpurities .and. search_start .eq. 0) then
    call stop_with_message(stderr, 'Impurities group not found')
  else if (algo%lImpurities .and. search_start .eq. -1) then
    call stop_with_message(stderr, 'Impurities group empty')
  else if (imp%nimp == 0 .and. search_start .gt. 0) then
    call stop_with_message(stderr, 'Impurities found but no NImp supplied')
  endif

  if (search_start .gt. 0) then ! found
    do iimp=1,imp%nimp
      write(str_imp,'(A2,I1,A2)') '[[',iimp,']]'
      call subgroup_find(str_imp, search_start, search_end, subsearch_start, subsearch_end)
      if (subsearch_start .le. 0) then
        call stop_with_message(stderr, 'Impurity sub group not found')
      endif
      if ((iimp .eq. imp%nimp) .and. (subsearch_end .ne. search_end)) then
        call stop_with_message(stderr,'More Impurity descriptions than provided in NImp')
      endif

      call floatn_find('Absolute', impurityinfo, subsearch_start, subsearch_end, found)
      if (.not. found) then
        call floatn_find('Valence', impurityinfo, subsearch_start, subsearch_end, found)
        if (.not. found) then
          call floatn_find('Conduction', impurityinfo, subsearch_start, subsearch_end, found)
          if (.not. found) then
            call floatn_find('Percentage', impurityinfo, subsearch_start, subsearch_end, found)
            if (.not. found) then
              call stop_with_message(stderr, 'Impurity Description not found')
            else
              imp%inputtype(iimp) = 3
            endif
          else
            imp%inputtype(iimp) = 2
          endif
        else
          imp%inputtype(iimp) = 1
        endif
      else
        imp%inputtype(iimp) = 0
      endif

      nshape = shape(impurityinfo)
      if (imp%inputtype(iimp) == 0) then
        if (nshape(1) /= 4) then
          call stop_with_message(stderr, 'Absolute impurity description has exactly 4 parameters')
        else
          imp%inputspin(iimp) = 1 ! default to spin up
          ! we don't really need a spin descprtion because this is an absolute energy level
        endif
      else
        if (nshape(1) == 4) then
          imp%inputspin(iimp) = 1 ! default to spin up
        else if (nshape(1) == 5) then
          imp%inputspin(iimp) = impurityinfo(5)
        else
          call stop_with_message(stderr, 'Relative impurity description have 4 or 5 parameters')
        endif
      endif

      if (abs(impurityinfo(1)) == 1.d0) then
        imp%Dopant(iimp)     = impurityinfo(1) ! +1 (donor) ; -1 (acceptor)
      else
        call stop_with_message(stderr, 'Dopant description is either +1 or -1')
      endif
      imp%Density(iimp)    = impurityinfo(2) ! density / unit cell  < 1
      imp%Energy(iimp)     = impurityinfo(3) ! energylevel [eV]
      imp%Degeneracy(iimp) = impurityinfo(4) ! degeneracy g
      deallocate(impurityinfo)
    enddo
  endif

  ! note here: the adjustment for the energy level
  ! will be after we read in the gap information

  deallocate(file_save)
end subroutine read_config


subroutine check_config(algo)
  implicit none
  type(algorithm) :: algo
  logical         :: there

  inquire (file=trim(adjustl(algo%input_energies)), exist=there)
  if (.not. there) then
    call stop_with_message(stderr, "Can not find the EnergyFile")
  endif

  if (algo%lScatteringFile) then
    inquire (file=trim(adjustl(algo%input_scattering)), exist=there)
    if (.not. there) then
      call stop_with_message(stderr, "Can not find the ScatteringFile")
    endif
  endif

end subroutine check_config

end module Mconfig
