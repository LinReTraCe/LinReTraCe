module Mconfig
  use Mlookup
  use Mtypes
  use Mparams
  use Maux
  implicit none

contains

subroutine read_config(algo, edisp, sct, temp, pot, imp)
  implicit none
  type(algorithm)   :: algo
  type(kpointmesh)  :: kmesh
  type(energydisp)  :: edisp
  type(scattering)  :: sct
  type(temperature) :: temp
  type(potential)   :: pot
  type(impurity)    :: imp

  character(20)        :: date, time, zone
  integer,dimension(8) :: time_date_values

  character(len=256) :: config_file, output
  character(len=256) :: str_temp, str_imp
  integer :: i,j,k,l,stat,iimp

  integer :: search_start, search_end
  integer :: subsearch_start, subsearch_end
  integer :: subsubsearch_start, subsubsearch_end
  integer :: pst, empty
  integer :: er

  real(8), allocatable :: impurityinfo(:)
  character(len=256), allocatable :: rootmethod(:)
  character(len=256), allocatable :: impdescription(:)
  character(len=256), allocatable :: imptype(:)
  character(len=256), allocatable :: dictionary(:)
  character(len=256) :: erstr
  integer :: nshape(1)
  real(8) :: floattemp

  logical :: found

  if (iargc() .ne. 1) then
    call stop_with_message(stderr, 'The program has to be executed with exactly one argument. (Name of config file)')
  end if

  call getarg(1,config_file)

  open(unit=10,file=trim(config_file),action='read',iostat=stat)
  if (stat .ne. 0) then
    call stop_with_message(stderr, 'Config file could not be opened') ! send to stderr
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
  algo%lTMODE         = .false.
  algo%lMUMODE        = .false.

  algo%output_file           = ''
  algo%input_energies        = ''
  algo%input_scattering_hdf5 = ''
  algo%input_scattering_text = ''
  algo%input_mu_text = ''
  algo%old_output_file       = ''
  algo%lBField               = .false.
  algo%rootMethod            = 2     ! 0 -> secant; 1 -> linint; 2 -> ridders; 3 -> bisection
  algo%muFermi               = .false. ! we evaluate the occupation with the digamma function
  algo%lQuad                 = .true.

  algo%lScatteringFile = .false.
  algo%lScatteringText = .false.

  algo%lInterbandQuantities = .false.
  algo%lIntrabandQuantities = .true.

  algo%lRedoMudft     = .false.

  algo%lEnergyOutput  = .false.
  algo%lBoltzmann     = .true.
  algo%lBoltzmannFermi= .true.       ! deprecated flag that switches between boltzmann ans psi1 approx
  algo%lFullOutput    = .false.
  sct%gamimp          = 0.d0
  imp%nimp            = 0

  pot%nMu             = 100

  temp%tlogarithmic   = .false.
  pot%mlogarithmic    = .false.

  edisp%doping        = 0.d0
  algo%ldoping        = .false.

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
  allocate(dictionary(13))
  dictionary(1) = 'EnergyFile'
  dictionary(2) = 'OutputFile'
  dictionary(3) = 'RunMode'
  dictionary(4) = 'Bandgap'
  dictionary(5) = 'ElectronOccupation'
  dictionary(6) = 'FullOutput'
  dictionary(7) = 'EnergyOutput'
  dictionary(8) = 'QuadResponse'
  dictionary(9) = 'Interband'
  dictionary(10) = 'Intraband'
  dictionary(11) = 'Boltzmann'
  dictionary(12) = 'BFieldMode'
  dictionary(13) = 'FermiOccupation'
  ! dictionary(14) = 'RootMethod'
  call spell_check(search_start,search_end, '[General]', dictionary, er, erstr)
  deallocate(dictionary)
  if (er /= 0) call stop_with_message(stdout, erstr)
  !--------------------------------------------------------------------------------

  call string_find('EnergyFile', algo%input_energies, search_start, search_end, found)
  call string_find('RunMode', str_temp, search_start, search_end, found)
  if (.not. found) then
    call stop_with_message(stderr, 'RunMode was not provided in [General]')
  else
    if (index('mu',to_lower(trim(str_temp))) .ne. 0) then
      algo%lMUMODE = .true.
      algo%muFermi = .true. ! change default value in mu-mode ... important for carrier concentration
    else if (index('temp',to_lower(trim(str_temp))) .ne. 0) then
      algo%lTMODE = .true.
    endif
  endif

  if ((algo%lMUMODE .and. algo%lTMODE) .or. (.not. (algo%lMUMODE .or. algo%lTMODE)))then
    call stop_with_message(stderr, 'Program has to be executed with RunMode = mu or temp.')
  endif

  call bool_find('BFieldMode', algo%lBfield, search_start, search_end, found)
  call string_find('OutputFile', algo%output_file, search_start, search_end, found)

  if (.not. found) then
    call date_and_time(date,time,zone,time_date_values)
    algo%output_file = 'lrtc-'//trim(date)//'-'//trim(time)//'-output.hdf5'
  endif

  call bool_find('FermiOccupation', algo%muFermi, search_start, search_end, found)

  call bool_find('FullOutput', algo%lFullOutput, search_start, search_end, found)
  call bool_find('EnergyOutput', algo%lEnergyOutput, search_start, search_end, found)
  call bool_find('Boltzmann', algo%lBoltzmann, search_start, search_end, found)
  ! call bool_find('BoltzFermi',algo%lBoltzmannFermi, search_start, search_end, found)
  call bool_find('Interband', algo%lInterbandQuantities, search_start, search_end, found)
  call bool_find('Intraband', algo%lIntrabandQuantities, search_start, search_end, found)
  call bool_find('QuadResponse', algo%lQuad, search_start, search_end, found)

  call floatn_find('Bandgap', edisp%scissors, search_start, search_end, found)
  if (found) then
    algo%lScissors  = .true.
    algo%lRedoMudft = .true.
  else
    algo%lScissors = .false.
  endif

  call float_find('ElectronOccupation', edisp%nelect_config, search_start, search_end, found)
  if (.not. found) then ! set to -1 if not found or wrong value
    edisp%nelect_config = -1.d0
  else
    if (edisp%nelect_config <= 0.d0) then
      edisp%nelect_config = -1.d0
      call stop_with_message(stderr, 'Error: ElectronOccupation must be > 0.0')
    else
      algo%lRedoMudft = .true.
    endif
  endif

  ! dont let user dictate root method
  ! allocate(rootmethod(0:3))
  ! rootmethod(0) = 'secant'
  ! rootmethod(1) = 'linint'
  ! rootmethod(2) = 'ridders'
  ! rootmethod(3) = 'bisection'
  ! call string_find('RootMethod', str_temp, search_start, search_end, found)
  ! if (found) then
  !   algo%rootMethod = -1
  !   do i=0,3
  !     if (index(trim(rootmethod(i)),to_lower(trim(str_temp))) .ne. 0) then
  !       algo%rootMethod = i
  !     endif
  !   enddo
  !   if (algo%rootMethod == -1) then
  !     call stop_with_message(stderr, 'RootMethod Description not available')
  !   endif
  ! endif


  if (algo%lMUMODE) then

    call group_find('[MuMode]', search_start, search_end)
    if (search_start .eq. 0) then
      call stop_with_message(stderr, 'MuMode Group not found')
    else if (search_start .eq. -1) then
      call stop_with_message(stderr, 'MuMode Group empty')
    endif
    call subgroup_find('[[Scattering]]', search_start, search_end, subsearch_start, subsearch_end)
    if (subsearch_start .le. 0) then
      call stop_with_message(stderr, 'Scattering group not found')
    endif

    !--------------------------------------------------------------------------------
    allocate(dictionary(9))
    dictionary(1) = 'ScatteringFile'
    dictionary(2) = 'ScatteringImpurity'
    dictionary(3) = 'Temperature'
    dictionary(4) = 'MuMinimum'
    dictionary(5) = 'MuMaximum'
    dictionary(6) = 'MuPoints'
    dictionary(7) = 'MLogarithmic'
    dictionary(8) = 'ScatteringRate'
    dictionary(9) = 'QuasiParticleWeight'
    call spell_check(subsearch_start,subsearch_end, '[MuMode] [[Scattering]]', dictionary, er, erstr)
    deallocate(dictionary)
    if (er /= 0) call stop_with_message(stdout, erstr)
    !--------------------------------------------------------------------------------

    call float_find('ScatteringImpurity', sct%gamimp, search_start, search_end, found)
    call string_find('ScatteringFile', algo%input_scattering_hdf5, subsearch_start, subsearch_end, found)
    if (found) then
      algo%lScatteringFile = .true.
    else
      algo%lScatteringFile = .false.
    endif

    if (.not. algo%lScatteringFile) then
      call int_find('MuPoints', pot%nMu, subsearch_start, subsearch_end, found)
      if (.not. found) call stop_with_message(stderr, 'MuPoints in MuMode group not found')
      call bool_find('MLogarithmic', pot%mlogarithmic, search_start, search_end, found)
      call float_find('Temperature', temp%temp, subsearch_start, subsearch_end, found)
      if (.not. found) call stop_with_message(stderr, 'Temperature in MuMode group not found')
      call float_find('MuMinimum', pot%MuMin, subsearch_start, subsearch_end, found)
      if (.not. found) call stop_with_message(stderr, 'MuMinimum in MuMode group not found')
      call float_find('MuMaximum', pot%MuMax, subsearch_start, subsearch_end, found)
      if (.not. found) call stop_with_message(stderr, 'MuMaximum in MuMode group not found')

      if (pot%MuMin > pot%MuMax) then
        call stop_with_message(stderr, 'MuMinimum must be smaller than MuMaximum')
      endif
      ! with respect to Fermi level at given Temperature
      ! shift afterwards

      call float_find('ScatteringRate', floattemp, subsearch_start, subsearch_end, found)
      if (.not. found) call stop_with_message(stderr, 'ScatteringRate in MuMode group not found')
      allocate(sct%gamcoeff(1))
      sct%gamcoeff(1) = floattemp
      call float_find('QuasiParticleWeight', floattemp, subsearch_start, subsearch_end, found)
      if (.not. found) call stop_with_message(stderr, 'QuasiParticleWeight in MuMode group not found')
      allocate(sct%zqpcoeff(1))
      sct%zqpcoeff(1) = floattemp

      edisp%lBandShift = .false. ! only with scattering HDF5 File where we have full control
    endif
  endif

  if (algo%lTMODE) then

    call group_find('[TempMode]', search_start, search_end)
    if (search_start .le. 0) then
      call stop_with_message(stderr, 'TempMode group not found')
    endif

    !--------------------------------------------------------------------------------
    allocate(dictionary(5))
    dictionary(1) = 'ChemicalPotential'
    dictionary(2) = 'OldOutput'
    dictionary(3) = 'OldOutputText'
    dictionary(4) = 'NImp'
    dictionary(5) = 'Doping'
    call spell_check(search_start,search_end, '[TempMode]', dictionary, er, erstr)
    deallocate(dictionary)
    if (er /= 0) call stop_with_message(stdout, erstr)
    !--------------------------------------------------------------------------------

    call float_find('ChemicalPotential', pot%mu_config, search_start, search_end, found)
    if (found) then
      algo%muSearch = .false.
    else
      algo%muSearch = .true.
    endif

    call string_find('OldOutput', algo%old_output_file, search_start, search_end, found)
    if (found) then
      algo%lOldmu   = .true.
      algo%muSearch = .false. !overwrite the previous option
    else
      algo%lOldmu = .false.
    endif

    if (.not. algo%lOldmu) then
      call string_find('OldOutputText', algo%input_mu_text, search_start, search_end, found)
      if (found) then
        algo%lOldmuText   = .true.
        algo%muSearch = .false. !overwrite the previous option
      else
        algo%lOldmuText = .false.
      endif
    endif

    call float_find('Doping', edisp%doping, search_start, search_end, found)
    if (found) then
      if (edisp%doping /= 0.d0) then
        algo%lDoping = .true.
      endif
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
        allocate(imp%Bandwidth(imp%nimp))
        allocate(imp%Bandtype(imp%nimp))
        allocate(imp%Bandcutoff(imp%nimp))
        allocate(imp%Band(imp%nimp))
      else if (imp%nimp == 0) then
        algo%lImpurities = .false.
      else
        call stop_with_message(stderr, 'Error: NImp must be > 0')
      endif
    else
      algo%lImpurities = .false.
    endif

    call subgroup_find('[[Scattering]]', search_start, search_end, subsearch_start, subsearch_end)
    if (subsearch_start .le. 0) then
      call stop_with_message(stderr, 'Scattering group not found')
    endif

    !--------------------------------------------------------------------------------
    allocate(dictionary(9))
    dictionary(1) = 'ScatteringFile'
    dictionary(2) = 'ScatteringText'
    dictionary(3) = 'ScatteringImpurity'
    dictionary(4) = 'TMinimum'
    dictionary(5) = 'TMaximum'
    dictionary(6) = 'TPoints'
    dictionary(7) = 'TLogarithmic'
    dictionary(8) = 'ScatteringCoefficients'
    dictionary(9) = 'QuasiParticleCoefficients'
    call spell_check(subsearch_start,subsearch_end, '[TempMode] [[Scattering]]', dictionary, er, erstr)
    deallocate(dictionary)
    if (er /= 0) call stop_with_message(stdout, erstr)
    !--------------------------------------------------------------------------------
    call float_find('ScatteringImpurity', sct%gamimp, search_start, search_end, found)
    call string_find('ScatteringFile', algo%input_scattering_hdf5, subsearch_start, subsearch_end, found)
    if (found) then
      algo%lScatteringFile = .true.
    else
      algo%lScatteringFile = .false.
    endif

    if (.not. algo%lScatteringFile) then
      call string_find('ScatteringText', algo%input_scattering_text, subsearch_start, subsearch_end, found)
      if (found) then
        algo%lScatteringText = .true.
      else
        algo%lScatteringText = .false.
      endif

      if (.not. algo%lScatteringText) then
        call float_find('TMinimum', temp%Tmin, search_start, search_end, found)
        if (.not. found) call stop_with_message(stderr, 'TMinimum in Scattering group not found')
        call float_find('TMaximum', temp%Tmax, search_start, search_end, found)
        if (.not. found) call stop_with_message(stderr, 'TMaximum in Scattering group not found')
        call int_find('TPoints', temp%nT, search_start, search_end, found)
        if (.not. found) call stop_with_message(stderr, 'TPoints in Scattering group not found')
        call bool_find('TLogarithmic', temp%tlogarithmic, search_start, search_end, found)
        call floatn_find('ScatteringCoefficients', sct%gamcoeff, search_start, search_end, found)
        if (.not. found) call stop_with_message(stderr, 'ScatteringCoefficients in Scattering group not found')
        call floatn_find('QuasiParticleCoefficients', sct%zqpcoeff, search_start, search_end, found)
        if (.not. found) call stop_with_message(stderr, 'QuasiParticleCoefficients in Scattering group not found')
      endif

      edisp%lBandShift = .false. ! only with scattering HDF5 File where we have full control
    endif



    if (algo%lImpurities) then
      call subgroup_find('[[Impurities]]', search_start, search_end, subsearch_start, subsearch_end)
      if (subsearch_start .le. 0) then
        call stop_with_message(stderr, 'Impurities group not found')
      endif

      allocate(impdescription(0:3))
      impdescription(0) = 'Absolute'
      impdescription(1) = 'Valence'
      impdescription(2) = 'Conduction'
      impdescription(3) = 'Percentage'
      allocate(imptype(0:3))
      imptype(0) = 'box'
      imptype(1) = 'lorentzian'
      imptype(2) = 'gaussian'
      imptype(3) = 'halfcircle'

      do iimp=1,imp%nimp
        write(str_imp,'(A2,I1,A2)') '[[[',iimp,']]]'
        call subsubgroup_find(str_imp, subsearch_start, subsearch_end, subsubsearch_start, subsubsearch_end)
        if (subsubsearch_start .le. 0) then
          call stop_with_message(stderr, 'Impurity sub group not found: '//str_imp)
        endif

        !--------------------------------------------------------------------------------
        allocate(dictionary(7))
        dictionary(1) = 'Absolute'
        dictionary(2) = 'Valence'
        dictionary(3) = 'Conduction'
        dictionary(4) = 'Percentage'
        dictionary(5) = 'Bandtype'
        dictionary(6) = 'Bandwidth'
        dictionary(7) = 'Bandcutoff'
        call spell_check(subsubsearch_start,subsubsearch_end, '[TempMode] [[Impurities]] '//str_imp, dictionary, er, erstr)
        deallocate(dictionary)
        if (er /= 0) call stop_with_message(stdout, erstr)
        !--------------------------------------------------------------------------------

        do i=0,3
          call floatn_find(impdescription(i), impurityinfo, subsubsearch_start, subsubsearch_end, found)
          if (.not. found) then
            cycle
          endif
          imp%inputtype(iimp) = i
          exit
        enddo

        ! if we didn't find an identifier
        if (.not. found) then
          call stop_with_message(stderr, 'Valid impurity Description not found in '//str_imp)
        endif

        call float_find('Bandwidth', imp%Bandwidth(iimp), subsubsearch_start, subsubsearch_end, found)
        if (.not. found) then
          imp%Bandwidth(iimp) = 0.d0
          imp%Band(iimp) = .false.
        else
          if (imp%Bandwidth(iimp) <= 0.d0) then
            imp%Bandwidth(iimp) = 0.d0
            imp%Band(iimp) = .false.
          else
            imp%Band(iimp) = .true.
          endif
        endif

        call string_find('Bandtype', str_temp, subsubsearch_start, subsubsearch_end, found)
        if (.not. found) then
          imp%Bandtype(iimp) = 0 ! box
        else
          imp%Bandtype(iimp) = -1
          do i=0,3
            if (index(trim(imptype(i)),to_lower(trim(str_temp))) .ne. 0) then
              imp%Bandtype(iimp) = i
            endif
          enddo
          if (imp%Bandtype(iimp) == -1) then
            call stop_with_message(stderr, 'Bandtype Description not available in '//str_imp)
          endif
        endif

        call float_find('Bandcutoff', imp%Bandcutoff(iimp), subsubsearch_start, subsubsearch_end, found)
        if (.not. found) then
          imp%Bandcutoff(iimp) = 1.d0 ! one standard deviation
        else
          if (imp%Bandcutoff(iimp) <= 0.d0) then
            imp%Bandcutoff(iimp) = 1.d0
          endif
        endif

        ! the saved information
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
      deallocate(impdescription)
      deallocate(imptype)
    endif !algo %lImpurities
  endif ! algo%TMODE

  ! note here: the adjustment for the energy level
  ! will be after we read in the gap information

  algo%lDebug = .false.
  algo%dbgstr = ''
  call group_find('[Debug]', search_start, search_end)
  if (search_start .ge. 1) then
     algo%lDebug = .true.
     algo%dbgstr = file_save(search_start)
  endif

  deallocate(file_save)
end subroutine read_config


subroutine check_files(algo)
  implicit none
  type(algorithm) :: algo
  logical         :: there

  inquire (file=trim(adjustl(algo%input_energies)), exist=there)
  if (.not. there) then
    call stop_with_message(stderr, "Can not find the EnergyFile")
  endif

  if (algo%lScatteringFile) then
    inquire (file=trim(adjustl(algo%input_scattering_hdf5)), exist=there)
    if (.not. there) then
      call stop_with_message(stderr, "Can not find the ScatteringFile")
    endif
  endif

  if (algo%lScatteringText) then
    inquire (file=trim(adjustl(algo%input_scattering_text)), exist=there)
    if (.not. there) then
      call stop_with_message(stderr, "Can not find the ScatteringText")
    endif
  endif

  if (algo%lOldmu) then
    inquire (file=trim(adjustl(algo%old_output_file)), exist=there)
    if (.not. there) then
      call stop_with_message(stderr, "Can not find the OldOutput file")
    endif
  endif

end subroutine check_files

end module Mconfig
