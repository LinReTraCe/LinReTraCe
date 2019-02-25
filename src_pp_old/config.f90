module Mconfig
  use Mlookup
  use Mparams
  implicit none

contains

subroutine read_config(kmesh, edisp, sct, outfile, er, erstr)
  implicit none
  type(kpointmesh) :: kmesh
  type(energydisp) :: edisp
  type(scatrate)   :: sct
  character(len=150), intent(out) :: outfile
  integer, intent(out)            :: er
  character(len=150), intent(out) :: erstr

  character(len=150) :: config_file, output
  character(len=150) :: str_temp
  integer :: i,j,k,l,stat

  integer :: search_start, search_end
  integer :: subsearch_start, subsearch_end
  integer :: pst, empty

  logical :: verbose, debug

  er = 0
  erstr = ''

  ! Config File checks
  if (iargc() .ne. 1) then
    er = 1
    erstr = 'The program has to be executed with exactly one argument. (Name of config file)'
    return
  end if
  call getarg(1,config_file)

  open(unit=10,file=trim(config_file),action='read',iostat=stat)
  if (stat .ne. 0) then
    close(10)
    erstr = 'Input file cannot be opened'; er = 2
    return
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
  algo%lw2k         = .false.
  algo%lvasp        = .false.
  algo%lderivatives = .false.
  outfile           = 'preproc_data.hdf5'

  ! search for General stuff + Allocation of values
  !--------------------------------------------------------------------------------
  call string_find('System', algo%mysyst, 1, lines)
  call string_find('SystemType',  str_temp, 1, lines)
  select case(trim(adjustl(str_temp)))
    case ('wien2k')
      algo%lw2k = .true.
    case ('Wien2k')
      algo%lw2k = .true.
    case ('WIEN2k')
      algo%lw2k = .true.
    case ('WIEN2K')
      algo%lw2k = .true.
    case ('vasp')
      algo%lvasp = .true.
    case ('Vasp')
      algo%lvasp = .true.
    case ('VASP')
      algo%lvasp = .true.
    case default
      erstr = 'Unknown input'; er = 5
      return
  end select

  call string_find('Derivatives', algo%myderivatives, 1, lines)
  call bool_find('OpticalElements', algo%loptic, 1, lines)
  call bool_find('DMFTData', algo%ldmft, 1, lines)
  call string_find('OutputFile', outfile, 1, lines)

  !--------------------------------------------------------------------------------
  !verbose = .false.
  !verbstr = ''
  !call group_find('[Verbose]', search_start, search_end)
  !if (search_start .ge. 1) then
  !   verbose = .true.
  !   verbstr = file_save(search_start)
  !endif

  !--------------------------------------------------------------------------------
  !debug = .false.
  !dbgstr = ''
  !call group_find('[Debug]', search_start, search_end)
  !if (search_start .ge. 1) then
  !   debug = .true.
  !   dbgstr = file_save(search_start)
  !endif

  deallocate(file_save)
  return
end subroutine read_config


! subroutine init_config(kmesh, edisp)
!   implicit none
!   type(kpointmesh) :: kmesh
!   type(energydisp) :: edisp
!   integer :: iband

!   ! fill out the rest of the datatype
!   kmesh%kred = kmesh%kx*kmesh%ky*kmesh%kz
!   kmesh%kful = (kmesh%kx+1)*(kmesh%ky+1)*(kmesh%kz+1)

!   lat%nalpha = 3 ! we need all three directions for the responses

!   if (algo%ltbind) then
!      write(*,*)'READ_CONFIG: reading TB parameters'
!      allocate(edisp%E0(edisp%nband_max),edisp%t(edisp%nband_max, edisp%tmax))
!      allocate(edisp%a(edisp%tmax))
!      open(unit=10, file=trim(adjustl(algo%mysyst)), status='old')
!      read(10,*) edisp%a
!      do iband=1,edisp%nband_max
!         read(10,*)edisp%E0(iband), edisp%t(iband,:)
!      enddo
!   endif
! end subroutine init_config

subroutine check_config(er,erstr)
  implicit none
  integer, intent(out)          :: er
  character(len=*), intent(out) :: erstr
  logical :: there

  if (algo%lw2k) then
     inquire (file=trim(adjustl(algo%mysyst))//'.scf', exist=there)
     if (.not. there) then
       erstr = "Error: Can not find the case.scf file"; er = 1
       return
     endif
     inquire (file=trim(adjustl(algo%mysyst))//'.klist', exist=there)
     if (.not. there) then
       erstr = "Error: Can not find the case.klist file"; er = 2
       return
     endif
     inquire (file=trim(adjustl(algo%mysyst))//'.struct', exist=there)
     if (.not. there) then
       erstr = "Error: Can not find the case.struct file"; er = 3
       return
     endif
     inquire (file=trim(adjustl(algo%mysyst))//'.energy', exist=there)
     if (.not. there) then
       erstr = "Error: Can not find the case.energy file"; er = 4
       return
     endif
     if (algo%loptic) then
        inquire (file=trim(adjustl(algo%mysyst))//'.symmat', exist=there)
        if (.not. there) then
          erstr = "Error: Can not find the case.symmat file"; er = 5
          return
        endif
     endif
  endif
  if (algo%lvasp) then
     inquire(file=trim(adjustl(algo%mysyst))//'vasprun.xml', exist=there)
     if (.not. there) then
       erstr = "Error: Can not find vasprun.xml"; er = 6
       return
     endif
  endif

end subroutine check_config

end module Mconfig
