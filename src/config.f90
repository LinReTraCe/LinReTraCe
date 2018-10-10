! This file is part of the Abinitio Dynamical Vertex Approximation (ADGA)
! package. It is an electronic structure code which allows the inclusion of
! non-local correlations beyond DMFT and the calculation of momentum-dependent
! susceptibilities.
!
! The public repository can be found at
! https://github.com/AbinitioDGA/ADGA
!
! The arXiv publication can be found at
! https://arxiv.org/abs/1710.06651
!
! Copyright (C) <2017, 2018>
! <Anna Galler*, Patrick Thunström, Josef Kaufmann, Matthias Pickem, Jan M. Tomczak, Karsten Held>
! * Corresponding author. E-mail address: galler.anna@gmail.com
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.

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


  ! setting up switches
  algo%ltbind  = .false.
  algo%ltetra  = .false.
  algo%lw2k    = .false.
  algo%lvasp   = .false.
  algo%loptic  = .false.
  algo%lBfield = .false.
  algo%ldmft   = .false.
  algo%lsymm   = .true.

  ! search for General stuff + Allocation of values
  !--------------------------------------------------------------------------------
  call group_find('[General]', search_start, search_end)
  if (search_start .eq. 0) then ! group was not found
    erstr = 'General Group not found'; er = 3
    return
  endif
  if (search_start .eq. -1) then ! group was not found
    erstr = 'General Group empty'; er = 4
    return
  endif

  call string_find('System', algo%mysyst, search_start, search_end)
  call int3_find('K-grid', kmesh%kx, kmesh%ky, kmesh%kz, search_start, search_end)
  call float_find('Nelect', edisp%nelect, search_start, search_end)
  call string_find('Input-type', str_temp, search_start, search_end)
  select case(trim(adjustl(str_temp)))
    case ('wien2k')
      algo%lw2k = .true.
    case ('vasp')
      algo%lvasp = .true.
    case ('tightbinding')
      algo%ltbind = .true.
    case default
      erstr = 'Unknown input'; er = 5
      return
  end select

  call bool_find('Optical', algo%loptic, search_start, search_end)
  call bool_find('BField', algo%lBfield, search_start, search_end)
  call bool_find('Dmft', algo%ldmft, search_start, search_end)
  call bool_find('Irreducible', algo%lsymm, search_start, search_end)
  call float_find('Zqp', edisp%ztmp, search_start, search_end)

  ! search for Lattice stuff
  call group_find('[Lattice]', search_start, search_end)
  if (search_start .eq. 0) then ! group was not found
    erstr = 'Lattice Group not found'; er = 6
    return
  endif
  if (search_start .eq. -1) then
    erstr = 'Lattice Group empty'; er = 7
    return
  endif

  call float_find('Alat', lat%alat, search_start, search_end)
  call float_find('Vol', lat%vol, search_start, search_end)
  call float3_find('Lat-vectors', lat%a(1), lat%a(2), lat%a(3), search_start, search_end)

  call group_find('[Output]', search_start, search_end)
  if (search_start .eq. 0) then ! group was not found
    erstr = 'Lattice Group not found'; er = 8
    return
  endif
  if (search_start .eq. -1) then
    erstr = 'Lattice Group empty'; er = 9
    return
  endif

  call string_find('Outfile', outfile, search_start, search_end)

  !--------------------------------------------------------------------------------
  !verbose = .false.
  !verbstr = ''
  !call group_find('[Verbose]', search_start, search_end)
  !if (search_start .ge. 1) then
  !   verbose = .true.
  !   verbstr = file_save(search_start)
  !endif

  !!--------------------------------------------------------------------------------
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


subroutine init_config(kmesh)
  implicit none
  type(kpointmesh) :: kmesh

  ! fill out the rest of the datatype
  kmesh%kred = kmesh%kx*kmesh%ky*kmesh%kz
  kmesh%kful = (kmesh%kx+1)*(kmesh%ky+1)*(kmesh%kz+1)

  lat%lcubic = .false.
  lat%ltetragonal = .false.
  lat%lorthorhombic = .false.
  lat%nalpha = 3 ! we need all three directions for the responses

  ! these data points come directly from the input file
  ! i.e. we shouldnt be worrying about double precision comparisons.
  if (lat%a(1) .eq. lat%a(2) .and. lat%a(1) .eq. lat%a(3) &
      .and. lat%a(1) .eq. lat%a(3)) then
     lat%lcubic = .true.
     lat%nalpha = 1 ! truncates the number of polarizations in cubic systems
  elseif(lat%a(1) .eq. lat%a(2) .and. lat%a(1) .ne. lat%a(3)) then
     lat%ltetragonal = .true.
  elseif(lat%a(1) .eq. lat%a(3) .and. lat%a(1) .ne. lat%a(2)) then
     lat%ltetragonal = .true.
  elseif(lat%a(1) .ne. lat%a(2) .and. lat%a(2) .eq. lat%a(3)) then
     lat%ltetragonal = .true.
  elseif(lat%a(1) .ne. lat%a(2) .and. lat%a(1) .ne. lat%a(3) &
         .and. lat%a(1) .ne. lat%a(3)) then
     lat%lorthorhombic = .true.
  endif
end subroutine init_config

subroutine check_config(er,erstr)
  implicit none
  integer, intent(out)          :: er
  character(len=*), intent(out) :: erstr
  logical :: there

  if (algo%lw2k) then
     inquire (file=trim(adjustl(algo%mysyst))//'.weight', exist=there)
     if (.not. there) then
       erstr = "Error: Can not find the case.weight file"; er = 1
       return
     endif
     inquire (file=trim(adjustl(algo%mysyst))//'.struct', exist=there)
     if (.not. there) then
       erstr = "Error: Can not find the case.struct file"; er = 2
       return
     endif
     inquire (file=trim(adjustl(algo%mysyst))//'.energy', exist=there)
     if (.not. there) then
       erstr = "Error: Can not find the case.energy file"; er = 3
       return
     endif
     if (algo%loptic) then
        inquire (file=trim(adjustl(algo%mysyst))//'.symmat', exist=there)
        if (.not. there) then
          erstr = "Error: Can not find the case.symmat file"; er = 4
          return
        endif
     endif
  endif
  if (algo%lvasp) then
     ! do something
  endif

  return
end subroutine check_config

end module Mconfig
