module Mlookup
  implicit none

  ! config file auxiliary variables
  integer :: lines

  ! comment character within config file
  character(len=1), parameter     :: cmnt = '#'
  ! seperator character within config file
  character(len=1), parameter     :: separator = '='
  ! multivalue separator character within config file
  character(len=1), parameter     :: multseparator = ' ' ! space
  ! where we save the whole config file
  character(len=256), allocatable :: file_temp(:), file_save(:)
  ! auxiliary variables
  integer :: i,j,pst
  character(len=256) :: str_temp, str_split, str_test

  public  :: lines, cmnt, separator, multseparator, file_temp, file_save
  private :: i,j,pst, str_temp, str_split, str_test

  contains

  ! find entry that needs to be saved as string
  ! search_string gets searched between search_start <= ... <= search_end
  ! if found:  found is set to .true. and entry is saved in save_string
  ! if not found: found is set to .false. (save_string contains garbage)
  subroutine string_find(search_string, save_string, search_start, search_end, found)
    character(*), intent(in)  :: search_string
    character(len=256), intent(inout) :: save_string ! keep default string
    integer, intent(in) :: search_start, search_end
    logical, intent(out) :: found

    found = .false.
    do i=search_start,search_end
      if (index(trim(file_save(i)),trim(search_string)) .ne. 0) then ! rough scan
        str_temp=file_save(i)
        pst=scan(str_temp,separator)
        if (pst > 0) then ! found the separator
          str_test    = trim(adjustl(str_temp(:pst-1))) ! left side of separator
          save_string = trim(adjustl(str_temp(pst+1:))) ! right side of separator
          if (trim(str_test) == trim(search_string)) then ! detailed comparison
            found = .true.
            exit
          endif
        endif
      endif
    enddo
  end subroutine string_find

  ! find entry that needs to be saved as integer
  ! search_string gets searched between search_start <= ... <= search_end
  ! if found:  found is set to .true. and entry is saved in save_int
  ! if not found: found is set to .false. (save_string contains garbage)
  subroutine int_find(search_string, save_int, search_start, search_end, found)
    character(*), intent(in)  :: search_string
    integer, intent(inout) :: save_int ! keep default values
    integer, intent(in) :: search_start, search_end
    logical, intent(out) :: found

    found = .false.
    do i=search_start,search_end
      if (index(trim(file_save(i)),trim(search_string)) .ne. 0) then
        str_temp=file_save(i)
        pst=scan(str_temp,separator)
        if (pst > 0) then ! found the separator
          str_test = trim(adjustl(str_temp(:pst-1))) ! left side of separator
          str_temp = trim(adjustl(str_temp(pst+1:))) ! right side of separator
          if (trim(str_test) == trim(search_string)) then ! detailed comparison
            read(str_temp,*) save_int
            found = .true.
            exit
          endif
        endif
      endif
    enddo
  end subroutine int_find

  ! find entry that needs to be saved as MULTIPLE integer
  ! automatically detects how many integers there are and saves it into the provided, unallocated array
  ! search_string gets searched between search_start <= ... <= search_end
  ! if found:  found is set to .true. and entry is saved in save_int
  ! if not found: found is set to .false. (save_int contains garbage)
  subroutine intn_find(search_string, int_array, search_start, search_end, found)
    character(*), intent(in) :: search_string
    integer, intent(in) :: search_start, search_end
    integer, intent(inout), allocatable :: int_array(:)
    logical, intent(out) :: found

    integer :: cnt
    character(len=256) :: str_original
    character(len=256) :: str_2save

    found = .false.
    do i=search_start,search_end
      if (index(trim(file_save(i)),trim(search_string)) .ne. 0) then
        str_temp=file_save(i)
        pst=scan(str_temp,separator)

        if (pst > 0) then ! found the separator
          str_test = trim(adjustl(str_temp(:pst-1))) ! left side of separator
          str_temp = trim(adjustl(str_temp(pst+1:))) ! right side of separator
          if (trim(str_test) == trim(search_string)) then ! detailed comparison
            str_original = str_temp

            ! we count the number of coefficients
            ! this behaves quite weirdly because of the empty space scan of an empty string
            ! definitely works as intended
            cnt = 0
            do
              pst=scan(str_temp,multseparator)
              if (pst == 1) then ! we find the empty space in an empty string at the first position
                exit
              else
                cnt = cnt + 1
                str_temp = trim(adjustl(str_temp(pst+1:)))
              endif
            enddo

            allocate(int_array(cnt))

            ! now that we know the number of floats we want to find
            ! we simply reset the string to its original content
            ! and do the same scan
            str_temp = str_original ! reset
            cnt = 1
            do
              pst=scan(str_temp,multseparator)
              if (pst == 1) then
                exit
              else
                str_2save = trim(adjustl(str_temp(:pst-1)))
                read(str_2save,*) int_array(cnt)
                str_temp = trim(adjustl(str_temp(pst+1:)))
                cnt = cnt + 1
              endif
            enddo

            found = .true.
          endif
        endif
      endif
    enddo
  end subroutine intn_find

  ! find entry that needs to be saved as floating number
  ! search_string gets searched between search_start <= ... <= search_end
  ! if found:  found is set to .true. and entry is saved in save_float
  ! if not found: found is set to .false. (save_float contains garbage)
  subroutine float_find(search_string, save_float, search_start, search_end, found)
    character(*), intent(in)  :: search_string
    real(8), intent(inout) :: save_float ! keep default values
    integer, intent(in) :: search_start, search_end
    logical, intent(out) :: found

    found = .false.
    do i=search_start,search_end
      if (index(trim(file_save(i)),trim(search_string)) .ne. 0) then
        str_temp=file_save(i)
        pst=scan(str_temp,separator)
        if (pst > 0) then ! found the separator
          str_test = trim(adjustl(str_temp(:pst-1))) ! left side of separator
          str_temp = trim(adjustl(str_temp(pst+1:))) ! right side of separator
          if (trim(str_test) == trim(search_string)) then ! detailed comparison
            read(str_temp,*) save_float
            found = .true.
            exit
          endif
        endif
      endif
    enddo
  end subroutine float_find

  ! find entry that needs to be saved as multiple floats
  ! automatically detects how many integers there are and saves it into the provided, unallocated array
  ! search_string gets searched between search_start <= ... <= search_end
  ! if found:  found is set to .true. and entry is saved in float_array
  ! if not found: found is set to .false. (float_array is not allocated)
  subroutine floatn_find(search_string, float_array, search_start, search_end, found)
    character(*), intent(in) :: search_string
    integer, intent(in) :: search_start, search_end
    real(8), intent(inout), allocatable :: float_array(:)
    logical, intent(out) :: found

    integer :: cnt
    character(len=256) :: str_original
    character(len=256) :: str_2save

    found = .false.
    do i=search_start,search_end
      if (index(trim(file_save(i)),trim(search_string)) .ne. 0) then
        str_temp=file_save(i)
        pst=scan(str_temp,separator)

        if (pst > 0) then ! found the separator
          str_test = trim(adjustl(str_temp(:pst-1))) ! left side of separator
          str_temp = trim(adjustl(str_temp(pst+1:))) ! right side of separator
          if (trim(str_test) == trim(search_string)) then ! detailed comparison
            str_original = str_temp

            ! we count the number of coefficients
            ! this behaves quite weirdly because of the empty space scan of an empty string
            ! definitely works as intended
            cnt = 0
            do
              pst=scan(str_temp,multseparator)
              if (pst == 1) then ! we find the empty space in an empty string at the first position
                exit
              else
                cnt = cnt + 1
                str_temp = trim(adjustl(str_temp(pst+1:)))
              endif
            enddo

            allocate(float_array(cnt))

            ! now that we know the number of floats we want to find
            ! we simply reset the string to its original content
            ! and do the same scan
            str_temp = str_original ! reset
            cnt = 1
            do
              pst=scan(str_temp,multseparator)
              if (pst == 1) then
                exit
              else
                str_2save = trim(adjustl(str_temp(:pst-1)))
                read(str_2save,*) float_array(cnt)
                str_temp = trim(adjustl(str_temp(pst+1:)))
                cnt = cnt + 1
              endif
            enddo

            found = .true.
          endif
        endif
      endif
    enddo
  end subroutine floatn_find

  ! find entry that needs to be saved as boolean (logical)
  ! search_string gets searched between search_start <= ... <= search_end
  ! if found:  found is set to .true. and entry is saved in save_bool
  ! if not found: found is set to .false. (save_bool is not allocated)
  subroutine bool_find(search_string, save_bool, search_start, search_end, found)
    character(*), intent(in)  :: search_string
    logical, intent(inout) :: save_bool
    integer, intent(in) :: search_start, search_end
    logical, intent(out) :: found

    found = .false.
    do i=search_start,search_end
      if (index(trim(file_save(i)),trim(search_string)) .ne. 0) then
        str_temp=file_save(i)
        pst=scan(str_temp,separator)
        if (pst > 0) then ! found the separator
          str_test = trim(adjustl(str_temp(:pst-1))) ! left side of separator
          str_temp = trim(adjustl(str_temp(pst+1:))) ! right side of separator
          if (trim(str_test) == trim(search_string)) then ! detailed comparison
            read(str_temp,*) save_bool
            found = .true.
            exit
          endif
        endif
      endif
    enddo
  end subroutine bool_find

  ! find Group ([...]) in config file
  ! save the start and end line numbers
  ! of the lines below until the next Group in save_start, save_end
  ! these numbers are then used for the above subroutines to look for parameters
  subroutine group_find(search_string, save_start, save_end)
    character(*), intent(in) :: search_string
    integer, intent(out) :: save_start, save_end
    save_start=0
    save_end=0

    do i=1,lines
      if (index(trim(file_save(i)),trim(search_string)) .ne. 0) then
        save_start=i+1
        exit
      endif
    enddo

    if (save_start .ge. 1) then ! group was found
      do i=save_start, lines
        if (index(trim(file_save(i)),'[') .eq. 1) then
          if (index(trim(file_save(i)),'[[') .eq. 1) then ! skip subgroups
            cycle
          endif
          save_end=i-1 ! one above the next session
          exit
        endif
      enddo

      if(save_end .eq. 0) then
        save_end = lines ! if nothing else is found, until the end of the file
      endif

      if (save_start .gt. save_end) then ! group found, but no content
        save_start = -1
      endif
    endif
    return
    ! save_start -> 0: not found; -1: found, but empty
  end subroutine group_find

  ! find SubGroup ([[...]]) in config file within Group range: search_start <= .. <= search_end
  ! save the start and end line numbers
  ! of the lines below until the next (Sub) Group in save_start, save_end
  ! these numbers are then used to look for parameters
  subroutine subgroup_find(search_string, search_start, search_end, save_start, save_end)
    character(*), intent(in) :: search_string
    integer, intent(in) :: search_start, search_end
    integer, intent(out) :: save_start, save_end
    save_start=0
    save_end=0

    do i=search_start, search_end
      if (index(trim(file_save(i)),trim(search_string)) .ne. 0) then
        save_start=i+1
        exit
      endif
    enddo

    if (save_start .ge. 1) then ! subgroup found
      do i=save_start, search_end
        if (index(trim(file_save(i)),'[') .eq. 1) then
          if (index(trim(file_save(i)),'[[[') .eq. 1) then ! skip subsubgroups
            cycle
          endif
          save_end=i-1 ! one above the next session
          exit
        endif
      enddo

      if(save_end .eq. 0) then
        save_end = search_end ! if nothing else is found, until the end of the group
                              ! whose size was already determined by group_find
      endif

      if (save_start .gt. save_end) then ! subgroup found, but no content
        save_start = -1
      endif
    endif
    return
    ! save_start -> 0: not found; -1: found, but empty
  end subroutine subgroup_find

  ! find SubSubGroup ([[[...]]]) in config file within SubGroup range: search_start <= .. <= search_end
  ! save the start and end line numbers
  ! of the lines below until the next (Sub) (Sub) Group in save_start, save_end
  ! these numbers are then used to look for parameters
  subroutine subsubgroup_find(search_string, search_start, search_end, save_start, save_end)
    character(*), intent(in) :: search_string
    integer, intent(in) :: search_start, search_end
    integer, intent(out) :: save_start, save_end
    save_start=0
    save_end=0

    do i=search_start, search_end
      if (index(trim(file_save(i)),trim(search_string)) .ne. 0) then
        save_start=i+1
        exit
      endif
    enddo

    if (save_start .ge. 1) then ! subgroup found
      do i=save_start, search_end
        if (index(trim(file_save(i)),'[') .eq. 1) then
          save_end=i-1 ! one above the next session
          exit
        endif
      enddo

      if(save_end .eq. 0) then
        save_end = search_end ! if nothing else is found, until the end of the group
                              ! whose size was already determined by group_find
      endif

      if (save_start .gt. save_end) then ! subgroup found, but no content
        save_start = -1
      endif
    endif
    return
    ! save_start -> 0: not found; -1: found, but empty
  end subroutine subsubgroup_find

  ! https://github.com/abinitiodga/adga
  ! + some adjustment regarding boundaries
  ! GPLv3 code
  ! spell check within provided lines search_start <= ... search_end
  ! with the help of a dictionary where all allowed keywords are documented
  ! if any inconsistency are found, return an error code + appropriate message
  subroutine spell_check(search_start, search_end, grname, dictionary, er, erstr)
    character(*), intent(in) :: grname
    character(*), intent(in) :: dictionary(:)
    integer, intent(in) :: search_start, search_end
    integer, intent(out) :: er
    character(len=256), intent(out) :: erstr
    character(len=1) :: boundary = "["

    do i=search_start,search_end
      str_temp = file_save(i)

      pst=scan(str_temp,boundary)
      if (pst .eq. 1) then ! found boundary
        er = 0
        exit
      endif

      pst=scan(str_temp,separator)
      if (pst .eq. 0) then
        er = 100
        erstr = 'Variable in '//trim(grname)//' group without argument: '//str_temp
        return
      endif

      str_temp=trim(adjustl(str_temp(:(pst-1))))
      er = 101
      do j = 1,size(dictionary)
        if (str_temp == dictionary(j)) then
          er = 0
          exit
        endif
      enddo

      if (er .ne. 0) then
        erstr = 'Spelling error or unknown variable in '//trim(grname)//' group: '//str_temp
        return
      endif

    enddo
  end subroutine spell_check

end module Mlookup
