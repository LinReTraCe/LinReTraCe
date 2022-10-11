program precisiontest
  implicit none

  write(*,*) 'begin quadruple precision tests'
  call double()
  call quad()
  write(*,*) 'end quadruple precision tests'
  write(*,*) '-----------'
  write(*,*) 'SUCCESSFUL.'
  write(*,*) '-----------'
  write(*,*)

  contains

  subroutine double()
    implicit none
    integer, parameter :: dp = selected_real_kind(15, 307)

    write(*,'(A)', advance='no') '      double precision:'
    if (dp .eq. 8) then
      write(*,*) 'success'
    else
      write(*,*) 'failed'
      write(*,*) '-------------------------------------------'
      write(*,*) 'Compiler does not support double precision.'
      write(*,*) '-------------------------------------------'
      stop
    endif
  end subroutine

  subroutine quad()
    implicit none
    integer, parameter :: qp = selected_real_kind(33, 4931)

    write(*,'(A)', advance='no') '      quad precision precision:'
    if (qp .eq. 16) then
      write(*,*) 'success'
    else
      write(*,*) 'failed'
      write(*,*) '-----------------------------------------'
      write(*,*) 'Compiler does not support quad precision.'
      write(*,*) '-----------------------------------------'
      stop
    endif
  end subroutine

end program
