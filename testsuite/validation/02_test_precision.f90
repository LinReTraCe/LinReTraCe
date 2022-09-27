program quadtest
  implicit none

  write(*,*) 'begin quadruple precision tests'
  call quad_real()
  call quad_complex()
  write(*,*) 'end quadruple precision tests'

  contains

  subroutine quad_real()
    implicit none
    real(16) :: quad
    quad = 1.0q0
    write(*,*) '  1.q0 =', quad
  end subroutine

  subroutine quad_complex()
    implicit none
    complex(16) :: quad
    quad = (0,1.q0)
    write(*,*) '  (0,1.q0) =', quad
  end subroutine

end program
