program precisiontest
  implicit none

  write(*,*) 'begin precision tests'
  call double_real()
  call double_complex()
  call quad_real()
  call quad_complex()
  write(*,*) 'end precision tests'

  contains

  subroutine double_real()
    implicit none
    real(8) :: double
    double = 1.0d0
    write(*,*) '  1d0 =', double
  end subroutine

  subroutine quad_real()
    implicit none
    real(16) :: quad
    quad = 1.0q0
    write(*,*) '  1q0 =', quad
  end subroutine

  subroutine double_complex()
    implicit none
    complex(8) :: double
    double = (0,1.d0)
    write(*,*) '  (0,1d0) =', double
  end subroutine

  subroutine quad_complex()
    implicit none
    complex(16) :: quad
    quad = (0,1.q0)
    write(*,*) '  (0,1q0) =', quad
  end subroutine

end program
