module Maux
  use Mmpi_org
#ifdef MPI
  use mpi
#endif
  implicit none
  public

contains

subroutine main_greeting(ounit)
  implicit none
  integer, intent(in) :: ounit
  write(ounit,*)
  write(ounit,*)'#####################################################'
  write(ounit,*)'#  Lin-ReTraCe -- Linear Response Transport Centre  #'
  write(ounit,*)'#####################################################'
  write(ounit,*)'#       M. Pickem, E. Maggio and J.M. Tomczak       #'
  write(ounit,*)'#####################################################'
  write(ounit,*)
end subroutine main_greeting

subroutine log_master(ounit, string)
  implicit none
  integer, intent(in)          :: ounit
  character(len=*), intent(in) :: string
  if (myid .eq. master) then
    write(ounit,*) trim(string)
  endif
end subroutine

subroutine stop_with_message(ounit, erstr, er)
  implicit none
  integer, intent(in)           :: ounit
  integer, intent(in), optional :: er
  character(len=*), intent(in)  :: erstr
  logical                       :: unitopened
#ifdef MPI
  integer :: mpierr
#endif

  write(ounit,*)
  if (present(er)) then
     write(ounit,*) "Error code: ", er
  endif
  write(ounit,*) "Error msg:  ", trim(adjustl(erstr))
  write(ounit,*)
#ifdef MPI
   call sleep(1)
   call MPI_Abort(MPI_COMM_WORLD,-1,mpierr)
#endif
   stop
end subroutine stop_with_message

function to_upper(strIn) result(strOut)
! Adapted from http://www.star.le.ac.uk/~cgp/fortran.html (25 May 2012)
! Original author: Clive Page
  implicit none
  character(len=*), intent(in) :: strIn
  character(len=len(strIn)) :: strOut
  integer :: i,j

  do i = 1, len(strIn)
    j = iachar(strIn(i:i))
    if (j>= iachar("a") .and. j<=iachar("z") ) then
      strOut(i:i) = achar(iachar(strIn(i:i))-32)
    else
      strOut(i:i) = strIn(i:i)
    end if
  end do
end function to_upper

function to_lower(strIn) result(strOut)
! Adapted from http://www.star.le.ac.uk/~cgp/fortran.html (25 May 2012)
! Original author: Clive Page
  implicit none
  character(len=*), intent(in) :: strIn
  character(len=len(strIn)) :: strOut
  integer :: i,j

  do i = 1, len(strIn)
    j = iachar(strIn(i:i))
    if (j>= iachar("a") .and. j<=iachar("z") ) then
      strOut(i:i) = strIn(i:i)
    else
      strOut(i:i) = achar(iachar(strIn(i:i))+32)
    end if
  end do
end function to_lower

end module
