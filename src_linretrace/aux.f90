module Maux
#ifdef MPI
  use mpi
#endif
  implicit none
  public

contains

subroutine greeting(ounit)
  implicit none
  integer, intent(in) :: ounit
  write(ounit,*)
  write(ounit,*)'#####################################################'
  write(ounit,*)'#  Lin-ReTraCe -- Linear Response Transport Centre  #'
  write(ounit,*)'#####################################################'
  write(ounit,*)'#       E. Maggio, M. Pickem and J.M. Tomczak       #'
  write(ounit,*)'#####################################################'
  write(ounit,*)
end subroutine main_greeting

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

end module
