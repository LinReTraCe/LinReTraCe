module Maux
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
   write(ounit,*)'#        Preprocessing band structure data          #'
   write(ounit,*)'#####################################################'
   write(ounit,*)'#        E. Maggio, M. Pickem, J.M. Tomczak         #'
   write(ounit,*)'#####################################################'
   write(ounit,*)
end subroutine preproc_greeting

subroutine stop_with_message(ounit, erstr, er)
   implicit none
   integer, intent(in)           :: ounit
   integer, intent(in), optional :: er
   character(len=*), intent(in)  :: erstr
   write(ounit,*)
   if (present(er)) then
      write(ounit,*) "Error code: ", er
   endif
   write(ounit,*) "Error msg:  ", trim(adjustl(erstr))
   write(ounit,*)
   stop
end subroutine stop_with_message

end module
