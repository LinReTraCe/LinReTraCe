module Mmymath
  use Mparams
  use nrtype
  use nrutil!, only : swap
  implicit none

  contains
   !***************************
   ! SUBROUTINE DERRICH
   !***************************
   ! the derivative is approximated using the symmetric 2-point expression:
   ! f'(x) ~ A0(h) = 1/2h (f(x+h)-f(x-h)) + O(h^2)
   ! where the truncation error is O(h^2). To improve on this error the 
   ! Richardson extrapolation is evaluated to second order. By doing so
   ! the truncation error is reduced by up to 5 orders of magnitude for
   ! elementary functions for h in the interval (10^-2, 10^-5) rounding
   ! error becomes leading for smaller values of h.
   subroutine derrich(x, y, dy)
      implicit none
      !passed variables
      double precision :: x(:), y(:)
      !output 
      double precision :: dy(:)
      !local variables
      integer :: i
      double precision :: xx         !spacing =h/4
      double precision :: a0, a1, a2 !Richardson's extrapolation of 0, 1st, 2nd order 
      double precision :: b2, b4     !Richardson's extrapolation of 0 order computed at h/2 and h/4
      double precision :: c2         !Richardson's extrapolation if 1st order computed at h/2

      dy(:)=0.0d0
      !around the endpoints use a 4-point evaluation w/o extrapolation
      do i=1,4
         xx = x(i+1)-x(i)
         dy(i) = (-2.5d1*y(i)+4.8d1*y(i+1)-3.6d1*y(i+2)+1.6d1*y(i+3)-3.0d0*y(i+4))/(12.0d0*xx)
      enddo
      do i=size(y)-3,size(y)
         xx = x(i)-x(i-1)
         dy(i) = (2.5d1*y(i)-4.8d1*y(i-1)+3.6d1*y(i-2)-1.6d1*y(i-3)+3.0d0*y(i-4))/(12.0d0*xx)
      enddo
      !for points far from the boundaries use the extrapolation 
      do i=5,size(y)-4
         xx = x(i)-x(i-1)
         !2-point expression for the derivatives for different step sizes h:
         a0 = (y(i+4)-y(i-4))/(8.0d0*xx)  !h=4xx
         b2 = (y(i+2)-y(i-2))/(4.0d0*xx)  !h=2xx
         b4 = (y(i+1)-y(i-1))/(2.0d0*xx)  !h=xx

         a1 = 2.0d0*b2 - a0
         c2 = 2.0d0*b4 - b2

         a2 = (4.0d0*c2 - a1)/3.0d0
         dy(i) = a2
      enddo

   end subroutine

   !***************************
   ! SUBROUTINE QUICKSORT_RV
   !***************************
   ! Numerical Recipes in F90 page 1169
   ! Sorts an input array arr into ascending order 
   ! by using the quicksort algorithm. Input 
   ! array is overwritten upon exit
   ! NN is the size of subarrays sorted by straight insertion
   ! NSTACK is the required auxiliary storage
   subroutine quicksort_rv (arr)
      implicit none
      !declarations
      real(4), intent(inout) :: arr(:)
      integer(I4B), parameter :: NN=15, NSTACK=50
      real(4) :: a
      integer(I4B) :: i, j, k, n, l, r, jstack
      integer(I4B), dimension(NSTACK) :: istack

      n=size(arr) 
      jstack=0
      l=1; r=n
      do
         if(r-l < NN) then
         !insertion sort when the subarray is small enough
            do j=l+1,r
               a=arr(j)
               do i=j-1,l,-1
                  if (arr(i) <= a ) exit
                  arr(i+1)=arr(i)
               enddo
               arr(i+1)=a !this is correct for a fortran90 compiler that keeps the promise of
                          !leaving the counters at the last value they had dunring the cycle +1
                          !not sure if that's the case with us 
            enddo    
            if (jstack==0) return
            ! pop stack and begin new round of partitioning
            r=istack(jstack)
            l=istack(jstack-1)
            jstack=jstack-2
         else
            ! Choose median of left, centre and right elements as partitioning element a
            k=(l+r)/2
            ! arrange so that a(l)<= a(l+1) <= a(r)
            call swap(arr(k), arr(l+1))
            call swap(arr(l), arr(r), arr(l)>arr(r))
            call swap(arr(l+1), arr(r), arr(l+1)>arr(r))
            call swap(arr(l), arr(l+1), arr(l)>arr(l+1))
            ! initialise pointers for partitioning
            i = l+1
            j = r
            a = arr(l+1)   !partitioning element
            do
               !scan up to find element >= a
               do
                  i = i+1
                  if (arr(i) >= a) exit
               enddo
               !scan down to find element <= a
               do
                  j = j-1
                  if (arr(j) <= a) exit
               enddo
               if (j<i) exit !POINTERS CROSSED, EXIT WITH PARTITIONING COMPLETE
               call swap(arr(i),arr(j)) !exchange elements
            enddo
            arr(l+1) = arr(j) !insert partitioning element
            arr(j) = a
            jstack = jstack+2
            ! Push pointers to larger subarray on stack; process smaller subarray immediately
            if (jstack > NSTACK) then
               write(*,*)'quicksort: NSTACK too small'
               STOP
            endif
            if (r-i+1 >= j-l) then
               istack(jstack)=r      
               istack(jstack-1)=i
               r=j-1      
            else
               istack(jstack)=j-1      
               istack(jstack-1)=l
               l=i      
            endif
         endif
      enddo

   end subroutine

end module
