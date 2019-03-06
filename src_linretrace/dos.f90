module Mdos
  use Mtypes
  use Mparams

contains

subroutine gendosel(kmesh, edisp, dos)
  implicit none

  type(kpointmesh) :: kmesh
  type(energydisp) :: edisp
  type(dosgrid)    :: dos
  !local variables
  integer :: i, ik, is, ikk, nb
  integer :: iband
  double precision :: br, de !broadening, energy spacing
  double precision :: maxenergy, minenergy, midenergy

  ! find the maximum energy
  maxenergy=edisp%band(edisp%nband_max,1,1)
  minenergy=edisp%band(1,1,1)

  ! scan the highest band for maxenergy
  ! scan the lowest band for minenergy
  do is = 1,edisp%ispin
    do ik = 2,kmesh%nkp
      if (maxenergy < edisp%band(edisp%nband_max,ik,is)) then
         maxenergy = edisp%band(edisp%nband_max,ik,is)
      endif
      if (minenergy > edisp%band(1,ik,is)) then
         maxenergy = edisp%band(1,ik,is)
      endif
    enddo
  enddo

  midenergy = (maxenergy-minenergy)/2.d0

  ! generate density of state boundaries
  ! by symmetrically increasing the energy window
  dos%emax= midenergy + (maxenergy-midenergy)*1.3d0
  dos%emin= midenergy + (minenergy-midenergy)*1.3d0
  dos%nnrg= 10001

  ! allocate the arrays and set to 0
  allocate (dos%enrg(dos%nnrg),dos%dos(dos%nnrg),dos%nos(dos%nnrg))
  dos%enrg= 0.d0
  dos%dos = 0.d0
  dos%nos = 0.d0

  ! get the energy increment along the window fixed in the dos datatype
  de=(dos%emax-dos%emin)/(real(dos%nnrg-1))
  do i=1,dos%nnrg
     dos%enrg(i)=dos%emin+(de*dble(i-1))
  enddo

  ! set artificial laurenthian broadening
  br = 3.0d0*de
  !lorentian bandshape
  do i=1,size(dos%enrg)
    do is=1,edisp%ispin
      do ik=1,kmesh%nkp
        do nb=1,edisp%nband_max
          dos%dos(i)=dos%dos(i)+((br/pi)*(1.0d0/(((dos%enrg(i)-edisp%band(nb,ik,is))**2)+(br**2)))) * &
                     kmesh%weight(ik)
          dos%nos(i)=dos%nos(i)+(0.5d0 + ((1.0d0/pi)*atan((dos%enrg(i)-edisp%band(nb,ik,is))/br))) * &
                     kmesh%weight(ik)
        enddo
      enddo
    enddo
  enddo

  ! spin prefactors are already taken care of with the weight array

end subroutine !GENDOSEL

subroutine findef(dos, edisp)
  implicit none
  type(dosgrid)    :: dos
  type(energydisp) :: edisp
  !local variables
  integer :: i,j
  integer :: pos
  real(8) :: ntol

  pos = 0
  do i=1,dos%nnrg
     if ( (dos%nos(i) - edisp%nelect) .ge. 0.d0 ) then ! sign changed
       pos = i
       exit
     endif
  enddo

  if (pos > 0) then ! found a changing sign
     if ( abs(dos%nos(pos) - edisp%nelect) .le. abs(dos%nos(pos-1) - edisp%nelect) ) then
        edisp%efer = dos%enrg(pos)
        i = pos
     else
        edisp%efer = dos%enrg(pos-1)
        i = pos-1
     endif
  else
     write(*,*) 'FINDEF: No root found for chemical potential, nelect = ', edisp%nelect
     edisp%efer = 0
     return
  endif

  ! find band gap and valence band maximum, conduction band minimum
  ! with the help of the number of states (nos)
  j = i
  do while(abs(dos%nos(j) - dos%nos(i)) < ntol)
     j = j-1
  enddo
  dos%vbm=dos%enrg(j)

  j = i
  do while(abs(dos%nos(j) - dos%nos(i)) < ntol)
     j = j+1
  enddo
  dos%cbm=dos%enrg(j)

  dos%gap=dos%cbm - dos%vbm
  if (dos%gap < 2.0d-2) then
    dos%gap=0.0d0
  endif

end subroutine

end module Mdos
