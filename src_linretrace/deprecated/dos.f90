module Mdos
  use Mtypes
  use Mparams
  use Maux

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
  maxenergy=edisp%band_original(edisp%nband_max,1,1)
  minenergy=edisp%band_original(1,1,1)

  ! scan the highest band for maxenergy
  ! scan the lowest band for minenergy
  do is = 1,edisp%ispin
    do ik = 2,kmesh%nkp
      if (maxenergy < edisp%band_original(edisp%nband_max,ik,is)) then
         maxenergy = edisp%band_original(edisp%nband_max,ik,is)
      endif
      if (minenergy > edisp%band_original(1,ik,is)) then
         minenergy = edisp%band_original(1,ik,is)
      endif
    enddo
  enddo

  midenergy = (maxenergy+minenergy)/2.d0

  ! generate density of state boundaries
  ! by symmetrically increasing the energy window
  dos%emax= midenergy + (maxenergy-midenergy)*1.3d0
  dos%emin= midenergy + (minenergy-midenergy)*1.3d0
  dos%nnrg= 10001

  ! allocate the arrays and set to 0
  allocate (dos%enrg(dos%nnrg),dos%dos(dos%nnrg, edisp%ispin),dos%nos(dos%nnrg, edisp%ispin))
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
          dos%dos(i,is)=dos%dos(i,is)+((br/pi)*(1.0d0/(((dos%enrg(i)-edisp%band_original(nb,ik,is))**2)+(br**2)))) * &
                     kmesh%weight(ik)
          dos%nos(i,is)=dos%nos(i,is)+(0.5d0 + ((1.0d0/pi)*atan((dos%enrg(i)-edisp%band_original(nb,ik,is))/br))) * &
                     kmesh%weight(ik)
        enddo
      enddo
    enddo
  enddo

  open(unit=11, file='dosnos.dat')
  do i=1,size(dos%enrg)
    write(11,*) dos%enrg(i), (dos%dos(i,is), is=1,edisp%ispin), (dos%nos(i,is), is=1,edisp%ispin)
  enddo
  close(11)

  ! spin prefactors are already taken care of with the weight array

end subroutine !GENDOSEL

subroutine findef(dos, edisp)
  implicit none
  type(dosgrid)    :: dos
  type(energydisp) :: edisp
  !local variables
  integer :: i,j,is
  integer :: pos
  real(8) :: nossum
  real(8) :: ntol


  pos = 0
  do i=1,dos%nnrg
     if ( (sum(dos%nos(i,:)) - edisp%nelect) .ge. 0.d0 ) then ! sign changed
       pos = i
       exit
     endif
  enddo

  if (pos > 0) then ! found a changing sign
     if ( abs(sum(dos%nos(pos,:)) - edisp%nelect) .le. abs(sum(dos%nos(pos-1,:)) - edisp%nelect) ) then
        edisp%efer = dos%enrg(pos)
        i = pos
     else
        edisp%efer = dos%enrg(pos-1)
        i = pos-1
     endif
     write(*,*) 'FINDEF: chemical potential found at: ', dos%enrg(i)
  else
     call stop_with_message(stderr, 'FINDEF: Could not find DFT Fermi level.')
  endif

  ! find band gap and valence band maximum, conduction band minimum
  ! with the help of the number of states (nos)

  allocate(dos%vbm(edisp%ispin))
  allocate(dos%cbm(edisp%ispin))
  allocate(dos%gap(edisp%ispin))

  ! we start with the fermi energy from above (saved in i)
  do is=1,edisp%ispin
    j = i
    do while(abs(dos%nos(j,is) - dos%nos(i,is)) < ntol)
       j = j-1
    enddo
    dos%vbm(is)=dos%enrg(j)

    j = i
    do while(abs(dos%nos(j,is) - dos%nos(i,is)) < ntol)
       j = j+1
    enddo
    dos%cbm(is)=dos%enrg(j)

    dos%gap(is)=dos%cbm(is) - dos%vbm(is)
    if (dos%gap(is) < 2.0d-2) then
      dos%gap(is)=0.0d0
    endif
  enddo

end subroutine

end module Mdos
