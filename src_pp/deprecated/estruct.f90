  subroutine gentbstr(kmesh, edisp)
    implicit none

    type(kpointmesh) :: kmesh
    type(energydisp) :: edisp
    !local variables
    integer :: kx, ky, kz, ktot, nband
    integer :: ik,ikx,iky,ikz, nk,nkx,nky,nkz, iband,idir,idir2 !counters
    double precision :: k(3)

    if (algo%ltetra) then
       ! in order to produce [0,1]
       kx = kmesh%kx+1
       ky = kmesh%ky+1
       kz = kmesh%kz+1
    else
       ! for the standard [0,1) interval
       kx = kmesh%kx
       ky = kmesh%ky
       kz = kmesh%kz
    endif
    ktot = kx*ky*kz
    nband=edisp%nband_max

    allocate(kmesh%k_id(kx, ky, kz))
    allocate(kmesh%k_coord(3,ktot))
    allocate(edisp%band(ktot, nband))
    allocate(edisp%Mopt(3, ktot, nband, nband))
    if (algo%lBfield) allocate(edisp%M2(3, 3, ktot, nband))

    if (allocated(kmesh%k_coord)) kmesh%k_coord=0.d0
    if (allocated(kmesh%k_id))    kmesh%k_id=0
    if (allocated(edisp%band))       edisp%band=0.d0
    if (allocated(edisp%Mopt))       edisp%Mopt=0.d0
    if (allocated(edisp%M2))         edisp%M2  =0.d0

    if(algo%ldebug) then !open files to write
       open(10,file='kcoords')
       open(11,file='ek')
       open(12,file='vk')
       if(algo%lBfield) open(13,file='vkk')
    endif

    ! generate e(k), e'(k) and possible e''(k)
    ik=0
    do ikx=1,kx
       k(1)=1.d0/real(kmesh%kx,kind=8)*real(ikx-1,kind=8)
       do iky=1,ky
          k(2)=1.d0/real(kmesh%ky,kind=8)*real(iky-1,kind=8)
          do ikz=1,kz
             k(3)=1.d0/real(kmesh%kz,kind=8)*real(ikz-1,kind=8)
             ik=ik+1
             kmesh%k_id(ikx,iky,ikz)=ik
             kmesh%k_coord(:,ik)=k
             do iband =1,nband
                edisp%band(ik,iband)=ek_sc(k,iband,edisp)
                do idir=1,3
                   edisp%Mopt(idir,ik,iband,iband)=vk_sc(idir,k,iband,edisp)
                enddo
             enddo
             if (algo%lBfield) then
                do iband =1,nband
                   do idir=1,3
                      edisp%M2(idir,idir,ik,iband)=vkk_sc(idir,idir,k,iband,edisp)
                   enddo
                enddo
             endif

             if(algo%ldebug) then !write to file
                write(10,'(1I14,3F12.6)')ik,k
                write(11,'(1I14,100F16.10)')ik,(edisp%band(ik,iband),iband=1,nband)
                write(12,'(1I14,500F16.10)')ik,((edisp%Mopt(idir,ik,iband,iband),idir=1,3),iband=1,nband)
                if (algo%lBfield) then
                   write(13,'(1I14,1000F16.10)')ik,(((edisp%M2(idir,idir2,ik,iband),idir=1,3),idir2=1,3),iband=1,nband)
                endif
             endif
          enddo !ikz
       enddo    !iky
    enddo       !ikx

    if(algo%ldebug) then !close files written on
       close(10)
       close(11)
       close(12)
       if(algo%lBfield) close(13)
    endif

    if (.not. allocated(edisp%Z)) then
       allocate(edisp%Z(kmesh%ktot,edisp%nband_max))
       edisp%Z=edisp%ztmp
    endif
    if (.not. allocated(edisp%Im)) then
       allocate(edisp%Im(kmesh%ktot,edisp%nband_max))
       edisp%Im=0.0d0
    endif

    ! if we are passing only 1 k-point chances are that we want to study a
    ! flat-band model; in this case the Fermi velocity gets overwritten

    if (kmesh%kx==1 .and. kmesh%ky==1 .and. kmesh%kz==1) then
       if (abs(kmesh%k_coord(1,1)) < 1d-4 .and. abs(kmesh%k_coord(2,1)) < 1.d-4 &
           .and. abs(kmesh%k_coord(3,1)) < 1d-4) then
          write(*,*) 'GENTBSTR: Gamma-point only calculation; assuming flatband model'
          do iband =1,nband
             do idir=1,3
                edisp%Mopt(idir,1,iband,iband) = 1.0d0
             enddo
          enddo
       endif
    endif

  end subroutine !gentbstr

   subroutine genelstr (kmesh, edisp)
     implicit none
     type(kpointmesh) :: kmesh
     type(energydisp) :: edisp

     integer          :: ik,ikx,iky,ikz,nb,nb1

     call getdisp (kmesh, edisp)
     if (algo%ldmft) then
        call getdmft(kmesh, edisp)
     endif
     ! read in the optical matrix elements
     call getopt (kmesh, edisp)
     ! get the rotations and translations from the appropriate w2k files
     ! call getsymop (kmesh, edisp)

     ! if (algo%lsymm) then
        ! generate a reducible kmesh from the set of symmetry operations and the irrek-mesh
        ! also save the rotations required in symm for the optical elements later
        ! call genredk  (kmesh)
        ! generate all the reducible data in one shot
        ! otherwise we calculate it on the fly
        ! TODO: MOVE THIS TO A DEBUGGING FLAG
        ! if (algo%lgenred) then
        !    call genredopt(kmesh, edisp)
        ! endif
     ! else !reducible BZ is already provided by W2k
        ! assign unique identifier to each k-point
        ! if (.not. allocated(kmesh%k_id))     allocate(kmesh%k_id(kmesh%kx, kmesh%ky, kmesh%kz)) ! kx, ky, kz always reducible
        ! if (.not. allocated(symm%symop_id))  allocate(symm%symop_id(2,kmesh%kred))
        ! ik=0
        ! do ikx=1,kmesh%kx
        !    do iky=1,kmesh%ky
        !       do ikz=1,kmesh%kz
        !          ik=ik+1
        !          kmesh%k_id(ikx, iky, ikz)=ik
        !          symm%symop_id(1,ik) = ik  ! one to one mapping
        !          symm%symop_id(2,ik) = 0   ! no operations necessary for already reducible points
        !       enddo
        !    enddo
        ! enddo
     ! endif !lsymm

     ! now we have everything in the reducible datatypes
   end subroutine !genelstr


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! GETSYMOP
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine extracts from the Wien2k files
! information about the symmetry operations that
! acting on the irreducible kmesh generate the reducible
! counterpart (see genredk below)
! I'm not sure about this... they act on the real space optical matrices...
!
  subroutine getsymop(kmesh, ek)
    implicit none

    type(kpointmesh) :: kmesh
    type(energydisp)      :: ek
    integer :: n, i, j, k, rest
    !for the additional detection of roto-translation in non-symmorphic groups
    character (len=80) :: line
    character (len=80) :: substring
    character (len=6)  :: ccrap
    integer :: ix, icrap
    !if centering is required additional temporary vectors are required
    !double precision, allocatable :: Mtmp(:,:,:), Ttmp(:,:)
    double precision, allocatable :: k_coord(:,:)
    double precision, allocatable :: band(:,:)

    substring="NUMBER OF SYMMETRY OPERATIONS"
    ix=0
    open(10,file=trim(adjustl(algo%mysyst))//'.struct',status='old')
    do i=1, 200
       read(10,'(A)') line
       ix = index(line, trim(substring))
       if (ix /= 0) exit
    enddo

    if (ix == 0) then
       write(*,*) 'GETSYMOP: error reading file *.struct'
       STOP
    else
       substring=trim(adjustl(line(:9)))
       read(substring,*) symm%rnsym
       write(*,*) 'GETSYMOP: total number of symmetry operations: ', symm%rnsym
       if (.not. allocated(symm%Msym)) allocate(symm%Msym(3,3,symm%rnsym))
       if (.not. allocated(symm%Tras)) allocate(symm%Tras(3,symm%rnsym))
       do n=1,symm%rnsym
          do j=1,3
             ! this is a bit nasty, but it 100% works
             read(10,'(A)') line
             substring=trim(adjustl(line(:2)))
             read(substring,*) symm%Msym(j,1,n)
             substring=trim(adjustl(line(3:4)))
             read(substring,*) symm%Msym(j,2,n)
             substring=trim(adjustl(line(5:6)))
             read(substring,*) symm%Msym(j,3,n)
             substring=trim(adjustl(line(7:)))
             read(substring,*) symm%Tras(j,n)
          enddo
          read(10,*)
          !write(*,'(I3, I3, 6A, 3F10.8)') icrap, n, ccrap, symm%Tras(:,n)
          !write(*,*) icrap, n, ccrap, symm%Tras(:,n)
       enddo
    endif
    close(10)

    100  FORMAT (I6)
    110  FORMAT (3(3f8.5/))
    ! set the descriptor for nonsymmorphic space groups
    do n=1,symm%rnsym
       if ((symm%Tras(1,n)/=0.0d0).or.(symm%Tras(2,n)/=0.0d0).or.(symm%Tras(3,n)/=0.0d0)) then
          symm%lnsymmr = .true.
       endif
    enddo

    if (symm%lnsymmr) then
       write(*,*) 'GETSYMOP: detected non-symmorphic space group'
    else
       write(*,*) 'GETSYMOP: detected symmorphic space group'
    endif


    if (symm%lnsymmr .and. symm%cntr=="B  ") then
      symm%knsym = symm%rnsym*2
    else
      symm%knsym = symm%rnsym
    endif

    if (algo%lsymm) then
       ! we get the irreducible k-list
       allocate(symm%Msym_reciprocal(3,3,symm%knsym))

       ! now we get the reciprocal space transformations
       substring="SYMMETRY MATRIX"
       ix=0
       open(10,file=trim(adjustl(algo%mysyst))//'.outputkgen',status='old')
       do i=1,200
          read(10,'(A)') line
          ix = index(line, trim(substring))
          if (ix /= 0) exit
       enddo

       if (ix == 0) then
          write(*,*) 'GETSYMOP: error reading file *.outputkgen'
          STOP
       else
          do i=1,symm%knsym/4 ! integer division
             read(10,*) ((symm%Msym_reciprocal(1,j,(i-1)*4+k), j=1,3), k=1,4)
             read(10,*) ((symm%Msym_reciprocal(2,j,(i-1)*4+k), j=1,3), k=1,4)
             read(10,*) ((symm%Msym_reciprocal(3,j,(i-1)*4+k), j=1,3), k=1,4)
             read(10,*)
          enddo
          rest = mod(symm%knsym, 4)
          if (rest>0) then
             ! read in the rest
             read(10,*) ((symm%Msym_reciprocal(1,j,(i-1)*4+k), j=1,3), k=1,rest)
             read(10,*) ((symm%Msym_reciprocal(2,j,(i-1)*4+k), j=1,3), k=1,rest)
             read(10,*) ((symm%Msym_reciprocal(3,j,(i-1)*4+k), j=1,3), k=1,rest)
          endif
      endif

      close(10)

    else
       ! we get the reducible k-list
       allocate(symm%Msym_reciprocal(3,3,1))
       symm%Msym_reciprocal(1,1,1) = 1.d0
       symm%Msym_reciprocal(2,2,1) = 1.d0
       symm%Msym_reciprocal(3,3,1) = 1.d0
    endif


    ! do i=1,symm%nsym
    !   symm%Msym_reciprocal(:,:,i) = matmul(transpose(symm%Msym_reciprocal(:,:,i)), &
    !      matmul(symm%Msym(:,:,i),symm%Msym_reciprocal(:,:,i)))
    ! enddo

    ! write(*,*) symm%Msym_reciprocal

  end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! GENREDK
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine  generates a reducible, uniformily spaced
! k-mesh with elements in the interval (-pi/a, +pi/a) of the BZ
! starting from the irredicible kmesh.
! This grid is then traslated into the interval (0, 2pi/a)
! in the following subroutine
!
  subroutine genredk (kmesh)
    implicit none
    !passed variables
    type(kpointmesh) :: kmesh  ! irreducible k-mesh generated by Wien2k

    !internal variables
    integer :: icrap, itmp,ntmp, i,j,k, ik, isym, isym2, nnsym, ineq
    integer :: g,h,l,m
    integer :: iexist, itest, ibackfold, ioutside
    real(8) :: G0(3,7)
    real(8) :: Gshift(3,27)
    real    :: rcrap
    double precision, allocatable :: tmpkall(:,:), tmpkall2(:,:)
    integer, allocatable          :: tmpoall(:,:), tmpoall2(:,:)
    double precision              :: tmpk(3), tmpk2(3), tmpk3(3)
    integer, allocatable :: counter(:)

    ! number of k-points generated by the symmetry operations acting on the irreducible mesh
    ! that is ... without the translations afterwards if we have a nonsymmorphic
    ! crystal structure
    itmp=kmesh%ktot*symm%knsym

    data G0/ &
      &      0,0,1, 0,1,0, 1,0,0, 0,1,1, 1,1,0, 1,0,1, 1,1,1 /

    ! activate the internal loop for non-symmorphic groups
    if (symm%lnsymmr) then
       nnsym = symm%rnsym
    else
       nnsym = 1
    endif

    allocate(counter(kmesh%ktot))
    counter = 0

    ineq = 0
    ! counter for the unique (inequivalent) k-points
    ! temporary arrays which gets extended on the fly
    ! if we find another unique k-point
    ! essentially we simulate here std:vector<int> from c++
    ! only much worse


    do ik=1,kmesh%ktot  !loop over irredk
       do isym=1,symm%knsym ! refers to the reciprocal matrices
          do j=1,3
             tmpk(j) = ((kmesh%k_coord(1,ik)*symm%Msym_reciprocal(j,1,isym))) &
                & + ((kmesh%k_coord(2,ik)*symm%Msym_reciprocal(j,2,isym))) &
                & + ((kmesh%k_coord(3,ik)*symm%Msym_reciprocal(j,3,isym)))
          enddo

          ! set numerically zero values to absolute zero
          ! + backfolding
          do j=1,3
             if(abs(tmpk(j)) .le. 1.d-5) tmpk(j) = 0.d0
             do while (tmpk(j) < 0.d0)
               tmpk(j) = tmpk(j) + 1.d0
             enddo
             do while (tmpk(j) >= 1.d0)
               tmpk(j) = tmpk(j) - 1.d0
             enddo
          enddo


          iexist = 0
          do j=1,ineq
             if ( (abs(  tmpk(1) - tmpkall(1,j) ) < 1.d-1/real(itmp)) &
             & .and. (abs(  tmpk(2) - tmpkall(2,j) ) < 1.d-1/real(itmp)) &
             & .and. (abs(  tmpk(3) - tmpkall(3,j) ) < 1.d-1/real(itmp)) ) then
                iexist=1
                exit
             endif
          enddo

          ! sanity check
          ! initialize the arrays
          if ((iexist==0) .and. (ineq==0)) then
             counter(ik) = counter(ik) + 1
             allocate(tmpkall(3,1), tmpoall(2,1))
             tmpkall(:,1) = tmpk
             tmpoall(1,1) = ik; tmpoall(2,1) = isym
             ineq = ineq+1
          ! add the new vector to the end of the existing arrays
          else if (iexist==0) then
             counter(ik) =  counter(ik) + 1
             ! k-point extension
             allocate(tmpkall2(3,ineq))
             allocate(tmpoall2(2,ineq))
             tmpkall2  = tmpkall ! save temporarily
             tmpoall2  = tmpoall
             ineq = ineq+1
             deallocate(tmpkall, tmpoall)
             allocate(tmpkall(3,ineq))  ! allocate a buffe which is of size +1
             allocate(tmpoall(2,ineq))
             tmpkall(:,:(ineq-1)) = tmpkall2 ! save back
             tmpoall(:,:(ineq-1)) = tmpoall2
             tmpkall(:,ineq) = tmpk ! add at the end
             tmpoall(1,ineq) = ik; tmpoall(2,ineq) = isym ! we only save the rotation number
             deallocate(tmpkall2, tmpoall2)
          endif
       enddo
    enddo


    write(*,*)
    do i=1,kmesh%ktot
      write(*,*) i, " : ", counter(i)
    enddo
    write(*,*)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! What we've got right now is a k-mesh with equivalent endpoints in a given
    ! direction if that direction is sampled by an even number of k-points
    ! either one has to assign weight 1/2 to those or remove the double counting
    ! when the BZ is translated (in trnredk below)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! mP : ???
    ! save the number of redkp found into the data structure
    ! redkm%kx = kmesh%kx+(1-mod(kmesh%kx,2))
    ! redkm%ky = kmesh%ky+(1-mod(kmesh%ky,2))
    ! redkm%kz = kmesh%kz+(1-mod(kmesh%kz,2))
    ! redkm%ktot = max(redkm%kx*redkm%ky*redkm%kz,ik)

    ! if ((ineq .ne. kmesh%kred) .and. (ineq*2 .ne. kmesh%kred)) then
    !    write(*,*) 'GENREDK: the number of k-points generated by symmetry is inconsistent',&
    !       ineq,kmesh%kred,kmesh%kx,kmesh%ky,kmesh%kz
    !    write(*,*) 'Constructed k-points:'
    !    write(*,*)
    !    do i=1,size(tmpkall,2)
    !       write(*,*) tmpkall(:,i), tmpoall(:,i)
    !    enddo
    !    write(*,*)
    !    STOP
    ! else if (ineq*2 .eq. kmesh%kred) then
    !    write(*,*) 'GENREDK: Modulo ( kpoints, generated points) = ', mod(kmesh%kred,ineq)
    !    if (algo%ltetra) then
    !       write(*,*) 'GENREDK: Cannot create tetrahedrons from half-filled Brillouin zone'
    !       stop
    !    endif
    ! endif
    write(*,*) 'GENREDK: total number of constructed reducible points: ', ineq, ', ', kmesh%kred

    write(*,*) 'Constructed k-points:'
    write(*,*)
    do i=1,size(tmpkall,2)
       write(*,*) tmpkall(:,i), tmpoall(:,i)
    enddo

    deallocate(kmesh%k_coord)
    allocate(kmesh%k_coord(3,kmesh%kred))
    allocate(kmesh%k_id(kmesh%kx, kmesh%ky, kmesh%kz))
    allocate(symm%symop_id(2,kmesh%kred))

    ! save the new coordinates into the data structure
    kmesh%k_coord = tmpkall
    symm%symop_id = tmpoall

    deallocate(tmpkall, tmpoall)

    ! do the mapping
    kmesh%k_id = 0
    do ik = 1,kmesh%kred
       kmesh%k_id(nint(kmesh%kx*kmesh%k_coord(1,ik))+1, &
                  nint(kmesh%ky*kmesh%k_coord(2,ik))+1, &
                  nint(kmesh%kz*kmesh%k_coord(3,ik))+1) = ik
    enddo

  end subroutine !GENREDK

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! GENREDOPT
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine generates the optical matrix elements
! starting from those read in on the irreducible k-mesh
!
  subroutine genredopt (kmesh, edisp)
    implicit none
    !passed variables
    type(kpointmesh) :: kmesh
    type(energydisp)      :: edisp
    !internal variables
    integer :: i, j, l, k, ik, isym, iks, nb, nb2, ikk
    real(8), allocatable :: Mopttmp(:,:,:,:)
    real(8), allocatable :: bandtmp(:,:)
    real(8), allocatable :: Ztmp(:,:)
    real(8), allocatable :: Imtmp(:,:)
    real(8), allocatable :: Mtmp(:,:,:,:)


    if (lat%lortho) then
       allocate(Mopttmp(3,kmesh%kred,edisp%nbopt_min:edisp%nbopt_max,edisp%nbopt_min:edisp%nbopt_max))
    else
       allocate(Mopttmp(6,kmesh%kred,edisp%nbopt_min:edisp%nbopt_max,edisp%nbopt_min:edisp%nbopt_max))
    endif

    if (lat%lortho) then
       ! create the cubic optical elements in Mopttmp
       ! after that save them in the original edisp%Mopt array.
       Mopttmp = 0.d0
       do ik=1,kmesh%kred
          ikk = symm%symop_id(1,ik)
          isym= symm%symop_id(2,ik)
          ! if (isym > symm%rnsym) then
          !   isym = isym - symm%rnsym
          ! endif
          ! check eq. 13.16 (pg 479) in "Symmetry and Condensed Matter Physics A Computational Approach" by M. El-Batanouny, F. Wooten, CUP
                ! Mopttmp(1,ik,nb,nb2) = Mopttmp(1,ik,nb,nb2) + edisp%Mopt(j,ikk,nb,nb2)*symm%Msym(j,1,isym)*symm%Msym(j,1,isym)
                ! Mopttmp(2,ik,nb,nb2) = Mopttmp(2,ik,nb,nb2) + edisp%Mopt(j,ikk,nb,nb2)*symm%Msym(j,2,isym)*symm%Msym(j,2,isym)
                ! Mopttmp(3,ik,nb,nb2) = Mopttmp(3,ik,nb,nb2) + edisp%Mopt(j,ikk,nb,nb2)*symm%Msym(j,3,isym)*symm%Msym(j,3,isym)

                ! if its unitary this is correct
                ! otherwise not so much
                ! Mopttmp(1,ik,nb,nb2) = Mopttmp(1,ik,nb,nb2) + &
                !         edisp%Mopt(j,ikk,nb,nb2)*symm%Msym_reciprocal(1,j,isym)*symm%Msym_reciprocal(j,1,isym)
                ! Mopttmp(2,ik,nb,nb2) = Mopttmp(2,ik,nb,nb2) + &
                !         edisp%Mopt(j,ikk,nb,nb2)*symm%Msym_reciprocal(2,j,isym)*symm%Msym_reciprocal(j,2,isym)
                ! Mopttmp(3,ik,nb,nb2) = Mopttmp(3,ik,nb,nb2) + &
                !         edisp%Mopt(j,ikk,nb,nb2)*symm%Msym_reciprocal(3,j,isym)*symm%Msym_reciprocal(j,3,isym)

          do nb2=edisp%nbopt_min,edisp%nbopt_max
          do nb=edisp%nbopt_min,edisp%nbopt_max
             ! this is wrong
             ! Mopttmp(1,ik,nb,nb2) = Mopttmp(1,ik,nb,nb2) + &
             !         edisp%Mopt(j,ikk,nb,nb2)*symm%Msym(j,1,isym)*symm%Msym(j,1,isym)
             ! Mopttmp(2,ik,nb,nb2) = Mopttmp(2,ik,nb,nb2) + &
             !         edisp%Mopt(j,ikk,nb,nb2)*symm%Msym(j,2,isym)*symm%Msym(j,2,isym)
             ! Mopttmp(3,ik,nb,nb2) = Mopttmp(3,ik,nb,nb2) + &
             !         edisp%Mopt(j,ikk,nb,nb2)*symm%Msym(j,3,isym)*symm%Msym(j,3,isym)

             ! also wrong
             ! Mopttmp(1,ik,nb,nb2) = Mopttmp(1,ik,nb,nb2) + &
             !         edisp%Mopt(j,ikk,nb,nb2)*symm%Msym(1,j,isym)*symm%Msym(1,j,isym)
             ! Mopttmp(2,ik,nb,nb2) = Mopttmp(2,ik,nb,nb2) + &
             !         edisp%Mopt(j,ikk,nb,nb2)*symm%Msym(2,j,isym)*symm%Msym(2,j,isym)
             ! Mopttmp(3,ik,nb,nb2) = Mopttmp(3,ik,nb,nb2) + &
             !         edisp%Mopt(j,ikk,nb,nb2)*symm%Msym(3,j,isym)*symm%Msym(3,j,isym)

             ! Mopttmp(1,ik,nb,nb2) = Mopttmp(1,ik,nb,nb2) + &
             !         edisp%Mopt(j,ikk,nb,nb2)*symm%Msym_reciprocal(1,j,isym)*symm%Msym_reciprocal(1,j,isym)
             ! Mopttmp(2,ik,nb,nb2) = Mopttmp(2,ik,nb,nb2) + &
             !         edisp%Mopt(j,ikk,nb,nb2)*symm%Msym_reciprocal(2,j,isym)*symm%Msym_reciprocal(2,j,isym)
             ! Mopttmp(3,ik,nb,nb2) = Mopttmp(3,ik,nb,nb2) + &
             !         edisp%Mopt(j,ikk,nb,nb2)*symm%Msym_reciprocal(3,j,isym)*symm%Msym_reciprocal(3,j,isym)

             ! Mopttmp(1,ik,nb,nb2) = Mopttmp(1,ik,nb,nb2) + &
             !         edisp%Mopt(j,ikk,nb,nb2)*symm%Msym_reciprocal(j,1,isym)*symm%Msym_reciprocal(j,1,isym)
             ! Mopttmp(2,ik,nb,nb2) = Mopttmp(2,ik,nb,nb2) + &
             !         edisp%Mopt(j,ikk,nb,nb2)*symm%Msym_reciprocal(j,2,isym)*symm%Msym_reciprocal(j,2,isym)
             ! Mopttmp(3,ik,nb,nb2) = Mopttmp(3,ik,nb,nb2) + &
             !         edisp%Mopt(j,ikk,nb,nb2)*symm%Msym_reciprocal(j,3,isym)*symm%Msym_reciprocal(j,3,isym)

             Mopttmp(1,ik,nb,nb2) = edisp%Mopt(1,ikk,nb,nb2)
             Mopttmp(2,ik,nb,nb2) = edisp%Mopt(2,ikk,nb,nb2)
             Mopttmp(3,ik,nb,nb2) = edisp%Mopt(3,ikk,nb,nb2)
          enddo
          enddo
       enddo
       deallocate(edisp%Mopt)
       allocate(edisp%Mopt(3,kmesh%kred,edisp%nbopt_min:edisp%nbopt_max,edisp%nbopt_min:edisp%nbopt_max))
       edisp%Mopt = Mopttmp
       deallocate(Mopttmp)

    ! The new matrix elements are generated by a matrix product of the type:
    !
    !  M'_ab = sum (R_aj M_jl R_bl)
    !          j,l
    ! along the cartesian indices. The Mopt matrix is indexed following the
    ! convention in the optic output of Wien2k:
    !  1  4  5
    !     2  6
    !        3
    ! whereas the symmetry matrices have rows and columns indexed independently.
    ! To treat the two objects in a consistent way the optical matrix elements are
    ! copied over to temporary matrices
    else
       allocate(Mtmp(3,3,edisp%nbopt_min:edisp%nbopt_max,edisp%nbopt_min:edisp%nbopt_max))
       Mopttmp = 0.d0
       do ik=1,kmesh%kred
          ikk = symm%symop_id(1,ik)
          isym= symm%symop_id(2,ik)

          ! copy over the optical matrix elements
          Mtmp(1,1,:,:) = edisp%Mopt(1,ikk,:,:)
          Mtmp(2,2,:,:) = edisp%Mopt(2,ikk,:,:)
          Mtmp(3,3,:,:) = edisp%Mopt(3,ikk,:,:)
          Mtmp(1,2,:,:) = edisp%Mopt(4,ikk,:,:)
          Mtmp(1,3,:,:) = edisp%Mopt(5,ikk,:,:)
          Mtmp(2,3,:,:) = edisp%Mopt(6,ikk,:,:)
          Mtmp(2,1,:,:) = Mtmp(1,2,:,:)
          Mtmp(3,1,:,:) = Mtmp(1,3,:,:)
          Mtmp(3,2,:,:) = Mtmp(2,3,:,:)

          ! do the contraction with the symmetry matrices
          ! check eq. 13.16 (pg 479) in "Symmetry and Condensed Matter Physics A Computational Approach" by M. El-Batanouny, F. Wooten, CUP
          do j = 1,3
             do l = 1,3
                Mopttmp(1,ik,:,:) = Mopttmp(1,ik,:,:) + symm%Msym(j,1,isym)*Mtmp(j,l,:,:)*symm%Msym(l,1,isym)
                Mopttmp(2,ik,:,:) = Mopttmp(2,ik,:,:) + symm%Msym(j,2,isym)*Mtmp(j,l,:,:)*symm%Msym(l,2,isym)
                Mopttmp(3,ik,:,:) = Mopttmp(3,ik,:,:) + symm%Msym(j,3,isym)*Mtmp(j,l,:,:)*symm%Msym(l,3,isym)
                Mopttmp(4,ik,:,:) = Mopttmp(4,ik,:,:) + symm%Msym(j,1,isym)*Mtmp(j,l,:,:)*symm%Msym(l,2,isym)
                Mopttmp(5,ik,:,:) = Mopttmp(5,ik,:,:) + symm%Msym(j,1,isym)*Mtmp(j,l,:,:)*symm%Msym(l,3,isym)
                Mopttmp(6,ik,:,:) = Mopttmp(6,ik,:,:) + symm%Msym(j,2,isym)*Mtmp(j,l,:,:)*symm%Msym(l,3,isym)
             enddo
          enddo
       enddo
       deallocate (Mtmp)
       deallocate(edisp%Mopt)
       allocate(edisp%Mopt(6,kmesh%kred,edisp%nbopt_min:edisp%nbopt_max,edisp%nbopt_min:edisp%nbopt_max))
       edisp%Mopt = Mopttmp
       deallocate(Mopttmp)
    endif

    ! save data
    allocate(Ztmp(kmesh%kred, edisp%nband_max))
    allocate(Imtmp(kmesh%kred, edisp%nband_max))
    allocate(bandtmp(kmesh%kred, edisp%nband_max))
    Ztmp    = edisp%Z
    Imtmp   = edisp%Im
    bandtmp = edisp%band

    ! increase size of arrays
    deallocate(edisp%Z, edisp%Im, edisp%band)
    allocate(edisp%band(kmesh%kred,edisp%nband_max))
    allocate(edisp%Z(kmesh%kred,edisp%nband_max))
    allocate(edisp%Im(kmesh%kred,edisp%nband_max))

    !now map the dispersion energies from the old grid to the new one (the energies do not change with the symmetry operation)
    do ik=1,size(symm%symop_id,2)
       ikk = symm%symop_id(1,ik) ! this counter runs over the irrekpoints
       isym= symm%symop_id(2,ik) ! not really necessary here

       edisp%band(ik,:) = bandtmp(ikk,:)
       edisp%Z(ik,:)    = Ztmp(ikk,:)
       edisp%Im(ik,:)   = Imtmp(ikk,:)
    enddo

    ! once generated, we have to set the symmetry operations
    ! to a default reducible grid
    ! necessary to access the right data from here on after
    do i=1,size(symm%symop_id,2)
       symm%symop_id(1,i) = i
    enddo
    symm%symop_id(2,:) = 0 ! this is what we use as signifier that we don't have to apply rotations

  end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! TRNREDK
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine translates the reducible
! k-mesh with elements in the interval (-pi/a, +pi/a)
! into the interval (0, 2pi/a), taking care also of
! the band structure and of the optical matrix elements
!

!
! mP: deprecated for the time being
!
  subroutine trnredk (irrkm, redkm, eredk)
    implicit none
    !passed variables
    type(kpointmesh) :: irrkm
    type(kpointmesh) :: redkm
    type(energydisp) :: eredk

    !local variables
    integer :: i, j, ik, ikx, iky, ikz, ibn, ibn2
    integer :: nk, nkx, nky, nkz
    integer :: ntmp, itest, iexist
    integer :: offdia !off-diagonal terms in the Mopt matrix?
    double precision :: dk(3), tmp1, tmp2
    double precision :: maxkx, maxky, maxkz
    double precision, allocatable :: cktmp(:,:), cktmp2(:,:) ! temporary k-point coordinates arrays
    double precision, allocatable :: bstmp(:,:), bstmp2(:,:) ! temporary bandstructure arrays
    double precision, allocatable :: Imtmp(:,:), Imtmp2(:,:) ! temporary bandstructure arrays
    double precision, allocatable :: Ztmp(:,:),  Ztmp2(:,:)  ! temporary bandstructure arrays
    double precision, allocatable :: Motmp(:,:,:,:), Motmp2(:,:,:,:) ! temporary optical transition matrices

    if (lat%lortho) then
       offdia=0  !no off-diagonal terms for cubic systems
    else
       offdia=1  !off-diagonal terms for non-cubic systems
    endif
    !allocation temporary arrays
    allocate(bstmp(redkm%ktot,eredk%nband_max)); allocate(bstmp2(redkm%ktot,eredk%nband_max))
    allocate(Imtmp(redkm%ktot,eredk%nband_max)); allocate(Imtmp2(redkm%ktot,eredk%nband_max))
    allocate(Ztmp(redkm%ktot,eredk%nband_max));  allocate(Ztmp2(redkm%ktot,eredk%nband_max))
    allocate(cktmp(3,redkm%ktot)); allocate(cktmp2(3,redkm%ktot))
    if (lat%lortho) then
       allocate(Motmp (3,redkm%ktot,eredk%nbopt_min:eredk%nbopt_max,eredk%nbopt_min:eredk%nbopt_max))
       allocate(Motmp2(3,redkm%ktot,eredk%nbopt_min:eredk%nbopt_max,eredk%nbopt_min:eredk%nbopt_max))
    else
       allocate(Motmp (6,redkm%ktot,eredk%nbopt_min:eredk%nbopt_max,eredk%nbopt_min:eredk%nbopt_max))
       allocate(Motmp2(6,redkm%ktot,eredk%nbopt_min:eredk%nbopt_max,eredk%nbopt_min:eredk%nbopt_max))
    endif
    cktmp = 0.0d0; cktmp2 = 0.0d0
    bstmp = 0.0d0; bstmp2 = 0.0d0
    Imtmp = 0.0d0; Imtmp2 = 0.0d0
    Ztmp  = 0.0d0; Ztmp2  = 0.0d0
    Motmp = 0.0d0; Motmp2 = 0.0d0

    do i=1,3
       tmp1 = 0.0d0
       do ik=2,redkm%ktot
          tmp2 = abs(redkm%k_coord(i,1) - redkm%k_coord(i,ik))
          if ((tmp2>tmp1)) tmp1=tmp2
       enddo
       dk(i)=tmp1
    enddo

    if (dk(1)==0 .or. dk(2)==0 .or. dk(3)==0) then
       write(*,*) 'TRNREDK: can not translate the k-mesh', dk(:)
       STOP
    endif

    !translate the points in the negative BZ into the +'ve interval
    do ik=1,redkm%ktot
       if ((abs(redkm%k_coord(1,ik)) > 1.0d0) .or. (abs(redkm%k_coord(2,ik)) > 1.0d0) .or.(abs(redkm%k_coord(3,ik)) > 1.0d0)) then
          STOP 'TRNREDK: something is seriously wrong here (e.g. improper traslation bigger than lattice constant)'
       endif
       do i = 1,3
          if (redkm%k_coord(i,ik)<0.0d0) then
             cktmp(i,ik)= 1.0d0+redkm%k_coord(i,ik)
          else
             cktmp(i,ik)= redkm%k_coord(i,ik)
          endif
       enddo
    enddo

    ! storing everything in the first temporary arrays
    do ibn =1,eredk%nband_max
       do ik=1,redkm%ktot
          bstmp(ik,ibn)= eredk%band(ik,ibn)
          Imtmp(ik,ibn)= eredk%Im(ik,ibn)
          Ztmp(ik,ibn) = eredk%Z(ik,ibn)
       enddo
    enddo

    !copy the optical transition matrix elements
    do ibn2=eredk%nbopt_min,eredk%nbopt_max
       do ibn =eredk%nbopt_min,eredk%nbopt_max
          do ik=1,redkm%ktot
             do i = 1,3+(offdia*3)
                Motmp(i,ik,ibn,ibn2) = eredk%Mopt(i,ik,ibn,ibn2)
             enddo
          enddo
       enddo
    enddo

    !!!!!!!!!!!!!!!!TEST
    !do ik=1,redkm%ktot
    !  write(776,'(A,I6,3f8.4)')'KP ',ik, cktmp(1,ik),cktmp(2,ik),cktmp(3,ik)
    !  do ibn=eredk%nbopt_min,eredk%nbopt_max
    !    do ibn2=ibn,eredk%nbopt_max
    !       write(776,170)ibn,ibn2,Motmp(1,ik,ibn,ibn2),Motmp(2,ik,ibn,ibn2),Motmp(3,ik,ibn,ibn2)
    !    enddo
    !  enddo
    !enddo
    !!!!!!!!!!!!!!!!TEST END

    !deallocate the old datastructure
    deallocate(redkm%k_coord)
    deallocate(redkm%k_id)
    deallocate(eredk%band)
    deallocate(eredk%Im)
    deallocate(eredk%Z)
    deallocate(eredk%Mopt)

    !!!!!!!!! K-POINT CLEAN UP !!!!!!!!!!!!!!!!
    !check which points in the redBZ have been saved already
    ntmp=1   !include the gamma point
    ik=1
    do i = ntmp, redkm%ktot
       do itest=i+1, redkm%ktot
          if ( (abs(  cktmp(1,itest) - cktmp(1,i) ) <0.1d0/real(redkm%kx)) &
          & .and. (abs(  cktmp(2,itest) - cktmp(2,i) ) <0.1d0/real(redkm%ky)) &
          & .and. (abs(  cktmp(3,itest) - cktmp(3,i) ) <0.1d0/real(redkm%kz))  ) cycle
          ! now it has found the first value of itest that corresponds to a different k-point from the
          ! the starting one

          ! now need to copy over this k-point corresponding to itest into cktmp2 BUT
          ! making sure that it hasn't appeared there already
          iexist=0
          do j=1,ik
             if ( (abs(  cktmp(1,itest) - cktmp2(1,j) ) <0.1d0/real(redkm%kx)) &
             & .and. (abs(  cktmp(2,itest) - cktmp2(2,j) ) <0.1d0/real(redkm%ky)) &
             & .and. (abs(  cktmp(3,itest) - cktmp2(3,j) ) <0.1d0/real(redkm%kz)) ) iexist=1
          enddo
          if (iexist==0) then
             ik=ik+1
             cktmp2(1,ik) = cktmp(1,itest)
             cktmp2(2,ik) = cktmp(2,itest)
             cktmp2(3,ik) = cktmp(3,itest)
             do ibn=1,eredk%nband_max
                bstmp2(1,ibn) = bstmp(1,ibn) !the gamma point was left out (min(itest)=2)
                bstmp2(ik,ibn)= bstmp(itest,ibn)
                Imtmp2(1,ibn) = Imtmp(1,ibn)
                Imtmp2(ik,ibn)= Imtmp(itest,ibn)
                Ztmp2(1,ibn)  = Ztmp(1,ibn)
                Ztmp2(ik,ibn) = Ztmp(itest,ibn)
                do ibn2=eredk%nbopt_min,eredk%nbopt_max
                    if ((ibn>eredk%nbopt_max) .or. (ibn<eredk%nbopt_min)) cycle
                    Motmp2(1,1,ibn,ibn2) = Motmp(1,1,ibn,ibn2)
                    Motmp2(1,ik,ibn,ibn2)= Motmp(1,itest,ibn,ibn2)
                    Motmp2(2,1,ibn,ibn2) = Motmp(2,1,ibn,ibn2)
                    Motmp2(2,ik,ibn,ibn2)= Motmp(2,itest,ibn,ibn2)
                    Motmp2(3,1,ibn,ibn2) = Motmp(3,1,ibn,ibn2)
                    Motmp2(3,ik,ibn,ibn2)= Motmp(3,itest,ibn,ibn2)
                    if (offdia == 1) then
                       Motmp2(4,1,ibn,ibn2) = Motmp(4,1,ibn,ibn2)
                       Motmp2(4,ik,ibn,ibn2)= Motmp(4,itest,ibn,ibn2)
                       Motmp2(5,1,ibn,ibn2) = Motmp(5,1,ibn,ibn2)
                       Motmp2(5,ik,ibn,ibn2)= Motmp(5,itest,ibn,ibn2)
                       Motmp2(6,1,ibn,ibn2) = Motmp(6,1,ibn,ibn2)
                       Motmp2(6,ik,ibn,ibn2)= Motmp(6,itest,ibn,ibn2)
                    endif
                enddo !ibn2
             enddo    !ibn
             ntmp=itest
          endif
       enddo
    enddo
    nkx = irrkm%kx
    nky = irrkm%ky
    nkz = irrkm%kz
    nk  =  nkx*nky*nkz
    if (ik .ne. nk) then
       write(*,*) 'TRNREDK: the number of k-points after clean-up is inconsistent',ik, nk, nkx, nky, nkz
       !do i=1,nk
       !   write(777,*)i, cktmp2(1,i),cktmp2(2,i),cktmp2(3,i)
       !enddo
       !STOP
       nk=ik
       if (algo%ltetra) STOP 'TRNREDK: tetrahedron method can not be used'
    endif

    !allocate the final datastructure with the correct dimensions
    allocate(redkm%k_coord(3,nk))
    allocate(redkm%k_id(nkx,nky,nkz))
    allocate(eredk%band(nk,eredk%nband_max))
    allocate(eredk%Im(nk,eredk%nband_max))
    allocate(eredk%Z(nk,eredk%nband_max))
    if (lat%lortho) then
       allocate(eredk%Mopt(3,nk,eredk%nbopt_min:eredk%nbopt_max,eredk%nbopt_min:eredk%nbopt_max))
    else
       allocate(eredk%Mopt(6,nk,eredk%nbopt_min:eredk%nbopt_max,eredk%nbopt_min:eredk%nbopt_max))
    endif

    !initialise the variables
    redkm%k_id(:,:,:)  = 0
    redkm%k_coord(:,:) = 0.0d0
    eredk%band(:,:)    = 0.0d0
    eredk%Im(:,:)      = 0.0d0
    eredk%Z(:,:)       = 0.0d0
    eredk%Mopt(:,:,:,:)= 0.0d0


    !update the k-points of the final mesh
    redkm%kx=nkx; redkm%ky=nky; redkm%kz=nkz
    redkm%ktot=nk

    !save the new coordinates, band dispersion and optical matrix elements into the data structure
    do i=1,3
       do ik=1,nk
          redkm%k_coord(i,ik)=cktmp2(i,ik)
       enddo
    enddo
    do ibn=1,eredk%nband_max
       do ik=1,nk
          eredk%band(ik,ibn)= bstmp2(ik,ibn)
          eredk%Im(ik,ibn)  = Imtmp2(ik,ibn)
          eredk%Z(ik,ibn)   = Ztmp2(ik,ibn)
       enddo
    enddo

    do ibn2=eredk%nbopt_min,eredk%nbopt_max
       do ibn=eredk%nbopt_min,eredk%nbopt_max
          do ik=1,nk
             do i=1,3+(offdia*3)
                eredk%Mopt(i,ik,ibn,ibn2)=Motmp2(i,ik,ibn,ibn2)
             enddo
          enddo
       enddo
    enddo


    deallocate (cktmp); deallocate (cktmp2)
    deallocate (bstmp); deallocate (bstmp2)
    deallocate (Imtmp); deallocate (Imtmp2)
    deallocate (Ztmp) ; deallocate (Ztmp2)
    deallocate (Motmp); deallocate (Motmp2)

    if (algo%ltetra) then
    !the k-point coordinates are all scrambled together, one has to construct
    !some identifier, otherwise is is difficult to add the terminal k-point
    !(equivalent to the Gamma point if we were in 1d)
       allocate(cktmp(3,nk))
       maxkx=0.d0; maxky=0.d0; maxkz=0.d0;
       do ik=1,nk
          if(redkm%k_coord(1,ik) > maxkx) maxkx= redkm%k_coord(1,ik)
          if(redkm%k_coord(2,ik) > maxky) maxky= redkm%k_coord(2,ik)
          if(redkm%k_coord(3,ik) > maxkz) maxkz= redkm%k_coord(3,ik)
       enddo

       !construct a linearly spaced grid equivalent to the one given
       ik=0
       do ikx=nkx,1,-1
          do iky=nky,1,-1
             do ikz=nkz,1,-1
                ik=ik+1
                cktmp(1,ik)= maxkx-real(ikx-1)*(maxkx/real(nkx-1))
                cktmp(2,ik)= maxky-real(iky-1)*(maxky/real(nky-1))
                cktmp(3,ik)= maxkz-real(ikz-1)*(maxkz/real(nkz-1))
             enddo
          enddo
       enddo

       itest=0
       do ikx=1,nkx
          do iky=1,nky
             do ikz=1,nkz
                itest=itest+1
                do ik=1,nk
                   if ( ( abs(cktmp(1,itest) - redkm%k_coord(1,ik)) <1.0d-1/real(redkm%ktot) ) &
                   .and. (abs(cktmp(2,itest) - redkm%k_coord(2,ik)) <1.0d-1/real(redkm%ktot) ) &
                   .and. (abs(cktmp(3,itest) - redkm%k_coord(3,ik)) <1.0d-1/real(redkm%ktot) ) ) then
                      redkm%k_id(ikx,iky,ikz)=ik
                   endif
                enddo
             enddo
          enddo
       enddo

       deallocate (cktmp)
    endif

    !!!!!!!!!!!!!!!!TEST
    !do ikx=1,nkx
    !  do iky=1,nky
    !    do ikz=1,nkz
    !      ik=redkm%k_id(ikx,iky,ikz)
    !      write(777,170)'KP ',ik, ikx, iky, ikz, redkm%k_coord(1,ik),redkm%k_coord(2,ik),redkm%k_coord(3,ik)
    !      do ibn=eredk%nbopt_min,eredk%nbopt_max
    !        do ibn2=ibn,eredk%nbopt_max
    !           write(777,171)ibn,ibn2,eredk%Mopt(1,ik,ibn,ibn2),eredk%Mopt(2,ik,ibn,ibn2),eredk%Mopt(3,ik,ibn,ibn2)
    !        enddo
    !      enddo
    !    enddo
    !  enddo
    !enddo
    !170  FORMAT  (A,4(I4,X),3(E12.6,X))
    !171  FORMAT  (2(I4,X),3(F12.6,X))
    !!!!!!!!!!!!!!!!TEST END

  end subroutine ! trnredk

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! GENFULKM
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine generates a reducible BZ that
! extends the k-point interval from [0,1) to
! [0,1] in each direction (i.e. if we were in
! 1D this would mean to include the Gamma point
! twice). This is necessary if we are going
! to construct tetrahedra on this k-mesh.
!
  subroutine genfulkm(kmesh, edisp)
    implicit none
    ! passed variables
    type(kpointmesh) :: kmesh
    type(energydisp) :: edisp
    ! local variables
    integer :: i, ik, ikx, iky, ikz, ibn, ibn2
    integer :: offdia !off-diagonal terms in the Mopt matrix?
    double precision :: dk(3), tmp1, tmp2

    integer, allocatable :: symop_tmp(:,:)
    integer, allocatable :: kid_tmp(:,:,:)
    real(8), allocatable :: kcoord_tmp(:,:)


    if (lat%lortho) then
       offdia=0  !no off-diagonal terms for cubic systems
    else
       offdia=1  !off-diagonal terms for non-cubic systems
    endif

    ! increase the size of the kmesh identifiers
    allocate(kid_tmp(kmesh%kx,kmesh%ky,kmesh%kz))
    allocate(kcoord_tmp(3,kmesh%kred))
    kid_tmp = kmesh%k_id
    kcoord_tmp = kmesh%k_coord

    deallocate(kmesh%k_id, kmesh%k_coord)
    allocate(kmesh%k_id(kmesh%kx+1, kmesh%ky+1, kmesh%kz+1))
    allocate(kmesh%k_coord(3,kmesh%kful))

    kmesh%k_id(:kmesh%kx, :kmesh%ky, :kmesh%kz) = kid_tmp
    kmesh%k_coord(:,:kmesh%kred) = kcoord_tmp
    deallocate(kid_tmp, kcoord_tmp)

    ! increase the size of symm%symop_id
    allocate(symop_tmp(2,kmesh%kred))
    symop_tmp = symm%symop_id
    deallocate(symm%symop_id)
    allocate(symm%symop_id(2,kmesh%kful))
    symm%symop_id(:,:kmesh%kred) = symop_tmp
    deallocate(symop_tmp)

    !Generate the minimum displacement necessary to
    !obtain the full BZ by adding it to the existing mesh
    do i=1,3
       tmp1 = 1.0d0
       do ik=2,kmesh%kred
          tmp2 = abs(kmesh%k_coord(i,1) - kmesh%k_coord(i,ik))
          if ((tmp2>0.0d0) .and. (tmp2<tmp1)) tmp1=tmp2
       enddo
       dk(i)=tmp1
    enddo

    if (dk(1)==0 .or. dk(2)==0 .or. dk(3)==0) then
       write(*,*) 'GENFULKM: cannot extend the k-mesh', dk(:)
       STOP
    endif

    !add on the terminal point to the BZ
    ik=kmesh%kred
    !select one face (z=const)
    do ikx=1,kmesh%kx
       do iky=1,kmesh%ky
          ik = ik+1
          kmesh%k_id(ikx,iky,kmesh%kz+1)=ik
          kmesh%k_coord(1,ik)=kmesh%k_coord(1,kmesh%k_id(ikx,iky,kmesh%kz))
          kmesh%k_coord(2,ik)=kmesh%k_coord(2,kmesh%k_id(ikx,iky,kmesh%kz))
          kmesh%k_coord(3,ik)=kmesh%k_coord(3,kmesh%k_id(ikx,iky,kmesh%kz))+dk(3)

          ! get the symmetry identifier from the equivalent point
          symm%symop_id(1,ik)=symm%symop_id(1,kmesh%k_id(ikx,iky,1))
          symm%symop_id(2,ik)=symm%symop_id(2,symm%symop_id(1,ik))
       enddo
    enddo
    !select the second face (y=const)
    do ikx=1,kmesh%kx
       do ikz=1,kmesh%kz+1
          ik = ik+1
          kmesh%k_id(ikx,kmesh%ky+1,ikz)=ik
          kmesh%k_coord(1,ik)=kmesh%k_coord(1,kmesh%k_id(ikx,kmesh%ky,ikz))
          kmesh%k_coord(2,ik)=kmesh%k_coord(2,kmesh%k_id(ikx,kmesh%ky,ikz))+dk(2)
          kmesh%k_coord(3,ik)=kmesh%k_coord(3,kmesh%k_id(ikx,kmesh%ky,ikz))

          ! get the symmetry identifier from the equivalent point
          symm%symop_id(1,ik)=symm%symop_id(1,kmesh%k_id(ikx,1,ikz))
          symm%symop_id(2,ik)=symm%symop_id(2,symm%symop_id(1,ik))
       enddo
    enddo
    !select the third face (x=const)
    do iky=1,kmesh%ky+1
       do ikz=1,kmesh%kz+1
          ik = ik+1
          kmesh%k_id(kmesh%kx+1,iky,ikz)=ik
          kmesh%k_coord(1,ik)=kmesh%k_coord(1,kmesh%k_id(kmesh%kx,iky,ikz))+dk(1)
          kmesh%k_coord(2,ik)=kmesh%k_coord(2,kmesh%k_id(kmesh%kx,iky,ikz))
          kmesh%k_coord(3,ik)=kmesh%k_coord(3,kmesh%k_id(kmesh%kx,iky,ikz))

          ! get the symmetry identifier from the equivalent point
          symm%symop_id(1,ik)=symm%symop_id(1,kmesh%k_id(1,iky,ikz))
          symm%symop_id(2,ik)=symm%symop_id(2,symm%symop_id(1,ik))
       enddo
    enddo

    !TODO: no second derivatives of the energy because at the moment there is no way to compute them

  end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! GENTETRA
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine produces the cutting of a uniformily spaced
! k-mesh (mesh) into tetrahedra.
! algo%tbind mesh=kmesh this includes the terminal BZ points
! but the counters are not updated (so the same mesh can be used
! for searching the chemical potential)
! algo%w2k mesh=fulkm
! REFERENCE: PRB (1994) 49,16223-16233
!
subroutine gentetra (kmesh, thdr)
 implicit none
 type(kpointmesh) :: kmesh
 type(tetramesh)  :: thdr
!local variables
 double precision :: x, y, z, xx, edmin, edmax, edgmax, edgmin
 !double precision :: x1, y1, z1, x2, y2, z2, x3, y3, z3
 double precision :: bk(3,3)           ! cartesian coordinates for the reciprocal lattice
 double precision :: p(3,4)            ! positions of the tetrahedron's vertices in cartesian coordinates
 double precision :: vltet             ! tetrahedron's volume in reciprocal space
 double precision :: vlcub             ! cube's volume in reciprocal space
 integer :: ntetd                      ! maximum number of tetrahedra allowed (it shouldn't exceed 6*nk)
 integer :: ntet                       ! actual number of tetrahedra
 integer :: kcut0(3,4,6), kcut(3,4,6)  ! partitioning of the kmesh
 integer :: itet, ic                   ! counter over #of tetrahedra and # of corners
 integer :: lx, ly, lxx, lyy           ! switches over taking one tetrahedron corner or another
 integer :: i1, i2, i3, j1, j2, j3, k1, k2, k3, i, ii, n, j
 integer :: nkx, nky, nkz
 integer, allocatable :: idtet(:,:)    ! identifies the tetrahedron
 integer :: iq(4),imc(0:1,0:1,0:1)
 integer :: icube, iccor
 integer :: ik
 double precision :: tmp1, tmp2

  ! ntetd is the maximum number of tetrahedra for the given mesh
  ! given that there are 6 tetraedra in each cubic cell it is reasonable to set
  ntetd=6*kmesh%kful
  allocate (idtet(0:4,ntetd))
  nkx=kmesh%kx+1; nky=kmesh%ky+1; nkz=kmesh%kz+1
      data kcut0/ &
     &         0,0,0, 0,1,0, 1,1,0, 1,1,1,  0,0,0, 1,0,0, 1,1,0, 1,1,1, &
     &         0,0,0, 1,0,0, 1,0,1, 1,1,1,  0,0,0, 0,1,0, 0,1,1, 1,1,1, &
     &         0,0,0, 0,0,1, 0,1,1, 1,1,1,  0,0,0, 0,0,1, 1,0,1, 1,1,1 /

 ! need to generate the cartesian basis for the cell
 ! for a simple cubic lattice the reciprocal lattice is also cubic and shrunk by a factor 2pi/alat
 bk=0.d0
 !generalisation to tetragonal and orthorhombic cases:

 ! mP note: this is definitely wrong
 ! bk(1,1)=lat%a(2)*lat%a(3)*(2.d0*pi/lat%alat)
 ! bk(2,2)=lat%a(3)*lat%a(1)*(2.d0*pi/lat%alat)
 ! bk(3,3)=lat%a(1)*lat%a(2)*(2.d0*pi/lat%alat)

 bk(1,1)=2.d0*pi/(lat%a(1))
 bk(2,2)=2.d0*pi/(lat%a(2))
 bk(3,3)=2.d0*pi/(lat%a(3))

! setting up the tetrahedra will be done cutting a cell with eight
! corners into six tetrahedra. the edges of the tetrahedra are given by
! three edges, two face diagonals and one space diagonal of the cell.
! giving the space diagonal, the way how to choose the rest is uniquely
! determined...  but now there are four different possibilities how to
! choose the space diagonal! Prefer the choice which gives the shortest
! edges for all tetrahedra ('the most compact tetrahedra') - just to
! avoid very long 'interpolation distances'
      lxx=0
      lyy=0
      edgmax=1.d30
      edgmin=0.d0
      icube=1
! for the four choices ...
      do lx=0,1
       do ly=0,1
! ... we set up the 'trial division' of a given cell into 6 tetrahedra:
         do itet=1,6
          do ic=1,4
            kcut(1,ic,itet)=kcut0(1,ic,itet)
            kcut(2,ic,itet)=kcut0(2,ic,itet)
            kcut(3,ic,itet)=kcut0(3,ic,itet)
            if (lx==1) kcut(1,ic,itet)=1-kcut0(1,ic,itet)
            if (ly==1) kcut(2,ic,itet)=1-kcut0(2,ic,itet)
          enddo
         enddo
         edmin=1.d30
         edmax=0.d0
! for this trial setting, loop over all tetrahedra ...,
         do 4 itet=1,6
! ... set up the cartesian coordinates of the four corner points ...,
            do 2 ic=1,4
               p(1,ic)=kcut(1,ic,itet)*bk(1,1)+ &
     &                 kcut(2,ic,itet)*bk(1,2)+ &
     &                 kcut(3,ic,itet)*bk(1,3)
               p(2,ic)=kcut(1,ic,itet)*bk(2,1)+ &
     &                 kcut(2,ic,itet)*bk(2,2)+ &
     &                 kcut(3,ic,itet)*bk(2,3)
               p(3,ic)=kcut(1,ic,itet)*bk(3,1)+ &
     &                 kcut(2,ic,itet)*bk(3,2)+ &
     &                 kcut(3,ic,itet)*bk(3,3)
    2       continue
! ... and get the shortest and longest distance between two points in
! each tetrahedron (minimum/maximum taken over all tetrahedra ...):
            do i=1,3
             do j=i+1,4
               xx=anrm2(p(1,i)-p(1,j),p(2,i)-p(2,j),p(3,i)-p(3,j))
               edmax=max(edmax,xx)
               edmin=min(edmin,xx)
             enddo
            enddo
    4    continue
! now look at the global maximum: have we found a cut with smaller
! maximum distance between two points within one tetrahedron than
! before? if yes: store it
         if (edmax<edgmax) then
            lxx=lx
            lyy=ly
            edgmax=edmax
            edgmin=edmin
         end if
       enddo
      enddo
! now set up the 'correct' cut giving the most compact tetrahdra ... :
      do itet=1,6
       do ic=1,4
         kcut(1,ic,itet)=kcut0(1,ic,itet)
         kcut(2,ic,itet)=kcut0(2,ic,itet)
         kcut(3,ic,itet)=kcut0(3,ic,itet)
         if (lxx==1) kcut(1,ic,itet)=1-kcut0(1,ic,itet)
         if (lyy==1) kcut(2,ic,itet)=1-kcut0(2,ic,itet)
       enddo
      enddo
! now start searching the tetrahedra ... :
      ntet=0
! for all k-points
      do i3=1,nkz
       do i2=1,nky
        do i1=1,nkx
         iccor=0
! set up microcell of 8 k-points (= 8 corners of unit cell of k-mesh):
         do k1=0,1
          j1=i1+k1
          if (j1 > nkx) cycle
          do k2=0,1
           j2=i2+k2
           if (j2 > nky) cycle
           do k3=0,1
            j3=i3+k3
            if (j3 > nkz) cycle
            iccor=iccor+1
! get the identifiers (the k-point connected to i1,i2,i3):
            ! for that particular choice of k1, k2, k3
            imc(k1,k2,k3)=kmesh%k_id(j1,j2,j3)
            !write(35,*) 'cube',icube, iccor, mesh%k_id(j1,j2,j3),j1,j2,j3,mesh%k_coord(:,mesh%k_id(j1,j2,j3))
            if (iccor==8) then
             vlcub=(abs(kmesh%k_coord(3,kmesh%k_id(j1,j2,j3))-kmesh%k_coord(3,kmesh%k_id(j1,j2,j3-1))))**3
             ! write(36,*) 'cube',icube, vlcub
            endif
           enddo
          enddo
         enddo
         icube = icube + 1
! from this cell we can cut out six tetrahedra:
         do 13 itet=1,6
            if (iccor < 8) cycle
! set the 4 corners (identifiers) of the actual tetrahedron:
            do 8 ic=1,4
               k1=kcut(1,ic,itet)
               k2=kcut(2,ic,itet)
               k3=kcut(3,ic,itet)
               iq(ic)=imc(k1,k2,k3)
    8       continue
! order the identifiers of the corners ...
            do j=1,3
             do i=1,4-j
               if (iq(i)>iq(i+1)) then
                  ii=iq(i)
                  iq(i)=iq(i+1)
                  iq(i+1)=ii
               end if
             enddo
            enddo
! first tetrahedron
            if (ntet==0) goto 11
! now test all tetrahedra found previously:
            do 10 n=1,ntet
               if ((idtet(1,n)==iq(1)) &
                 & .and.(idtet(2,n)==iq(2)) &
                 & .and.(idtet(3,n)==iq(3)) &
                 & .and.(idtet(4,n)==iq(4))) then
! we have found the same combination previously, so increment the
! counter for this type of tetrahedron ...
                  idtet(0,n)=idtet(0,n)+1
! ... and go to the next tetrahedron:
                  goto 13
               end if
   10       continue
! new tetrahedron found if arriving here:
   11       continue
! count it, ...
            ntet=ntet+1
! ... store the corner coordiantes (identifier) ...
            do i=1,4
               idtet(i,ntet)=iq(i)
            enddo
! ... and initialize the counter for this new type of tetrahedron:
            idtet(0,ntet)=1
   13    continue
        enddo
       enddo
      enddo

      if (ntet > ntetd) then
         write(*,*) 'GENTETRA: more tetrahedra than the number allowed',ntet, ntetd
         STOP
      endif
! now tell us the result:
      write(*,15) ntet,kmesh%kful
   15 FORMAT(1x,'GENTETRA: found ',i6,' inequivalent tetrahedra from ',i8,' k-points' )

!!!!!!! COMPUTE THE TETRAHEDRON'S VOLUME !!!!!!!!!!
! The expression for the volume was retrieved on http://mathworld.wolfram.com/Tetrahedron.html
! and it require computing the determinant of the square matrix formed with the coordinates of the tetrahedron's
! vertices: the 4th column has all 1 entries. The determinant is multiplied by 1/3!
      if (.not. allocated(thdr%vltet)) allocate (thdr%vltet(ntet))
      do itet=1,ntet
         thdr%vltet(itet) = 0.0d0
         !access the coordinates of the vertices
         do ic=1,4
            do i=1,3
               p(i,ic)=kmesh%k_coord(i,idtet(ic,itet))
            enddo
         enddo
!!!!!!!!!!!!!!!!!!!!!
!Wolfram expression
         x=p(1,3)-p(1,4)
         y=p(2,3)-p(2,4)
         z=p(3,3)-p(3,4)
         vltet=(x*((p(2,1)*p(3,2))-(p(2,2)*p(3,1)))-y*((p(1,1)*p(3,2))-(p(1,2)*p(3,1)))+z*((p(1,1)*p(2,2))-(p(1,2)*p(2,1))))/6.0d0
         thdr%vltet(itet) = vltet
         !write(*,*) x,y,z,vltet

         x=p(1,1)-p(1,2)
         y=p(2,1)-p(2,2)
         z=p(3,1)-p(3,2)
         vltet=(x*((p(2,3)*p(3,4))-(p(2,4)*p(3,3)))-y*((p(1,3)*p(3,4))-(p(1,4)*p(3,3)))+z*((p(1,3)*p(2,4))-(p(1,4)*p(2,3))))/6.0d0
         thdr%vltet(itet) = abs(thdr%vltet(itet) + vltet)
!!!!!!!!!!!!!!!!!!!!!
!Wikipedia expression
!          x1=p(1,1)-p(1,4)
!          y1=p(2,1)-p(2,4)
!          z1=p(3,1)-p(3,4)
!          x2=p(1,2)-p(1,4)
!          y2=p(2,2)-p(2,4)
!          z2=p(3,2)-p(3,4)
!          x3=p(1,3)-p(1,4)
!          y3=p(2,3)-p(2,4)
!          z3=p(3,3)-p(3,4)
!          vltet=abs(x1*((y2*z3)-(y3*z2)) -x2*((y1*z3)-(y3*z1)) +x3*((y1*z2)-(y2*z1)) )/6.0d0
!          !thdr%vltet(itet) = vltet

      enddo

!store the result away in the datatype
      if (.not. allocated(thdr%idtet)) allocate (thdr%idtet(0:4, ntet))
      thdr%ntet = ntet
      do itet=1, ntet
         thdr%idtet(0,itet) = idtet(0,itet)
         thdr%idtet(1,itet) = idtet(1,itet)
         thdr%idtet(2,itet) = idtet(2,itet)
         thdr%idtet(3,itet) = idtet(3,itet)
         thdr%idtet(4,itet) = idtet(4,itet)
      enddo
      deallocate(idtet)


      if (algo%ldebug) then
        open(unit=33, file='thdr_1')
        open(unit=34, file='thdr_2')
        thdr%vltot=0.0d0
        tmp1=0.0d0
        do itet=1,thdr%ntet
           thdr%vltot=thdr%vltot+thdr%vltet(itet)
           do ik =1,4
              write(33,*) itet, kmesh%k_coord(:,thdr%idtet(ik,itet))
           enddo
           tmp1=tmp1+thdr%idtet(0,itet)
           if (algo%ldebug) write(34,*) itet, thdr%idtet(0,itet), thdr%vltet(itet)
        enddo
        tmp2=((2*pi)**3)/lat%vol
        write(*,*) 'GENTETRA: tetrahedra volume',thdr%vltot,thdr%vltot*((2*pi)**3)/lat%vol,tmp2
      endif

     return

end subroutine  !gentetra

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! INTETRA
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This routine computes the DOS using the tetrahedron method
! Also the integrated DOS (number of states, NOS) is computed;
! expressions are given in appendix C and A respectively of
! PRB (1994) 49, 16223-16233
!
subroutine intetra (kmesh, edisp, thdr, dos)
  implicit none

   type(kpointmesh) :: kmesh
   type(energydisp) :: edisp
   type(tetramesh)  :: thdr
   type(dosgrid)    :: dos

   !local variables
   integer :: i, j, i00, itet, nb, istart, istop, iband
   integer :: iq(4)
   double precision :: de, ec(4), ec1(4), es(4)
   double precision :: e1, e2, e3, e4
   double precision :: c0, c1, c2, c3, cc12, cc34  ! constants
   double precision :: wthdr   ! weight of the tetrahedron
   double precision :: eact, x ! free energy variables
   double precision :: adddos  ! accumulation variable for the dos
   double precision :: maxenergy

   ! SANITY CHECKS
   if (kmesh%kful<4) then
      write(*,*)'INTETRA: tetrahedron method fails (number of k-points < 4)',kmesh%kful
      STOP
   endif

   ! find the energy interval
   maxenergy=0.d0
   do iband=1,edisp%nband_max
      ! band_fill_value is large enough that we don't have to worry about it in the tb case
      if ((edisp%band(1,iband) < band_fill_value) .and. (maxenergy < abs(edisp%band(1,iband)))) then
         maxenergy = abs(edisp%band(1,iband))
      endif
   enddo
   dos%emax= 2.d0*maxenergy
   dos%emin=-dos%emax
   dos%nnrg= 5001

   ! initialize arrays for dos/number of states
   !write(*,*)'INTETRA: constructing energy mesh'
   if (.not. allocated(dos%enrg)) allocate (dos%enrg(dos%nnrg))
   if (.not. allocated(dos%dos )) allocate (dos%dos(dos%nnrg))
   if (.not. allocated(dos%nos )) allocate (dos%nos(dos%nnrg))
   dos%enrg= 0.d0
   dos%dos = 0.d0
   dos%nos = 0.d0
! get the energy increment along the window fixed in the dos datatype
   de=(dos%emax-dos%emin)/(real(dos%nnrg-1))
   do i=0,dos%nnrg-1
      dos%enrg(i+1)=dos%emin+(de*real(i))
   enddo
! loop over tetrahedra:
   do itet=1,thdr%ntet
! get the four corner points:
      iq(1) = symm%symop_id(1,thdr%idtet(1,itet))
      iq(2) = symm%symop_id(1,thdr%idtet(2,itet))
      iq(3) = symm%symop_id(1,thdr%idtet(3,itet))
      iq(4) = symm%symop_id(1,thdr%idtet(4,itet))
      wthdr = thdr%vltet(itet)

      do nb=1,edisp%nband_max
        if (edisp%band(iq(1),nb)> band_fill_value) cycle
        if (edisp%band(iq(2),nb)> band_fill_value) cycle
        if (edisp%band(iq(3),nb)> band_fill_value) cycle
        if (edisp%band(iq(4),nb)> band_fill_value) cycle
        ! get the band energy at each corner of the tetrahedron:
        ec(1) = edisp%band(iq(1),nb)
        ec(2) = edisp%band(iq(2),nb)
        ec(3) = edisp%band(iq(3),nb)
        ec(4) = edisp%band(iq(4),nb)

! sort the energies at the four corners (array ec) into array es
         do i=1,4
            ec1(i)=ec(i)
         enddo
         do 3 i=1,4
            i00=1
            do j=2,4
               if (ec1(j)<ec1(i00)) i00=j
            enddo
            es(i)=ec1(i00)
            ec1(i00)=1.d30
     3   continue

! lowest energy still above emax ---> no contributions to dos/nos:
         if (es(1)>=(dos%emax+0.00000001d0*de)) return
! highest energy still below emin ---> no contribution to dos and
! contribution of complete tetrahedron to nos (1*wthdr):
         if (es(4)<=(dos%emin-0.00000001d0*de)) then
            do 4 i=1,dos%nnrg
               dos%nos(i)=dos%nos(i)+wthdr
       4    continue
            return
         end if
! now the rest
         e1=es(1)
         e2=es(2)
         e3=es(3)
         e4=es(4)
         !write(*,*) 'itet',itet,'nband',nb, e1,e2,e3,e4
! now get the minimum and maximum index for the range we have to update
! dos(i) and nos(i) [so that emin>e(istart) and emax<e(istop)] ... :
         istart=max((int((e1-dos%emin)/de-0.00000001d0)),1)
         istart=min(istart,dos%nnrg)
         istop=min((int((e4-dos%emin)/de+0.00000001d0))+2,dos%nnrg)
         istop=max(istop,1)

! constants occuring in the integration formulas:
         if ((e3-e2)>0.d0) then
            c3= wthdr*(e1+e2-e3-e4)/((e3-e1)*(e4-e1)*(e3-e2)*(e4-e2))
            c2= wthdr*3.d0/((e3-e1)*(e4-e1))
         else
            c3=0.d0
            c2=0.d0
         endif
         c1= c2*(e2-e1)
         c0= c1*(e2-e1)/3.d0
         if ((e2-e1)>0.d0) then
            cc12= wthdr/((e2-e1)*(e3-e1)*(e4-e1))
         else
            cc12= 0.d0
         endif
         if ((e4-e3)>0.d0) then
            cc34= wthdr/((e3-e4)*(e2-e4)*(e1-e4))
         else
            cc34=0.d0
         endif

! LOOP OVER FREE ENERGY VARIABLE
         do 7 i=istart,istop
            eact=dos%emin+(de*real(i-1))
            adddos=0.d0
! case eact between e2,e3:
            if ((e2<eact).and.(eact<=e3)) then
               x=eact-e2
               adddos=c1+x*(2.d0*c2+3.d0*x*c3)
               dos%nos(i)=dos%nos(i)+c0+x*(c1+x*(c2+x*c3))
! case eact between e1,e2:
            else if ((e1<eact).and.(eact<=e2)) then
               x=eact-e1
               adddos=3.d0*cc12*x*x
               dos%nos(i)=dos%nos(i)+cc12*x*x*x
! case eact between e3,e4:
            else if ((e3<eact).and.(eact<=e4)) then
               x=eact-e4
               adddos=-3.d0*cc34*x*x
               dos%nos(i)=dos%nos(i)+wthdr-cc34*x*x*x
! case eact greater than e4 (might probably happen for i=istop):
            else if (e4<=eact) then
               dos%nos(i)=dos%nos(i)+wthdr
            end if
            dos%dos(i)=dos%dos(i)+adddos
       7 continue

! all energies higer than e(istop) give same contribution to nos as
! in the case eact greater than e4 above:
         if (istop<dos%nnrg) then
            do 10 i=istop+1,dos%nnrg
               dos%nos(i)=dos%nos(i)+wthdr
      10    continue
         end if

      enddo  ! over bands
   enddo     ! over tetrahedra

   ! spin multiplicity
   dos%dos = 2.d0 * dos%dos
   dos%nos = 2.d0 * dos%nos

 return

end subroutine   !INTETRA


  ! OLD FINDEF with RIDDER
  !subroutine findef(dos, ek)
  !  implicit none

  !  type(dosgrid) :: dos
  !  type(energydisp)   :: ek
  !  !local variables
  !  double precision :: F(4), P(4)
  !  double precision :: s
  !  double precision :: psave, ptol, ntol
  !  integer  :: I(4), iter, maxiter, itmp

  !  ! initialise the varibles
  !  I(1)= 1
  !  I(2)= dos%nnrg
  !  P(1)= dos%enrg(1)
  !  P(2)= dos%enrg(dos%nnrg)
  !  F(1)= dos%nos(1)-ek%nelect
  !  F(2)= dos%nos(dos%nnrg)-ek%nelect
  !  ptol   =  1.0d-16
  !  psave  = -1.1d30
  !  maxiter= 60

  !  do iter = 1, maxiter
  !     itmp = I(1)+I(2)
  !     I(3) = int(itmp/2)
  !     P(3) = dos%enrg(I(3))
  !     F(3) = dos%nos(I(3))-ek%nelect
  !     s = sqrt((F(3)**2)-(F(1)*F(2)))
  !     if (s==0.0d0) then
  !        write(*,*) 'Error in Ridders search for Fermi level'
  !        write(*,*) 'ITER', iter, 'x1', P(1),'  x2',P(2),'  x3', P(3)
  !        write(*,*) 'ITER', iter, 'F1', F(1),'  F2',F(2),'  F3', F(3)
  !        goto 400
  !     endif
  !     I(4) = I(3)+(I(3)-I(1))*int(sign(1.0d0,F(1)-F(2))*F(3)/s)
  !     P(4) = dos%enrg(I(4))

  !     if(abs(P(4)-psave)<=ptol) goto 400
  !     psave= P(4)
  !     F(4) = dos%nos(I(4))-ek%nelect
  !     if (F(4) ==0.0d0) goto 400
  !     if (sign(F(3), F(4)) /= F(3)) then
  !     !change of sign btw x3 and x4 then reduce search interval
  !        I(1)  = I(3)
  !        P(1)  = P(3)
  !        F(1)  = F(3)
  !        I(2)  = I(4)
  !        P(2)  = P(4)
  !        F(2)  = F(4)
  !     else if (sign(F(1), F(4)) /= F(1)) then
  !     !change of sign btw x1 and x4 then reduce search interval
  !        I(2)  = I(4)
  !        P(2)  = P(4)
  !        F(2)  = F(4)
  !     else if (sign(F(2), F(4)) /= F(2)) then
  !     !change of sign btw x2 and x4 then reduce search interval
  !        I(1)  = I(4)
  !        P(1)  = P(4)
  !        F(1)  = F(4)
  !     endif
  !     !condition for termination
  !     if (abs(P(2)-P(1)) <= ptol) goto 400
  !  enddo ! over number of iterations
  !  write (*,*) 'here 6'
  !  400 if (iter == maxiter) write(*,*) 'Ridders seach might not have converged'
  !  ek%efer=P(4)
  !  !find the band gap
  !  ntol=4.0d-2
  !  I(3)  = I(4)
  !  P(3)  = P(4)
  !  F(3)  = F(4)

  !  do while(abs(F(3)-F(4)) < ntol)
  !     I(3) = I(3)-1
  !     P(3) = dos%enrg(I(3))
  !     F(3) = dos%nos(I(3))-ek%nelect
  !  enddo
  !  dos%vbm=P(3)

  !  I(3)  = I(4)
  !  P(3)  = P(4)
  !  F(3)  = F(4)
  !  do while(abs(F(3)-F(4)) < ntol)
  !     I(3) = I(3)+1
  !     P(3) = dos%enrg(I(3))
  !     F(3) = dos%nos(I(3))-ek%nelect
  !  enddo
  !  dos%cbm=P(3)
  !  dos%gap=dos%cbm - dos%vbm
  !  if (dos%gap < 2.0d-2) dos%gap=0.0d0

  !end subroutine ! FINDEF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! INTERPTRA_RE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This routine linearly interpolates the product
! of the optical matrix elemets with the transport kernel
! defined on the vertices of a tetrahedron, by computing a
! weighted sum of these values according to eq.6 in
! PRB (1994) 49, 16223-16233. The expressions for the
! weights are given in app. B therein.
!
subroutine interptra_re (iT, itet, mu, lBoltz, mesh, ek, thdr, sct, resp, hpresp )
  implicit none

   class(dp_resp)   :: resp
   type(kpointmesh) :: mesh
   type(energydisp)      :: ek
   type(tetramesh)  :: thdr
   type(scatrate)   :: sct
   type(qp_resp),optional :: hpresp
   integer, intent(in)           :: iT     !temperature index
   integer, intent(in)           :: itet   !tetrahedron identifier
   double precision, intent (in) :: mu     !chemical potential
   logical :: lBoltz
!local variables
   integer :: ik, ikk, iktet
   integer :: i, j, i00, nb, nb1, iswitch
   integer :: ix, iy        !polarisation directions in the response functions
   integer :: iq(4), iqs(4) !contains the tetrahedron identifiers before and after the corners have been sorted for increasing energy
   double precision :: ef, ec(4), ec1(4), es(4)
   double precision :: c0, c1, c2, c3, cc1, cc4  !constants
   double precision :: wthdr          !weight of the tetrahedron
   double precision :: x1, x2, x3, x4 !energies zeroed w.r.t. Fermi energy: xj = es(j)-ef
   double precision :: w(4)           !weights of the interpolation formula
   integer :: itmp


   !*******************************************
   ! At the beginning this is almost a literal
   ! copy from the INTETRA routine
   !*******************************************
   ! SANITY CHECK
   if (mesh%ktot<4) then
      write(*,*)'INTERPTRA: tetrahedron method fails (number of k-points < 4)',mesh%ktot
      STOP
   endif

   !initialise the accumulation variables
   resp%s_tet(:,:) = 0.0d0 ; resp%a_tet(:,:) = 0.0d0
   resp%s(:,:,:) = 0.0d0 ; resp%a(:,:,:) = 0.0d0
   if (allocated(resp%sB)) then
      resp%sB_tet(:,:) = 0.0d0 ; resp%aB_tet(:,:) = 0.0d0
      resp%sB(:,:,:) = 0.0d0 ; resp%aB(:,:,:) = 0.0d0
   endif
   if (present(hpresp)) then
      hpresp%s_tet(:,:) = 0.0q0 ; hpresp%a_tet(:,:) = 0.0q0
      hpresp%s(:,:,:) = 0.0q0 ; hpresp%a(:,:,:) = 0.0q0
      if (allocated(hpresp%sB)) then
         hpresp%sB_tet(:,:) = 0.0q0 ; hpresp%aB_tet(:,:) = 0.0q0
         hpresp%sB(:,:,:) = 0.0q0 ; hpresp%aB(:,:,:) = 0.0q0
      endif
   endif

   ef = mu ! as the value of the chemical potential changes, so does the occupation of the tetrahedra
   ! get the four corner points:
   ! thdr -> redk -> data
   iq(1) = symm%symop_id(1,thdr%idtet(1,itet))
   iq(2) = symm%symop_id(1,thdr%idtet(2,itet))
   iq(3) = symm%symop_id(1,thdr%idtet(3,itet))
   iq(4) = symm%symop_id(1,thdr%idtet(4,itet))

   !wthdr = thdr%vltet(itet)
   wthdr = real(thdr%idtet(0,itet))*thdr%vltet(itet) !the volume of the individual tetrahedra is already in units of the reciprocal unit cell

   do nb=1,ek%nband_max !loop over bands

      if (nb < ek%nbopt_min) cycle !if there are no optical matrix elements there is nothing to interpolate
      if (nb > ek%nbopt_max) cycle
      ! get the band energy at each corner of the tetrahedron:
      ec(1) = ek%band(iq(1),nb)
      ec(2) = ek%band(iq(2),nb)
      ec(3) = ek%band(iq(3),nb)
      ec(4) = ek%band(iq(4),nb)

      ! sort the energies at the four corners (array ec) into array es
      do i=1,4
         ec1(i)=ec(i)
      enddo
      do 3 i=1,4
         i00=1
         do j=2,4
            if (ec1(j)<ec1(i00)) i00=j
         enddo
         es(i) = ec1(i00)
         iqs(i)= iq(i00)
         ec1(i00)=1.d30
      3   continue
!!!!!!!!!!!TEST
      !write(*,*) 'tetrahedron no.',itet
      !do i=1,4
      !   write(*,'(F8.3, 2I6)') ec(i),iq(i)
      !enddo
      !do i=1,4
      !   write(*,'(F8.3, 2I6)') es(i),iqs(i)
      !enddo
      !STOP
!!!!!!!!!!!TEST END

      ! define the constants required for later
      x1 = es(1)-ef
      x2 = es(2)-ef
      x3 = es(3)-ef
      x4 = es(4)-ef
      c0 = wthdr/4.d0

      if (es(1) > ef) then
         iswitch=1

      else if ((es(1) < ef) .and. (es(2) > ef)) then
         iswitch=2
         c1 = -c0*x1*x1*x1/((es(2)-es(1))*(es(3)-es(1))*(es(4)-es(1)))
         cc1= (1.d0/(es(2)-es(1))) + (1.d0/(es(3)-es(1))) + (1.d0/(es(4)-es(1)))

      else if ((es(2) < ef) .and. (es(3) > ef)) then
         iswitch=3
         c1 = c0*x1*x1/((es(4)-es(1))*(es(3)-es(1)))
         c2 = c0*x1*x2*x3/((es(4)-es(1))*(es(3)-es(2))*(es(3)-es(1)))
         c3 = c0*x2*x2*x4/((es(4)-es(2))*(es(3)-es(2))*(es(4)-es(1)))

      else if ((es(3) < ef) .and. (es(4) > ef)) then
         iswitch=4
         c1 = c0*x4*x4*x4/((es(4)-es(1))*(es(4)-es(2))*(es(4)-es(3)))
         cc4= (1.d0/(es(4)-es(1))) + (1.d0/(es(4)-es(2))) + (1.d0/(es(4)-es(3)))

      else if (es(4) < ef) then
         iswitch=5
      else
         write(*,*)'INTERPTRA_RE: the ordering of your thetrahedron vertices is not consistent', itet
         STOP
      endif
      !09.04.2018: the switch case above cuts out the
      !the contribution from unoccupied states, this is
      !erroneous since also states above the fermi level
      !contribute. The easiest fix is to set:
      iswitch=5

      select case (iswitch)
         case(1)
            w(1) = 0.d0
            w(2) = 0.d0
            w(3) = 0.d0
            w(4) = 0.d0

         case(2)
            w(1) = c1*(4+(cc1*x1)) !in eq. B2 there is a factor 4-(ef-e1), hence the opposite sign
            w(2) = c1*x1/(es(1)-es(2))
            w(3) = c1*x1/(es(1)-es(3))
            w(4) = c1*x1/(es(1)-es(4))

         case(3)
            w(1) = c1 + ((c1+c2)*x3/(es(3)-es(1))) + ((c1+c2+c3)*x4/(es(4)-es(1)))
            w(2) = c1+c2+c3 + ((c2+c3)*x3/(es(3)-es(2))) + (c3*x4/(es(4)-es(2)))
            w(3) = ((c1+c2)*x1/(es(1)-es(3))) + ((c2+c3)*x2/(es(2)-es(3)))
            w(4) = ((c1+c2+c3)*x1/(es(1)-es(4))) + (c3*x2/(es(2)-es(4)))

         case(4)
            w(1) = c0 - (c1*x4/(es(4)-es(1)))
            w(2) = c0 - (c1*x4/(es(4)-es(2)))
            w(3) = c0 - (c1*x4/(es(4)-es(3)))
            w(4) = c0 - (c1*(4-(cc4*x4)))

         case(5)
            w(1) = c0
            w(2) = c0
            w(3) = c0
            w(4) = c0

      end select
         ! now that the weights are set need to perform the integration within the tetrahedron for each band
         ! and trace them over

      do ix=1,lat%nalpha
         do iy=ix,lat%nalpha
            do iktet=1,4 !linear interpolation within the  tetrahedron
               ! set the value of the scattering rate for the specific temperature, band
               ik = thdr%idtet(iktet, itet) ! from thdr identifier -> k
               ikk = symm%symop_id(1,ik)    ! from k -> to the data we have saved
               if (allocated(sct%ykb)) then
                  resp%gamma=real(ek%z(ikk,nb)*(sct%gam(iT)+sct%ykb(iT,ikk,nb)),8)
                  if (present(hpresp)) hpresp%gamma=real(ek%z(ikk,nb)*(sct%gam(iT)+ &
                                                    sct%ykb(iT,ikk,nb)),16)
               else
                  resp%gamma=real(ek%z(ikk,nb)*sct%gam(iT),8)
                  if (present(hpresp)) hpresp%gamma=real(ek%z(ikk,nb)*sct%gam(iT),16)

               endif
               if (lBoltz) then
                  ! no interband transitions for Boltzmann response
                  resp%s(nb,ix,iy) = resp%s(nb,ix,iy) + (real(w(iktet),8)*resp%s_tmp(iktet,nb,ix,iy))/(resp%gamma)
                  resp%a(nb,ix,iy) = resp%a(nb,ix,iy) + (real(w(iktet),8)*resp%a_tmp(iktet,nb,ix,iy))/(resp%gamma)
                  if (allocated(resp%sB)) then
                     resp%sB(nb,ix,iy) = resp%sB(nb,ix,iy) + &
                                       (real(w(iktet),8)*resp%sB_tmp(iktet,nb,ix,iy))/(resp%gamma**2)
                     resp%aB(nb,ix,iy) = resp%aB(nb,ix,iy) + &
                                       (real(w(iktet),8)*resp%aB_tmp(iktet,nb,ix,iy))/(resp%gamma**2)
                  endif
               else
                  select type(resp)
                     type is(dp_resp)
                     ! intraband response
                     resp%s(nb,ix,iy) = resp%s(nb,ix,iy) + (real(w(iktet),8)*resp%s_tmp(iktet,nb,ix,iy))*beta/(resp%gamma )
                     resp%a(nb,ix,iy) = resp%a(nb,ix,iy) + (real(w(iktet),8)*resp%a_tmp(iktet,nb,ix,iy))*(beta**2)/(resp%gamma )
                     type is(dp_respinter)
                     ! interband response
                     resp%s(nb,ix,iy) = resp%s(nb,ix,iy) + (real(w(iktet),8)*resp%s_tmp(iktet,nb,ix,iy)) !all the prefactors have been included
                     resp%a(nb,ix,iy) = resp%a(nb,ix,iy) + (real(w(iktet),8)*resp%a_tmp(iktet,nb,ix,iy)) !in RESPINTERT or in RESDERTET
                  end select

                  if (allocated(resp%sB)) then
                     select type(resp)
                        type is(dp_resp)
                        resp%sB(nb,ix,iy) = resp%sB(nb,ix,iy) + &
                                          (real(w(iktet),8)*resp%sB_tmp(iktet,nb,ix,iy))*beta/(resp%gamma**2)
                        resp%aB(nb,ix,iy) = resp%aB(nb,ix,iy) + &
                                          (real(w(iktet),8)*resp%aB_tmp(iktet,nb,ix,iy))*beta/(resp%gamma**2)
                        type is(dp_respinter)
                        resp%sB(nb,ix,iy) = 0.0d0
                        resp%aB(nb,ix,iy) = 0.0d0
                     end select
                  endif
               endif
            enddo !over corners of the tetrahedron

            ! trace over the bands
            resp%s_tet(ix,iy) = resp%s_tet(ix,iy) + resp%s(nb,ix,iy)
            resp%a_tet(ix,iy) = resp%a_tet(ix,iy) + resp%a(nb,ix,iy)
            if (allocated(resp%sB)) then
              resp%sB_tet(ix,iy) = resp%sB_tet(ix,iy) + resp%sB(nb,ix,iy)
              resp%aB_tet(ix,iy) = resp%aB_tet(ix,iy) + resp%aB(nb,ix,iy)
            endif

            if (present(hpresp)) then
               do iktet=1,4 !linear interpolation within the  tetrahedron
                  hpresp%s(nb,ix,iy) = hpresp%s(nb,ix,iy) + &
                                     (real(w(iktet),16)*hpresp%s_tmp(iktet,nb,ix,iy))*betaQ/(hpresp%gamma)
                  hpresp%a(nb,ix,iy) = hpresp%a(nb,ix,iy) + &
                                     (real(w(iktet),16)*hpresp%a_tmp(iktet,nb,ix,iy))*(betaQ**2)/(hpresp%gamma)
                  if (allocated(hpresp%sB)) then
                     hpresp%sB(nb,ix,iy) = hpresp%sB(nb,ix,iy) + &
                                         (real(w(iktet),16)*hpresp%sB_tmp(iktet,nb,ix,iy))*betaQ/(hpresp%gamma**2)
                     hpresp%aB(nb,ix,iy) = hpresp%aB(nb,ix,iy) + &
                                         (real(w(iktet),16)*hpresp%aB_tmp(iktet,nb,ix,iy))*(betaQ/hpresp%gamma)**2
                  endif
               enddo

               ! trace over the bands
               hpresp%s_tet(ix,iy) = hpresp%s_tet(ix,iy) + hpresp%s(nb,ix,iy)
               hpresp%a_tet(ix,iy) = hpresp%a_tet(ix,iy) + hpresp%a(nb,ix,iy)
               if (allocated(hpresp%sB)) then
                 hpresp%sB_tet(ix,iy) = hpresp%sB_tet(ix,iy) + hpresp%sB(nb,ix,iy)
                 hpresp%aB_tet(ix,iy) = hpresp%aB_tet(ix,iy) + hpresp%aB(nb,ix,iy)
               endif

            endif !hpresp

         enddo !iy
      enddo !ix

   enddo !over bands

  !!!!!!!!!!!!!!!!!!!TEST
  ! if (mod(itet,1000) ==0) then
  !   write(*,*) 'resp%s  allocated?',allocated(resp%s)
  !   write(*,*) 'hpresp type present?',present(hpresp)
  !   write(itet+1,*) resp%s_tmp(1,:,1,1)
  !   write(itet+2,*) resp%s_tmp(2,:,1,1)
  !   write(itet+3,*) resp%s_tmp(3,:,1,1)
  !   write(itet+4,*) resp%s_tmp(4,:,1,1)
  ! endif
  !!!!!!!!!!!!!!!!!!!TEST END

  !!!!!!!!!!!!!!!!!!!TEST 11.06.2018
  ! if (itet ==15) then
  !   select type(resp)
  !      type is(dp_resp)
  !      do nb = 1, ek%nband_max
  !         write(itet,*) resp%s(nb,1,2)
  !      enddo
  !   end select
  ! endif
  !!!!!!!!!!!!!!!!!!!TEST 11.06.2018 END

  !!!!!!!!!!!!!!!!!!!SANITY CHECK
  ! for the TB model at small gamma values there are many!
  !select type(resp)
  !   type is(dp_resp)
  !   if (resp%s_tet(1,1) < 0.0d0) then
  !      write(*,*) 'WARNING: negative response in INTERPTRA_MU', itet, resp%s_tet(1,1)
  !   endif
  !end select


 return
end subroutine !INTERPTRA_RE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! INTERPTRA_MU
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
subroutine interptra_mu (vltet, occ_tet, occ_intp)
  implicit none

   double precision, intent (in) :: vltet      !tetrahedron volume
   double precision, intent (in) :: occ_tet(4) !occupation numbers at tetrahedra vertices
   double precision, intent (out):: occ_intp   !occupation numbers at tetrahedra vertices
!local variables
   integer :: i
   double precision :: c0
   double precision :: wthdr  !weight of the tetrahedron
   double precision :: w(4)   !weights of the interpolation formula

   !wthdr = real(thdr%idtet(0,itet))*thdr%vltet(itet) !the volume of the individual tetrahedra is already in units of the reciprocal unit cell
   wthdr = vltet !the volume of the individual tetrahedra is already in units of the reciprocal unit cell

   c0 = wthdr/4.d0

   w(1) = c0
   w(2) = c0
   w(3) = c0
   w(4) = c0

   occ_intp = 0.d0
   do i=1,4 !linear interpolation within the  tetrahedron
      occ_intp = occ_intp + (w(i)*occ_tet(i))
   enddo !over corners of the tetrahedron

 return
end subroutine !INTERPTRA_MU

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! INTERPTRA_MUQ
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
subroutine interptra_muQ (vltet, target_tet, target_intp)
  implicit none

   double precision, intent (in) :: vltet      !tetrahedron volume
   real(16), intent (in) :: target_tet(4) !occupation numbers at tetrahedra vertices
   real(16), intent (out):: target_intp   !occupation numbers at tetrahedra vertices
!local variables
   integer :: i
   real(16) :: c0
   real(16) :: wthdr  !weight of the tetrahedron
   real(16) :: w(4)   !weights of the interpolation formula

   wthdr = real(vltet,16) !the volume of the individual tetrahedra is already in units of the reciprocal unit cell

   c0 = wthdr/4.q0

   w(1) = c0
   w(2) = c0
   w(3) = c0
   w(4) = c0

   target_intp = 0.q0
   do i=1,4 !linear interpolation within the  tetrahedron
      target_intp = target_intp + (w(i)*target_tet(i))
   enddo !over corners of the tetrahedron

 return
end subroutine !INTERPTRA_MUQ

  ! tight binding functions for the nearest neighbor case
  ! ek_sc -> e(k)
  ! vk_sc -> e'(k)
  ! vkk_sc -> e''(k)
  double precision function ek_sc(k,iband,edisp)
    implicit none
    type(energydisp)      :: edisp
    double precision :: k(3),ek,bw
    integer          :: iband,i

    ek=edisp%E0(iband)
    bw=edisp%t(iband,1)
    do i=1,3
       ek=ek + 2.d0*bw*cos(2.d0*pi*k(i))
    enddo
    ek_sc=ek
    return
  end function ek_sc

  double precision function vk_sc(idir,k,iband,edisp)
    implicit none
    type(energydisp) :: edisp
    double precision k(3),bw
    integer iband,i,idir

    bw=edisp%t(iband,1)
    vk_sc=-2.d0*bw*sin(2.d0*pi*k(idir))*lat%a(idir)
    return
  end function vk_sc

  double precision function vkk_sc(idir,idir2,k,iband,edisp)
    implicit none
    type(energydisp) :: edisp
    double precision k(3),bw
    integer iband,idir,idir2

    bw =edisp%t(iband,1)
    if (idir.eq.idir2) then
       vkk_sc=-2.d0*bw*cos(2.d0*pi*k(idir))*(lat%a(idir))**2
    else
       vkk_sc=0.d0
    endif
    return
  end function vkk_sc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! This function evaluates the length of the tetrahedron edge
!! required by subroutine GENTETRA
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! elemental function -> you can put in here also 3 equally sized arrays
  ! and get out the result in form of an identically sized array
  ! the function is then applied element wise
  pure elemental double precision function anrm2(x,y,z)
    implicit none
    double precision, intent(in) :: x, y, z
    anrm2=x*x*1.00001d0+y*y*1.00002d0+z*z*1.00003d0 &
      &             -x*0.000004d0-y*0.000003d0-z*0.000002d0
  end function
