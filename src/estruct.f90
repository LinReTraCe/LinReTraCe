module Mestruct
  use Mparams
  use Mtypes
  use Mmymath
  implicit none

  real(8) :: gmax,gminall,vol ! gmax is maximal gamma of bands that determine [gap]
  contains

   subroutine estruct_init(algo, kmesh, redkm, fulkm, eirrk, eredk, efulk, thdr, dos, sct)
     implicit none

     type(algorithm)  :: algo
     type(kpointmesh) :: kmesh  ! contains k-point mesh specifiers and logical switches on how to get the mesh from
     type(kpointmesh) :: redkm  ! contains k-point mesh specifiers and logical switches on how to get the mesh from
     type(kpointmesh) :: fulkm  ! contains k-point mesh specifiers and logical switches on how to get the mesh from
     type(edisp)      :: eirrk  ! contains the band dispersion energy and the optical matrix elements (when which > 2) along the irr-k-mesh
     type(edisp)      :: eredk  ! contains the band dispersion energy and the optical matrix elements (when which > 2) along the red-k-mesh
     type(edisp)      :: efulk  ! contains the band dispersion energy and the optical matrix elements for the red-k-mesh including BZ endpoints
     type(symop)      :: symm   ! contains the symmetry operations
     type(tetramesh)  :: thdr   ! contains the tetrahedra (should you need them...)
     type(dosgrid)    :: dos
     type(scatrate)   :: sct
     integer          :: ierr, iband, i
     double precision :: maxtmp

!eM: ikstart and ikend are parameters set by a previous mpi call in main.f90
! I think that the parallelisation will have to take place after the construction of the
! tetrahedra

     !the total number of k-points has to be increased because we need also the endpoints in the BZ
     if (algo%ltetra) then
        kmesh%ktot=(kmesh%kx +1)*(kmesh%ky +1)*(kmesh%kz +1)
     else
        kmesh%ktot=kmesh%kx*kmesh%ky*kmesh%kz
     endif

     !these optical quantities are properly set only for algo%loptic==.true.
     !*************************************
     ! get electronic structure
     if (algo%ltbind) call gentbstr(algo, eirrk, kmesh)
     if (algo%lw2k .or. algo%loptic) call genelstr (algo, kmesh, redkm, eirrk, eredk, symm)
     !*************************************
     if (algo%ltbind) then
        eirrk%nbopt_max = eirrk%nband_max
        eirrk%nbopt_min = 1
     endif


     if (algo%ltetra) then
 ! evaluate the density of states (DOS) and the integrated DOS (NOS)
 ! the latter can be used to evaluate the Fermi level/chemical potential

        if (algo%ltbind) then
            ! generate the energy grid for the DOS
            maxtmp=0.d0
            do iband=1,eirrk%nband_max
               if (maxtmp < abs(eirrk%E0(iband))) maxtmp = abs(eirrk%E0(iband))
            enddo
            dos%emax= 2.0d0*maxtmp
            ! this is quite critical cutoff value, if the energy window is too small
            ! one loses DOS from tetrahedra with vertices whose energy is outside the window
            dos%emin=-dos%emax
            dos%nnrg=2.0d0*dos%emax/0.007d0
            !energy spacing of roughly 0.1 eV
            call gentetra(algo, kmesh, thdr)             ! generates the tetrahedra
            call intetra (kmesh, thdr, dos, eirrk%band, eirrk%nband_max)   ! computes the dos and the integrated dos
            !call intetra (kmesh, eirrk, thdr, dos)       ! computes the dos and the integrated dos
            !spin multiplicity
            do i=1,size(dos%enrg)
               dos%dos(i)=2.0d0*dos%dos(i)
               dos%nos(i)=2.0d0*dos%nos(i)
            enddo
            open(10,file='tetrados',status='unknown')
            do i=1,size(dos%enrg)
               write (10,*) dos%enrg(i), dos%dos(i), dos%nos(i)
            enddo
            close(10)
            if (algo%imurestart == 0) call findef(dos, eirrk)   ! finds the Fermi level
        else
        !using electronic structure from Wien2k
            maxtmp=0.d0
            do iband=1,eredk%nband_max
               if ((eredk%band(1,iband) < 1.0d2) .and. (maxtmp < abs(eredk%band(1,iband)))) &
                  maxtmp = abs(eredk%band(1,iband))
            enddo
            dos%emax= 20.0d0*maxtmp
            dos%emin=-dos%emax
            dos%nnrg=2.0d0*dos%emax/0.007d0
            !energy spacing of roughly 0.1 Ry
            write(*,*) 'estruct: Wien2k electronic structure tetrahedron method calculation'

            ! a new k-point mesh has to be generated from redkm to include also
            ! the terminal k-points on the BZ, let's call it fulkm
            call genfulkm(redkm, fulkm, eredk, efulk, symm)
            call gentetra(algo, fulkm, thdr)         ! generates the tetrahedra
            call intetra (fulkm, thdr, dos, efulk%band, efulk%nband_max)   ! computes the dos and the integrated dos
            !call intetra (fulkm, efulk, thdr, dos)   ! computes the dos and the integrated dos
            !spin multiplicity only for Wien2k calculations
            do i=1,size(dos%enrg)
               dos%dos(i)=2.0d0*dos%dos(i)
               dos%nos(i)=2.0d0*dos%nos(i)
            enddo
            open(10,file='tetrados',status='unknown')
            do i=1,size(dos%enrg)
               write (10,*) dos%enrg(i), dos%dos(i), dos%nos(i)
            enddo
            close(10)
            if (algo%imurestart == 0) call findef(dos, efulk)   ! finds the Fermi level
        endif
     else  !ltetra = F
     !generate the density of states for the linear grid
        if (algo%ltbind) then
           !generate energy grid exactly as before
           maxtmp=0.d0
           do iband=1,eirrk%nband_max
              if (maxtmp < abs(eirrk%E0(iband))) maxtmp = abs(eirrk%E0(iband))
           enddo
           dos%emax= 1.2d0*maxtmp
           dos%emin=-dos%emax
           dos%nnrg=5001
           !energy spacing of roughly 0.1 eV

           call gendosel(kmesh, eirrk, dos)
           !spin multiplicity
           dos%dos=2.0d0*dos%dos
           dos%nos=2.0d0*dos%nos

           open(10,file='linrados',status='unknown')
           do i=1,size(dos%enrg)
              write (10,*) dos%enrg(i), dos%dos(i), dos%nos(i)
           enddo
           close(10)
           if (algo%imurestart == 0) call findef(dos, eirrk)   ! finds the Fermi level

        else !Wien2k calculation
        !brute force summation for the DOS, with each
        !Dirac's delta replaced by a Lorentzian bandshape
           !generate energy grid exactly as before
           maxtmp=0.d0
           do iband=1,eredk%nband_max
              if ((eredk%band(1,iband) < band_fill_value) .and. (maxtmp < abs(eredk%band(1,iband)))) &
                 maxtmp = abs(eredk%band(1,iband))
           enddo
           dos%emax= 1.2d0*maxtmp
           dos%emin=-dos%emax
           dos%nnrg=5001
           !energy spacing of roughly 0.1 Ry
           write(*,*) 'estruct: Wien2k electronic structure linear k-mesh DOS calculation'

           call gendosel (redkm, eredk, dos)

           !spin multiplicity
           dos%dos=2.0d0*dos%dos
           dos%nos=2.0d0*dos%nos

           open(10,file='linrados',status='unknown')
           do i=1,size(dos%enrg)
              write (10,*) dos%enrg(i), dos%dos(i), dos%nos(i)
           enddo
           close(10)
           if (algo%imurestart == 0) call findef(dos, eredk)   ! finds the Fermi level
        endif

     endif !ltetra
   end subroutine !estruct_init

   subroutine read_config(algo, ek, kmesh, sct)
     implicit none
    type (algorithm) :: algo
    type(edisp)      :: ek     ! contains the band dispersion energy and the optical matrix elements (when which > 2) along the irr-k-mesh
    type(kpointmesh) :: kmesh  ! contains k-point mesh specifiers and logical switches on how to get the mesh from
    type(scatrate)   :: sct
  !local variables
    !double precision :: Z, ReS ! renormalisation factor, Re{Selfenergy} (a.k.a. quasi-particle shift)
    integer :: which           ! switches the type of algorithm (then stored in type(algorithm))
    integer :: which_tetra     ! =0 no tetrahedron integration =1 use tetrahedron
    integer :: which_grid      ! =0 no symetry operations used  =1 use symmetry to generate BZ
    integer :: tmax            ! if not equal to 1 it introduces anisotropic dispersion in the TB model
    integer :: iband
    character(len=100) :: ctmp
    integer :: itmp

    !read in lattice parameters...
    open(10,file='inp.only',status='old')
    read(10,*)ek%nband_max, tmax
    read(10,*)kmesh%alat,kmesh%a(:) ! ALAT in A, a(1:3) in ALAT
    read(10,*)ek%nelect             ! read the number of electrons, the fermi level will be computed, consistently with the LINRETRACE routines
    read(10,*)vol

! at this point one should choose in the input file what to do:
! case = 1 use the uniform k-mesh generated with tight binding
! case = 2 read the k-points and band dispersion from Wien2k
! case = 3 read also the optical matrix elements from Wien2k
! case = 4 read the renormalised qp energies, renormalisation factor and -Im{selfnrg}
    read(10,*)which, which_tetra, which_grid

    if (which==1) then
       algo%ltbind =.true.
       algo%lw2k   =.false.
       algo%loptic =.false.
       algo%ldmft  =.false.
    else if (which==2) then
       algo%ltbind =.false.
       algo%lw2k   =.true.
       algo%loptic =.false.
       algo%lBfield=.false.
       algo%ldmft  =.false.
    else if (which==3) then
       algo%ltbind =.false.
       algo%lw2k   =.false.
       algo%loptic =.true.
       algo%lBfield=.false.
       algo%ldmft  =.false.
    else if (which==4) then
       algo%ltbind =.false.
       algo%lw2k   =.false.
       algo%loptic =.true.
       algo%lBfield=.false.
       algo%ldmft  =.true.
    else
       write(*,*)'Not valid algorithm selection'
       STOP
    endif

    if (which_tetra==0) then
       algo%ltetra=.false.
       write(*,*)'Tetrahedron method not selected'
    else if (which_tetra==1) then
       algo%ltetra=.true.
       write(*,*)'Tetrahedron method selected'
    else
       write(*,*)'Integration method does not exist',which_tetra
       STOP
    endif

    if (which_grid==0) then
       algo%lsymm=.false.
       write(*,*)'Symmetry switched OFF, REDucible BZ required'
    elseif(which_grid==1) then
       algo%lsymm=.true.
       write(*,*)'Symmetry switched ON, IRReducible BZ required'
    else
       write(*,*)'You must either provide reducible BZ or switch on symmetry',which_grid
       STOP
    endif
    if (algo%ltbind) algo%lsymm=.false.

    ! this is unfortunately necessary if we want arbitrary paths
    ! this whole section will be replaced anyways ...
    read(10,'(A)') ctmp
    itmp = index(ctmp,'#')
    if (itmp .gt. 0) then
      algo%mysyst = trim(adjustl(ctmp(:itmp-1)))
    else
      algo%mysyst = ctmp
    endif

    read(10,*)kmesh%kx,kmesh%ky,kmesh%kz

    if (algo%ltbind) then
       write(*,*)'readinfile: reading TB parameters'
       allocate(ek%E0(ek%nband_max),ek%t(ek%nband_max, tmax))
       do iband=1,ek%nband_max
          read(10,*)ek%E0(iband),ek%t(iband,:)
       enddo
       read(10,*)algo%lBfield
    endif
    !now read the temperature variables
    read(10,*)sct%Tmin,sct%Tmax,sct%dT
    read(10,*)algo%imurestart, ek%efer
   !at this point I read in also the coefficients of gamma
    read(10,*)sct%ng
    allocate(sct%gc(0:sct%ng))
    read(10,*)sct%gc(:)
   !read in also the renormalisation factor and the quasi-particle shift
    read(10,*)ek%ztmp

    close(10)

   end subroutine

   subroutine gentbstr(algo, eirrk, kmesh)
    implicit none

    type(algorithm) :: algo
    type(edisp) :: eirrk                 ! contains the band dispersion energy and the optical matrix elements (when which > 2) along the irr-k-mesh
    type(kpointmesh) :: kmesh            ! contains k-point mesh specifiers and logical switches on how to get the mesh from
  !local variables
    integer :: ik,ikx,iky,ikz, nk,nkx,nky,nkz, iband,nband, idir,idir2 !counters
    double precision :: k(3)
    double precision, allocatable :: xtmp(:), ytmp(:), dytmp(:)

   if (algo%ltetra) then
      nkx=kmesh%kx+1; nky=kmesh%ky+1; nkz=kmesh%kz+1
   else
      nkx=kmesh%kx; nky=kmesh%ky; nkz=kmesh%kz
   endif
   nband=eirrk%nband_max

   allocate(kmesh%k_id(nkx,nky,nkz))
   allocate(kmesh%k_coord(3,kmesh%ktot))
   allocate(eirrk%band( kmesh%ktot, nband ))
   allocate(eirrk%Mopt(3, kmesh%ktot, nband, nband))
   if (algo%lBfield) allocate(eirrk%M2(3, 3, kmesh%ktot, nband))

   if (allocated(kmesh%k_coord)) kmesh%k_coord=0.d0
   if (allocated(kmesh%k_id)) kmesh%k_id=0
   if (allocated(eirrk%band)) eirrk%band=0.d0
   if (allocated(eirrk%Mopt)) eirrk%Mopt=0.d0
   if (allocated(eirrk%M2))   eirrk%M2  =0.d0

   if(algo%ldebug) then !open files to write
      open(10,file='kcoords')
      open(11,file='ek')
      open(12,file='vk')
      if(algo%lBfield) open(13,file='vkk')
   endif
   !! eM: for a tight binding calculations one needs to generate a
   !! separate k-mesh because the search for the chemical potential
   !! is carried out on a linear mesh anyway, hence one does not
   !! include the BZ endpoint in this linear mesh.
   !! See definition of nkx,nky,nkz above.
   ik=0
   do ikx=1,nkx
      k(1)=1.d0/real(kmesh%kx,kind=8)*real(ikx-1,kind=8)
      do iky=1,nky
         k(2)=1.d0/real(kmesh%ky,kind=8)*real(iky-1,kind=8)
         do ikz=1,nkz
            k(3)=1.d0/real(kmesh%kz,kind=8)*real(ikz-1,kind=8)
            ik=ik+1
            kmesh%k_id(ikx,iky,ikz)=ik
            kmesh%k_coord(1,ik)=k(1)
            kmesh%k_coord(2,ik)=k(2)
            kmesh%k_coord(3,ik)=k(3)
            do iband =1,nband
               eirrk%band(ik,iband)=ek_sc(k,iband,eirrk)
               do idir=1,3
                  eirrk%Mopt(idir,ik,iband,iband)=vk_sc(idir,k,iband,eirrk,kmesh)
               enddo
            enddo
            if (algo%lBfield) then
               do iband =1,nband
                  do idir=1,3
                     eirrk%M2(idir,idir,ik,iband)=vkk_sc(idir,idir,k,iband,eirrk,kmesh)
                  enddo
               enddo
            endif

            if(algo%ldebug) then !write to file
               write(10,'(1I14,3F12.6)')ik,k
               write(11,'(1I14,100F16.10)')ik,(eirrk%band(ik,iband),iband=1,nband)
               write(12,'(1I14,500F16.10)')ik,((eirrk%Mopt(idir,ik,iband,iband),idir=1,3),iband=1,nband)
               if (algo%lBfield) &
      write(13,'(1I14,1000F16.10)')ik,(((eirrk%M2(idir,idir2,ik,iband),idir=1,3),idir2=1,3),iband=1,nband)
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
   !!!!!!!!!!!!!!!TEST
   ! costruct the numerical derivative for the band dispersion
   ! and compare it with the analitycal estimate
   ! select the direction in which you want to estimate the derivative
   !allocate(xtmp(1:nkx))
   !allocate(ytmp(1:nkx))
   !allocate(dytmp(1:nkx))
   !iband=1
   !do ikx=1,nkx
   !   ik = kmesh%k_id(ikx,2,2)
   !   xtmp(ikx) = kmesh%k_coord(1,ik)
   !   ytmp(ikx) = eirrk%band(ik,iband)
   !   dytmp(ikx)= 0.0d0
   !enddo
   !
   !call derrich(xtmp,ytmp,dytmp)
   !do ikx=1,nkx
   !   write(70,*)xtmp(ikx),ytmp(ikx),dytmp(ikx)
   !   write(71,*)xtmp(ikx),dytmp(ikx),2.0d0*pi*eirrk%Mopt(1,kmesh%k_id(ikx,2,2),iband,iband)/kmesh%alat
   !enddo
   !deallocate(xtmp,ytmp,dytmp)
   !STOP
   !!!!!!!!!!!!!!!TEST END

   ! if we are passing only 1 k-point chances are that we want to study a
   ! flat-band model; in this case the Fermi velocity gets overwritten
   ! TODO: testing
   if (kmesh%kx==1 .and. kmesh%ky==1 .and. kmesh%kz==1) then
   !if (nkx==1 .and. nky==1 .and. nkz==1) then !equivalent to the line above because with 1 k-point you cannot select algo%tetra
      if (kmesh%k_coord(1,1)==0.0d0 .and. kmesh%k_coord(2,1)==0.0d0 .and. kmesh%k_coord(3,1)==0.0d0) then
         write(*,*)  'Gamma-point only calculation; assuming flatband model'
         do iband =1,nband
            do idir=1,3
               eirrk%Mopt(idir,1,iband,iband)=1.0d0
            enddo
         enddo
      endif
   endif

   end subroutine !gentbstr

   subroutine genelstr (algo, kmesh, redkm, eirrk, eredk, symm)
     implicit none
    type(algorithm) :: algo
    type(kpointmesh) :: kmesh
    type(kpointmesh) :: redkm
    type(edisp) :: eirrk
    type(edisp) :: eredk
    type(symop) :: symm
    integer :: ik,ikx,iky,ikz, nb, nb1

   if (algo%lsymm) then
      call getirrk  (algo%mysyst, symm%cntr, kmesh, eirrk)
      call getsymop (algo%mysyst, symm,  kmesh, eirrk)
      ! generate a reucible kmesh from the set of symmetry operations and the irrek-mesh
      call genredk  (algo%mysyst, symm, kmesh, redkm)
      ! copy also the number of electrons and the Fermi level from the old datastructure to the new one
      eredk%nband_max = eirrk%nband_max
      eredk%nelect = eirrk%nelect
      eredk%efer = eirrk%efer
      redkm%alat = kmesh%alat
      redkm%a    = kmesh%a
      ! read in the optical matrix elements on the irreducible grid
      ! and if algo%ldmft is true also the selfenergy file
      call getirropt (algo%mysyst, kmesh, eirrk, algo%loptic, symm%lcubic, algo%ldmft)
      ! generate the optical matrix elements evaluated on the new redkm grid
      call genredopt (redkm, symm, eirrk, eredk)
      ! translate the reducible k-mesh and take care of the bandstructure
      call trnredk (kmesh, redkm, eredk, symm, algo)

   else !reducible BZ is already provided by W2k
      call getirrk (algo%mysyst, symm%cntr, redkm, eredk)
      !the values below are stored in the irreducible mesh by readinfile
      eredk%nband_max = eirrk%nband_max
      eredk%nelect = eirrk%nelect
      eredk%efer = eirrk%efer
      redkm%alat = kmesh%alat
      redkm%a    = kmesh%a
      redkm%kx = kmesh%kx; redkm%ky = kmesh%ky; redkm%kz = kmesh%kz
      call getsymop (algo%mysyst, symm, redkm, eredk) !need to call this routine to set the logical switches in symmop type

      if (.not. allocated(redkm%k_id)) allocate(redkm%k_id(1:redkm%kx, 1:redkm%ky, 1:redkm%kz)) !if algo%lsymm=.false. redkm%k_id has not been allocated
      ik=0
      do ikx=1,redkm%kx
         do iky=1,redkm%ky
            do ikz=1,redkm%kz
               ik=ik+1
               redkm%k_id(ikx, iky, ikz)=ik
            enddo
         enddo
      enddo
      !!!!!!!!!!!!!!!!TEST
      !do ik=1,redkm%ktot
      !   write(667,*)'KP',ik,redkm%k_coord(1,ik),redkm%k_coord(2,ik),redkm%k_coord(3,ik)
      !   do nb=1,eredk%nband_max
      !      write(667,160)nb,eredk%band(ik,nb)
      !   enddo
      !enddo
      !160  FORMAT  (I4,X,F12.8,X)
      !!!!!!!!!!!!!!!!TEST END

      ! read in the optical matrix elements on the reducible grid
      ! and if algo%ldmft is true also the selfenergy file
      call getirropt ( algo%mysyst, redkm, eredk, algo%loptic, symm%lcubic, algo%ldmft)
   endif !lsymm

   if (.not. algo%loptic) then
      eredk%Mopt(:,:,:,:)=0.0d0
      write(*,*)  'GENELSTR: setting Mopt = identity matrix'
      write(*,*) 'size Mopt =',size(eredk%Mopt,1)

      do ik=1,redkm%ktot
         do nb=eredk%nbopt_min,eredk%nbopt_max
            eredk%Mopt(1,ik,nb,nb)=1.0d0
            eredk%Mopt(2,ik,nb,nb)=1.0d0
            eredk%Mopt(3,ik,nb,nb)=1.0d0
         enddo
      enddo
      !!!!!!!!!!!!!!!!TEST
      !do ik=1,redkm%ktot
      !   do nb=eredk%nbopt_min,eredk%nbopt_max
      !      do nb1=nb+1,eredk%nbopt_max
      !         eredk%Mopt(1,ik,nb,nb1)=0.0d0
      !         eredk%Mopt(2,ik,nb,nb1)=0.0d0
      !         eredk%Mopt(3,ik,nb,nb1)=0.0d0
      !      enddo
      !   enddo
      !   do nb=eredk%nbopt_min,eredk%nbopt_max
      !      do nb1=nb+1,eredk%nbopt_max
      !         eredk%Mopt(1,ik,nb1,nb)=0.0d0
      !         eredk%Mopt(2,ik,nb1,nb)=0.0d0
      !         eredk%Mopt(3,ik,nb1,nb)=0.0d0
      !      enddo
      !   enddo
      !enddo
      !!!!!!!!!!!!!!!!TEST END
   endif
   end subroutine !genelstr


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! GETIRRK
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine extracts from the Wien2k files
! information about the irreducible k-mesh, this is stored
! either in the kpointmesh or in edisp type
!

subroutine getirrk (mysyst, cntr, kmesh, edspk )
     implicit none

!passed variables
  character(len=*) :: mysyst ! file name to open
  character(3)  :: cntr     ! centering of the unit cell: P,F,B,H,C??
  type(kpointmesh) :: kmesh ! k-mesh generated by Wien2k
  type(edisp)      :: edspk ! energy dispersion and optical matrix elements over k-mesh generated by Wien2k
!internal variables
  integer :: icrap, itmp, nband, nband_loc, i, j, ik
  integer :: iatm !number of non-equivalent atoms in cell
  real :: rcrap
  character(100) :: ccrap ! this works for arbitrary size
  double precision :: dtmp
  double precision, allocatable :: band_tmp(:,:)    ! bigger datastructure that will contain all the energy dispersion curves

   !get the number of k-points
   open(10,file=trim(adjustl(mysyst))//'.weight',status='old')
   read(10,*) rcrap, rcrap, ccrap
   read(10,*) kmesh%ktot, ccrap
   write(*,*) 'total number of irred k-points read from W2k',kmesh%ktot
   close(10)

   allocate(kmesh%k_coord(3,kmesh%ktot))
   allocate(kmesh%k_weight(kmesh%ktot))
   allocate(kmesh%mult(kmesh%ktot))

   !get the number of nonequivalent atoms in cell
   open(10,file=trim(adjustl(mysyst))//'.struct',status='old')
   read(10,*)
   read(10,*) cntr, ccrap, ccrap, iatm
   !write(*,*) 'number of inequivalent atoms in cell',iatm
   close(10)

   !get k-points coord's and energy dispersion curves
   open(11,file=trim(adjustl(mysyst))//'.energy',status='old')
   do i=1,2*iatm
     read(11,*)
   enddo
   edspk%nband_max=0
   do ik=1,kmesh%ktot
      read(11,*)kmesh%k_coord(1,ik),kmesh%k_coord(2,ik),kmesh%k_coord(3,ik),icrap,icrap,nband_loc

      ! allocate and initialise the temporary array
      if (.not. allocated (band_tmp)) allocate(band_tmp(kmesh%ktot,2*nband_loc))
      band_tmp(ik,:) = band_fill_value  !set the band dispersion to a large value

      if (nband_loc .gt. edspk%nband_max) edspk%nband_max=nband_loc

      ! this check has to be already made here
      if (edspk%nband_max .gt. size(band_tmp,2)) then
         write(*,*) 'you are trying to access energy bands that have not been stored by Wien2k'
         STOP
      endif

      do nband=1,nband_loc
         read(11,*) icrap,band_tmp(ik,nband)
      enddo
   enddo
   !done reading the energy dispersion curves
   close(11)

   !now fill in the actual data structure by cropping band_tmp
   if (.not. allocated(edspk%band) ) allocate(edspk%band(kmesh%ktot,edspk%nband_max))
   edspk%band(:,:) = band_fill_value
   do ik=1,kmesh%ktot
      do nband=1,edspk%nband_max
         edspk%band(ik,nband) = band_tmp(ik,nband)
      enddo
   enddo

   deallocate(band_tmp)

   !get k-points weight
   !TODO: remove because we never use weighted k-point sums
   open(12,file=trim(adjustl(mysyst))//'.klist',status='old')
   do ik=1,kmesh%ktot
      read(12,*)icrap,icrap,icrap,icrap,icrap,dtmp
      kmesh%mult(ik)=dtmp
      kmesh%k_weight(ik)=1.d0/dtmp
   enddo
   close(12)

end subroutine !getirrk

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! GETIRROPT
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine extracts from the Wien2k files
! information about the optical matrix elements,
! this is stored in the edisp type variable Mopt
! the first column runs over the polarisation index
! with a convention consistent with the optic module
! in Wien2k, namely:
! 1...Re<x><x>
! 2...Re<y><y>
! 3...Re<z><z>
! 4...Re<x><y>
! 5...Re<x><z>
! 6...Re<y><z>
! columns 4-6 are allocated only for non-cubic systems.
! If DMFT reference state is used the one-particle energies
! are overwritten and self-energy data are read in.
!
subroutine getirropt (mysyst, kmesh, eirrk, loptic, lcubic, ldmft )
     implicit none

!passed variables
  character(len=*) :: mysyst   ! file name to open
  type(kpointmesh) :: kmesh ! irreducible k-mesh generated by Wien2k
  type(edisp)      :: eirrk ! energy dispersion and optical matrix elements over irr k-mesh generated by Wien2k
  logical :: loptic
  logical :: lcubic
  logical :: ldmft          ! read self-energy file if true
!internal variables
  integer :: icrap, itmp, i, j, ik
  integer :: ierr
  integer :: nband, nband_loc, min_nopt, max_nopt
  real :: rcrap
  character(len=6) :: ccrap
  double precision :: dtmp
  double precision, allocatable :: Mopt_tmp(:,:,:,:)  ! temporary matrices where the Wien2k optical matrices are stored


   open(10,file=trim(adjustl(mysyst))//'.symmat',status='old')
   read(10,*)   !there is a heading line specifying which component of M is printed on file
   ! allocate and initialise the temporary arrays
   if (lcubic) then
      allocate(Mopt_tmp(1:3,kmesh%ktot,1:eirrk%nband_max,1:eirrk%nband_max))
   else
      allocate(Mopt_tmp(1:6,kmesh%ktot,1:eirrk%nband_max,1:eirrk%nband_max))
   endif
   Mopt_tmp=0.d0

   eirrk%nbopt_min=1000
   eirrk%nbopt_max=0
   if (loptic) then
   do ik=1,kmesh%ktot
      read(10,*)   !there is an empty line
      read(10,*)ccrap,itmp,ccrap,ccrap,ccrap,min_nopt,max_nopt
      read(10,*)   !there is an empty line

      !sanity tests
      if ( itmp .ne. ik) then
         write(*,*) 'ERROR: there is a mismatch between k-points in case.energy and case.symmat'
         STOP
      endif
      if ( max_nopt .gt. eirrk%nband_max) then
         write(*,*) 'ERROR: there are more bands computed in the optical routine than we have energies for'
         STOP
      endif

      !identify the smallest min_nopt and the biggest max_nopt
      if (eirrk%nbopt_min .gt. min_nopt) eirrk%nbopt_min = min_nopt
      if (eirrk%nbopt_max .lt. max_nopt) eirrk%nbopt_max = max_nopt
      ! (basically I want the M matrices to have the same size for all the k-points considered)
      !read the matrix elements
      if (lcubic) then
         do i=min_nopt,max_nopt
            do j=i,max_nopt
               read(10,130)icrap,icrap,Mopt_tmp(1,ik,i,j),Mopt_tmp(2,ik,i,j),Mopt_tmp(3,ik,i,j)
            enddo
         enddo
      else
         do i=min_nopt,max_nopt
            do j=i,max_nopt
               read(10,160)icrap, icrap, Mopt_tmp(1,ik,i,j),Mopt_tmp(2,ik,i,j),Mopt_tmp(3,ik,i,j),&
                 Mopt_tmp(4,ik,i,j),Mopt_tmp(5,ik,i,j),Mopt_tmp(6,ik,i,j)
            enddo
         enddo
      endif

   enddo !over kpoints
   endif

   130  FORMAT (4X,I3,X,I3,3(X,E12.6))
   160  FORMAT (4X,I3,X,I3,6(X,E12.6))
   close(10)
!!!!TEST
   !write(*,*) eirrk%nbopt_min,eirrk%nbopt_max,size(eirrk%Mopt_tmp,3)
   !write(*,*) 'temporary structure'
   !write(*,*) icrap
   !write(*,*) Mopt_tmp(1,1,10,30)
   !write(*,*) Mopt_tmp(2,1,10,31)
   !write(*,*) Mopt_tmp(3,1,10,32)
   !write(*,*) Mopt_tmp(1,10,1,3)
   !write(*,*) Mopt_tmp(2,10,1,3)
   !write(*,*) Mopt_tmp(3,10,1,3)
   !STOP
!!!!TEST END
   !eirrk%nbopt_max=int(eirrk%nbopt_max/10)
   ! copy over the data to the permanent data structure (with the minimal number of entries consistent for all the k-point considered)
   if ((.not. allocated(eirrk%Mopt)) .and. (lcubic) ) &
     & allocate(eirrk%Mopt(1:3,kmesh%ktot,eirrk%nbopt_min:eirrk%nbopt_max,eirrk%nbopt_min:eirrk%nbopt_max))
   if ((.not. allocated(eirrk%Mopt)) .and. (.not. lcubic) ) &
     & allocate(eirrk%Mopt(1:6,kmesh%ktot,eirrk%nbopt_min:eirrk%nbopt_max,eirrk%nbopt_min:eirrk%nbopt_max))

   itmp=6
   if (lcubic) itmp=3

   do i=1,itmp
      do ik=1,kmesh%ktot
         eirrk%Mopt(i,ik,eirrk%nbopt_min:eirrk%nbopt_max,eirrk%nbopt_min:eirrk%nbopt_max)=&
           & Mopt_tmp(i,ik,eirrk%nbopt_min:eirrk%nbopt_max,eirrk%nbopt_min:eirrk%nbopt_max)
      enddo
   enddo

   deallocate(Mopt_tmp)

!!!!TEST
   !write(*,*) eirrk%nbopt_min,eirrk%nbopt_max,size(eirrk%Mopt,3)
   !write(*,*) 'permanent structure'
   !write(*,*) eirrk%Mopt(1,1,1,1)
   !write(*,*) eirrk%Mopt(2,1,1,1)
   !write(*,*) eirrk%Mopt(3,1,1,1)
   !write(*,*) eirrk%Mopt(4,1,1,1)
   !write(*,*) eirrk%Mopt(5,1,1,1)
   !write(*,*) eirrk%Mopt(6,1,1,1)
   !STOP
!!!!TEST END

   !allocation of renormalised bandstructure
   !the -Im{Sigma} read in here is added to the temperature dependent scattering rate
   !in the response module
   if (.not. allocated(eirrk%Z)) then
      !allocate(eirrk%Z(kmesh%ktot,eirrk%nbopt_min:eirrk%nbopt_max))
      allocate(eirrk%Z(kmesh%ktot,eirrk%nband_max))
      eirrk%Z=1.0d0
   endif
   if (.not. allocated(eirrk%Im)) then
      !allocate(eirrk%Im(kmesh%ktot,eirrk%nbopt_min:eirrk%nbopt_max))
      allocate(eirrk%Im(kmesh%ktot,eirrk%nband_max))
      eirrk%Im=0.0d0
   endif
   if (ldmft) then
      open(11,file=trim(adjustl(mysyst))//'.dmft',status='old')
      read(11,*) eirrk%efer
      do ik=1,kmesh%ktot
         read(11,*) i, min_nopt, max_nopt
         !sanity checks
         if (min_nopt < eirrk%nbopt_min) STOP 'Self-energy bands do not match optical matrix elements'
         if (max_nopt > eirrk%nbopt_max) STOP 'Self-energy bands do not match optical matrix elements'
         !overwrite the bandstructure with the renormalised one
         do j=min_nopt,max_nopt
            read(11,*) eirrk%band(ik,j), eirrk%Z(ik,j), eirrk%Im(ik,j)
         enddo
      enddo
      close(11)
   else ! if there are no dmft data copy over the values read from input file
      eirrk%Z(:,:)=eirrk%ztmp
   endif

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! GETSYMOP
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine extracts from the Wien2k files
! information about the symmetry operations that
! acting on the irreducible kmesh generate the reducible
! counterpart (see genredk below)
!
subroutine getsymop (mysyst, symm, kmesh, ek)
  implicit none

  character(len=*) :: mysyst
  type(symop) :: symm
  type(kpointmesh) :: kmesh ! irreducible k-mesh generated by Wien2k
  type(edisp) :: ek
  integer :: n, i, j
  !for the additional detection of roto-translation in non-symmorphic groups
  character (len=80) :: line
  character (len=30) :: substring
  character (len=6)  :: ccrap
  integer :: ix, icrap
  !if centering is required additional temporary vectors are required
  !double precision, allocatable :: Mtmp(:,:,:), Ttmp(:,:)
  double precision, allocatable :: k_coord(:,:)
  double precision, allocatable :: band(:,:)


! read the matrix transformations
  open(9,file=trim(adjustl(mysyst))//'.symop',status='old')
  read(9,100)symm%nsym
  if (.not. allocated(symm%Msym)) allocate(symm%Msym(3,3,symm%nsym))
  if (.not. allocated(symm%Tras)) allocate(symm%Tras(3,symm%nsym))
  do n=1,symm%nsym
     read(9,110) ((symm%Msym(j,i,n),i=1,3),j=1,3)
  enddo
  close(9)

! if the space group is non-symmorphic we have to read in also the
! translation vectors to generate the Brillouin zone
  substring="NUMBER OF SYMMETRY OPERATIONS"
  ix=0
  symm%Tras(:,:) = 0.0d0
  open(10,file=trim(adjustl(mysyst))//'.struct',status='old')
  do i=1, 100
     read(10,'(A)') line
     ix = index(line, substring)
     if (ix /= 0) exit
  enddo

  if (ix == 0) then
     write(*,*) 'GETSYMOP: error reading file *.struct'
     STOP
  else
     do n=1,symm%nsym
        do j=1,3
           !thanks to Matthias for pointing out this trick! ;-)
           read(10,'(A)') line
           substring=trim(adjustl(line(7:)))
           read(substring,*) symm%Tras(j,n)
        enddo
        read(10,*) icrap
        !write(*,'(I3, I3, 6A, 3F10.8)') icrap, n, ccrap, symm%Tras(:,n)
        !write(*,*) icrap, n, ccrap, symm%Tras(:,n)
     enddo
  endif
  close(10)

  ! if the unit cell is not primitive we need to include centering
  select case (symm%cntr)

     case('B  ') !body centered cell
        write(*,*) 'body centering unit cell'
        !allocate(Mtmp(3,3,2*symm%nsym))
        !allocate(Ttmp(3,2*symm%nsym))
        !do n=1,symm%nsym
        !   Mtmp(:,:,n)= symm%Msym(:,:,n)
        !   Ttmp(:,n)  = symm%Tras(:,n)
        !   Mtmp(:,:,n+symm%nsym)= symm%Msym(:,:,n)
        !   Ttmp(:,n+symm%nsym)  = modulo(symm%Tras(:,n)+0.5d0,1.0d0)
        !enddo
        !deallocate(symm%Msym); deallocate(symm%Tras)
        !symm%nsym=2*symm%nsym
        !allocate(symm%Msym(3,3,symm%nsym))
        !allocate(symm%Tras(3,symm%nsym))
        !do n=1,symm%nsym
        !   !do i=1,3 ; do j=1,3
        !   !symm%Msym(i,j,n)= Mtmp(j,i,n)  !do I take the direct or the inverse rotation??
        !   !enddo
        !   !symm%Tras(i,n)  = Ttmp(i,n)
        !   !enddo
        !   symm%Msym(:,:,n)= Mtmp(:,:,n)  !do I take the direct or the inverse rotation??
        !   symm%Tras(:,n)  = Ttmp(:,n)
        !enddo
        !deallocate(Mtmp); deallocate(Ttmp)
        !allocate(band(2*kmesh%ktot, ek%nband_max))
        !allocate(k_coord(3,2*kmesh%ktot))
        !do i =1, kmesh%ktot
        !   k_coord(:,i)=kmesh%k_coord(:,i)
        !   k_coord(:,kmesh%ktot+i)=modulo(kmesh%k_coord(:,i)+0.5d0,1.0d0)
        !   band(i,:)=ek%band(i,:)
        !   band(kmesh%ktot+i,:)=ek%band(i,:)
        !enddo
        !kmesh%ktot=2*kmesh%ktot
        !deallocate(ek%band,kmesh%k_coord)
        !allocate(ek%band(kmesh%ktot,ek%nband_max))
        !allocate(kmesh%k_coord(3,kmesh%ktot))
        !do i=1,kmesh%ktot
        !   kmesh%k_coord(:,i)=k_coord(:,i)
        !   ek%band(i,:)=band(i,:)
        !enddo
        !deallocate(k_coord,band)

     case('F  ') !face centered cell
        write(*,*) 'face centering unit cell'
        !STOP

     case('H  ') !hexagonal cell?
        write(*,*) 'hexagonal centering unit cell'
        !STOP

     case('CXY') !base centered along xy plane
        write(*,*) 'base centering unit cell'
        !STOP

     case('CXZ') !base centered along xz plane
        write(*,*) 'base centering unit cell'
        !STOP

     case('CYZ') !base centered along yz plane
        write(*,*) 'base centering unit cell'
        !STOP

  end select

  100  FORMAT (I6)
  110  FORMAT (3(3f8.5/))
!!!!!!!!!!!!!!!!!!!!!!!TEST
!  do n=1,symm%nsym
!     write(40,120) ((symm%Msym(j,i,n),i=1,3), symm%Tras(j,n),j=1,3)
!  enddo
!  120  FORMAT (3(4f8.5/))
!  STOP
!!!!!!!!!!!!!!!!!!!!!!!TEST END

! set up the symmetry group descriptors
  if ((kmesh%a(1).eq.kmesh%a(2)) .and. (kmesh%a(3).eq.kmesh%a(2))) then
     symm%lcubic=.true.
     symm%ltetrag=.false.
     symm%lorthor=.false.
  elseif((kmesh%a(1).eq.kmesh%a(2)) .and. (kmesh%a(3).ne.kmesh%a(2))) then
     symm%lcubic=.false.
     symm%ltetrag=.true.
     symm%lorthor=.false.
  elseif((kmesh%a(1).eq.kmesh%a(3)) .and. (kmesh%a(3).ne.kmesh%a(2))) then
     symm%lcubic=.false.
     symm%ltetrag=.true.
     symm%lorthor=.false.
  elseif((kmesh%a(2).eq.kmesh%a(3)) .and. (kmesh%a(1).ne.kmesh%a(2))) then
     symm%lcubic=.false.
     symm%ltetrag=.true.
     symm%lorthor=.false.
  elseif((kmesh%a(1).ne.kmesh%a(2)) .and. &
   (kmesh%a(2).ne.kmesh%a(3)) .and. (kmesh%a(3).ne.kmesh%a(1)) ) then
     symm%lcubic=.false.
     symm%ltetrag=.false.
     symm%lorthor=.true.
  else
     write(*,*) 'ERROR: could not work out symmetry group'
     STOP
  endif

! set the descriptor for nonsymmorphic space groups
  do n=1,symm%nsym
     if ((symm%Tras(1,n)/=0.0d0).or.(symm%Tras(2,n)/=0.0d0).or.(symm%Tras(3,n)/=0.0d0)) then
        symm%lnsymmr = .true.
     endif
  enddo

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
subroutine genredk (mysyst, symm, kmesh, redkm )
   implicit none
   !passed variables
   character(len=*) :: mysyst ! file name to open
   type(symop)      :: symm   ! contains symmetry operations
   type(kpointmesh) :: kmesh  ! irreducible k-mesh generated by Wien2k
   type(kpointmesh) :: redkm  ! reducible k-mesh generated here

   !internal variables
   integer :: icrap, itmp,itmp2,ntmp, i,j,k, ik, isym, isym2, nnsym, ineq
   integer :: iexist, itest
   integer :: G0(3,7)
   real    :: rcrap
   double precision, allocatable :: tmpkall(:,:), tmpkall2(:,:)
   integer, allocatable          :: tmpoall(:,:), tmpoall2(:,:)
   double precision              :: tmpk(3), tmpk2(3)

   itmp=kmesh%ktot*symm%nsym ! number of k-points generated by the symmetry operations acting on the irreducible mesh
   itmp2=1

   data G0/ &
     &      0,0,1, 0,1,0, 1,0,0, 0,1,1, 1,1,0, 1,0,1, 1,1,1 /

   if (symm%lnsymmr) then
      nnsym = symm%nsym
   else
      nnsym = 1
   endif

   ineq=0 ! counter for the unique k-points
   ! temporary arrays which gets extended on the fly
   ! if we find another unique k-point

   do ik=1,kmesh%ktot  !loop over irredk
      do isym=1,symm%nsym
         do isym2=1,nnsym
            !for symmorphic groups the factor system is equal to one for all k-points
            do j=1,3
            ! create the new k-vector
               tmpk(j) = ((kmesh%k_coord(1,ik)*symm%Msym(1,j,isym))) &
                  & + ((kmesh%k_coord(2,ik)*symm%Msym(2,j,isym))) &
                  & + ((kmesh%k_coord(3,ik)*symm%Msym(3,j,isym)))
            enddo

            !!!!!!!!!!!!!!!!!!!!!!!!
            !If the space-group is nonsymmorphic we have to
            !construct the Herring's group which is the direct
            !product of the representations of the point-group
            !and the non-unitary translations
            !!!!!!!!!!!!!!!!!!!!!!!!
            if (symm%lnsymmr) then
               do j=1,3
                  tmpk2(j)=tmpk(j)*exp(-2.0d0*pi*ci*dot_product(G0(:,6),symm%Tras(:,isym2)))
                  ! the choice of this reciprocal lattice vector (G0(:,6)) is a bit euristic
                  ! (tested on FeSi 3x3x3, 4x4x4, 5x5x5)
               enddo
            else
               tmpk2 = tmpk
            endif

            iexist = 0
            do j=1,ineq
               if ( (abs(  tmpk2(1) - tmpkall(1,j) ) < 1.d-1/real(itmp)) &
               & .and. (abs(  tmpk2(2) - tmpkall(2,j) ) < 1.d-1/real(itmp)) &
               & .and. (abs(  tmpk2(3) - tmpkall(3,j) ) < 1.d-1/real(itmp)) ) then
                  iexist=1
                  exit
               endif
            enddo
            if ((iexist==0) .and. (ineq==0)) then
               allocate(tmpkall(3,1), tmpoall(2,1))
               tmpkall(:,1) = tmpk2
               tmpoall(1,1) = ik; tmpoall(2,1) = isym
               ineq = ineq+1
            else if (iexist==0) then
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
               tmpkall(:,ineq) = tmpk2 ! add at the end
               tmpoall(1,ineq) = ik; tmpoall(2,ineq) = isym ! we only save the rotation number
               deallocate(tmpkall2, tmpoall2)
            endif
         enddo
      enddo
   enddo

   !!!!!!!!!!!!!!!!!!!!!!!TEST
   ! do k=1, size(tmpkall,2)
   !    write(*,'(A,I4,3f8.4)')'KP',k,tmpkall(1,k),tmpkall(2,k),tmpkall(3,k)
   ! enddo

   !!!!!!!!!!!!!!!!!!!!!!!TEST END


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! What we've got right now is a k-mesh with equivalent endpoints in a given
   ! direction if that direction is sampled by an even number of k-points
   ! either one has to assign weight 1/2 to those or remove the double counting
   ! when the BZ is translated (in trnredk below)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! save the number of redkp found into the data structure
   redkm%kx = kmesh%kx+(1-mod(kmesh%kx,2))
   redkm%ky = kmesh%ky+(1-mod(kmesh%ky,2))
   redkm%kz = kmesh%kz+(1-mod(kmesh%kz,2))
   redkm%ktot = max(redkm%kx*redkm%ky*redkm%kz,ik)

   ! if ((ineq .ne. redkm%ktot) .and. (mod(redkm%ktot,ineq) .ne. 0) ) then
   !    write(*,*) 'GENREDK: the number of k-points generated by symmetry is inconsistent',ineq,redkm%ktot,redkm%kx,redkm%ky,redkm%kz
   !    STOP
   ! else if (ineq .ne. redkm%ktot) then
   !    write(*,*) 'GENREDK: WARNING: number of generated k-points inconsistent', ineq, redkm%ktot, redkm%kx, redkm%ky, redkm%kz
   !    write(*,*) 'GENREDK: Modulo ( kpoints, generated points) = ', mod(redkm%ktot,ineq)
   ! endif
   if (.not. allocated(redkm%mult))     allocate(redkm%mult(1:redkm%ktot))
   if (.not. allocated(redkm%k_weight)) allocate(redkm%k_weight(1:redkm%ktot))
   if (.not. allocated(redkm%k_coord))  allocate(redkm%k_coord(3,redkm%ktot))
   if (.not. allocated(redkm%k_id))     allocate(redkm%k_id(redkm%kx, redkm%ky, redkm%kz))
   if (.not. allocated(symm%symop_id))  allocate(symm%symop_id(2,redkm%ktot))
   redkm%k_coord=0.0d0
   redkm%k_id=0
   symm%symop_id=0

   ! save the new coordinates into the data structure
   redkm%k_coord(:,:ineq) = tmpkall
   symm%symop_id(:,:ineq) = tmpoall
   deallocate(tmpkall, tmpoall)

end subroutine !GENREDK

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! GENREDOPT
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine generates the optical matrix elements
! starting from those read in on the irreducible k-mesh
!
subroutine genredopt (redkm, symm, eirrk, eredk )
   implicit none
!passed variables
   type(kpointmesh) :: redkm ! reducible k-mesh
   type(symop) :: symm       ! contains symmetry operations
   type(edisp) :: eirrk      ! energy dispersion and optical matrix elements over irr k-mesh generated by Wien2k
   type(edisp) :: eredk
!internal variables
   integer :: i, j, l, k, ik, isym, iks, nb, nb2
   double precision, allocatable :: osm1(:,:), osm2(:,:), osm3(:,:)
   double precision, allocatable :: Mtmp(:,:,:,:) !only necessary for non-cubic systems

   eredk%nband_max = eirrk%nband_max
   eredk%nbopt_max = eirrk%nbopt_max
   eredk%nbopt_min = eirrk%nbopt_min
   eredk%nelect = eirrk%nelect
   eredk%efer = eirrk%efer

   if (.not. allocated(eredk%Mopt) .and. (symm%lcubic)) &
     allocate(eredk%Mopt(1:3,redkm%ktot,eredk%nbopt_min:eredk%nbopt_max,eredk%nbopt_min:eredk%nbopt_max))
   if (.not. allocated(eredk%Mopt) .and. (.not.symm%lcubic)) &
     allocate(eredk%Mopt(1:6,redkm%ktot,eredk%nbopt_min:eredk%nbopt_max,eredk%nbopt_min:eredk%nbopt_max))
   if (.not. allocated(eredk%band)) &
     allocate(eredk%band(redkm%ktot,eredk%nband_max))
   if (.not. allocated(eredk%Z)) &
     allocate(eredk%Z(redkm%ktot,eredk%nband_max))
   if (.not. allocated(eredk%Im)) &
     allocate(eredk%Im(redkm%ktot,eredk%nband_max))
   eredk%Mopt(:,:,:,:)= 0.d0
   eredk%band(:,:)   = 0.d0
   eredk%Im(:,:)    = 0.d0
   eredk%Z(:,:)    = 0.d0
   if (symm%lcubic) then
      allocate(osm1(eredk%nbopt_min:eredk%nbopt_max,eredk%nbopt_min:eredk%nbopt_max))
      allocate(osm2(eredk%nbopt_min:eredk%nbopt_max,eredk%nbopt_min:eredk%nbopt_max))
      allocate(osm3(eredk%nbopt_min:eredk%nbopt_max,eredk%nbopt_min:eredk%nbopt_max))
   else
      allocate(Mtmp(3,3,eredk%nbopt_min:eredk%nbopt_max,eredk%nbopt_min:eredk%nbopt_max))
   endif

   ! scan the symm%symop_id to check what are the combinations of k-point and symmetry that produce non-redundant element in redkm
   if (symm%lcubic) then
      do i=1,size(symm%symop_id, 2) ! should be equivalent to redkm%ktot but currently the cerium compound misbehaves
         ik  = symm%symop_id(1,i) ! this counter runs over the irrekpoints
         isym= symm%symop_id(2,i)

         ! TODO: FIX THIS
         ! happens if we only recover a divisior of the kmesh
         if (isym .eq. 0) cycle   ! temporary fix
         !

         osm1=0.d0; osm2=0.d0; osm3=0.d0
         ! check eq. 13.16 (pg 479) in "Symmetry and Condensed Matter Physics A Computational Approach" by M. El-Batanouny, F. Wooten, CUP
         do j=1,3
            do nb=eredk%nbopt_min,eredk%nbopt_max ; do nb2=eredk%nbopt_min,eredk%nbopt_max
               osm1(nb,nb2) = osm1(nb,nb2) + eirrk%Mopt(j,ik,nb,nb2)*symm%Msym(j,1,isym)*symm%Msym(j,1,isym)
               osm2(nb,nb2) = osm2(nb,nb2) + eirrk%Mopt(j,ik,nb,nb2)*symm%Msym(j,2,isym)*symm%Msym(j,2,isym)
               osm3(nb,nb2) = osm3(nb,nb2) + eirrk%Mopt(j,ik,nb,nb2)*symm%Msym(j,3,isym)*symm%Msym(j,3,isym)
            enddo ; enddo
         enddo
         eredk%Mopt(1,i,:,:)=osm1(:,:)
         eredk%Mopt(2,i,:,:)=osm2(:,:)
         eredk%Mopt(3,i,:,:)=osm3(:,:)
      enddo
      deallocate(osm1,osm2,osm3)

      !!!!!!!!!!!!!!!!TEST
      !write(666,'(A,I6,3f8.4)')'KP ',k,redkm%k_coord(1,k),redkm%k_coord(2,k),redkm%k_coord(3,k)
      !do nb=eredk%nbopt_min,eredk%nbopt_max
      !   do nb2=nb,eredk%nbopt_max
      !      write(666, '(2(I4,X),3(E12.6,X))')nb,nb2,eredk%Mopt(1,k,nb,nb2),eredk%Mopt(2,k,nb,nb2),eredk%Mopt(3,k,nb,nb2)
      !   enddo
      !enddo
      !STOP
      !!!!!!!!!!!!!!!!TEST END

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
      do i=1,redkm%ktot
         ik  = symm%symop_id(1,i) ! this counter runs over the irrekpoints
         isym= symm%symop_id(2,i)
         ! copy over the optical matrix elements
         Mtmp(1,1,:,:) = eirrk%Mopt(1,ik,:,:)
         Mtmp(2,2,:,:) = eirrk%Mopt(2,ik,:,:)
         Mtmp(3,3,:,:) = eirrk%Mopt(3,ik,:,:)
         Mtmp(1,2,:,:) = eirrk%Mopt(4,ik,:,:)
         Mtmp(1,3,:,:) = eirrk%Mopt(5,ik,:,:)
         Mtmp(2,3,:,:) = eirrk%Mopt(6,ik,:,:)
         Mtmp(2,1,:,:) = Mtmp(1,2,:,:)
         Mtmp(3,1,:,:) = Mtmp(1,3,:,:)
         Mtmp(3,2,:,:) = Mtmp(2,3,:,:)

         ! do the contraction with the symmetry matrices
         ! check eq. 13.16 (pg 479) in "Symmetry and Condensed Matter Physics A Computational Approach" by M. El-Batanouny, F. Wooten, CUP
         do j = 1,3
            do l = 1,3
               eredk%Mopt(1,k,:,:) = eredk%Mopt(1,k,:,:) + symm%Msym(j,1,isym)*Mtmp(j,l,:,:)*symm%Msym(l,1,isym)
               eredk%Mopt(2,k,:,:) = eredk%Mopt(2,k,:,:) + symm%Msym(j,2,isym)*Mtmp(j,l,:,:)*symm%Msym(l,2,isym)
               eredk%Mopt(3,k,:,:) = eredk%Mopt(3,k,:,:) + symm%Msym(j,3,isym)*Mtmp(j,l,:,:)*symm%Msym(l,3,isym)
               eredk%Mopt(4,k,:,:) = eredk%Mopt(4,k,:,:) + symm%Msym(j,1,isym)*Mtmp(j,l,:,:)*symm%Msym(l,2,isym)
               eredk%Mopt(5,k,:,:) = eredk%Mopt(5,k,:,:) + symm%Msym(j,1,isym)*Mtmp(j,l,:,:)*symm%Msym(l,3,isym)
               eredk%Mopt(6,k,:,:) = eredk%Mopt(6,k,:,:) + symm%Msym(j,2,isym)*Mtmp(j,l,:,:)*symm%Msym(l,3,isym)
            enddo
         enddo
      enddo
      deallocate (Mtmp)
   endif


   !now map the dispersion energies from the old grid to the new one (the energies do not change with the symmetry operation)
   do i=1,size(symm%symop_id,2)
      ik  = symm%symop_id(1,i) ! this counter runs over the irrekpoints
      isym= symm%symop_id(2,i)

      ! FIXME
      if (isym .eq. 0 ) cycle

      do nb=1,eredk%nband_max
         eredk%band(i,nb) = eirrk%band(ik,nb)
         eredk%Im(i,nb) = eirrk%Im(ik,nb)
         eredk%Z(i,nb) = eirrk%Z(ik,nb)
      enddo
   enddo

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! TRNREDK
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine translates the reducible
! k-mesh with elements in the interval (-pi/a, +pi/a)
! into the interval (0, 2pi/a), taking care also of
! the band structure and of the optical matrix elements
!
subroutine trnredk (irrkm, redkm, eredk, symm, algo)
  implicit none
!passed variables
  type(kpointmesh) :: irrkm
  type(kpointmesh) :: redkm
  type(edisp) :: eredk
  type(symop) :: symm
  type(algorithm) :: algo

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

  if (symm%lcubic) then
     offdia=0  !no off-diagonal terms for cubic systems
  else
     offdia=1  !off-diagonal terms for non-cubic systems
  endif
  !allocation temporary arrays
  allocate(bstmp(redkm%ktot,eredk%nband_max)); allocate(bstmp2(redkm%ktot,eredk%nband_max))
  allocate(Imtmp(redkm%ktot,eredk%nband_max)); allocate(Imtmp2(redkm%ktot,eredk%nband_max))
  allocate(Ztmp(redkm%ktot,eredk%nband_max));  allocate(Ztmp2(redkm%ktot,eredk%nband_max))
  allocate(cktmp(3,redkm%ktot)); allocate(cktmp2(3,redkm%ktot))
  if (symm%lcubic) then
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
     if ((redkm%k_coord(1,ik) > 1.0d0) .or. (redkm%k_coord(2,ik) > 1.0d0) .or.(redkm%k_coord(3,ik) > 1.0d0)) &
        STOP 'TRNREDK: something is seriously wrong here (e.g. improper traslation bigger than lattice constant)'
     do i = 1,3
        if (redkm%k_coord(i,ik)<0.0d0) then
           cktmp(i,ik)=mod(1.0d0+redkm%k_coord(i,ik),1.0d0)
        else
           cktmp(i,ik)=mod(redkm%k_coord(i,ik),1.0d0)
        endif
     enddo
     do ibn =1,eredk%nband_max
        bstmp(ik,ibn)= eredk%band(ik,ibn)
        Imtmp(ik,ibn)= eredk%Im(ik,ibn)
        Ztmp(ik,ibn) = eredk%Z(ik,ibn)
     enddo
  enddo

  !copy the optical transition matrix elements
  do i = 1,3+(offdia*3)
     do ik=1,redkm%ktot
        do ibn =eredk%nbopt_min,eredk%nbopt_max
           do ibn2=eredk%nbopt_min,eredk%nbopt_max
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
  do i =ntmp, redkm%ktot
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
  if (symm%lcubic) then
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
  do ik=1,nk
     do i=1,3
        redkm%k_coord(i,ik)=cktmp2(i,ik)
     enddo
     do ibn=1,eredk%nband_max
        eredk%band(ik,ibn)= bstmp2(ik,ibn)
        eredk%Im(ik,ibn)  = Imtmp2(ik,ibn)
        eredk%Z(ik,ibn)   = Ztmp2(ik,ibn)
     enddo
  enddo

  do i=1,3+(offdia*3)
     do ik=1,nk
        do ibn=eredk%nbopt_min,eredk%nbopt_max
           do ibn2=eredk%nbopt_min,eredk%nbopt_max
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
subroutine genfulkm( redkm, fulkm, eredk, efulk, symm)
 implicit none
!passed variables
 type(kpointmesh) :: redkm
 type(kpointmesh) :: fulkm
 type(edisp) :: eredk
 type(edisp) :: efulk
 type(symop) :: symm
!local variables
 integer :: i, ik, ikx, iky, ikz, ibn, ibn2
 integer :: nk, nkx, nky, nkz, nband
 integer :: offdia !off-diagonal terms in the Mopt matrix?
 double precision :: dk(3), tmp1, tmp2

!initialisation new datatype
  efulk%nband_max = eredk%nband_max
  efulk%nbopt_max = eredk%nbopt_max
  efulk%nbopt_min = eredk%nbopt_min
  efulk%nelect = eredk%nelect
  efulk%efer = eredk%efer

  if (symm%lcubic) then
     offdia=0  !no off-diagonal terms for cubic systems
  else
     offdia=1  !off-diagonal terms for non-cubic systems
  endif

  nband = eredk%nband_max
  nk=redkm%ktot; nkx=redkm%kx; nky=redkm%ky; nkz=redkm%kz

!add the missing k-point to the grid, so it includes the terminal BZ points
  fulkm%kx = nkx+1; fulkm%ky = nky+1; fulkm%kz = nkz+1
  fulkm%ktot = fulkm%kx * fulkm%ky * fulkm%kz
  fulkm%a    = redkm%a
  fulkm%alat = redkm%alat
  !if (.not. allocated(fulkm%mult))     allocate(fulkm%mult(1:fulkm%ktot))      *** TODO: REDUNDANT BECAUSE NEVER USED,
  !if (.not. allocated(fulkm%k_weight)) allocate(fulkm%k_weight(1:fulkm%ktot))      HENCE NOT ALLOCATED, CONSIDER REMOVING FROM DATATYPE ***
  if (.not. allocated(fulkm%k_id))    allocate(fulkm%k_id(1:fulkm%kx, 1:fulkm%ky, 1:fulkm%kz))
  if (.not. allocated(fulkm%k_coord)) allocate(fulkm%k_coord(3,fulkm%ktot))
  !bandstructure allocation
  if (.not. allocated(efulk%band)) allocate(efulk%band(fulkm%ktot, nband))
  if (.not. allocated(efulk%Im))   allocate(efulk%Im(fulkm%ktot, nband))
  if (.not. allocated(efulk%Z))    allocate(efulk%Z(fulkm%ktot, nband))
  if (.not. allocated(efulk%Mopt) .and. (offdia==0)) &
    allocate(efulk%Mopt(3,fulkm%ktot,efulk%nbopt_min:efulk%nbopt_max,efulk%nbopt_min:efulk%nbopt_max))
  if (.not. allocated(efulk%Mopt) .and. (offdia==1)) &
    allocate(efulk%Mopt(6,fulkm%ktot,efulk%nbopt_min:efulk%nbopt_max,efulk%nbopt_min:efulk%nbopt_max))

   fulkm%k_id(:,:,:) = 0
   fulkm%k_coord(:,:)= 0.d0
   efulk%band(:,:) = 0.d0
   efulk%Im(:,:)  = 0.d0
   efulk%Z(:,:)  = 0.d0
   efulk%Mopt(:,:,:,:)=0.d0

  do ik=1,nk
     do i=1,3
        fulkm%k_coord(i,ik)=redkm%k_coord(i,ik)
     enddo
     do ibn=1,nband
        efulk%band(ik,ibn)= eredk%band(ik,ibn)
        efulk%Im(ik,ibn)  = eredk%Im(ik,ibn)
        efulk%Z(ik,ibn)   = eredk%Z(ik,ibn)
     enddo
  enddo

  do i=1,3+(offdia*3)
     do ik=1, nk
        do ibn=efulk%nbopt_min,efulk%nbopt_max
           do ibn2=efulk%nbopt_min,efulk%nbopt_max
             efulk%Mopt(i,ik,ibn,ibn2)=eredk%Mopt(i,ik,ibn,ibn2)
           enddo
        enddo
     enddo
  enddo

  do ikx=1,nkx
     do iky=1,nky
        do ikz=1,nkz
           fulkm%k_id(ikx,iky,ikz)=redkm%k_id(ikx,iky,ikz)
        enddo
     enddo
  enddo

  !Now that the k-point identifier has been copied over,
  !generate the minimum displacement necessary to
  !obtain the full BZ by adding it to the existing mesh
  do i=1,3
     tmp1 = 1.0d0
     do ik=2,nk
        tmp2 = abs(fulkm%k_coord(i,1) - fulkm%k_coord(i,ik))
        if ((tmp2>0.0d0) .and. (tmp2<tmp1)) tmp1=tmp2
     enddo
     dk(i)=tmp1
  enddo

  if (dk(1)==0 .or. dk(2)==0 .or. dk(3)==0) then
     write(*,*) 'GENFULKM can not extend the k-mesh', dk(:)
     STOP
  endif

  !add on the terminal point to the BZ
  ik=nk
  !select one face (z=const)
  do ikx=1,nkx
     do iky=1,nky
        ik = ik+1
        fulkm%k_id(ikx,iky,fulkm%kz)=ik
        fulkm%k_coord(1,ik)=fulkm%k_coord(1,fulkm%k_id(ikx,iky,nkz))
        fulkm%k_coord(2,ik)=fulkm%k_coord(2,fulkm%k_id(ikx,iky,nkz))
        fulkm%k_coord(3,ik)=fulkm%k_coord(3,fulkm%k_id(ikx,iky,nkz))+dk(3)
        do ibn=1,nband
           efulk%band(ik,ibn)= efulk%band(fulkm%k_id(ikx,iky,1),ibn)
           efulk%Im(ik,ibn)  = efulk%Im(fulkm%k_id(ikx,iky,1),ibn)
           efulk%Z(ik,ibn)   = efulk%Z(fulkm%k_id(ikx,iky,1),ibn)
           if ((ibn>efulk%nbopt_max) .or. (ibn<efulk%nbopt_min)) cycle
           do ibn2=efulk%nbopt_min,efulk%nbopt_max
             efulk%Mopt(1,ik,ibn,ibn2)=efulk%Mopt(1,fulkm%k_id(ikx,iky,1),ibn,ibn2)
             efulk%Mopt(2,ik,ibn,ibn2)=efulk%Mopt(2,fulkm%k_id(ikx,iky,1),ibn,ibn2)
             efulk%Mopt(3,ik,ibn,ibn2)=efulk%Mopt(3,fulkm%k_id(ikx,iky,1),ibn,ibn2)
             if (offdia == 1) then
                efulk%Mopt(4,ik,ibn,ibn2)=efulk%Mopt(4,fulkm%k_id(ikx,iky,1),ibn,ibn2)
                efulk%Mopt(5,ik,ibn,ibn2)=efulk%Mopt(5,fulkm%k_id(ikx,iky,1),ibn,ibn2)
                efulk%Mopt(6,ik,ibn,ibn2)=efulk%Mopt(6,fulkm%k_id(ikx,iky,1),ibn,ibn2)
             endif
           enddo
        enddo
     enddo
  enddo
  !select the second face (y=const)
  do ikx=1,nkx
     do ikz=1,fulkm%kz
        ik = ik+1
        fulkm%k_id(ikx,fulkm%ky,ikz)=ik
        fulkm%k_coord(1,ik)=fulkm%k_coord(1,fulkm%k_id(ikx,nky,ikz))
        fulkm%k_coord(2,ik)=fulkm%k_coord(2,fulkm%k_id(ikx,nky,ikz))+dk(2)
        fulkm%k_coord(3,ik)=fulkm%k_coord(3,fulkm%k_id(ikx,nky,ikz))
        do ibn=1,nband
           efulk%band(ik,ibn)= efulk%band(fulkm%k_id(ikx,1,ikz),ibn)
           efulk%Im(ik,ibn)  = efulk%Im(fulkm%k_id(ikx,1,ikz),ibn)
           efulk%Z(ik,ibn)   = efulk%Z(fulkm%k_id(ikx,1,ikz),ibn)
           if ((ibn>efulk%nbopt_max) .or. (ibn<efulk%nbopt_min)) cycle
           do ibn2=efulk%nbopt_min,efulk%nbopt_max
             efulk%Mopt(1,ik,ibn,ibn2)=efulk%Mopt(1,fulkm%k_id(ikx,1,ikz),ibn,ibn2)
             efulk%Mopt(2,ik,ibn,ibn2)=efulk%Mopt(2,fulkm%k_id(ikx,1,ikz),ibn,ibn2)
             efulk%Mopt(3,ik,ibn,ibn2)=efulk%Mopt(3,fulkm%k_id(ikx,1,ikz),ibn,ibn2)
             if (offdia == 1) then
                efulk%Mopt(4,ik,ibn,ibn2)=efulk%Mopt(4,fulkm%k_id(ikx,1,ikz),ibn,ibn2)
                efulk%Mopt(5,ik,ibn,ibn2)=efulk%Mopt(5,fulkm%k_id(ikx,1,ikz),ibn,ibn2)
                efulk%Mopt(6,ik,ibn,ibn2)=efulk%Mopt(6,fulkm%k_id(ikx,1,ikz),ibn,ibn2)
             endif
           enddo
        enddo
     enddo
  enddo
  !select the third face (x=const)
  do iky=1,fulkm%ky
     do ikz=1,fulkm%kz
        ik = ik+1
        fulkm%k_id(fulkm%kx,iky,ikz)=ik
        fulkm%k_coord(1,ik)=fulkm%k_coord(1,fulkm%k_id(nkx,iky,ikz))+dk(1)
        fulkm%k_coord(2,ik)=fulkm%k_coord(2,fulkm%k_id(nkx,iky,ikz))
        fulkm%k_coord(3,ik)=fulkm%k_coord(3,fulkm%k_id(nkx,iky,ikz))
        do ibn=1,nband
           efulk%band(ik,ibn)= efulk%band(fulkm%k_id(1,iky,ikz),ibn)
           efulk%Im(ik,ibn)  = efulk%Im(fulkm%k_id(1,iky,ikz),ibn)
           efulk%Z(ik,ibn)   = efulk%Z(fulkm%k_id(1,iky,ikz),ibn)
           if ((ibn>efulk%nbopt_max) .or. (ibn<efulk%nbopt_min)) cycle
           do ibn2=efulk%nbopt_min,efulk%nbopt_max
             efulk%Mopt(1,ik,ibn,ibn2)=efulk%Mopt(1,fulkm%k_id(1,iky,ikz),ibn,ibn2)
             efulk%Mopt(2,ik,ibn,ibn2)=efulk%Mopt(2,fulkm%k_id(1,iky,ikz),ibn,ibn2)
             efulk%Mopt(3,ik,ibn,ibn2)=efulk%Mopt(3,fulkm%k_id(1,iky,ikz),ibn,ibn2)
             if (offdia == 1) then
                efulk%Mopt(4,ik,ibn,ibn2)=efulk%Mopt(4,fulkm%k_id(1,iky,ikz),ibn,ibn2)
                efulk%Mopt(5,ik,ibn,ibn2)=efulk%Mopt(5,fulkm%k_id(1,iky,ikz),ibn,ibn2)
                efulk%Mopt(6,ik,ibn,ibn2)=efulk%Mopt(6,fulkm%k_id(1,iky,ikz),ibn,ibn2)
             endif
           enddo
        enddo
     enddo
  enddo

  !!!!!!!!!!!!!!!!TEST
  !do ik=1,fulkm%ktot
  !   write(778,171)'KP ',ik,fulkm%k_coord(1,ik),fulkm%k_coord(2,ik),fulkm%k_coord(3,ik)
  !   do ibn=efulk%nbopt_min,efulk%nbopt_max
  !      do ibn2=ibn,efulk%nbopt_max
  !         write(778,170)ibn,ibn2,efulk%Mopt(1,ik,ibn,ibn2),efulk%Mopt(2,ik,ibn,ibn2),efulk%Mopt(3,ik,ibn,ibn2)
  !      enddo
  !   enddo
  !enddo
  !170  FORMAT  (2(I4,X),3(E12.6,X))
  !171  FORMAT  (A,I4,X,3(F12.6,X))
  !!!!!!!!!!!!!!!!TEST END

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
subroutine gentetra (algo, mesh, thdr )
 implicit none
 type(algorithm)  :: algo
 type(kpointmesh) :: mesh
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

  ! ntetd is the maximum number of tetrahedra for the given mesh
  ! given that there are 6 tetraedra in each cubic cell it is reasonable to set
  ntetd=6*(mesh%kx)*(mesh%ky)*(mesh%kz)
  allocate (idtet(0:4,ntetd))
  if (algo%ltbind) then
     nkx=mesh%kx+1; nky=mesh%ky+1; nkz=mesh%kz+1
  else
     nkx=mesh%kx; nky=mesh%ky; nkz=mesh%kz
  endif
      data kcut0/ &
     &         0,0,0, 0,1,0, 1,1,0, 1,1,1,  0,0,0, 1,0,0, 1,1,0, 1,1,1, &
     &         0,0,0, 1,0,0, 1,0,1, 1,1,1,  0,0,0, 0,1,0, 0,1,1, 1,1,1, &
     &         0,0,0, 0,0,1, 0,1,1, 1,1,1,  0,0,0, 0,0,1, 1,0,1, 1,1,1 /

 ! need to generate the cartesian basis for the cell
 ! for a simple cubic lattice the reciprocal lattice is also cubic and shrunk by a factor 2pi/alat
 bk=0.d0
 !generalisation to tetragonal and orthorhombic cases:
 bk(1,1)=mesh%a(2)*mesh%a(3)*(2.d0*pi/mesh%alat)
 bk(2,2)=mesh%a(3)*mesh%a(1)*(2.d0*pi/mesh%alat)
 bk(3,3)=mesh%a(1)*mesh%a(2)*(2.d0*pi/mesh%alat)

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
            imc(k1,k2,k3)=mesh%k_id(j1,j2,j3)
            !write(35,*) 'cube',icube, iccor, mesh%k_id(j1,j2,j3),j1,j2,j3,mesh%k_coord(:,mesh%k_id(j1,j2,j3))
            if (iccor==8) then
             vlcub=(abs(mesh%k_coord(3,mesh%k_id(j1,j2,j3))-mesh%k_coord(3,mesh%k_id(j1,j2,j3-1))))**3
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
      write(*,15) ntet,mesh%kx*mesh%ky*mesh%kz
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
               p(i,ic)=mesh%k_coord(i,idtet(ic,itet))
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
subroutine intetra (mesh, thdr, dos, band, nband_max)
  implicit none

   type(kpointmesh) :: mesh
   type(edisp)      :: ek
   type(tetramesh)  :: thdr
   type(dosgrid)    :: dos
   real(8),intent(in):: band(:,:)
   integer,intent(in):: nband_max

!local variables
   integer :: i, j, i00, itet, nb, istart, istop
   integer :: iq(4)
   double precision :: de, ec(4), ec1(4), es(4)
   double precision :: e1, e2, e3, e4
   double precision :: c0, c1, c2, c3, cc12, cc34  ! constants
   double precision :: wthdr   ! weight of the tetrahedron
   double precision :: eact, x ! free energy variables
   double precision :: adddos  ! accumulation variable for the dos

! SANITY CHECK
   if (mesh%ktot<4) then
      write(*,*)'INTETRA: tetrahedron method fails (number of k-points < 4)',mesh%ktot
      STOP
   endif
! initialize arrays for dos/number of states
   !write(*,*)'INTETRA: constructing energy mesh'
   if (.not. allocated(dos%enrg)) allocate (dos%enrg(1:dos%nnrg))
   if (.not. allocated(dos%dos )) allocate (dos%dos(1:dos%nnrg))
   if (.not. allocated(dos%nos )) allocate (dos%nos(1:dos%nnrg))
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
      iq(1) = thdr%idtet(1,itet)
      iq(2) = thdr%idtet(2,itet)
      iq(3) = thdr%idtet(3,itet)
      iq(4) = thdr%idtet(4,itet)
      wthdr = thdr%vltet(itet)

      do nb=1,nband_max
        if (band(iq(1),nb)> band_fill_value) cycle
        if (band(iq(2),nb)> band_fill_value) cycle
        if (band(iq(3),nb)> band_fill_value) cycle
        if (band(iq(4),nb)> band_fill_value) cycle
! get the band energy at each corner of the tetrahedron:
        ec(1) = band(iq(1),nb)
        ec(2) = band(iq(2),nb)
        ec(3) = band(iq(3),nb)
        ec(4) = band(iq(4),nb)

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

 return

end subroutine   !INTETRA

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! FINDEF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This routine finds the Fermi level using the
! dos%nos variable and the number of electrons
! provided as an input
!
subroutine findef(dos, ek)
  implicit none

   type(dosgrid) :: dos
   type(edisp)   :: ek
!local variables
   double precision :: F(4), P(4)
   double precision :: s
   double precision :: psave, ptol, ntol
   integer  :: I(4), iter, maxiter, itmp

! initialise the varibles
    I(1)= 1
    I(2)= dos%nnrg
    P(1)= dos%enrg(1)
    P(2)= dos%enrg(dos%nnrg)
    F(1)= dos%nos(1)-ek%nelect
    F(2)= dos%nos(dos%nnrg)-ek%nelect
    ptol   =  1.0d-16
    psave  = -1.1d30
    maxiter= 60

    do iter = 1, maxiter
       itmp = I(1)+I(2)
       I(3) = int(itmp/2)
       P(3) = dos%enrg(I(3))
       F(3) = dos%nos(I(3))-ek%nelect
       s = sqrt((F(3)**2)-(F(1)*F(2)))
       if (s==0.0d0) then
          write(*,*) 'Error in Ridders search for Fermi level'
          write(*,*) 'ITER', iter, 'x1', P(1),'  x2',P(2),'  x3', P(3)
          write(*,*) 'ITER', iter, 'F1', F(1),'  F2',F(2),'  F3', F(3)
          goto 400
       endif
       I(4) = I(3)+(I(3)-I(1))*int(sign(1.0d0,F(1)-F(2))*F(3)/s)
       P(4) = dos%enrg(I(4))

       if(abs(P(4)-psave)<=ptol) goto 400
       psave= P(4)
       F(4) = dos%nos(I(4))-ek%nelect
       if (F(4) ==0.0d0) goto 400
       if (sign(F(3), F(4)) /= F(3)) then
       !change of sign btw x3 and x4 then reduce search interval
          I(1)  = I(3)
          P(1)  = P(3)
          F(1)  = F(3)
          I(2)  = I(4)
          P(2)  = P(4)
          F(2)  = F(4)
       else if (sign(F(1), F(4)) /= F(1)) then
       !change of sign btw x1 and x4 then reduce search interval
          I(2)  = I(4)
          P(2)  = P(4)
          F(2)  = F(4)
       else if (sign(F(2), F(4)) /= F(2)) then
       !change of sign btw x2 and x4 then reduce search interval
          I(1)  = I(4)
          P(1)  = P(4)
          F(1)  = F(4)
       endif
       !condition for termination
       if (abs(P(2)-P(1)) <= ptol) goto 400
    enddo ! over number of iterations
    write (*,*) 'here 6'
 400 if (iter == maxiter) write(*,*) 'Ridders seach might not have converged'
    ek%efer=P(4)
    !find the band gap
    ntol=4.0d-2
    I(3)  = I(4)
    P(3)  = P(4)
    F(3)  = F(4)

    do while(abs(F(3)-F(4)) < ntol)
       I(3) = I(3)-1
       P(3) = dos%enrg(I(3))
       F(3) = dos%nos(I(3))-ek%nelect
    enddo
    dos%vbm=P(3)

    I(3)  = I(4)
    P(3)  = P(4)
    F(3)  = F(4)
    do while(abs(F(3)-F(4)) < ntol)
       I(3) = I(3)+1
       P(3) = dos%enrg(I(3))
       F(3) = dos%nos(I(3))-ek%nelect
    enddo
    dos%cbm=P(3)
    dos%gap=dos%cbm - dos%vbm
    if (dos%gap < 2.0d-2) dos%gap=0.0d0

end subroutine ! FINDEF

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
subroutine interptra_re (iT, itet, mu, nalpha, lBoltz, mesh, ek, thdr, sct, resp, hpresp )
  implicit none

   class(dp_resp)   :: resp
   type(kpointmesh) :: mesh
   type(edisp)      :: ek
   type(tetramesh)  :: thdr
   type(scatrate)   :: sct
   type(qp_resp),optional :: hpresp
   integer, intent(in)           :: iT     !temperature index
   integer, intent(in)           :: itet   !tetrahedron identifier
   double precision, intent (in) :: mu     !chemical potential
   integer, intent(in)           :: nalpha !truncates the number of polarisation directions required in a cubic system
   logical :: lBoltz
!local variables
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
   iq(1) = thdr%idtet(1,itet)
   iq(2) = thdr%idtet(2,itet)
   iq(3) = thdr%idtet(3,itet)
   iq(4) = thdr%idtet(4,itet)

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

      do ix=1,nalpha
         do iy=ix,nalpha
            do i=1,4 !linear interpolation within the  tetrahedron
               ! set the value of the scattering rate for the specific temperature, band
               if (allocated(sct%ykb)) then
                  resp%gamma=real(ek%z(thdr%idtet(i,itet),nb)*(sct%gam(iT)+sct%ykb(iT,thdr%idtet(i,itet),nb)),8)
                  if (present(hpresp)) hpresp%gamma=real(ek%z(thdr%idtet(i,itet),nb)*(sct%gam(iT)+ &
                                                    sct%ykb(iT,thdr%idtet(i,itet),nb)),16)
               else
                  resp%gamma=real(ek%z(thdr%idtet(i,itet),nb)*sct%gam(iT),8)
                  if (present(hpresp)) hpresp%gamma=real(ek%z(thdr%idtet(i,itet),nb)*sct%gam(iT),16)

               endif
               if (lBoltz) then
                  ! no interband transitions for Boltzmann response
                  resp%s(nb,ix,iy) = resp%s(nb,ix,iy) + (real(w(i),8)*resp%s_tmp(i,nb,ix,iy))/(resp%gamma)
                  resp%a(nb,ix,iy) = resp%a(nb,ix,iy) + (real(w(i),8)*resp%a_tmp(i,nb,ix,iy))/(resp%gamma)
                  if (allocated(resp%sB)) then
                     resp%sB(nb,ix,iy) = resp%sB(nb,ix,iy) + &
                                       (real(w(i),8)*resp%sB_tmp(i,nb,ix,iy))/(resp%gamma**2)
                     resp%aB(nb,ix,iy) = resp%aB(nb,ix,iy) + &
                                       (real(w(i),8)*resp%aB_tmp(i,nb,ix,iy))/(resp%gamma**2)
                  endif
               else
                  select type(resp)
                     type is(dp_resp)
                     ! intraband response
                     resp%s(nb,ix,iy) = resp%s(nb,ix,iy) + (real(w(i),8)*resp%s_tmp(i,nb,ix,iy))*beta/(resp%gamma )
                     resp%a(nb,ix,iy) = resp%a(nb,ix,iy) + (real(w(i),8)*resp%a_tmp(i,nb,ix,iy))*(beta**2)/(resp%gamma )
                     type is(dp_respinter)
                     ! interband response
                     resp%s(nb,ix,iy) = resp%s(nb,ix,iy) + (real(w(i),8)*resp%s_tmp(i,nb,ix,iy)) !all the prefactors have been included
                     resp%a(nb,ix,iy) = resp%a(nb,ix,iy) + (real(w(i),8)*resp%a_tmp(i,nb,ix,iy)) !in RESPINTERT or in RESDERTET
                  end select

                  if (allocated(resp%sB)) then
                     select type(resp)
                        type is(dp_resp)
                        resp%sB(nb,ix,iy) = resp%sB(nb,ix,iy) + &
                                          (real(w(i),8)*resp%sB_tmp(i,nb,ix,iy))*beta/(resp%gamma**2)
                        resp%aB(nb,ix,iy) = resp%aB(nb,ix,iy) + &
                                          (real(w(i),8)*resp%aB_tmp(i,nb,ix,iy))*beta/(resp%gamma**2)
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
               do i=1,4 !linear interpolation within the  tetrahedron
                  hpresp%s(nb,ix,iy) = hpresp%s(nb,ix,iy) + &
                                     (real(w(i),16)*hpresp%s_tmp(i,nb,ix,iy))*betaQ/(hpresp%gamma)
                  hpresp%a(nb,ix,iy) = hpresp%a(nb,ix,iy) + &
                                     (real(w(i),16)*hpresp%a_tmp(i,nb,ix,iy))*(betaQ**2)/(hpresp%gamma)
                  if (allocated(hpresp%sB)) then
                     hpresp%sB(nb,ix,iy) = hpresp%sB(nb,ix,iy) + &
                                         (real(w(i),16)*hpresp%sB_tmp(i,nb,ix,iy))*betaQ/(hpresp%gamma**2)
                     hpresp%aB(nb,ix,iy) = hpresp%aB(nb,ix,iy) + &
                                         (real(w(i),16)*hpresp%aB_tmp(i,nb,ix,iy))*(betaQ/hpresp%gamma)**2
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

   do i=1,4 !linear interpolation within the  tetrahedron
      target_intp = target_intp + (w(i)*target_tet(i))
   enddo !over corners of the tetrahedron

 return
end subroutine !INTERPTRA_MUQ

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!SUBROUTINE GENDOSEL
! Generates the Density of States starting
! from the Wien2k eigenvalues
! by replacing the Dirac's delta with Lorentzians
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine gendosel(mesh, ek, dos)
  implicit none

   type(kpointmesh) :: mesh
   type(edisp)      :: ek
   type(dosgrid)    :: dos
!local variables
   integer :: i, ik, nb !energy, k-point, band counters
   double precision :: br, de !broadening, energy spacing

   allocate (dos%enrg(dos%nnrg),dos%dos(dos%nnrg),dos%nos(dos%nnrg))
   dos%enrg= 0.d0
   dos%dos = 0.d0
   dos%nos = 0.d0
! get the energy increment along the window fixed in the dos datatype
   de=(dos%emax-dos%emin)/(real(dos%nnrg-1))
   do i=0,dos%nnrg-1
      dos%enrg(i+1)=dos%emin+(de*real(i))
   enddo
   !broadening twice the energy spacing
   br = 2.0d0*de
   !lorentian bandshape
   do i =1,size(dos%enrg)
      do nb=1,ek%nband_max
         do ik=1,mesh%ktot
            if (ek%band(ik,nb) > band_fill_value) cycle !necessary because big eig'vals
                                               !have been introduced to make matrices square
            dos%dos(i)=dos%dos(i)+((br/pi)*(1.0d0/(((dos%enrg(i)-ek%band(ik,nb))**2)+(br**2))))
            dos%nos(i)=dos%nos(i)+(0.5d0 + ((1.0d0/pi)*atan((dos%enrg(i)-ek%band(ik,nb))/br)))
         enddo
      enddo
   enddo

   dos%dos = dos%dos / mesh%ktot ! normalizing
   dos%nos = dos%nos / mesh%ktot

end subroutine !GENDOSEL

  double precision function ek_sc(k,iband,eirrk)
   implicit none
   type(edisp) :: eirrk
   double precision k(3),ek,bw
   integer iband,i

   ek=eirrk%E0(iband)
   bw=eirrk%t(iband,1)
   do i=1,3
      ek=ek + 2.d0*bw*cos(2.d0*pi*k(i))
   enddo
   ek_sc=ek

  return
  end function ek_sc

  double precision function vk_sc(idir,k,iband,eirrk,kmesh)
   implicit none
   type(edisp) :: eirrk
   type(kpointmesh) :: kmesh
   double precision k(3),bw
   integer iband,i,idir

   bw=eirrk%t(iband,1)
   vk_sc=-2.d0*bw*sin(2.d0*pi*k(idir))*kmesh%a(idir)*kmesh%alat

  return
  end function vk_sc

  double precision function vkk_sc(idir,idir2,k,iband,eirrk,kmesh)
   implicit none
   type(edisp) :: eirrk
   type(kpointmesh) :: kmesh
   double precision k(3),bw
   integer iband,idir,idir2

   bw =eirrk%t(iband,1)
   if (idir.eq.idir2) then
      vkk_sc=-2.d0*bw*cos(2.d0*pi*k(idir))*(kmesh%a(idir)*kmesh%alat)**2
   else
      vkk_sc=0.d0
   endif

  return
  end function vkk_sc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! This function evaluates the length of the tetrahedron edge
!! required by subroutine GENTETRA
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  double precision  function anrm2(x,y,z)
   implicit none
   double precision :: x, y, z
   anrm2=x*x*1.00001d0+y*y*1.00002d0+z*z*1.00003d0 &
     &             -x*0.000004d0-y*0.000003d0-z*0.000002d0
  end function

end module Mestruct
