

! somehow mpif90 on hclm doesn like the more elegant iso way...    
! #define qp 16
! #define dp 8

module types

type algorithm
   logical :: ldebug       ! evaluate quad precision data?
   logical :: ltbind       ! tight binding lattice?
   logical :: ltetra       ! use tetrahedron method?
   logical :: lw2k         ! use Wien2k input  
   logical :: loptic       ! use also optic matrix elements
   logical :: lBfield      ! calculations in the presence of a magnetic field?
   logical :: lsymm        ! reads in the full BZ or generates it using symmetry operations?
   integer :: imurestart   ! restart with a privided value of the chemical potential?
   character(10) :: mysyst ! label that is used to open the Wien2k files 
end type

type kpointmesh
   integer, allocatable :: k_id(:,:,:)           ! counter idenfifying a specific k-point through values of ikx, iky, ikz
   double precision, allocatable :: k_coord(:,:) ! coordinates in reciprocal space associated to a k-point (exclunding the BZ endpoints)
   double precision, allocatable :: k_weight(:)  ! k-point weights  
   double precision, allocatable :: mult(:)      ! k-point multiplicity  
   integer :: kx                                 ! number of k-points in each cartesian direction
   integer :: ky
   integer :: kz
   integer :: ktot                               ! total number of k-points (which varies depending on the sampling method for the BZ)
   double precision :: alat                      ! lattice constant
   double precision :: a(3)                      ! direct lattice vectors in units of ALAT. needed for derivatives --> units of vk's... ???
   integer :: nkeq                               ! number of direction with the same number of k-point (nkeq=3 uniform k-mesh)
end type

type edisp
   double precision, allocatable :: E0(:)            ! TIGHT BINDING PARAMETER: band energy at Gamma point (nbands)
   double precision, allocatable :: t(:,:)           ! TIGHT BINDING PARAMETER: band width (nbands,tmax), tmax/=1 is used to set anisotropic dispersion
   double precision, allocatable :: M2(:,:,:,:)      ! TIGHT BINDING PARAMETER: 2nd derivative of the energy dispersion (idir,idir2,nk,nbands), not used if ltbind=F
   integer :: nband_max                              
   integer :: nbopt_min                              ! number of bands (interval) included in the optical matrix elements
   integer :: nbopt_max                              ! number of bands (interval) included in the optical matrix elements
   double precision :: efer                          ! Fermi energy 
   double precision :: nelect                        ! number of electrons (provided from inputfile at first iteration)
   double precision :: occ_tot                       ! number of electrons computed for the updated value of chemical potential (non-eq value) 
   double precision, allocatable :: band(:,:)        ! band(k_id,nband) contains nband_max dispersion curves, if for some k-points
                                                     ! fewer bands are computed, then the missing elements are set to 1000 (arbitrarily large value)
   double precision, allocatable :: occ(:,:)         ! occ(k_id,nband) contains the occupation numbers  
   double precision, allocatable :: Mopt(:,:,:,:)    ! M(x,k,n',n)= | <n',k|p.e_x|n,k> |^2
   double precision, allocatable :: Mopt_tetra(:,:,:,:)    ! M(x,t,n',n)= sum_k(j) {w(j)*|<n',k(j)|p.e_x|n,k(j)>|^2 } with w(j) the tetrahedron weights
end type

type symop
   integer, allocatable :: symop_id(:,:)        ! symop_id(3,irrkp*nsym) this counter tells me if for a given k-point (2nd entry) a given symmetry operation (3rd entry) 
                                                ! produces a new element or a redundant one (1 or 0 in the 1st entry) 
   integer :: nsym                              ! number of symmetry operations 
   double precision, allocatable :: Msym(:,:,:) ! Msym(3,3,nsym) matrix containing the 3x3 symmetry transformations 
   double precision, allocatable :: Tras(:,:)   ! Tras(3,nsym) matrix containing additional lattice traslations for non-symmorphic groups 
   logical :: lcubic
   logical :: ltetrag
   logical :: lorthor
   logical :: lnsymmr                           ! .true. for nonsymmorphic space groups
   logical :: lBZpnt 
end type

type tetramesh
   integer :: ntet                            ! number of inequivalent tetrahedra found starting from an equally spaced k-grid      
   integer, allocatable :: idtet(:,:)         ! dimension (0:4,ntet); idtet(0,ntet) is the tetrahedron's multiplicity, 1:4 are the identifiers of the vertices
   double precision :: vltot                  ! vltot is the total tetrahedron's volume (in reciprocal space)
   double precision, allocatable :: vltet(:)  ! dimension (ntet); vltet(itet) is the tetrahedron's volume (in reciprocal space)
   double precision, allocatable :: nocc(:)   ! dimension (nbands); number of occupied tetrahedra at a given temperature      
end type

type dosgrid
   integer :: nnrg                         ! number of points in the energy window
   double precision :: emin                ! bottom of the energy window
   double precision :: emax                ! top of the energy window
   double precision :: vbm                 ! valence band maximum    
   double precision :: cbm                 ! conduction band minimum 
   double precision :: gap                 ! band gap 
   double precision, allocatable :: enrg(:)! energy grid 
   double precision, allocatable :: dos(:) ! density of states (as computed in PRB,49,16223, appx C ) 
   double precision, allocatable :: nos(:) ! number  of states (as computed in PRB,49,16223, appx A ) 
end type

type scatrate
   !temperature grid
   integer :: nT                           ! number of points in the temperature window
   double precision :: Tmin                ! bottom of the temperature window
   double precision :: Tmax                ! top of the temperature window
   double precision :: dT                  ! temperature interval
   double precision, allocatable :: TT(:)  ! temperature grid 
   double precision, allocatable :: mu(:)  ! chemical potential (temperature dependent)
   !double precision, allocatable :: drhodT(:) ! 1st derivative of rho w.r.t. T 
   double precision, allocatable :: d1(:)  ! square of the 1st derivative of sigma 
   double precision, allocatable :: d2(:)  ! product of the 2nd derivative times sigma
   double precision, allocatable :: d0(:)  ! linear combination of the two above, whose zero corresponds to T*
   double precision :: Tstar               ! temperature for which (d^2 rho)/(d beta^2)=0
                                           ! in practice it is the T for which (d^2 sigma)/(d beta^2) changes sign
   double precision :: Tflat               ! T for which (d sigma)/(d beta) changes sign (onset of saturation)
   double precision :: z                   ! renormalisation factor (read from file) 
   !gamma 
   integer :: ng                           ! degree of the polynomial temperature dependence of gamma
   double precision, allocatable :: gc(:)  ! coefficient of the polynomial temperature dependence of gamma
   double precision, allocatable :: gam(:,:)  ! temparature and band dependent scattering rate
end type

end module types

module estruct
  real(8) :: gmax,gminall,vol ! gmax is maximal gamma of bands that determine [gap] 
  contains

   subroutine estruct_init(algo, kmesh, redkm, fulkm, eirrk, eredk, efulk, thdr, dos, sct)
     !use mpi_org
     use params
     use types
     implicit none

     type (algorithm) :: algo
     type(kpointmesh) :: kmesh  ! contains k-point mesh specifiers and logical switches on how to get the mesh from
     type(kpointmesh) :: redkm  ! contains k-point mesh specifiers and logical switches on how to get the mesh from
     type(kpointmesh) :: fulkm  ! contains k-point mesh specifiers and logical switches on how to get the mesh from
     type(edisp) :: eirrk       ! contains the band dispersion energy and the optical matrix elements (when which > 2) along the irr-k-mesh
     type(edisp) :: eredk       ! contains the band dispersion energy and the optical matrix elements (when which > 2) along the red-k-mesh
     type(edisp) :: efulk       ! contains the band dispersion energy and the optical matrix elements for the red-k-mesh including BZ endpoints
     type(symop) :: symm        ! contains the symmetry operations
     type(tetramesh) :: thdr    ! contains the tetrahedra (should you need them...)
     type(dosgrid) :: dos 
     type(scatrate):: sct 
     integer :: ierr, iband, i
     double precision :: maxtmp

!eM: ikstart and ikend are parameters set by a previous mpi call in main.f90
! I think that the parallelisation will have to take place after the construction of the 
! tetrahedra
     call readinfile(algo, eirrk, kmesh, sct)
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
            call intetra (kmesh, eirrk, thdr, dos)       ! computes the dos and the integrated dos
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
            call intetra (fulkm, efulk, thdr, dos)   ! computes the dos and the integrated dos
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
           dos%emax= 2.0d0*maxtmp
           dos%emin=-dos%emax
           dos%nnrg=2.0d0*dos%emax/0.007d0
           !energy spacing of roughly 0.1 eV
           
           call gendosel(kmesh, eirrk, dos)
           !spin multiplicity 
           do i=1,size(dos%enrg)
              dos%dos(i)=2.0d0*dos%dos(i)
              dos%nos(i)=2.0d0*dos%nos(i)
           enddo
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
              if ((eredk%band(1,iband) < 1.0d2) .and. (maxtmp < abs(eredk%band(1,iband)))) &
                 maxtmp = abs(eredk%band(1,iband))
           enddo             
           dos%emax= 20.0d0*maxtmp
           dos%emin=-dos%emax
           dos%nnrg=2.0d0*dos%emax/0.007d0
           !energy spacing of roughly 0.1 Ry
           write(*,*) 'estruct: Wien2k electronic structure linear k-mesh DOS calculation'

           call gendosel (redkm, eredk, dos)

           !spin multiplicity 
           do i=1,size(dos%enrg)
              dos%dos(i)=2.0d0*dos%dos(i)
              dos%nos(i)=2.0d0*dos%nos(i)
           enddo
           open(10,file='linrados',status='unknown')
           do i=1,size(dos%enrg)
              write (10,*) dos%enrg(i), dos%dos(i), dos%nos(i)
           enddo
           close(10)
           if (algo%imurestart == 0) call findef(dos, eredk)   ! finds the Fermi level
        endif

     endif !ltetra
   end subroutine !estruct_init

   subroutine readinfile(algo, ek, kmesh, sct)
     use types
     implicit none
    
    type (algorithm) :: algo
    type(edisp) :: ek          ! contains the band dispersion energy and the optical matrix elements (when which > 2) along the irr-k-mesh
    type(kpointmesh) :: kmesh  ! contains k-point mesh specifiers and logical switches on how to get the mesh from
    type(scatrate) :: sct
  !local variables
    integer :: which           ! switches the type of algorithm (then stored in type(algorithm))
    integer :: which_tetra     ! =0 no tetrahedron integration =1 use tetrahedron  
    integer :: which_grid      ! =0 no symetry operations used  =1 use symmetry to generate BZ   
    integer :: tmax            ! if not equal to 1 it introduces anisotropic dispersion in the TB model
    integer :: iband

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
    read(10,*)which, which_tetra, which_grid

    if (which==1) then
       algo%ltbind =.true.
       algo%lw2k   =.false.
       algo%loptic =.false.
    else if (which==2) then
       algo%ltbind =.false.
       algo%lw2k   =.true.
       algo%loptic =.false.
       algo%lBfield=.false.
    else if (which==3) then
       algo%ltbind =.false.
       algo%lw2k   =.false.
       algo%loptic =.true.
       algo%lBfield=.false.
    else 
       write(*,*)'No k-mesh selected'
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

    read(10,'(A)')algo%mysyst
    read(10,*)kmesh%kx,kmesh%ky,kmesh%kz
    !
    !TODO
    !consider removing from datatype
    if ((kmesh%kx == kmesh%ky) .AND. (kmesh%ky == kmesh%kz) ) then
       kmesh%nkeq=3
    elseif ((kmesh%kx == kmesh%ky) .AND. (kmesh%ky /= kmesh%kz)) then
       kmesh%nkeq=2
    elseif ((kmesh%kx == kmesh%kz) .AND. (kmesh%ky /= kmesh%kz)) then
       kmesh%nkeq=2
    else
       kmesh%nkeq=1
    endif
 
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
   !read in also the renormalisation factor 
    read(10,*)sct%z
   !originally there was a loop over bands that I'm skipping 
   !out of the complication to generate a long input file for 40 bands or so... 

    close(10)

   end subroutine !readinfile

   subroutine gentbstr(algo, eirrk, kmesh)
     use types
     implicit none

    type(algorithm) :: algo
    type(edisp) :: eirrk                 ! contains the band dispersion energy and the optical matrix elements (when which > 2) along the irr-k-mesh
    type(kpointmesh) :: kmesh            ! contains k-point mesh specifiers and logical switches on how to get the mesh from
  !local variables
    integer :: ik,ikx,iky,ikz, nk,nkx,nky,nkz, iband,nband, idir,idir2 !counters
    double precision :: k(3)
    double precision, external :: ek_sc, vk_sc, vkk_sc

   if (algo%ltetra) then
      nkx=kmesh%kx+1; nky=kmesh%ky+1; nkz=kmesh%kz+1
   else
      nkx=kmesh%kx; nky=kmesh%ky; nkz=kmesh%kz
   endif
   nband=eirrk%nband_max

   allocate(kmesh%k_id(1:nkx, 1:nky, 1:nkz))
   allocate(kmesh%k_coord( 3, 1:kmesh%ktot))
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
     use types
     implicit none
    type(algorithm) :: algo
    type(kpointmesh) :: kmesh  
    type(kpointmesh) :: redkm  
    type(edisp) :: eirrk       
    type(edisp) :: eredk       
    type(symop) :: symm
    integer :: ik,ikx,iky,ikz, nb, nb1

   if (algo%lsymm) then
      call getirrk  (algo%mysyst, kmesh, eirrk) 
      call getsymop (algo%mysyst, symm, kmesh)
      ! generate a reucible kmesh from the set of symmetry operations and the irrek-mesh
      call genredk  (algo%mysyst, symm, kmesh, redkm)
      ! copy also the number of electrons and the Fermi level from the old datastructure to the new one
      eredk%nband_max = eirrk%nband_max
      eredk%nelect = eirrk%nelect
      eredk%efer = eirrk%efer
      redkm%alat = kmesh%alat
      redkm%a    = kmesh%a
      ! read in the optical matrix elements on the irreducible grid 
      call getirropt (algo%mysyst, kmesh, eirrk, symm%lcubic)
      ! generate the optical matrix elements evaluated on the new redkm grid
      call genredopt (redkm, symm, eirrk, eredk)
      ! translate the reducible k-mesh and take care of the bandstructure
      call trnredk (kmesh, redkm, eredk, symm)

   else !reducible BZ is already provided by W2k
      call getirrk (algo%mysyst, redkm, eredk) 
      !the values below are stored in the irreducible mesh by readinfile
      eredk%nband_max = eirrk%nband_max
      eredk%nelect = eirrk%nelect
      eredk%efer = eirrk%efer
      redkm%alat = kmesh%alat
      redkm%a    = kmesh%a
      redkm%kx = kmesh%kx; redkm%ky = kmesh%ky; redkm%kz = kmesh%kz
      call getsymop (algo%mysyst, symm, redkm) !need to call this routine to set the logical switches in symmop type

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
      call getirropt ( algo%mysyst, redkm, eredk, symm%lcubic)
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

subroutine getirrk (mysyst, kmesh, edspk )
     use types
     implicit none

!passed variables
  character(10) :: mysyst ! file name to open
  type(kpointmesh) :: kmesh ! k-mesh generated by Wien2k
  type(edisp)      :: edspk ! energy dispersion and optical matrix elements over k-mesh generated by Wien2k
!internal variables
  integer :: icrap, itmp, nband, nband_loc, i, j, ik  
  integer :: iatm !number of non-equivalent atoms in cell
  real :: rcrap 
  character(27) :: ccrap
  double precision :: dtmp
  double precision, allocatable :: band_tmp(:,:)    ! bigger datastructure that will contain all the energy dispersion curves 

   !get the number of k-points 
   open(10,file=trim(adjustl(mysyst))//'.weight',status='old')
   read(10,*) rcrap, rcrap, ccrap
   read(10,*) kmesh%ktot, ccrap
   write(*,*) 'total number of irred k-points read from W2k',kmesh%ktot
   close(10)
 
   allocate(kmesh%k_coord( 3, 1:kmesh%ktot))
   allocate(kmesh%k_weight( 1:kmesh%ktot))
   allocate(kmesh%mult( 1:kmesh%ktot))
   
   !get the number of nonequivalent atoms in cell
   open(10,file=trim(adjustl(mysyst))//'.struct',status='old')
   read(10,*)
   read(10,*) ccrap, ccrap, ccrap, iatm
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
      band_tmp(ik,:) = 1.0d2  !set the band dispersion to a large value

      if (nband_loc .gt. edspk%nband_max) edspk%nband_max=nband_loc

      do nband=1,nband_loc
         read(11,*) icrap,band_tmp(ik,nband)
      enddo
   enddo
   !done reading the energy dispersion curves
   close(11)

   !now fill in the actual data structure by cropping band_tmp
   if (.not. allocated(edspk%band) ) allocate(edspk%band(kmesh%ktot,edspk%nband_max))
   edspk%band(:,:) = 1.0d2
   do ik=1,kmesh%ktot
     if (edspk%nband_max .le. size(band_tmp,2)) then
        do nband=1,edspk%nband_max
           edspk%band(ik,nband) = band_tmp(ik,nband)
        enddo
     else
         write(*,*) 'you are trying to access energy bands that have not been stored by Wien2k'
         STOP
     endif
   enddo
 
   deallocate(band_tmp)

   !get k-points weight 
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
! columns 4-6 are allocated only for non-cubic systems
!
subroutine getirropt (mysyst, kmesh, eirrk, lcubic )
     use types
     implicit none

!passed variables
  character(10) :: mysyst   ! file name to open
  type(kpointmesh) :: kmesh ! irreducible k-mesh generated by Wien2k
  type(edisp)      :: eirrk ! energy dispersion and optical matrix elements over irr k-mesh generated by Wien2k
  logical :: lcubic 
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

   enddo

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

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! GETSYMOP
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine extracts from the Wien2k files 
! information about the symmetry operations that 
! acting on the irreducible kmesh generate the reducible 
! counterpart (see genredk below) 
!
subroutine getsymop (mysyst, symm, kmesh)
     use types
  implicit none
   
  character(10) :: mysyst
  type(symop) :: symm
  type(kpointmesh) :: kmesh ! irreducible k-mesh generated by Wien2k
  integer :: n, i, j
  !for the additional detection of roto-translation in non-symmorphic groups
  character (len=80) :: line
  character (len=30) :: substring
  character (len=6)  :: ccrap
  integer :: ix, icrap
  !check for the presence of inversion symmetry and introduce it
  integer :: iexist, nsymold
  double precision, allocatable :: Mtmp(:,:,:)
  double precision, allocatable :: Ttmp(:,:)
 
  
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
           !read(10,'(6A,X,f10.8)') ccrap, symm%Tras(j,n) 
           !thanks to Mathias for pointing out this trick! ;-)
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

  100  FORMAT (I6)
  110  FORMAT (3(3f8.5/))

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
     use types
     use params
  implicit none 
!passed variables
  character(10) :: mysyst   ! file name to open
  type(symop) :: symm       ! contains symmetry operations
  type(kpointmesh) :: kmesh ! irreducible k-mesh generated by Wien2k
  type(kpointmesh) :: redkm ! reducible k-mesh generated here 
!internal variables
  integer :: icrap, itmp,itmp2,ntmp, i,j,k, ik,isym 
  integer :: iexist, itest
  integer :: redkpW2k
  integer :: G0(3,7)
  real :: rcrap
  double precision, allocatable :: cktmp(:,:),cktmp2(:,:),cktmp3(:,:) ! temporary k-point coordinates arrays
  double precision :: a, b, c, res
  complex :: iu

  
! READ THE NUMBER OF REDUCIBLE K-POINTS
   !fundamentally redundant, the value is never used and for non-uniform grids it is 0 
   open(12,file=trim(adjustl(mysyst))//'.klist',status='old')
   read(12,*)icrap,icrap,icrap,icrap,icrap, rcrap,rcrap,rcrap,redkpW2k  !not used at the moment
   close(12)

   itmp=kmesh%ktot*symm%nsym ! number of k-points generated by the symmetry operations acting on the irreducible mesh
   itmp2=1
   allocate(cktmp(3,itmp))
   if (symm%lnsymmr) then
      itmp2=kmesh%ktot*(symm%nsym**2)
      allocate(cktmp2(3,itmp2))
      allocate(cktmp3(3,itmp2))
      if (.not. allocated(symm%symop_id)) allocate(symm%symop_id(4,1:itmp2))
   else
      allocate(cktmp2(3,itmp))
      if (.not. allocated(symm%symop_id)) allocate(symm%symop_id(3,1:itmp))
   endif

      data G0/ &
     &         0,0,1, 0,1,0, 1,0,0, 0,1,1, 1,1,0, 1,0,1, 1,1,1 / 
   if (symm%lnsymmr) then
   ! Check if any of the given k-points lies on the Bragg hyperplanes
      do ik=1,kmesh%ktot  !loop over irredk
         do i =1,7
            a =real(G0(1,i))*(kmesh%k_coord(1,ik)-0.5d0*real(G0(1,i)))
            b =real(G0(2,i))*(kmesh%k_coord(2,ik)-0.5d0*real(G0(2,i)))
            c =real(G0(3,i))*(kmesh%k_coord(3,ik)-0.5d0*real(G0(3,i)))
            res =a+b+c
            if (abs(res)<1.0d-4) then
               write(*,'(A,3f8.5,X,A,3I3)') 'k-point', kmesh%k_coord(:,ik), 'lies on the Bragg plane normal to',G0(:,i)
               symm%lBZpnt=.true.
            endif
         enddo
      enddo
   endif
   
   cktmp = 0.0d0   
   cktmp2= 0.0d0
   symm%symop_id = 0
   i=0
   iu=cmplx(0.0d0,1.0d0)
   do ik=1,kmesh%ktot  !loop over irredk
      do isym=1,symm%nsym
         i=i+1
         symm%symop_id(2,i) = ik
         symm%symop_id(3,i) = isym
         !for symmorphic groups the factor system is equal to one for all k-points
         do j=1,3
            cktmp(j,i)=((kmesh%k_coord(1,ik)*symm%Msym(1,j,isym))) &
              & + ((kmesh%k_coord(2,ik)*symm%Msym(2,j,isym))) &
              & + ((kmesh%k_coord(3,ik)*symm%Msym(3,j,isym))) 
         enddo
      enddo
   enddo
   !!!!!!!!!!!!!!!!!!!!!!!TEST  
   !rotations only
   !write(*,*) 'write 555',size(cktmp,2)
   !do k=1, size(cktmp,2)
   !   write(555,'(A,I4,3f8.4)')'KP',k,cktmp(1,k),cktmp(2,k),cktmp(3,k) 
   !enddo
   !!!!!!!!!!!!!!!!!!!!!!!TEST END 
  
   !!!!!!!!!!!!!!!!!!!!!!!!
   !If the space-group is nonsymmorphic we have to 
   !construct the Herring's group which is the direct
   !product of the representations of the point-group
   !and the non-unitary translations 
   !!!!!!!!!!!!!!!!!!!!!!!!
   if (symm%lnsymmr) then
      do i=1,(itmp2/itmp)-1
         symm%symop_id(2,(i*itmp)+1:(i+1)*itmp)=symm%symop_id(2,1:itmp)
         symm%symop_id(3,(i*itmp)+1:(i+1)*itmp)=symm%symop_id(3,1:itmp)
      enddo
      i=0
      do isym=1,symm%nsym
         do k=1, size(cktmp,2)
            i=i+1
            symm%symop_id(4,i) = isym
            do j=1,3
               cktmp3(j,i)=cktmp(j,k)*exp(-2.0d0*pi*iu*dot_product(G0(:,6),symm%Tras(:,isym)))
               ! the choice of this reciprocal lattice vector is a bit euristic 
               ! (tested on FeSi 3x3x3, 4x4x4, 5x5x5)
            enddo
         enddo
      enddo
      !allocate the variable with the new size so the clean-up can proceed unchanged 
      !w.r.t. the symmorphic case
      deallocate(cktmp)
      allocate(cktmp(3,itmp2))
      do i=1, size(cktmp,2)
         do j=1,3
            cktmp(j,i)=cktmp3(j,i)
         enddo
      enddo
      deallocate(cktmp3)
      
      !!!!!!!!!!!!!!!!!!!!!!!TEST  
      !write(*,*) 'write 556',size(cktmp,2)
      !do k=1, size(cktmp,2)
      !   write(556,'(A,I4,3f8.4)')'KP',k,cktmp(1,k),cktmp(2,k),cktmp(3,k) 
      !enddo
      !!!!!!!!!!!!!!!!!!!!!!!TEST END  
   endif
   

   !!!!!!!!!!!!!!!!!!!!!!!!
   !Clean up 
   !check which points in the redBZ have been saved already
   !!!!!!!!!!!!!!!!!!!!!!!!
   ntmp=2   !exclude the gamma point 
   ik=1     
   do i =ntmp,size(cktmp,2)
      do itest=i+1,size(cktmp,2)
         if ( (abs(  cktmp(1,itest) - cktmp(1,i) ) <1.d-1/real(itmp)) & 
         & .and. (abs(  cktmp(2,itest) - cktmp(2,i) ) <1.d-1/real(itmp)) &  
         & .and. (abs(  cktmp(3,itest) - cktmp(3,i) ) <1.d-1/real(itmp))  ) cycle
         ! now it has found the first value of itest that corresponds to a different k-point from the 
         ! the starting one

         ! now need to copy over this k-point corresponding to itest into cktmp2 BUT
         ! making sure that it hasn't appeared there already
         iexist=0
         do j=1,ik
            if ( (abs(  cktmp(1,itest) - cktmp2(1,j) ) <1.d-1/real(itmp)) & 
            & .and. (abs(  cktmp(2,itest) - cktmp2(2,j) ) <1.d-1/real(itmp)) &  
            & .and. (abs(  cktmp(3,itest) - cktmp2(3,j) ) <1.d-1/real(itmp)) ) iexist=1
         enddo     
         if (iexist==0) then 
            ik=ik+1 
            cktmp2(1,ik) = cktmp(1,itest)
            cktmp2(2,ik) = cktmp(2,itest)
            cktmp2(3,ik) = cktmp(3,itest)
            symm%symop_id(1,itest)=1  !the combination of kpoint and isym in the compound index itest produces a non-redundant result
            ntmp=itest
         endif 
      enddo 
   enddo
  !the loop above leaves out the Gamma point; here I'm fixing MANUALLY this by saying that the only non-redundant symmetry operation 
  !for the first k-point (the Gamma point) is the identity:
  symm%symop_id(1,1) = 1
   
  deallocate (cktmp)
  !!!!!!!!!!!!!!!!!!!!!!!TEST  
  !write(*,*) 'write 666', ik
  !do k=1, ik
  !   write(666,'(A,I4,3f8.4)')'KP',k,cktmp2(1,k),cktmp2(2,k),cktmp2(3,k) 
  !enddo
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
   redkm%nkeq = kmesh%nkeq

   !if ((ik .ne. redkm%ktot) .and. (.not. symm%lnsymmr) ) then
      !write(*,*) 'GENREDK: the number of k-points generated by symmetry is: ',ik,redkm%ktot,redkm%kx,redkm%ky,redkm%kz
      !STOP
   !endif
   if (.not. allocated(redkm%mult))     allocate(redkm%mult(1:redkm%ktot))
   if (.not. allocated(redkm%k_weight)) allocate(redkm%k_weight(1:redkm%ktot))
   if (.not. allocated(redkm%k_coord))  allocate(redkm%k_coord(3,redkm%ktot))
   if (.not. allocated(redkm%k_id))     allocate(redkm%k_id(1:redkm%kx, 1:redkm%ky, 1:redkm%kz))
   redkm%k_coord=0.0d0
   redkm%k_id=0
! save the new coordinates into the data structure
   do i=1,redkm%ktot
      do j=1,3   
         redkm%k_coord(j,i)=cktmp2(j,i)
      enddo
   enddo
   deallocate(cktmp2)

end subroutine !GENREDK

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! GENREDOPT
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine generates the optical matrix elements
! starting from those read in on the irreducible k-mesh 
!
subroutine genredopt (redkm, symm, eirrk, eredk )
   use types
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
  
   iks = size(symm%symop_id,2) !composite index = irrk*nsym (symmorphic groups) or irrk*nsym^2 (nonsymmorphic groups)
   !since the number of bands doesn't change with the symmetry operations we may set those parameters in eredk here
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
   eredk%Mopt=0.d0
   eredk%band=0.d0
   if (symm%lcubic) then
      allocate(osm1(eredk%nbopt_min:eredk%nbopt_max,eredk%nbopt_min:eredk%nbopt_max))
      allocate(osm2(eredk%nbopt_min:eredk%nbopt_max,eredk%nbopt_min:eredk%nbopt_max))
      allocate(osm3(eredk%nbopt_min:eredk%nbopt_max,eredk%nbopt_min:eredk%nbopt_max))
   else
      allocate(Mtmp(3,3,eredk%nbopt_min:eredk%nbopt_max,eredk%nbopt_min:eredk%nbopt_max))
   endif

   ! scan the symm%symop_id to check what are the combinations of k-point and symmetry that produce non-redundant element in redkm
   k=0
   if (symm%lcubic) then
      do i=1,iks
         if(symm%symop_id(1,i) == 0) cycle 
         !this summation runs in practice over the new reducible k-points k (because only for those symop_id(1) /=0 )
         k=k+1
         ik  = symm%symop_id(2,i) ! this counter runs over the irrekpoints 
         isym= symm%symop_id(3,i)
         osm1=0.d0; osm2=0.d0; osm3=0.d0
         ! check eq. 13.16 (pg 479) in "Symmetry and Condensed Matter Physics A Computational Approach" by M. El-Batanouny, F. Wooten, CUP
         do j=1,3
            do nb=eredk%nbopt_min,eredk%nbopt_max ; do nb2=eredk%nbopt_min,eredk%nbopt_max 
               osm1(nb,nb2) = osm1(nb,nb2) + eirrk%Mopt(j,ik,nb,nb2)*symm%Msym(j,1,isym)*symm%Msym(j,1,isym)
               osm2(nb,nb2) = osm2(nb,nb2) + eirrk%Mopt(j,ik,nb,nb2)*symm%Msym(j,2,isym)*symm%Msym(j,2,isym)
               osm3(nb,nb2) = osm3(nb,nb2) + eirrk%Mopt(j,ik,nb,nb2)*symm%Msym(j,3,isym)*symm%Msym(j,3,isym)
            enddo ; enddo
         enddo
         eredk%Mopt(1,k,:,:)=osm1(:,:)
         eredk%Mopt(2,k,:,:)=osm2(:,:)
         eredk%Mopt(3,k,:,:)=osm3(:,:)
      enddo
      deallocate(osm1,osm2,osm3) 
  
      !!!!!!!!!!!!!!!!TEST
      !write(666,*)'KP',k,redkm%k_coord(1,k),redkm%k_coord(2,k),redkm%k_coord(3,k) 
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
      do i=1,iks
         if(symm%symop_id(1,i) == 0) cycle 
         !this summation runs in practice over the new reducible k-points k (because only for those symop_id(1) /=0 )
         k=k+1
         ik  = symm%symop_id(2,i) ! this counter runs over the irrekpoints 
         isym= symm%symop_id(3,i)
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
   k=0
   do i=1,iks
      if(symm%symop_id(1,i) == 0) cycle 
      !this summation runs in practice over the new reducible k-points k (because only for those symop_id(1) /=0 )
      k=k+1
      ik  = symm%symop_id(2,i) ! this counter runs over the irrekpoints 
      isym= symm%symop_id(3,i)
      do nb=1,eredk%nband_max
         eredk%band(k,nb) = eirrk%band(ik,nb)
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
subroutine trnredk (irrkm, redkm, eredk, symm)
  use types
  implicit none 
!passed variables
  type(kpointmesh) :: irrkm
  type(kpointmesh) :: redkm
  type(edisp) :: eredk
  type(symop) :: symm

!local variables
  integer :: i, j, ik, ikx, iky, ikz, ibn, ibn2 
  integer :: nk, nkx, nky, nkz
  integer :: ntmp, itest, iexist
  integer :: offdia !off-diagonal terms in the Mopt matrix?
  double precision :: dk(3), tmp1, tmp2
  double precision :: maxkx, maxky, maxkz
  double precision, allocatable :: cktmp(:,:),cktmp2(:,:) ! temporary k-point coordinates arrays
  double precision, allocatable :: bstmp(:,:),bstmp2(:,:) ! temporary bandstructure arrays
  double precision, allocatable :: Motmp(:,:,:,:),Motmp2(:,:,:,:) ! temporary optical transition matrices

  if (symm%lcubic) then
     offdia=0  !no off-diagonal terms for cubic systems   
  else
     offdia=1  !off-diagonal terms for non-cubic systems
  endif
  !allocation temporary arrays
  allocate(bstmp(redkm%ktot,eredk%nband_max)); allocate(bstmp2(redkm%ktot,eredk%nband_max))
  allocate(cktmp(3,redkm%ktot)); allocate(cktmp2(3,redkm%ktot))
  if (symm%lcubic) then 
     allocate(Motmp (3,redkm%ktot,eredk%nbopt_min:eredk%nbopt_max,eredk%nbopt_min:eredk%nbopt_max))
     allocate(Motmp2(3,redkm%ktot,eredk%nbopt_min:eredk%nbopt_max,eredk%nbopt_min:eredk%nbopt_max))
  else
     allocate(Motmp (6,redkm%ktot,eredk%nbopt_min:eredk%nbopt_max,eredk%nbopt_min:eredk%nbopt_max))
     allocate(Motmp2(6,redkm%ktot,eredk%nbopt_min:eredk%nbopt_max,eredk%nbopt_min:eredk%nbopt_max))
  endif
  cktmp = 0.0d0   
  cktmp2= 0.0d0
  bstmp = 0.0d0   
  bstmp2= 0.0d0
  Motmp = 0.0d0   
  Motmp2= 0.0d0
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
  do i = 1,3 
     do ik=1,redkm%ktot
        if ((redkm%k_coord(1,ik) > 1.0d0) .or. (redkm%k_coord(2,ik) > 1.0d0) .or.(redkm%k_coord(3,ik) > 1.0d0)) &
           STOP 'TRNREDK: something is seriously wrong here (e.g. improper traslation bigger than lattice constant)' 
        if (redkm%k_coord(i,ik)<0.0d0) then
           cktmp(i,ik)=mod(1.0d0+redkm%k_coord(i,ik),1.0d0)
           do ibn =1,eredk%nband_max
              bstmp(ik,ibn)=eredk%band(ik,ibn)
           enddo
        else 
           cktmp(i,ik)=mod(redkm%k_coord(i,ik),1.0d0)
           do ibn =1,eredk%nband_max
              bstmp(ik,ibn)=eredk%band(ik,ibn)
           enddo
        endif
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
  !  write(776,'(A,I4,3f8.4)')'KP',ik, cktmp(1,ik),cktmp(2,ik),cktmp(3,ik) 
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
  deallocate(eredk%Mopt)

  !!!!!!!!! k-point clean up !!!!!!!!!!!!!!!!
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
              bstmp2(1,ibn)  = bstmp(1,ibn) !the gamma point was left out (min(itest)=2) 
              bstmp2(ik,ibn) = bstmp(itest,ibn)
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
     STOP
  endif

  !allocate the final datastructure with the correct dimensions
  allocate(redkm%k_coord(3,nk))
  allocate(redkm%k_id(nkx,nky,nkz))
  allocate(eredk%band(nk,eredk%nband_max))
  if (symm%lcubic) then 
     allocate(eredk%Mopt(3,nk,eredk%nbopt_min:eredk%nbopt_max,eredk%nbopt_min:eredk%nbopt_max))
  else
     allocate(eredk%Mopt(6,nk,eredk%nbopt_min:eredk%nbopt_max,eredk%nbopt_min:eredk%nbopt_max))
  endif
  
  !initialise the variables
  eredk%band(:,:)    = 0.0d0
  redkm%k_coord(:,:) = 0.0d0
  redkm%k_id(:,:,:)  = 0
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
        eredk%band(ik,ibn)=bstmp2(ik,ibn)
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
  deallocate (Motmp); deallocate (Motmp2) 

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

  !!!!!!!!!!!!!!!!TEST
  !do ikx=1,nkx   
  !  do iky=1,nky   
  !    do ikz=1,nkz   
  !      ik=redkm%k_id(ikx,iky,ikz)
  !      write(777,*)'KP',ik, ikx, iky, ikz, redkm%k_coord(1,ik),redkm%k_coord(2,ik),redkm%k_coord(3,ik) 
  !      do ibn=eredk%nbopt_min,eredk%nbopt_max
  !        do ibn2=ibn,eredk%nbopt_max
  !           write(777,170)ibn,ibn2,eredk%Mopt(1,ik,ibn,ibn2),eredk%Mopt(2,ik,ibn,ibn2),eredk%Mopt(3,ik,ibn,ibn2)
  !        enddo
  !      enddo
  !    enddo
  !  enddo
  !enddo
  !170  FORMAT  (2(I4,X),3(E12.6,X))
  !!!!!!!!!!!!!!!!TEST END

end subroutine ! trnredk

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! GENFULKM   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine translates the BZ from 0 to -kmax
! to the interval [1-kmax,1] and takes care of the 
! bandstructure as well
!
subroutine genfulkm( redkm, fulkm, eredk, efulk, symm)
 use types 
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
  if (.not. allocated(fulkm%k_coord)) allocate(fulkm%k_coord(3,fulkm%ktot))
  if (.not. allocated(fulkm%k_id))    allocate(fulkm%k_id(1:fulkm%kx, 1:fulkm%ky, 1:fulkm%kz))
  !bandstructure allocation 
  if (.not. allocated(efulk%band)) allocate(efulk%band( fulkm%ktot, nband ))
  if (.not. allocated(efulk%Mopt) .and. (offdia==0)) &
    allocate(efulk%Mopt(3,fulkm%ktot,efulk%nbopt_min:efulk%nbopt_max,efulk%nbopt_min:efulk%nbopt_max))
  if (.not. allocated(efulk%Mopt) .and. (offdia==1)) &
    allocate(efulk%Mopt(6,fulkm%ktot,efulk%nbopt_min:efulk%nbopt_max,efulk%nbopt_min:efulk%nbopt_max))

  if (allocated(fulkm%k_coord)) fulkm%k_coord=0.d0
  if (allocated(fulkm%k_id)) fulkm%k_id=0
  if (allocated(efulk%band)) efulk%band(:,:)=0.d0
  if (allocated(efulk%Mopt)) efulk%Mopt(:,:,:,:)=0.d0

  do ik=1,nk
     do i=1,3
        fulkm%k_coord(i,ik)=redkm%k_coord(i,ik)
     enddo
     do ibn=1,nband
        efulk%band(ik,ibn)=eredk%band(ik,ibn)
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
           efulk%band(ik,ibn)=efulk%band(fulkm%k_id(ikx,iky,1),ibn)
           do ibn2=efulk%nbopt_min,efulk%nbopt_max
             if ((ibn>efulk%nbopt_max) .or. (ibn<efulk%nbopt_min)) cycle
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
           efulk%band(ik,ibn)=efulk%band(fulkm%k_id(ikx,1,ikz),ibn)
           do ibn2=efulk%nbopt_min,efulk%nbopt_max
             if ((ibn>efulk%nbopt_max) .or. (ibn<efulk%nbopt_min)) cycle
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
           efulk%band(ik,ibn)=efulk%band(fulkm%k_id(1,iky,ikz),ibn)
           do ibn2=efulk%nbopt_min,efulk%nbopt_max
             if ((ibn>efulk%nbopt_max) .or. (ibn<efulk%nbopt_min)) cycle
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
  !   write(778,*)'KP',ik,fulkm%k_coord(1,ik),fulkm%k_coord(2,ik),fulkm%k_coord(3,ik) 
  !   do ibn=efulk%nbopt_min,efulk%nbopt_max
  !      do ibn2=ibn,efulk%nbopt_max
  !         write(778,170)ibn,ibn2,efulk%Mopt(1,ik,ibn,ibn2),efulk%Mopt(2,ik,ibn,ibn2),efulk%Mopt(3,ik,ibn,ibn2)
  !      enddo
  !   enddo
  !enddo
  !170  FORMAT  (2(I4,X),3(E12.6,X))
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
 use types 
 use params
 implicit none
 type(algorithm)  :: algo
 type(kpointmesh) :: mesh
 type(tetramesh)  :: thdr
!local variables
 double precision :: x, y, z, xx, edmin, edmax, edgmax, edgmin
 !double precision :: x1, y1, z1, x2, y2, z2, x3, y3, z3
 double precision :: bk(3,3)           ! cartesian coordinates for the reciprocal lattice
 double precision, external :: anrm2   ! computes the distance between the tetrahedron's vertices 
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
   15 format(1x,'GENTETRA: found ',i6,' inequivalent tetrahedra from ',i8,' k-points' )

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
subroutine intetra (mesh, ek, thdr, dos)
  use types
  implicit none

   type(kpointmesh) :: mesh
   type(edisp)      :: ek 
   type(tetramesh)  :: thdr
   type(dosgrid)    :: dos
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
   allocate (dos%enrg(1:dos%nnrg), dos%dos(1:dos%nnrg),dos%nos(1:dos%nnrg))
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

      do nb=1,ek%nband_max
        if (ek%band(iq(1),nb)>99.0d0) cycle 
        if (ek%band(iq(2),nb)>99.0d0) cycle 
        if (ek%band(iq(3),nb)>99.0d0) cycle 
        if (ek%band(iq(4),nb)>99.0d0) cycle 
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
  use types
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
! INTERPTRA 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This routine linearly interpolates the optical matrix elements
! defined on the vertices of a tetrahedron, by computing a
! weighted sum of these values according to eq.6 in
! PRB (1994) 49, 16223-16233. The expressions for the 
! weights are given in app. B therein.
!
subroutine interptra (mesh, ek, thdr )
  use types
  implicit none

   type(kpointmesh) :: mesh
   type(edisp)      :: ek 
   type(tetramesh)  :: thdr
!local variables
   integer :: i, j, i00, itet, nb, nb1, iswitch
   integer :: iq(4), iqs(4) !contains the tetrahedron identifiers before and after the corners have been sorted for increasing energy
   double precision :: ef, ec(4), ec1(4), es(4)
   double precision :: c0, c1, c2, c3, cc1, cc4  !constants
   double precision :: wthdr          !weight of the tetrahedron
   double precision :: x1, x2, x3, x4 !energies zeroed w.r.t. Fermi energy: xj = es(j)-ef  
   double precision :: w(4)           !weights of the interpolation formula
   double precision, allocatable :: Mopt_tetra(:,:,:) !(x, n, n') interpolates the optical matrix element in a tetrahedron 

   !********************
   ! At the beginning this is almost a literal
   ! copy from the INTETRA routine
   !********************
! SANITY CHECK
   if (mesh%ktot<4) then
      write(*,*)'INTERPTRA: tetrahedron method fails (number of k-points < 4)',mesh%ktot
      STOP
   endif
   ef = ek%efer !local variable for the Fermi energy
   allocate(Mopt_tetra(3,ek%nbopt_min:ek%nbopt_max,ek%nbopt_min:ek%nbopt_max))
   if (.not. allocated(ek%Mopt_tetra)) &
     &  allocate(ek%Mopt_tetra(3, 1:thdr%ntet, ek%nbopt_min:ek%nbopt_max, ek%nbopt_min:ek%nbopt_max))
! loop over tetrahedra:
   do itet=1,thdr%ntet
      Mopt_tetra=0.d0
! get the four corner points:
      iq(1) = thdr%idtet(1,itet)
      iq(2) = thdr%idtet(2,itet)
      iq(3) = thdr%idtet(3,itet)
      iq(4) = thdr%idtet(4,itet)
      wthdr = thdr%vltet(itet)

      do nb=1,ek%nband_max
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
!!!!!!!!!!!test
         !write(*,*) 'tetrahedron no.',itet
         !do i=1,4
         !   write(*,'(F8.3, 2I6)') ec(i),iq(i)
         !enddo
         !do i=1,4
         !   write(*,'(F8.3, 2I6)') es(i),iqs(i)
         !enddo
         !STOP 
!!!!!!!!!!!test end

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
            c1 = c0*x1*x1*x1/((es(1)-es(2))*(es(3)-es(1))*(es(4)-es(1)))
            cc1= (1.d0/(es(2)-es(1))) + (1.d0/(es(3)-es(1))) + (1.d0/(es(4)-es(1)))

         else if ((es(2) < ef) .and. (es(3) > ef)) then
            iswitch=3 
            c1 = c0*x1*x1/((es(4)-es(1))*(es(3)-es(1)))
            c2 = c0*x1*x2*x3/((es(4)-es(1))*(es(3)-es(2))*(es(3)-es(1)))
            c3 = c0*x2*x2*x4/((es(4)-es(2))*(es(3)-es(2))*(es(4)-es(1))) 

         else if ((es(3) < ef) .and. (es(4) > ef)) then
            iswitch=4 
            c1 = c0*x4*x4*x4/((es(4)-es(1))*(es(4)-es(2))*(es(4)-es(3)))
            cc4= (1.d0/(es(1)-es(4))) + (1.d0/(es(2)-es(4))) + (1.d0/(es(3)-es(4)))

         else if (es(4) < ef) then
            iswitch=5 
         else 
            write(*,*)'INTERPTRA: the ordering of your thetrahedron vertices is not consistent', itet
            !STOP
         endif

         select case (iswitch)
            case(1)
               w(1) = 0.d0
               w(2) = 0.d0
               w(3) = 0.d0
               w(4) = 0.d0

            case(2)
               w(1) = c1*(4+(cc1*x1))
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
               w(4) = c0 - (c1*(4+(cc4*x4))) 

            case(5)
               w(1) = c0 
               w(2) = c0 
               w(3) = c0 
               w(4) = c0 

         end select
 ! now that the weights are set need to perform the integration within the tetrahedron
         do nb1=nb,ek%nbopt_max
            do j=1,3
               do i=1,4
                  Mopt_tetra(j,nb,nb1)=Mopt_tetra(j,nb,nb1) + (ek%Mopt(j,iqs(i),nb,nb1 )*w(i))
               enddo
            enddo 
         enddo !over nb1
      enddo    ! over nb band
      ek%Mopt_tetra(1,itet,:,: )=Mopt_tetra(1,:,:)
      ek%Mopt_tetra(2,itet,:,: )=Mopt_tetra(2,:,:)
      ek%Mopt_tetra(3,itet,:,: )=Mopt_tetra(3,:,:)
   enddo       ! over tetrahedra

   deallocate(Mopt_tetra)

 return
end subroutine !INTERPTRA 


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! INTERPTRA_mu 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This routine linearly interpolates the product 
! of the optical matrix elemets with the transport kernel
! defined on the vertices of a tetrahedron, by computing a
! weighted sum of these values according to eq.6 in
! PRB (1994) 49, 16223-16233. The expressions for the 
! weights are given in app. B therein.
!
subroutine interptra_mu (iT, itet, mu, nalpha, lBoltz, mesh, ek, thdr, sct, resp, hpresp )
  use types
  use response
  use params

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


 !!!!!!!!!!!!!!!!test
  ! select type(resp)
  !    type is(dp_resp)
  !    write (*,*) 'dp_resp type'
  !    type is(dp_respinter)
  !    write (*,*) 'dp_respinter type'
  ! end select
  ! STOP
 !!!!!!!!!!!!!!!!test end

   !********************
   ! At the beginning this is almost a literal
   ! copy from the INTETRA routine
   !********************
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
      ! set the value of the scattering rate for the specific temperature, band
      resp%gamma=real(sct%gam(iT,nb)*sct%z,8)
      if (present(hpresp)) hpresp%gamma=real(sct%gam(iT,nb)*sct%z,16)

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
         write(*,*)'INTERPTRA_MU: the ordering of your thetrahedron vertices is not consistent', itet
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
                                         (real(w(i),16)*hpresp%aB_tmp(i,nb,ix,iy))*betaQ/(hpresp%gamma**2)
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
end subroutine !INTERPTRA_mu 


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!SUBROUTINE GENDOSEL
! Generates the Density of States starting 
! from the Wien2k eigenvalues  
! by replacing the Dirac's delta with Lorentzians
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine gendosel(mesh, ek, dos)
  use types
  use params
  implicit none

   type(kpointmesh) :: mesh
   type(edisp)      :: ek 
   type(dosgrid)    :: dos
!local variables
   integer :: i, ik, nb !energy, k-point, band counters
   double precision :: br, de !broadening, energy spacing

   allocate (dos%enrg(1:dos%nnrg),dos%dos(1:dos%nnrg),dos%nos(1:dos%nnrg))
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
      do ik=1,mesh%ktot
         do nb=1,ek%nband_max
            if (ek%band(ik,nb) > 99.0d0) cycle !necessary because big eig'vals 
                                               !have been introduced to make matrices square 
            dos%dos(i)=dos%dos(i)+((br/pi)*(1.0d0/(((dos%enrg(i)-ek%band(ik,nb))**2)+(br**2))))/mesh%ktot 
            dos%nos(i)=dos%nos(i)+(0.5d0 + ((1.0d0/pi)*atan((dos%enrg(i)-ek%band(ik,nb))/br)))/mesh%ktot   
         enddo
      enddo
   enddo

end subroutine !GENDOSEL

end module estruct

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  double precision function ek_sc(k,iband,eirrk)
   use types
   use params
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
   use types
   use params
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
   use types
   use params
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

 
