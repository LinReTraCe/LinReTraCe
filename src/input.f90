module Minput
  use Mparams
  use Mtypes
  use hdf5
  use hdf5_wrapper
  implicit none

  contains

  subroutine read_config(kmesh, edisp, sct)
    implicit none
    type(kpointmesh) :: kmesh  ! contains k-point mesh specifiers and logical switches on how to get the mesh from
    type(energydisp) :: edisp     ! contains the band dispersion energy and the optical matrix elements (when which > 2) along the irr-k-mesh
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
    read(10,*) edisp%nband_max, tmax
    read(10,*) lat%alat, lat%a(:) ! ALAT in A, a(1:3) in ALAT
    read(10,*) edisp%nelect              ! read the number of electrons, the fermi level will be computed, consistently with the LINRETRACE routines
    read(10,*) lat%vol

    ! at this point one should choose in the input file what to do:
    ! case = 1 use the uniform k-mesh generated with tight binding
    ! case = 2 read the k-points and band dispersion from Wien2k
    ! case = 3 read also the optical matrix elements from Wien2k
    ! case = 4 read the renormalised qp energies, renormalisation factor and -Im{selfnrg}
    read(10,*) which, which_tetra, which_grid

    algo%ltbind  = .false.
    algo%lw2k    = .false.
    algo%loptic  = .false.
    algo%lBfield = .false.
    algo%ldmft   = .false.

    select case(which)
       case(1)
          algo%ltbind =.true.
       case(2)
          algo%lw2k   =.true.
       case(3)
          algo%loptic =.true.
       case(4)
          algo%loptic =.true.
          algo%ldmft  =.true.
       case default
          write(*,*)'Not a valid algorithm selection'
          STOP
    end select

    select case(which_tetra)
       case(0)
          algo%ltetra=.false.
          write(*,*)'Tetrahedron method not selected'
       case(1)
          algo%ltetra=.true.
          write(*,*)'Tetrahedron method selected'
       case default
          write(*,*)'Integration method does not exist',which_tetra
          STOP
    end select

    select case(which_grid)
       case(0)
          algo%lsymm=.false.
          write(*,*)'Symmetry switched OFF, REDucible BZ required'
       case(1)
          algo%lsymm=.true.
          write(*,*)'Symmetry switched ON, IRReducible BZ required'
       case default
          write(*,*)'You must either provide reducible BZ or switch on symmetry',which_grid
          STOP
    end select

    if (algo%ltbind) algo%lsymm = .false.

    ! this is unfortunately necessary if we want arbitrary paths
    ! this whole section will be replaced anyways ...
    read(10,'(A)') ctmp
    itmp = index(ctmp,'#')
    if (itmp .gt. 0) then
       algo%mysyst = trim(adjustl(ctmp(:itmp-1)))
    else
       algo%mysyst = trim(adjustl(ctmp))
    endif

    read(10,*) kmesh%kx,kmesh%ky,kmesh%kz

    if (algo%ltbind) then
       write(*,*)'READ_CONFIG: reading TB parameters'
       allocate(edisp%E0(edisp%nband_max),edisp%t(edisp%nband_max, tmax))
       do iband=1,edisp%nband_max
          read(10,*)edisp%E0(iband),edisp%t(iband,:)
       enddo
       read(10,*)algo%lBfield
    endif
    !now read the temperature variables
    read(10,*)sct%Tmin,sct%Tmax,sct%dT
    read(10,*)algo%imurestart, edisp%efer
    !at this point I read in also the coefficients of gamma
    read(10,*)sct%ng
    allocate(sct%gc(0:sct%ng))
    read(10,*)sct%gc(:)
    !read in also the renormalisation factor and the quasi-particle shift
    read(10,*)edisp%ztmp
    close(10)


    ! fill out the rest of the datatype
    kmesh%kred = kmesh%kx*kmesh%ky*kmesh%kz
    kmesh%kful = (kmesh%kx+1)*(kmesh%ky+1)*(kmesh%kz+1)

  end subroutine


  subroutine read_preproc_data(fname, kmesh, edisp, thdr, dos)
     implicit none
     character(len=*), intent(in) :: fname
     type(kpointmesh) :: kmesh
     type(energydisp) :: edisp
     type(tetramesh)  :: thdr
     type(dosgrid)    :: dos

     integer(hid_t)       :: ifile
     integer              :: loccubic, kpoints, i
     character(len=6)   :: nmbstring
     real(8), allocatable :: rank1arr(:)
     real(8), allocatable :: rank3arr(:,:,:)

     call hdf5_init()
     call hdf5_open_file(trim(adjustl(fname)), ifile, rdonly=.true.)

     ! mesh
     call hdf5_read_data(ifile, "/.kmesh/k_coord", kmesh%k_coord)
     call hdf5_read_data(ifile, "/.kmesh/kx",      kmesh%kx)
     call hdf5_read_data(ifile, "/.kmesh/ky",      kmesh%ky)
     call hdf5_read_data(ifile, "/.kmesh/kz",      kmesh%kz)
     call hdf5_read_data(ifile, "/.kmesh/ktot",    kmesh%ktot) ! number of k-points from w2k
     call hdf5_read_data(ifile, "/.kmesh/kred",    kmesh%kred)
     call hdf5_read_data(ifile, "/.kmesh/kful",    kmesh%kful)

     ! symmetry
     call hdf5_read_data(ifile, "/.symmetry/nsym",      symm%knsym)
     call hdf5_read_data(ifile, "/.symmetry/rotations", symm%Msym_reciprocal)
     call hdf5_read_data(ifile, "/.symmetry/mapping",   symm%symop_id)

     ! crystal
     call hdf5_read_data(ifile, "/.crystal/alat",   lat%alat)
     call hdf5_read_data(ifile, "/.crystal/a",      rank1arr)
     call hdf5_read_data(ifile, "/.crystal/vol",    lat%vol)
     call hdf5_read_data(ifile, "/.crystal/nalpha", lat%nalpha)
     lat%a = rank1arr
     deallocate(rank1arr)

     ! bands
     call hdf5_read_data(ifile, "/.bands/band_max",         edisp%nband_max)
     call hdf5_read_data(ifile, "/.bands/optical_band_min", edisp%nbopt_min)
     call hdf5_read_data(ifile, "/.bands/optical_band_max", edisp%nbopt_max)

     ! DOS
     call hdf5_read_data(ifile, '/.dos/grid', dos%enrg)
     call hdf5_read_data(ifile, '/.dos/dos',  dos%dos)
     call hdf5_read_data(ifile, '/.dos/nos',  dos%nos)
     dos%nnrg = size(dos%enrg, 1)
     dos%emin = dos%enrg(1)
     dos%emax = dos%enrg(dos%nnrg)

     call hdf5_read_data(ifile, '/.dos/mu', edisp%efer)

     ! number of saved k-points
     kpoints = hdf5_get_number_groups(ifile, "/kpoint")
     write(*,*) "MAIN: found ", kpoints, " kpoints in the preprocessed file."

     ! k-point information
     allocate(edisp%band(kpoints, edisp%nband_max))
     allocate(edisp%Z(kpoints, edisp%nband_max))
     if (lat%lcubic) then
        allocate(edisp%Mopt(3, kpoints, edisp%nbopt_min:edisp%nbopt_max, edisp%nbopt_min:edisp%nbopt_max))
     else
        allocate(edisp%Mopt(6, kpoints, edisp%nbopt_min:edisp%nbopt_max, edisp%nbopt_min:edisp%nbopt_max))
     endif
     do i=1, kpoints
        write(nmbstring,'(I6.6)') i
        call hdf5_read_data(ifile, "/kpoint/"//nmbstring//"/energies", rank1arr)
        edisp%band(i,:)     = rank1arr
        deallocate(rank1arr)
        call hdf5_read_data(ifile, "/kpoint/"//nmbstring//"/zqp",      rank1arr)
        edisp%Z(i,:)        = rank1arr
        deallocate(rank1arr)
        call hdf5_read_data(ifile, "/kpoint/"//nmbstring//"/optical",  rank3arr)
        edisp%Mopt(:,i,:,:) = rank3arr
        deallocate(rank3arr)
     enddo

     if (algo%ltetra) then
        call hdf5_read_data(ifile, '/.tetrahedrons/ntet',     thdr%ntet)
        call hdf5_read_data(ifile, '/.tetrahedrons/thdr_id',  thdr%idtet)
        call hdf5_read_data(ifile, '/.tetrahedrons/thdr_vol', thdr%vltet)
     endif

  end subroutine

end module Minput
