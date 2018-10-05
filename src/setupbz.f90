program setupbz
  use Mparams
  use Mtypes
  use Mestruct
  use hdf5_wrapper
  use hdf5
  implicit none

  type(kpointmesh)         :: irrkm
  type(kpointmesh), target :: redkm
  type(kpointmesh), target :: fulkm
  type(edisp)              :: eirrk
  type(edisp), target      :: eredk
  type(edisp), target      :: efulk
  type(tetramesh)          :: thdr
  type(dosgrid)            :: dos
  type(scatrate)           :: sct

  integer            :: i,j,k
  integer            :: kk
  integer(hid_t)     :: ifile
  character(len=100) :: string

  class(kpointmesh), pointer :: kpointer
  class(edisp), pointer      :: epointer

  write(*,*)
  write(*,*)'#####################################################'
  write(*,*)'#  Lin-ReTraCe -- Linear Response Transport Centre  #'
  write(*,*)'#####################################################'
  write(*,*)'#  Preprocessing irreducible Brillouin Zone         #'
  write(*,*)'#####################################################'
  write(*,*)'#  J.M. Tomczak, E. Maggio, M. Pickem               #'
  write(*,*)'#####################################################'
  write(*,*)

  algo%ldebug = .true.
  algo%lgenred = .false.
  algo%lpreproc = .true.

  call read_config(irrkm, eirrk, sct)
  call setup_algo(irrkm, redkm, fulkm, eirrk, eredk, efulk)
  call estruct_init(irrkm, redkm, fulkm, eirrk, eredk, efulk, thdr, dos, sct)

  if (algo%ltetra) then
     STOP 'does not work for tetrahedrons'
     ! kpointer => fulkm
     ! epointer => efulk
  else
     kpointer => redkm
     epointer => eredk
  endif

  ! file setup
  call hdf5_init()
  call hdf5_create_file('test.hdf5')
  call hdf5_open_file('test.hdf5', ifile)

  call hdf5_write_data(ifile, '/.kmesh/k_coord', kpointer%k_coord)
  call hdf5_write_data(ifile, '/.kmesh/kx', kpointer%kx)
  call hdf5_write_data(ifile, '/.kmesh/ky', kpointer%ky)
  call hdf5_write_data(ifile, '/.kmesh/kz', kpointer%kz)
  call hdf5_write_data(ifile, '/.kmesh/ktot', irrkm%ktot)

  ! symmetry information
  call hdf5_write_data(ifile, '/.symmetry/nsym', symm%nsym)
  call hdf5_write_data(ifile, '/.symmetry/rotations', symm%Msym)
  call hdf5_write_data(ifile, '/.symmetry/mapping', symm%symop_id)

  ! lattice information
  call hdf5_write_data(ifile, '/.crystal/alat', lat%alat)
  call hdf5_write_data(ifile, '/.crystal/a', lat%a)
  call hdf5_write_data(ifile, '/.crystal/vol', lat%vol)
  call hdf5_write_data(ifile, '/.crystal/nalpha', lat%nalpha)
  if (lat%lcubic) then
     call hdf5_write_data(ifile, '/.crystal/lcubic', 1)
  else
     call hdf5_write_data(ifile, '/.crystal/lcubic', 0)
  endif

  do i=1,size(epointer%band, 1)
     write(string,"('/kpoint/',I6.6, '/energies')") i
     call hdf5_write_data(ifile, trim(string), epointer%band(i,:))
  enddo

  do i=1,size(epointer%Z, 1)
     write(string,"('/kpoint/',I6.6, '/zqp')") i
     call hdf5_write_data(ifile, trim(string), epointer%Z(i,:))
  enddo

  ! M(x,k,n',n)= | <n',k|p.e_x|n,k> |^2
  do i=1,size(epointer%Mopt, 2)
     write(string, "('/kpoint/',I6.6, '/optical/')") i
     call hdf5_write_data(ifile, trim(string), epointer%Mopt(:,i,:,:))
  enddo

  call hdf5_write_data(ifile, '/.bands/band_max',  epointer%nband_max)
  call hdf5_write_data(ifile, '/.bands/optical_band_min', epointer%nbopt_min)
  call hdf5_write_data(ifile, '/.bands/optical_band_max', epointer%nbopt_max)
  call hdf5_write_data(ifile, '/.mu_lda', epointer%efer)

  !also save all the tetrahedron information
  if (algo%ltetra) then
     call hdf5_write_data(ifile, '/.thdr_id', thdr%idtet)
     call hdf5_write_data(ifile, '/.thdr_vol', thdr%vltet)
  endif

  call hdf5_close_file(ifile)
  call hdf5_finalize()

end program setupbz
