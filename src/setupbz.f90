program setupbz
  use Mparams
  use Mtypes
  use Mestruct
  use hdf5_wrapper
  use hdf5
  implicit none

  type(algorithm)          :: algo
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

  call read_config(algo,  irrkm, eirrk, sct)
  call setup_meshes(algo, irrkm, redkm, fulkm, eirrk, eredk, efulk)
  call estruct_init(algo, irrkm, redkm, fulkm, eirrk, eredk, efulk, thdr, dos, sct)

  ! if we work with tetrahedrons the reducible BZ got extended to full
  ! i.e. we have some double counting which was necessary for the
  ! construction of the tetrahedrons
  if (algo%ltetra) then
     kpointer => fulkm
     epointer => efulk
  else
     kpointer => redkm
     epointer => eredk
  endif

  call hdf5_init()
  call hdf5_create_file('test.hdf5')
  call hdf5_open_file('test.hdf5', ifile)

  call hdf5_write_data(ifile, '/.kmesh', kpointer%k_coord)

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

  call hdf5_write_data(ifile, '/.optical_band_min', epointer%nbopt_min)
  call hdf5_write_data(ifile, '/.optical_band_max', epointer%nbopt_max)
  call hdf5_write_data(ifile, '/.mu_lda', epointer%efer)

  !also save all the tetrahedron information
  if (algo%ltetra) then
     call hdf5_write_data(ifile, '/.thdr', thdr%idtet)
     call hdf5_write_data(ifile, '/.thdr_vol', thdr%vltet)
  endif

  call hdf5_close_file(ifile)
  call hdf5_finalize()

end program setupbz
