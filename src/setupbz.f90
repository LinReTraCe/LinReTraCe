program setupbz
  use Mparams
  use Mtypes
  use Mestruct
  use Mlookup
  use Mconfig
  use hdf5_wrapper
  use hdf5
  implicit none

  type(kpointmesh) :: kmesh
  type(energydisp) :: edisp
  type(tetramesh)  :: thdr
  type(dosgrid)    :: dos
  type(scatrate)   :: sct

  integer            :: i
  integer(hid_t)     :: ifile
  character(len=150) :: outfile
  character(len=150) :: string

  integer            :: er
  character(len=150) :: erstr

  write(*,*)
  write(*,*)'#####################################################'
  write(*,*)'#  Lin-ReTraCe -- Linear Response Transport Centre  #'
  write(*,*)'#####################################################'
  write(*,*)'#  Preprocessing Band structure data                #'
  write(*,*)'#####################################################'
  write(*,*)'#  J.M. Tomczak, E. Maggio, M. Pickem               #'
  write(*,*)'#####################################################'
  write(*,*)

  algo%ldebug = .true.
  algo%lgenred = .false.

  call read_config(kmesh, edisp, sct, outfile, er, erstr)
  if (er /= 0) then
     write(*,*) erstr
     stop
  endif
  write(*,*) 'SETUPBZ: writing processed data to: ', adjustl(trim(outfile))
  call init_config(kmesh)
  call check_config(er,erstr)
  if (er /= 0) then
     write(*,*) erstr
     stop
  endif
  call estruct_init(kmesh, edisp, thdr, dos, sct)

  ! file setup
  call hdf5_init()
  call hdf5_create_file(trim(adjustl(outfile)))
  call hdf5_open_file(trim(adjustl(outfile)), ifile)

  call hdf5_write_data(ifile, '/.kmesh/k_coord', kmesh%k_coord)
  call hdf5_write_data(ifile, '/.kmesh/kx',      kmesh%kx)
  call hdf5_write_data(ifile, '/.kmesh/ky',      kmesh%ky)
  call hdf5_write_data(ifile, '/.kmesh/kz',      kmesh%kz)
  call hdf5_write_data(ifile, '/.kmesh/ktot',    kmesh%ktot)
  call hdf5_write_data(ifile, '/.kmesh/kred',    kmesh%kred)
  call hdf5_write_data(ifile, '/.kmesh/kful',    kmesh%kful)

  ! symmetry information
  call hdf5_write_data(ifile, '/.symmetry/nsym',      symm%nsym)
  call hdf5_write_data(ifile, '/.symmetry/rotations', symm%Msym)
  call hdf5_write_data(ifile, '/.symmetry/mapping',   symm%symop_id)

  ! lattice information
  call hdf5_write_data(ifile, '/.crystal/alat',      lat%alat)
  call hdf5_write_data(ifile, '/.crystal/a',         lat%a)
  call hdf5_write_data(ifile, '/.crystal/vol',       lat%vol)
  call hdf5_write_data(ifile, '/.crystal/nalpha',    lat%nalpha)
  if (lat%lcubic) then
     call hdf5_write_data(ifile, '/.crystal/lcubic', 1)
  else
     call hdf5_write_data(ifile, '/.crystal/lcubic', 0)
  endif

  ! band information
  call hdf5_write_data(ifile, '/.bands/band_max',         edisp%nband_max)
  call hdf5_write_data(ifile, '/.bands/optical_band_min', edisp%nbopt_min)
  call hdf5_write_data(ifile, '/.bands/optical_band_max', edisp%nbopt_max)

  ! dispersion
  do i=1,size(edisp%band, 1)
     write(string,"('/kpoint/',I6.6, '/energies')") i
     call hdf5_write_data(ifile, trim(string), edisp%band(i,:))
  enddo

  ! quasiparticle weights
  do i=1,size(edisp%Z, 1)
     write(string,"('/kpoint/',I6.6, '/zqp')") i
     call hdf5_write_data(ifile, trim(string), edisp%Z(i,:))
  enddo

  ! DOS
  call hdf5_write_data(ifile, '/.dos/grid', dos%enrg)
  call hdf5_write_data(ifile, '/.dos/dos',  dos%dos)
  call hdf5_write_data(ifile, '/.dos/nos',  dos%nos)
  call hdf5_write_data(ifile, '/.dos/mu',   edisp%efer)

  ! optical elements
  ! M(x,k,n',n)= | <n',k|p.e_x|n,k> |^2
  do i=1,size(edisp%Mopt, 2)
     write(string, "('/kpoint/',I6.6, '/optical/')") i
     call hdf5_write_data(ifile, trim(string), edisp%Mopt(:,i,:,:))
  enddo

  !also save all the tetrahedron information
  if (algo%ltetra) then
     call hdf5_write_data(ifile, '/.tetrahedrons/ntet',     thdr%ntet)
     call hdf5_write_data(ifile, '/.tetrahedrons/thdr_id',  thdr%idtet)
     call hdf5_write_data(ifile, '/.tetrahedrons/thdr_vol', thdr%vltet)
  endif

  call hdf5_close_file(ifile)
  call hdf5_finalize()

end program setupbz
