program setupbz
#ifdef iso_fortran_env
   use, intrinsic :: iso_fortran_env, only : stdin =>input_unit, &
                                             stdout=>output_unit, &
                                             stderr=>error_unit
#else
#define stdin  5
#define stdout 6
#define stderr 0
#endif

  use Mparams
  use Maux
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
  integer            :: filesize
  integer(hid_t)     :: ifile
  character(len=256) :: outfile
  character(len=256) :: string

  integer            :: er
  character(len=256) :: erstr

  call greeting(stdout)

  ! Read the Input
  call read_config(kmesh, edisp, sct, outfile, er, erstr)
  if (er /= 0) call stop_with_message(stderr, erstr, er)

  ! Check for existance of files
  call check_config(er,erstr)
  if (er /= 0) call stop_with_message(stderr, erstr, er)

  ! Create the energy dispersion, kmesh, etc
  call estruct_init(kmesh, edisp, thdr, dos, sct)
  write(*,*) 'SETUPBZ: writing processed data to: ', adjustl(trim(outfile))

  ! file setup
  call hdf5_init()
  call hdf5_create_file(trim(adjustl(outfile)))
  call hdf5_open_file(trim(adjustl(outfile)), ifile)

  ! call hdf5_write_data(ifile, '/.kmesh/k_coord', kmesh%k_coord)
  call hdf5_write_data(ifile, '/.kmesh/kx',           kmesh%kx)
  call hdf5_write_data(ifile, '/.kmesh/ky',           kmesh%ky)
  call hdf5_write_data(ifile, '/.kmesh/kz',           kmesh%kz)
  call hdf5_write_data(ifile, '/.kmesh/ktot',         kmesh%ktot)
  call hdf5_write_data(ifile, '/.kmesh/kred',         kmesh%kred)
  call hdf5_write_data(ifile, '/.kmesh/multiplicity', kmesh%multiplicity)
  call hdf5_write_data(ifile, '/.kmesh/weight',       kmesh%weight)

  ! symmetry information
  ! call hdf5_write_data(ifile, '/.symmetry/knsym',     symm%knsym)
  ! call hdf5_write_data(ifile, '/.symmetry/rotations', symm%Msym_reciprocal)
  ! call hdf5_write_data(ifile, '/.symmetry/mapping',   symm%symop_id)

  call hdf5_write_data(ifile, '/.edisp/nelect', edisp%nelect)

  ! lattice information
  call hdf5_write_data(ifile, '/.crystal/a',         lat%a)
  call hdf5_write_data(ifile, '/.crystal/vol',       lat%vol)
  call hdf5_write_data(ifile, '/.crystal/nalpha',    lat%nalpha)
  if (lat%lortho) then
     call hdf5_write_data(ifile, '/.crystal/lortho', 1)
  else
     call hdf5_write_data(ifile, '/.crystal/lortho', 0)
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
  ! if (algo%ltetra) then
  !    call hdf5_write_data(ifile, '/.tetrahedrons/ntet',     thdr%ntet)
  !    call hdf5_write_data(ifile, '/.tetrahedrons/thdr_id',  thdr%idtet)
  !    call hdf5_write_data(ifile, '/.tetrahedrons/thdr_vol', thdr%vltet)
  ! endif

  call hdf5_close_file(ifile)
  call hdf5_finalize()

  inquire(file=trim(adjustl(outfile)), SIZE=filesize)
  if (filesize /= -1) then
    write(*,*) "Preprocessed file: ", trim(adjustl(outfile)), " successfully created."
    write(*,"(A,F8.2,A)") "File size: ", dble(filesize)/1000000, " MB."
  else
    write(*,*) "Error in creation of preprocessed file."
    stop
  endif

end program setupbz
