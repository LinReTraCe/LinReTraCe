program setupbz
  use Mparams
  use Mtypes
  use Mestruct
  use hdf5_wrapper
  use hdf5
  implicit none

  type(algorithm)  :: algo
  type(kpointmesh) :: kmesh
  type(kpointmesh) :: redkm
  type(kpointmesh) :: fulkm
  type(edisp)      :: eirrk
  type(edisp)      :: eredk
  type(edisp)      :: efulk
  type(tetramesh)  :: thdr
  type(dosgrid)    :: dos
  type(scatrate)   :: sct

  integer :: i,j,k
  integer :: kk
  integer(hid_t) :: ifile
  character(len=100) :: string

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

  call read_config(algo, eirrk, kmesh, sct)
  call estruct_init(algo, kmesh, redkm, fulkm, eirrk, eredk, efulk, thdr, dos, sct)

  call hdf5_init()
  call hdf5_create_file('test.hdf5')
  call hdf5_open_file('test.hdf5', ifile)

  ! here we take the kmesh variable
  if (algo%ltbind) then
     call hdf5_write_data(ifile, '/.kmesh', kmesh%k_coord)

     do i=1,size(eirrk%band, 1)
        write(string,"(I5.5, '/energies')") i
        call hdf5_write_data(ifile, trim(string), eirrk%band(i,:))
     enddo

     ! ! Z values are not defined in the construction
     ! ! of the tight-binding case
     ! do i=1,size(eirrk%Z, 1)
     !   write(string,"(I5.5, '/zqp', I5.5)") i
     !   call hdf5_write_data(ifile, trim(string), eirrk%Z(i,:))
     ! enddo

     ! M(x,k,n',n)= | <n',k|p.e_x|n,k> |^2
     do i=1,size(eirrk%Mopt, 2)
        write(string, "(I5.5, '/optical/')") i
        call hdf5_write_data(ifile, trim(string), eirrk%Mopt(:,i,:,:))
     enddo

     call hdf5_write_data(ifile, '/.optical_band_min', eirrk%nbopt_min)
     call hdf5_write_data(ifile, '/.optical_band_max', eirrk%nbopt_max)
  ! here we simply take the eredk variable
  else
     call hdf5_write_data(ifile, '/.kmesh', redkm%k_coord)
     do i=1,size(eredk%band, 1)
        write(string,"(I5.5, '/energies')") i
        call hdf5_write_data(ifile, trim(string), eredk%band(i,:))
     enddo

     do i=1,size(eredk%Z, 1)
        write(string,"(I5.5, '/zqp', I5.5)") i
        call hdf5_write_data(ifile, trim(string), eredk%Z(i,:))
     enddo

     ! M(x,k,n',n)= | <n',k|p.e_x|n,k> |^2
     do i=1,size(eredk%Mopt, 2)
        write(string, "(I5.5, '/optical/')") i
        call hdf5_write_data(ifile, trim(string), eredk%Mopt(:,i,:,:))
     enddo

     call hdf5_write_data(ifile, '/.optical_band_min', eredk%nbopt_min)
     call hdf5_write_data(ifile, '/.optical_band_max', eredk%nbopt_max)
  endif


  call hdf5_close_file(ifile)
  call hdf5_finalize()

end program setupbz
