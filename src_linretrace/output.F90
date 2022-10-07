module Moutput
  use Mmpi_org
  use Mtypes
  use Mparams
  use hdf5
  use hdf5_wrapper

  contains

! output all auxiliary data to the LRTC HDF5 output file
! includes, configuration details, input structure details, scattering rate details
subroutine output_auxiliary(algo, info, pot, temp, kmesh, edisp, sct, imp)
  implicit none
  type(algorithm)   :: algo
  type(runinfo)     :: info
  type(potential)   :: pot
  type(temperature) :: temp
  type(kpointmesh)  :: kmesh
  type(energydisp)  :: edisp
  type(scattering)  :: sct
  type(impurity)    :: imp

  character(len=256)   :: string
  integer(hid_t)       :: ifile
  integer              :: iimp
  real(8), allocatable :: weights(:)

  call hdf5_open_file(algo%output_file, ifile)
  call hdf5_write_attribute(ifile, '/', 'identifier', 'LRTCoutput')

  call hdf5_create_group(ifile, '.config')
  call hdf5_write_attribute(ifile, '.config', 'tmode', algo%lTMODE)
  call hdf5_write_attribute(ifile, '.config', 'mumode',algo%lMUMODE)
  call hdf5_write_attribute(ifile, '.config', 'doping',algo%ldoping)
  call hdf5_write_attribute(ifile, '.config', 'quad',algo%lQuad)
  call hdf5_write_attribute(ifile, '.config', 'debug', algo%lDebug)
  call hdf5_write_attribute(ifile, '.config', 'bfield', algo%lBfield)
  call hdf5_write_attribute(ifile, '.config', 'rootmethod', algo%rootMethod)
  call hdf5_write_attribute(ifile, '.config', 'musearch', algo%muSearch)
  call hdf5_write_attribute(ifile, '.config', 'mufermi', algo%muFermi)
  call hdf5_write_attribute(ifile, '.config', 'oldmuhdf5', algo%lOldmuHdf5)
  call hdf5_write_attribute(ifile, '.config', 'oldmutext', algo%lOldmuText)
  call hdf5_write_attribute(ifile, '.config', 'scatteringfile', algo%lScatteringFile)
  call hdf5_write_attribute(ifile, '.config', 'scatteringtext', algo%lScatteringText)
  call hdf5_write_attribute(ifile, '.config', 'interbandquantities', algo%lInterBandQuantities)
  call hdf5_write_attribute(ifile, '.config', 'intrabandquantities', algo%lIntraBandQuantities)
  call hdf5_write_attribute(ifile, '.config', 'fulloutput', algo%fullOutput)
  ! call hdf5_write_attribute(ifile, '.config', 'energyOutput', algo%lEnergyOutput)
  call hdf5_write_attribute(ifile, '.config', 'boltzmann', algo%lBoltzmann)
  call hdf5_write_attribute(ifile, '.config', 'scissors', algo%lScissors)
  call hdf5_write_attribute(ifile, '.config', 'steps', algo%steps)
  call hdf5_write_attribute(ifile, '.config', 'step_dir', algo%step_dir)
  call hdf5_write_attribute(ifile, '.config', 'impurities', algo%lImpurities)
  call hdf5_write_attribute(ifile, '.config', 'input_energies', trim(algo%input_energies))
  if (len(trim(algo%input_scattering_hdf5)) == 0) then
    call hdf5_write_attribute(ifile, '.config', 'input_scattering_hdf5', "-")
  else
    call hdf5_write_attribute(ifile, '.config', 'input_scattering_hdf5', trim(algo%input_scattering_hdf5))
  endif
  if (len(trim(algo%input_scattering_text)) == 0) then
    call hdf5_write_attribute(ifile, '.config', 'input_scattering_text', "-")
  else
    call hdf5_write_attribute(ifile, '.config', 'input_scattering_text', trim(algo%input_scattering_text))
  endif
  if (len(trim(algo%input_mu_hdf5)) == 0) then
    call hdf5_write_attribute(ifile, '.config', 'input_mu_hdf5', "-")
  else
    call hdf5_write_attribute(ifile, '.config', 'input_mu_hdf5', trim(algo%input_mu_hdf5))
  endif
  if (len(trim(algo%input_mu_text)) == 0) then
    call hdf5_write_attribute(ifile, '.config', 'input_mu_text', "-")
  else
    call hdf5_write_attribute(ifile, '.config', 'input_mu_text', trim(algo%input_mu_text))
  endif
  if (len(trim(algo%dbgstr)) == 0) then
    call hdf5_write_attribute(ifile, '.config', 'dbgstr', "-")
  else
    call hdf5_write_attribute(ifile, '.config', 'dbgstr', trim(algo%dbgstr))
  endif

  call hdf5_create_group(ifile, '.scattering')
  if (allocated(sct%gamcoeff)) then
    call hdf5_write_data(ifile, '.scattering/gamcoeff', sct%gamcoeff)
  endif
  if (allocated(sct%zqpcoeff)) then
    call hdf5_write_data(ifile, '.scattering/zqpcoeff', sct%zqpcoeff)
  endif
  if (allocated(sct%gamtext)) then
    call hdf5_write_data(ifile, '.scattering/gamtext', sct%gamtext)
  endif
  if (allocated(sct%zqptext)) then
    call hdf5_write_data(ifile, '.scattering/zqptext', sct%zqptext)
  endif
  if (algo%lScatteringFile .or. algo%lScatteringText) then
    call hdf5_write_data(ifile, '.scattering/gamimp', sct%gamimp)
  endif

  call hdf5_write_data(ifile, ".structure/ispin",   edisp%ispin)
  call hdf5_write_data(ifile, ".structure/charge",  edisp%nelect) ! this might have been changed by config
  call hdf5_write_data(ifile, '.structure/mudft',   pot%mu_dft)    ! this also might have changed
  call hdf5_write_data(ifile, '.structure/nkp',     kmesh%nkp)

  allocate(weights(kmesh%nkp))
  weights = kmesh%multiplicity / real(kmesh%nkp,8) ! since we only have them on the MPI range saved
  call hdf5_write_data(ifile, '.structure/weights', weights)
  deallocate(weights)

  ! remove this for the time being
  ! if (edisp%ispin == 1) then
  !   call hdf5_write_data(ifile, "/.structure/energies", edisp%band(:,:,1))
  ! else
  !   call hdf5_write_data(ifile, "/.structure/energies/up", edisp%band(:,:,1))
  !   call hdf5_write_data(ifile, "/.structure/energies/dn", edisp%band(:,:,2))
  ! endif

  ! output bandgap information to have access to it
  ! in the general output file
  if (edisp%ispin == 1) then
    call hdf5_write_data(ifile, "/.structure/bandgap/gapped", edisp%gapped(1))
    if (edisp%gapped(1)) then
      call hdf5_write_data(ifile, "/.structure/bandgap/gapsize", edisp%gap(1))
      call hdf5_write_data(ifile, "/.structure/bandgap/ene_vband", edisp%ene_valenceBand(1))
      call hdf5_write_data(ifile, "/.structure/bandgap/ene_cband", edisp%ene_conductionBand(1))
      call hdf5_write_data(ifile, "/.structure/bandgap/vband", edisp%valenceBand(1))
      call hdf5_write_data(ifile, "/.structure/bandgap/cband", edisp%conductionBand(1))
    endif
  else
    call hdf5_write_data(ifile, "/.structure/bandgap/up/gapped", edisp%gapped(1))
    if (edisp%gapped(1)) then
      call hdf5_write_data(ifile, "/.structure/bandgap/up/gapsize", edisp%gap(1))
      call hdf5_write_data(ifile, "/.structure/bandgap/up/ene_vband", edisp%ene_valenceBand(1))
      call hdf5_write_data(ifile, "/.structure/bandgap/up/ene_cband", edisp%ene_conductionBand(1))
      call hdf5_write_data(ifile, "/.structure/bandgap/up/vband", edisp%valenceBand(1))
      call hdf5_write_data(ifile, "/.structure/bandgap/up/cband", edisp%conductionBand(1))
    endif

    call hdf5_write_data(ifile, "/.structure/bandgap/dn/gapped", edisp%gapped(2))
    if (edisp%gapped(2)) then
      call hdf5_write_data(ifile, "/.structure/bandgap/dn/gapsize", edisp%gap(2))
      call hdf5_write_data(ifile, "/.structure/bandgap/dn/ene_vband", edisp%ene_valenceBand(2))
      call hdf5_write_data(ifile, "/.structure/bandgap/dn/ene_cband", edisp%ene_conductionBand(2))
      call hdf5_write_data(ifile, "/.structure/bandgap/dn/vband", edisp%valenceBand(1))
      call hdf5_write_data(ifile, "/.structure/bandgap/dn/cband", edisp%conductionBand(1))
    endif
  endif

  call hdf5_write_data(ifile, ".quantities/doping", edisp%doping)
  call hdf5_write_data(ifile, '.quantities/tempAxis', temp%TT)
  call hdf5_write_data(ifile, '.quantities/betaAxis', temp%BB)

  if (algo%lTMODE) then
    call hdf5_write_attribute(ifile, '.quantities', 'mode', 'temp')
  else if (algo%lMUMODE) then
    call hdf5_write_attribute(ifile, '.quantities', 'mode', 'mu')
  endif

  call hdf5_write_data(ifile, "/.quantities/impurities/nimp", imp%nimp)
  if (algo%lImpurities) then
    do iimp = 1, imp%nimp
      write(string,'("/.quantities/impurities/imp-",I3.3,"/energy")') iimp
      call hdf5_write_data(ifile, string, imp%Energy(iimp))
      write(string,'("/.quantities/impurities/imp-",I3.3,"/density")') iimp
      call hdf5_write_data(ifile, string, imp%Density(iimp))
      write(string,'("/.quantities/impurities/imp-",I3.3,"/degeneracy")') iimp
      call hdf5_write_data(ifile, string, imp%Degeneracy(iimp))
      write(string,'("/.quantities/impurities/imp-",I3.3,"/dopant")') iimp
      call hdf5_write_data(ifile, string, imp%Dopant(iimp))
      write(string,'("/.quantities/impurities/imp-",I3.3,"/width")') iimp
      call hdf5_write_data(ifile, string, imp%Bandwidth(iimp))
    enddo
  endif

  call hdf5_write_data(ifile, "/.unitcell/dims", kmesh%dims)
  call hdf5_write_data(ifile, "/.unitcell/ndim", kmesh%ndim)
  call hdf5_write_data(ifile, "/.unitcell/vol",  kmesh%vol)

  call hdf5_close_file(ifile)

end subroutine

! subroutine output_energies(algo, edisp, kmesh, sct, info)
!   implicit none
!   type(algorithm)  :: algo
!   type(energydisp) :: edisp
!   type(kpointmesh) :: kmesh
!   type(scattering) :: sct
!   type(runinfo)    :: info

!   integer(hid_t)     :: ifile
!   real(8), allocatable :: enrgy(:,:,:)
!   character(len=128) :: string

!   allocate(enrgy(edisp%nband_max,kmesh%nkp,edisp%ispin))

!   enrgy = sct%zqp * (edisp%band - info%mu)

!   write(string,'(I6.6,"/energies")') info%iStep
!   call hdf5_open_file(algo%output_file, ifile)
!   call hdf5_write_data(ifile, string, enrgy)

!   write(string,'(I6.6)') info%iStep
!   call hdf5_write_attribute(ifile, string, "temperature", info%temp)
!   call hdf5_write_attribute(ifile, string, "invtemperature", info%beta)
!   call hdf5_write_attribute(ifile, string, "chemicalpotential", info%mu)

!   call hdf5_close_file(ifile)

!   deallocate(enrgy)

! end subroutine



! ######################################
! main routine for the output of double precision response data
! this subroutine gets called at each step of the temperature/chemical_potential loop
!
! perform partial (band, momentum) summation, if required
! perform full (band, momentum) summation
! apply pre, post- factors to achieve SI units
!
! if a full / partial output is used, write it to the HDF5 file
! save fully summed quantities in appropriate array and write it at the final step of the temperature
! / chemical potential loop
subroutine output_response_D(resp, gname, edisp, algo, info, temp, kmesh, lBfield)
  implicit none
  type (response_dp)  :: resp
  character(len=*)    :: gname

  type(energydisp)    :: edisp
  type(algorithm)     :: algo
  type(runinfo)       :: info
  type(kpointmesh)    :: kmesh
  type(temperature)   :: temp

  logical, optional, intent(in) :: lBfield
  logical :: lBoutput

  character(len=128) :: string
  character(len=5)   :: sumstyle
  integer(hid_t)     :: ifile

  integer :: iband, ik

  ! partial sums
  complex(8), allocatable :: s_partial_sum(:,:,:,:,:)
  complex(8), allocatable :: sB_partial_sum(:,:,:,:,:,:)
  complex(8), allocatable :: a_partial_sum(:,:,:,:,:)
  complex(8), allocatable :: aB_partial_sum(:,:,:,:,:,:)
  complex(8), allocatable :: x_partial_sum(:,:,:,:,:)
  complex(8), allocatable :: xB_partial_sum(:,:,:,:,:,:)

  ! gather arrays for MPI
  complex(8), allocatable :: s_gather(:,:,:,:,:)
  complex(8), allocatable :: sB_gather(:,:,:,:,:,:)
  complex(8), allocatable :: a_gather(:,:,:,:,:)
  complex(8), allocatable :: aB_gather(:,:,:,:,:,:)
  complex(8), allocatable :: x_gather(:,:,:,:,:)
  complex(8), allocatable :: xB_gather(:,:,:,:,:,:)

  if (algo%lBfield) then
    if(present(lBfield)) then
      if (lBfield) then
        lBoutput = .true.
      else
        lBoutput = .false.
      endif
    else
      lBoutput = .true.
    endif
  else
    lBoutput = .false.
  endif

  if (myid.eq.master) then
    call hdf5_open_file(algo%output_file, ifile)
  endif

  ! conductivity and seebeck coefficient without B-field
  if (algo%fullOutput > 0) then

    ! L11
    select case (algo%fullOutput)
      case (1) ! full
        sumstyle = '/full'
        allocate(s_partial_sum(3,3,edisp%nband_max,edisp%ispin,ikstr:ikend))
      case (2) ! ksum
        sumstyle = '/ksum'
        allocate(s_partial_sum(3,3,edisp%nband_max,edisp%ispin,1))
      case (3) ! bsum
        sumstyle = '/bsum'
        allocate(s_partial_sum(3,3,1,edisp%ispin,ikstr:ikend))
    end select
    s_partial_sum = 0.d0
    select case (algo%fullOutput)
      case (1) ! shift the optical range into the full energy band range
        s_partial_sum(:,:,edisp%nbopt_min:edisp%nbopt_max,:,:) = resp%s_full
      case (2) ! shift optical range + momentum sum
        do ik=ikstr,ikend
          s_partial_sum(:,:,edisp%nbopt_min:edisp%nbopt_max,:,1) = &
                s_partial_sum(:,:,edisp%nbopt_min:edisp%nbopt_max,:,1) + &
                resp%s_full(:,:,edisp%nbopt_min:edisp%nbopt_max,:,ik) * kmesh%weight(ik)
        enddo
      case (3) ! band sum
        do iband=edisp%nbopt_min,edisp%nbopt_max
          s_partial_sum(:,:,1,:,:) = &
                s_partial_sum(:,:,1,:,:) + resp%s_full(:,:,iband,:,:)
        enddo
    end select

    if (myid .eq. master) then
      select case (algo%fullOutput)
        case (1)
          allocate(s_gather(3,3,edisp%nband_max,edisp%ispin,kmesh%nkp))
        case (2)
          allocate(s_gather(3,3,edisp%nband_max,edisp%ispin,1))
        case (3)
          allocate(s_gather(3,3,1,edisp%ispin,kmesh%nkp))
      end select
    else
      allocate(s_gather(1,1,1,1,1))
    endif
#ifdef MPI
    select case (algo%fullOutput)
      case (1)
        call MPI_gatherv(s_partial_sum,(ikend-ikstr+1)*9*edisp%nband_max*edisp%ispin, &
                         MPI_DOUBLE_COMPLEX, s_gather, rcounts*9*edisp%nband_max*edisp%ispin, &
                         displs*9*edisp%nband_max*edisp%ispin, MPI_DOUBLE_COMPLEX, master, &
                         MPI_COMM_WORLD, mpierr)
      case (2)
        if (myid.eq.master) then
          call MPI_REDUCE(MPI_IN_PLACE, s_partial_sum, &
               9*edisp%ispin*edisp%nband_max, MPI_DOUBLE_COMPLEX, MPI_SUM, master, MPI_COMM_WORLD, mpierr)
          s_gather = s_partial_sum
        else
          call MPI_REDUCE(s_partial_sum, s_partial_sum, &
               9*edisp%ispin*edisp%nband_max, MPI_DOUBLE_COMPLEX, MPI_SUM, master, MPI_COMM_WORLD, mpierr)
        endif
      case (3)
        call MPI_gatherv(s_partial_sum,(ikend-ikstr+1)*9*edisp%ispin, &
                         MPI_DOUBLE_COMPLEX, s_gather, rcounts*9*edisp%ispin, &
                         displs*9*edisp%ispin, MPI_DOUBLE_COMPLEX, master, &
                         MPI_COMM_WORLD, mpierr)
    end select
#else
    s_gather = s_partial_sum
#endif
    deallocate(s_partial_sum)

    ! L12
    select case (algo%fullOutput)
      case (1) ! full
        allocate(a_partial_sum(3,3,edisp%nband_max,edisp%ispin,ikstr:ikend))
      case (2) ! ksum
        allocate(a_partial_sum(3,3,edisp%nband_max,edisp%ispin,1))
      case (3) ! bsum
        allocate(a_partial_sum(3,3,1,edisp%ispin,ikstr:ikend))
    end select
    a_partial_sum = 0.d0
    select case (algo%fullOutput)
      case (1) ! shift the optical range into the full energy band range
        a_partial_sum(:,:,edisp%nbopt_min:edisp%nbopt_max,:,:) = resp%a_full
      case (2) ! shift optical range + momentum sum
        do ik=ikstr,ikend
          a_partial_sum(:,:,edisp%nbopt_min:edisp%nbopt_max,:,1) = &
                a_partial_sum(:,:,edisp%nbopt_min:edisp%nbopt_max,:,1) + &
                resp%a_full(:,:,edisp%nbopt_min:edisp%nbopt_max,:,ik) * kmesh%weight(ik)
        enddo
      case (3) ! band sum
        do iband=edisp%nbopt_min,edisp%nbopt_max
          a_partial_sum(:,:,1,:,:) = &
                a_partial_sum(:,:,1,:,:) + resp%a_full(:,:,iband,:,:)
        enddo
    end select

    if (myid .eq. master) then
      select case (algo%fullOutput)
        case (1)
          allocate(a_gather(3,3,edisp%nband_max,edisp%ispin,kmesh%nkp))
        case (2)
          allocate(a_gather(3,3,edisp%nband_max,edisp%ispin,1))
        case (3)
          allocate(a_gather(3,3,1,edisp%ispin,kmesh%nkp))
      end select
    else
      allocate(a_gather(1,1,1,1,1))
    endif
#ifdef MPI
    select case (algo%fullOutput)
      case (1)
        call MPI_gatherv(a_partial_sum,(ikend-ikstr+1)*9*edisp%nband_max*edisp%ispin, &
                         MPI_DOUBLE_COMPLEX, a_gather, rcounts*9*edisp%nband_max*edisp%ispin, &
                         displs*9*edisp%nband_max*edisp%ispin, MPI_DOUBLE_COMPLEX, master, &
                         MPI_COMM_WORLD, mpierr)
      case (2)
        if (myid.eq.master) then
          call MPI_REDUCE(MPI_IN_PLACE, a_partial_sum, &
               9*edisp%ispin*edisp%nband_max, MPI_DOUBLE_COMPLEX, MPI_SUM, master, MPI_COMM_WORLD, mpierr)
          a_gather = a_partial_sum
        else
          call MPI_REDUCE(a_partial_sum, a_partial_sum, &
               9*edisp%ispin*edisp%nband_max, MPI_DOUBLE_COMPLEX, MPI_SUM, master, MPI_COMM_WORLD, mpierr)
        endif
      case (3)
        call MPI_gatherv(a_partial_sum,(ikend-ikstr+1)*9*edisp%ispin, &
                         MPI_DOUBLE_COMPLEX, a_gather, rcounts*9*edisp%ispin, &
                         displs*9*edisp%ispin, MPI_DOUBLE_COMPLEX, master, &
                         MPI_COMM_WORLD, mpierr)
    end select
#else
    a_gather = a_partial_sum
#endif
    deallocate(a_partial_sum)

    ! L22
    select case (algo%fullOutput)
      case (1) ! full
        allocate(x_partial_sum(3,3,edisp%nband_max,edisp%ispin,ikstr:ikend))
      case (2) ! ksum
        allocate(x_partial_sum(3,3,edisp%nband_max,edisp%ispin,1))
      case (3) ! bsum
        allocate(x_partial_sum(3,3,1,edisp%ispin,ikstr:ikend))
    end select
    x_partial_sum = 0.d0
    select case (algo%fullOutput)
      case (1) ! shift the optical range into the full energy band range
        x_partial_sum(:,:,edisp%nbopt_min:edisp%nbopt_max,:,:) = resp%x_full
      case (2) ! shift optical range + momentum sum
        do ik=ikstr,ikend
          x_partial_sum(:,:,edisp%nbopt_min:edisp%nbopt_max,:,1) = &
                x_partial_sum(:,:,edisp%nbopt_min:edisp%nbopt_max,:,1) + &
                resp%x_full(:,:,edisp%nbopt_min:edisp%nbopt_max,:,ik) * kmesh%weight(ik)
        enddo
      case (3) ! band sum
        do iband=edisp%nbopt_min,edisp%nbopt_max
          x_partial_sum(:,:,1,:,:) = &
                x_partial_sum(:,:,1,:,:) + resp%x_full(:,:,iband,:,:)
        enddo
    end select

    if (myid .eq. master) then
      select case (algo%fullOutput)
        case (1)
          allocate(x_gather(3,3,edisp%nband_max,edisp%ispin,kmesh%nkp))
        case (2)
          allocate(x_gather(3,3,edisp%nband_max,edisp%ispin,1))
        case (3)
          allocate(x_gather(3,3,1,edisp%ispin,kmesh%nkp))
      end select
    else
      allocate(x_gather(1,1,1,1,1))
    endif
#ifdef MPI
    select case (algo%fullOutput)
      case (1)
        call MPI_gatherv(x_partial_sum,(ikend-ikstr+1)*9*edisp%nband_max*edisp%ispin, &
                         MPI_DOUBLE_COMPLEX, x_gather, rcounts*9*edisp%nband_max*edisp%ispin, &
                         displs*9*edisp%nband_max*edisp%ispin, MPI_DOUBLE_COMPLEX, master, &
                         MPI_COMM_WORLD, mpierr)
      case (2)
        if (myid.eq.master) then
          call MPI_REDUCE(MPI_IN_PLACE, x_partial_sum, &
               9*edisp%ispin*edisp%nband_max, MPI_DOUBLE_COMPLEX, MPI_SUM, master, MPI_COMM_WORLD, mpierr)
          x_gather = x_partial_sum
        else
          call MPI_REDUCE(x_partial_sum, x_partial_sum, &
               9*edisp%ispin*edisp%nband_max, MPI_DOUBLE_COMPLEX, MPI_SUM, master, MPI_COMM_WORLD, mpierr)
        endif
      case (3)
        call MPI_gatherv(x_partial_sum,(ikend-ikstr+1)*9*edisp%ispin, &
                         MPI_DOUBLE_COMPLEX, x_gather, rcounts*9*edisp%ispin, &
                         displs*9*edisp%ispin, MPI_DOUBLE_COMPLEX, master, &
                         MPI_COMM_WORLD, mpierr)
    end select
#else
    x_gather = x_partial_sum
#endif
    deallocate(x_partial_sum)

    s_gather = s_gather * pi * ( echarge / (kmesh%vol*hbarevs)) * 1.d10 ! -> 1/(Ohm*m) = A / (V * m)
    a_gather = a_gather * pi * ( echarge / (kmesh%vol*hbarevs)) * 1.d10 ! -> A / m
    x_gather = x_gather * pi * ( echarge / (kmesh%vol*hbarevs)) * 1.d10 ! -> VA / m

    if (myid .eq. master) then
      write(string,'(I6.6)') info%iStep
      string = 'step/' // trim(string) // "/L11/" // trim(adjustl(gname)) // sumstyle
      call hdf5_write_data(ifile, string, s_gather)

      write(string,'(I6.6)') info%iStep
      string = 'step/' // trim(string) // "/L12/" // trim(adjustl(gname)) // sumstyle
      call hdf5_write_data(ifile, string, a_gather)

      write(string,'(I6.6)') info%iStep
      string = 'step/' // trim(string) // "/L22/" // trim(adjustl(gname)) // sumstyle
      call hdf5_write_data(ifile, string, x_gather)
    endif

    deallocate(s_gather)
    deallocate(a_gather)
    deallocate(x_gather)
  ! full output end ... this is always for each T-point
  ! arrays would be too large for this
  endif

  ! always do the full summation!
  ! perform a local summation
  resp%s_sum = 0
  resp%a_sum = 0
  resp%x_sum = 0
  do ik = ikstr,ikend
    do iband = edisp%nbopt_min,edisp%nbopt_max
      resp%s_sum(:,:,:) = resp%s_sum(:,:,:) + resp%s_full(:,:,iband,:,ik) * kmesh%weight(ik)
      resp%a_sum(:,:,:) = resp%a_sum(:,:,:) + resp%a_full(:,:,iband,:,ik) * kmesh%weight(ik)
      resp%x_sum(:,:,:) = resp%x_sum(:,:,:) + resp%x_full(:,:,iband,:,ik) * kmesh%weight(ik)
    enddo
  enddo

  ! perform MPI summation
#ifdef MPI
  if (myid.eq.master) then
    call MPI_REDUCE(MPI_IN_PLACE, resp%s_sum, &
         9*edisp%ispin, MPI_DOUBLE_COMPLEX, MPI_SUM, master, MPI_COMM_WORLD, mpierr)
  else
    call MPI_REDUCE(resp%s_sum, resp%s_sum, &
         9*edisp%ispin, MPI_DOUBLE_COMPLEX, MPI_SUM, master, MPI_COMM_WORLD, mpierr)
  endif

  if (myid.eq.master) then
    call MPI_REDUCE(MPI_IN_PLACE, resp%a_sum, &
         9*edisp%ispin, MPI_DOUBLE_COMPLEX, MPI_SUM, master, MPI_COMM_WORLD, mpierr)
  else
    call MPI_REDUCE(resp%a_sum, resp%a_sum, &
         9*edisp%ispin, MPI_DOUBLE_COMPLEX, MPI_SUM, master, MPI_COMM_WORLD, mpierr)
  endif
  if (myid.eq.master) then
    call MPI_REDUCE(MPI_IN_PLACE, resp%x_sum, &
         9*edisp%ispin, MPI_DOUBLE_COMPLEX, MPI_SUM, master, MPI_COMM_WORLD, mpierr)
  else
    call MPI_REDUCE(resp%x_sum, resp%x_sum, &
         9*edisp%ispin, MPI_DOUBLE_COMPLEX, MPI_SUM, master, MPI_COMM_WORLD, mpierr)
  endif
#endif

  resp%s_sum = resp%s_sum * pi * ( echarge / (kmesh%vol*hbarevs)) * 1.d10 ! e**2 / hbar[Js] -> e / hbar[eVs]
  resp%a_sum = resp%a_sum * pi * ( echarge / (kmesh%vol*hbarevs)) * 1.d10
  resp%x_sum = resp%x_sum * pi * ( echarge / (kmesh%vol*hbarevs)) * 1.d10

  if (myid .eq. master) then
    ! gather the data in the arrays
    resp%s_sum_range(:,:,:,info%iStep) = resp%s_sum
    resp%a_sum_range(:,:,:,info%iStep) = resp%a_sum
    resp%x_sum_range(:,:,:,info%iStep) = resp%x_sum

    ! output at the last temperature step
    if ((algo%step_dir==1 .and. info%iStep==algo%steps) .or. (algo%step_dir==-1 .and. info%iStep==1)) then
      string = "/L11/" // trim(adjustl(gname)) // "/sum"
      call hdf5_write_data(ifile, string, resp%s_sum_range)
      string = "/L12/" // trim(adjustl(gname)) // "/sum"
      call hdf5_write_data(ifile, string, resp%a_sum_range)
      string = "/L22/" // trim(adjustl(gname)) // "/sum"
      call hdf5_write_data(ifile, string, resp%x_sum_range)
    endif
  endif


  if (lBoutput) then
    if (algo%fullOutput > 0) then

      ! L11B
      select case (algo%fullOutput)
        case (1) ! full
          allocate(sB_partial_sum(3,3,3,edisp%nband_max,edisp%ispin,ikstr:ikend))
        case (2) ! ksum
          allocate(sB_partial_sum(3,3,3,edisp%nband_max,edisp%ispin,1))
        case (3) ! bsum
          allocate(sB_partial_sum(3,3,3,1,edisp%ispin,ikstr:ikend))
      end select
      sB_partial_sum = 0.d0
      select case (algo%fullOutput)
        case (1) ! shift the optical range into the full energy band range
          sB_partial_sum(:,:,:,edisp%nbopt_min:edisp%nbopt_max,:,:) = resp%sB_full
        case (2) ! shift optical range + momentum sum
          do ik=ikstr,ikend
            sB_partial_sum(:,:,:,edisp%nbopt_min:edisp%nbopt_max,:,1) = &
                  sB_partial_sum(:,:,:,edisp%nbopt_min:edisp%nbopt_max,:,1) + &
                  resp%sB_full(:,:,:,edisp%nbopt_min:edisp%nbopt_max,:,ik) * kmesh%weight(ik)
          enddo
        case (3) ! band sum
          do iband=edisp%nbopt_min,edisp%nbopt_max
            sB_partial_sum(:,:,:,1,:,:) = &
                  sB_partial_sum(:,:,:,1,:,:) + resp%sB_full(:,:,:,iband,:,:)
          enddo
      end select

      if (myid .eq. master) then
        select case (algo%fullOutput)
          case (1)
            allocate(sB_gather(3,3,3,edisp%nband_max,edisp%ispin,kmesh%nkp))
          case (2)
            allocate(sB_gather(3,3,3,edisp%nband_max,edisp%ispin,1))
          case (3)
            allocate(sB_gather(3,3,3,1,edisp%ispin,kmesh%nkp))
        end select
      else
        allocate(sB_gather(1,1,1,1,1,1))
      endif
#ifdef MPI
      select case (algo%fullOutput)
        case (1)
          call MPI_gatherv(sB_partial_sum,(ikend-ikstr+1)*27*edisp%nband_max*edisp%ispin, &
                           MPI_DOUBLE_COMPLEX, sB_gather, rcounts*27*edisp%nband_max*edisp%ispin, &
                           displs*27*edisp%nband_max*edisp%ispin, MPI_DOUBLE_COMPLEX, master, &
                           MPI_COMM_WORLD, mpierr)
        case (2)
          if (myid.eq.master) then
            call MPI_REDUCE(MPI_IN_PLACE, sB_partial_sum, &
                 27*edisp%ispin*edisp%nband_max, MPI_DOUBLE_COMPLEX, MPI_SUM, master, MPI_COMM_WORLD, mpierr)
            sB_gather = sB_partial_sum
          else
            call MPI_REDUCE(sB_partial_sum, sB_partial_sum, &
                 27*edisp%ispin*edisp%nband_max, MPI_DOUBLE_COMPLEX, MPI_SUM, master, MPI_COMM_WORLD, mpierr)
          endif
        case (3)
          call MPI_gatherv(sB_partial_sum,(ikend-ikstr+1)*27*edisp%ispin, &
                           MPI_DOUBLE_COMPLEX, sB_gather, rcounts*27*edisp%ispin, &
                           displs*27*edisp%ispin, MPI_DOUBLE_COMPLEX, master, &
                           MPI_COMM_WORLD, mpierr)
      end select
#else
      sB_gather = sB_partial_sum
#endif
      deallocate(sB_partial_sum)

      ! L12B
      select case (algo%fullOutput)
        case (1) ! full
          allocate(aB_partial_sum(3,3,3,edisp%nband_max,edisp%ispin,ikstr:ikend))
        case (2) ! ksum
          allocate(aB_partial_sum(3,3,3,edisp%nband_max,edisp%ispin,1))
        case (3) ! bsum
          allocate(aB_partial_sum(3,3,3,1,edisp%ispin,ikstr:ikend))
      end select
      aB_partial_sum = 0.d0
      select case (algo%fullOutput)
        case (1) ! shift the optical range into the full energy band range
          aB_partial_sum(:,:,:,edisp%nbopt_min:edisp%nbopt_max,:,:) = resp%aB_full
        case (2) ! shift optical range + momentum sum
          do ik=ikstr,ikend
            aB_partial_sum(:,:,:,edisp%nbopt_min:edisp%nbopt_max,:,1) = &
                  aB_partial_sum(:,:,:,edisp%nbopt_min:edisp%nbopt_max,:,1) + &
                  resp%aB_full(:,:,:,edisp%nbopt_min:edisp%nbopt_max,:,ik) * kmesh%weight(ik)
          enddo
        case (3) ! band sum
          do iband=edisp%nbopt_min,edisp%nbopt_max
            aB_partial_sum(:,:,:,1,:,:) = &
                  aB_partial_sum(:,:,:,1,:,:) + resp%aB_full(:,:,:,iband,:,:)
          enddo
      end select

      if (myid .eq. master) then
        select case (algo%fullOutput)
          case (1)
            allocate(aB_gather(3,3,3,edisp%nband_max,edisp%ispin,kmesh%nkp))
          case (2)
            allocate(aB_gather(3,3,3,edisp%nband_max,edisp%ispin,1))
          case (3)
            allocate(aB_gather(3,3,3,1,edisp%ispin,kmesh%nkp))
        end select
      else
        allocate(aB_gather(1,1,1,1,1,1))
      endif
#ifdef MPI
      select case (algo%fullOutput)
        case (1)
          call MPI_gatherv(aB_partial_sum,(ikend-ikstr+1)*27*edisp%nband_max*edisp%ispin, &
                           MPI_DOUBLE_COMPLEX, aB_gather, rcounts*27*edisp%nband_max*edisp%ispin, &
                           displs*27*edisp%nband_max*edisp%ispin, MPI_DOUBLE_COMPLEX, master, &
                           MPI_COMM_WORLD, mpierr)
        case (2)
          if (myid.eq.master) then
            call MPI_REDUCE(MPI_IN_PLACE, aB_partial_sum, &
                 27*edisp%ispin*edisp%nband_max, MPI_DOUBLE_COMPLEX, MPI_SUM, master, MPI_COMM_WORLD, mpierr)
            aB_gather = aB_partial_sum
          else
            call MPI_REDUCE(aB_partial_sum, aB_partial_sum, &
                 27*edisp%ispin*edisp%nband_max, MPI_DOUBLE_COMPLEX, MPI_SUM, master, MPI_COMM_WORLD, mpierr)
          endif
        case (3)
          call MPI_gatherv(aB_partial_sum,(ikend-ikstr+1)*27*edisp%ispin, &
                           MPI_DOUBLE_COMPLEX, aB_gather, rcounts*27*edisp%ispin, &
                           displs*27*edisp%ispin, MPI_DOUBLE_COMPLEX, master, &
                           MPI_COMM_WORLD, mpierr)
      end select
#else
      aB_gather = aB_partial_sum
#endif
      deallocate(aB_partial_sum)

      ! L22B
      select case (algo%fullOutput)
        case (1) ! full
          allocate(xB_partial_sum(3,3,3,edisp%nband_max,edisp%ispin,ikstr:ikend))
        case (2) ! ksum
          allocate(xB_partial_sum(3,3,3,edisp%nband_max,edisp%ispin,1))
        case (3) ! bsum
          allocate(xB_partial_sum(3,3,3,1,edisp%ispin,ikstr:ikend))
      end select
      xB_partial_sum = 0.d0
      select case (algo%fullOutput)
        case (1) ! shift the optical range into the full energy band range
          xB_partial_sum(:,:,:,edisp%nbopt_min:edisp%nbopt_max,:,:) = resp%xB_full
        case (2) ! shift optical range + momentum sum
          do ik=ikstr,ikend
            xB_partial_sum(:,:,:,edisp%nbopt_min:edisp%nbopt_max,:,1) = &
                  xB_partial_sum(:,:,:,edisp%nbopt_min:edisp%nbopt_max,:,1) + &
                  resp%xB_full(:,:,:,edisp%nbopt_min:edisp%nbopt_max,:,ik) * kmesh%weight(ik)
          enddo
        case (3) ! band sum
          do iband=edisp%nbopt_min,edisp%nbopt_max
            xB_partial_sum(:,:,:,1,:,:) = &
                  xB_partial_sum(:,:,:,1,:,:) + resp%xB_full(:,:,:,iband,:,:)
          enddo
      end select

      if (myid .eq. master) then
        select case (algo%fullOutput)
          case (1)
            allocate(xB_gather(3,3,3,edisp%nband_max,edisp%ispin,kmesh%nkp))
          case (2)
            allocate(xB_gather(3,3,3,edisp%nband_max,edisp%ispin,1))
          case (3)
            allocate(xB_gather(3,3,3,1,edisp%ispin,kmesh%nkp))
        end select
      else
        allocate(xB_gather(1,1,1,1,1,1))
      endif
#ifdef MPI
      select case (algo%fullOutput)
        case (1)
          call MPI_gatherv(xB_partial_sum,(ikend-ikstr+1)*27*edisp%nband_max*edisp%ispin, &
                           MPI_DOUBLE_COMPLEX, xB_gather, rcounts*27*edisp%nband_max*edisp%ispin, &
                           displs*27*edisp%nband_max*edisp%ispin, MPI_DOUBLE_COMPLEX, master, &
                           MPI_COMM_WORLD, mpierr)
        case (2)
          if (myid.eq.master) then
            call MPI_REDUCE(MPI_IN_PLACE, xB_partial_sum, &
                 27*edisp%ispin*edisp%nband_max, MPI_DOUBLE_COMPLEX, MPI_SUM, master, MPI_COMM_WORLD, mpierr)
            xB_gather = xB_partial_sum
          else
            call MPI_REDUCE(xB_partial_sum, xB_partial_sum, &
                 27*edisp%ispin*edisp%nband_max, MPI_DOUBLE_COMPLEX, MPI_SUM, master, MPI_COMM_WORLD, mpierr)
          endif
        case (3)
          call MPI_gatherv(xB_partial_sum,(ikend-ikstr+1)*27*edisp%ispin, &
                           MPI_DOUBLE_COMPLEX, xB_gather, rcounts*27*edisp%ispin, &
                           displs*27*edisp%ispin, MPI_DOUBLE_COMPLEX, master, &
                           MPI_COMM_WORLD, mpierr)
      end select
#else
      xB_gather = xB_partial_sum
#endif
      deallocate(xB_partial_sum)

      sB_gather = sB_gather * 4.d0 / 3.d0 * pi**2 * ( echarge / (kmesh%vol*hbarevs)) * (1.d-10 / hbarevs) ! -> A * m / (V**2 * s)
      aB_gather = aB_gather * 4.d0 / 3.d0 * pi**2 * ( echarge / (kmesh%vol*hbarevs)) * (1.d-10 / hbarevs) ! -> A**2 * m / V
      xB_gather = xB_gather * 4.d0 / 3.d0 * pi**2 * ( echarge / (kmesh%vol*hbarevs)) * (1.d-10 / hbarevs) ! -> A**3 * m * s

      if (myid .eq. master) then
        write(string,'(I6.6)') info%iStep
        string = 'step/' // trim(string) // "/L11B/" // trim(adjustl(gname)) // sumstyle
        call hdf5_write_data(ifile, string, sB_gather)

        write(string,'(I6.6)') info%iStep
        string = 'step/' // trim(string) // "/L12B/" // trim(adjustl(gname)) // sumstyle
        call hdf5_write_data(ifile, string, aB_gather)

        write(string,'(I6.6)') info%iStep
        string = 'step/' // trim(string) // "/L22B/" // trim(adjustl(gname)) // sumstyle
        call hdf5_write_data(ifile, string, xB_gather)

      endif

      deallocate(sB_gather)
      deallocate(aB_gather)
      deallocate(xB_gather)
    endif ! full output

    ! perform a local summation
    ! these are already initialized to 0
    do ik = ikstr,ikend
      do iband = edisp%nbopt_min,edisp%nbopt_max
        resp%sB_sum(:,:,:,:) = resp%sB_sum(:,:,:,:) + resp%sB_full(:,:,:,iband,:,ik) * kmesh%weight(ik)
        resp%aB_sum(:,:,:,:) = resp%aB_sum(:,:,:,:) + resp%aB_full(:,:,:,iband,:,ik) * kmesh%weight(ik)
        resp%xB_sum(:,:,:,:) = resp%xB_sum(:,:,:,:) + resp%xB_full(:,:,:,iband,:,ik) * kmesh%weight(ik)
      enddo
    enddo

  ! perform MPI summation
#ifdef MPI
  if (myid.eq.master) then
    call MPI_REDUCE(MPI_IN_PLACE, resp%sB_sum, &
         27*edisp%ispin, MPI_DOUBLE_COMPLEX, MPI_SUM, master, MPI_COMM_WORLD, mpierr)
  else
    call MPI_REDUCE(resp%sB_sum, resp%sB_sum, &
         27*edisp%ispin, MPI_DOUBLE_COMPLEX, MPI_SUM, master, MPI_COMM_WORLD, mpierr)
  endif

  if (myid.eq.master) then
    call MPI_REDUCE(MPI_IN_PLACE, resp%aB_sum, &
         27*edisp%ispin, MPI_DOUBLE_COMPLEX, MPI_SUM, master, MPI_COMM_WORLD, mpierr)
  else
    call MPI_REDUCE(resp%aB_sum, resp%aB_sum, &
         27*edisp%ispin, MPI_DOUBLE_COMPLEX, MPI_SUM, master, MPI_COMM_WORLD, mpierr)
  endif

  if (myid.eq.master) then
    call MPI_REDUCE(MPI_IN_PLACE, resp%xB_sum, &
         27*edisp%ispin, MPI_DOUBLE_COMPLEX, MPI_SUM, master, MPI_COMM_WORLD, mpierr)
  else
    call MPI_REDUCE(resp%xB_sum, resp%xB_sum, &
         27*edisp%ispin, MPI_DOUBLE_COMPLEX, MPI_SUM, master, MPI_COMM_WORLD, mpierr)
  endif
#endif

    resp%sB_sum = resp%sB_sum * 4.d0 / 3.d0 * pi**2 * ( echarge / (kmesh%vol*hbarevs)) * (1.d-10 / hbarevs)
    resp%aB_sum = resp%aB_sum * 4.d0 / 3.d0 * pi**2 * ( echarge / (kmesh%vol*hbarevs)) * (1.d-10 / hbarevs)
    resp%xB_sum = resp%xB_sum * 4.d0 / 3.d0 * pi**2 * ( echarge / (kmesh%vol*hbarevs)) * (1.d-10 / hbarevs)

    if (myid .eq. master) then
      ! gather the data in the arrays
      resp%sB_sum_range(:,:,:,:,info%iStep) = resp%sB_sum
      resp%aB_sum_range(:,:,:,:,info%iStep) = resp%aB_sum
      resp%xB_sum_range(:,:,:,:,info%iStep) = resp%xB_sum

      ! output at the last temperature step
      if ((algo%step_dir==1 .and. info%iStep==algo%steps) .or. (algo%step_dir==-1 .and. info%iStep==1)) then
        string = "/L11B/" // trim(adjustl(gname)) // "/sum"
        call hdf5_write_data(ifile, string, resp%sB_sum_range)
        string = "/L12B/" // trim(adjustl(gname)) // "/sum"
        call hdf5_write_data(ifile, string, resp%aB_sum_range)
        string = "/L22B/" // trim(adjustl(gname)) // "/sum"
        call hdf5_write_data(ifile, string, resp%xB_sum_range)
      endif
    endif
  endif ! Boutput

  if (myid.eq.master) then
    call hdf5_close_file(ifile)
  endif

end subroutine


! main routine for the output of quad precision response data
! this subroutine gets called at each step of the temperature/chemical_potential loop
!
! perform partial (band, momentum) summation, if required
! perform full (band, momentum) summation
! apply pre, post- factors to achieve SI units
!
! if a full / partial output is used, write it to the HDF5 file
! save fully summed quantities in appropriate array and write it at the final step of the temperature
! / chemical potential loop
subroutine output_response_Q(resp, gname, edisp, algo, info, temp, kmesh, lBfield)
  ! for the quad precision we don't have a full output
  ! since the implemented MPI routines only support double-precision
  implicit none
  type (response_qp)  :: resp
  character(len=*)    :: gname

  type(energydisp)    :: edisp
  type(algorithm)     :: algo
  type(runinfo)       :: info
  type(kpointmesh)    :: kmesh
  type(temperature)   :: temp

  logical, optional, intent(in) :: lBfield
  logical :: lBoutput

  character(len=128) :: string
  character(len=5)   :: sumstyle
  integer(hid_t)     :: ifile

  ! partial sums
  complex(16), allocatable :: s_partial_sum(:,:,:,:,:)
  complex(16), allocatable :: sB_partial_sum(:,:,:,:,:,:)
  complex(16), allocatable :: a_partial_sum(:,:,:,:,:)
  complex(16), allocatable :: aB_partial_sum(:,:,:,:,:,:)
  complex(16), allocatable :: x_partial_sum(:,:,:,:,:)
  complex(16), allocatable :: xB_partial_sum(:,:,:,:,:,:)

  ! partial sums
  complex(8), allocatable :: s_partial_sum_dp(:,:,:,:,:)
  complex(8), allocatable :: sB_partial_sum_dp(:,:,:,:,:,:)
  complex(8), allocatable :: a_partial_sum_dp(:,:,:,:,:)
  complex(8), allocatable :: aB_partial_sum_dp(:,:,:,:,:,:)
  complex(8), allocatable :: x_partial_sum_dp(:,:,:,:,:)
  complex(8), allocatable :: xB_partial_sum_dp(:,:,:,:,:,:)

  ! gather arrays for MPI
  complex(8), allocatable :: s_gather(:,:,:,:,:)
  complex(8), allocatable :: sB_gather(:,:,:,:,:,:)
  complex(8), allocatable :: a_gather(:,:,:,:,:)
  complex(8), allocatable :: aB_gather(:,:,:,:,:,:)
  complex(8), allocatable :: x_gather(:,:,:,:,:)
  complex(8), allocatable :: xB_gather(:,:,:,:,:,:)

  ! quadruple response s/a/x array
  real(16), allocatable :: qrsarr(:,:,:) ! to collect the data
  real(16), allocatable :: qraarr(:,:,:)
  real(16), allocatable :: qrxarr(:,:,:)
  real(16), allocatable :: qisarr(:,:,:)
  real(16), allocatable :: qiaarr(:,:,:)
  real(16), allocatable :: qixarr(:,:,:)

  ! magnetic quadruple response s/a/x array
  real(16), allocatable :: qrsarrB(:,:,:,:) ! to collect the data
  real(16), allocatable :: qraarrB(:,:,:,:)
  real(16), allocatable :: qrxarrB(:,:,:,:)
  real(16), allocatable :: qisarrB(:,:,:,:)
  real(16), allocatable :: qiaarrB(:,:,:,:)
  real(16), allocatable :: qixarrB(:,:,:,:)

  ! quad gather for partial momentum summation
  real(16), allocatable :: qrksumarr(:,:,:,:,:)
  real(16), allocatable :: qiksumarr(:,:,:,:,:)
  real(16), allocatable :: qrksumarrB(:,:,:,:,:,:)
  real(16), allocatable :: qiksumarrB(:,:,:,:,:,:)

  integer :: iband, ik, is
  integer :: ii, ij

  allocate(qrsarr(3,3,edisp%ispin))
  allocate(qraarr(3,3,edisp%ispin))
  allocate(qrxarr(3,3,edisp%ispin))
  allocate(qisarr(3,3,edisp%ispin))
  allocate(qiaarr(3,3,edisp%ispin))
  allocate(qixarr(3,3,edisp%ispin))

  ! we need these switches since the interband quantities
  ! do not have these arrays
  if (algo%lBfield) then
    if(present(lBfield)) then
      if (lBfield) then
        lBoutput = .true.
      else
        lBoutput = .false.
      endif
    else
      lBoutput = .true.
    endif
  else
    lBoutput = .false.
  endif

  if (lBoutput) then
    allocate(qrsarrB(3,3,3,edisp%ispin))
    allocate(qraarrB(3,3,3,edisp%ispin))
    allocate(qrxarrB(3,3,3,edisp%ispin))
    allocate(qisarrB(3,3,3,edisp%ispin))
    allocate(qiaarrB(3,3,3,edisp%ispin))
    allocate(qixarrB(3,3,3,edisp%ispin))
  endif

  if (myid.eq.master) then
    call hdf5_open_file(algo%output_file, ifile)
  endif

  ! conductivity and seebeck coefficient without B-field
  if (algo%fullOutput > 0) then
    ! L11
    select case (algo%fullOutput)
      case (1) ! full
        sumstyle = '/full'
        allocate(s_partial_sum_dp(3,3,edisp%nband_max,edisp%ispin,ikstr:ikend))
      case (2) ! ksum
        sumstyle = '/ksum'
        allocate(s_partial_sum(3,3,edisp%nband_max,edisp%ispin,1))
        allocate(s_partial_sum_dp(3,3,edisp%nband_max,edisp%ispin,1))
        s_partial_sum = 0.q0
      case (3) ! bsum
        sumstyle = '/bsum'
        allocate(s_partial_sum(3,3,1,edisp%ispin,ikstr:ikend))
        allocate(s_partial_sum_dp(3,3,1,edisp%ispin,ikstr:ikend))
        s_partial_sum = 0.q0
    end select
    s_partial_sum_dp = 0.d0
    select case (algo%fullOutput)
      case (1) ! shift the optical range into the full energy band range
        s_partial_sum_dp(:,:,edisp%nbopt_min:edisp%nbopt_max,:,:) = &
            resp%s_full * piQ * (echarge / (kmesh%vol * hbarevs)) * 1.q10
      case (2) ! shift optical range + momentum sum
        do ik=ikstr,ikend
          s_partial_sum(:,:,edisp%nbopt_min:edisp%nbopt_max,:,1) = &
                s_partial_sum(:,:,edisp%nbopt_min:edisp%nbopt_max,:,1) + &
                resp%s_full(:,:,edisp%nbopt_min:edisp%nbopt_max,:,ik) * kmesh%weightQ(ik)
        enddo
        s_partial_sum = s_partial_sum * piQ * ( echarge / (kmesh%vol * hbarevs)) * 1.q10
      case (3) ! band sum
        do iband=edisp%nbopt_min,edisp%nbopt_max
          s_partial_sum(:,:,1,:,:) = &
                s_partial_sum(:,:,1,:,:) + resp%s_full(:,:,iband,:,:)
        enddo
        s_partial_sum_dp = s_partial_sum * piQ * (echarge / (kmesh%vol * hbarevs)) * 1.q10
        deallocate(s_partial_sum)
    end select
    ! note that I put the unit factors to get the most out of the accuracy

    if (myid .eq. master) then
      select case (algo%fullOutput)
        case (1)
          allocate(s_gather(3,3,edisp%nband_max,edisp%ispin,kmesh%nkp))
        case (2)
          allocate(s_gather(3,3,edisp%nband_max,edisp%ispin,1))
        case (3)
          allocate(s_gather(3,3,1,edisp%ispin,kmesh%nkp))
      end select
    else
      allocate(s_gather(1,1,1,1,1))
    endif
#ifdef MPI
    select case (algo%fullOutput)
      case (1)
        call MPI_gatherv(s_partial_sum_dp,(ikend-ikstr+1)*9*edisp%nband_max*edisp%ispin, &
                         MPI_DOUBLE_COMPLEX, s_gather, rcounts*9*edisp%nband_max*edisp%ispin, &
                         displs*9*edisp%nband_max*edisp%ispin, MPI_DOUBLE_COMPLEX, master, &
                         MPI_COMM_WORLD, mpierr)
      case (2)
        allocate(qrksumarr(3,3,edisp%nband_max,edisp%ispin,1))
        allocate(qiksumarr(3,3,edisp%nband_max,edisp%ispin,1))
        qrksumarr = 0.q0
        qiksumarr = 0.q0
        ! we want to keep the k-sum accuracy over all cores
        do ii=1,3
          do ij=1,3
            do iband=1,edisp%nband_max
              do is=1,edisp%ispin
                call mpi_reduce_quad(real(s_partial_sum(ii,ij,iband,is,1)),qrksumarr(ii,ij,iband,is,1))
                call mpi_reduce_quad(aimag(s_partial_sum(ii,ij,iband,is,1)),qiksumarr(ii,ij,iband,is,1))
              enddo
            enddo
          enddo
        enddo
        s_gather = cmplx(real(qrksumarr,8),real(qiksumarr,8))
        deallocate(qrksumarr)
        deallocate(qiksumarr)
      case (3)
        call MPI_gatherv(s_partial_sum_dp,(ikend-ikstr+1)*9*edisp%ispin, &
                         MPI_DOUBLE_COMPLEX, s_gather, rcounts*9*edisp%ispin, &
                         displs*9*edisp%ispin, MPI_DOUBLE_COMPLEX, master, &
                         MPI_COMM_WORLD, mpierr)
    end select
#else
    if (algo%fullOutput == 2) then
      s_gather = s_partial_sum
    else
      s_gather = s_partial_sum_dp
    endif
#endif
    if (allocated(s_partial_sum)) deallocate(s_partial_sum)
    if (allocated(s_partial_sum_dp)) deallocate(s_partial_sum_dp)

    ! L12
    select case (algo%fullOutput)
      case (1) ! full
        allocate(a_partial_sum_dp(3,3,edisp%nband_max,edisp%ispin,ikstr:ikend))
      case (2) ! ksum
        allocate(a_partial_sum(3,3,edisp%nband_max,edisp%ispin,1))
        allocate(a_partial_sum_dp(3,3,edisp%nband_max,edisp%ispin,1))
        a_partial_sum = 0.q0
      case (3) ! bsum
        allocate(a_partial_sum(3,3,1,edisp%ispin,ikstr:ikend))
        allocate(a_partial_sum_dp(3,3,1,edisp%ispin,ikstr:ikend))
        a_partial_sum = 0.q0
    end select
    a_partial_sum_dp = 0.d0
    select case (algo%fullOutput)
      case (1) ! shift the optical range into the full energy band range
        a_partial_sum_dp(:,:,edisp%nbopt_min:edisp%nbopt_max,:,:) = &
            resp%a_full * piQ * (echarge / (kmesh%vol * hbarevs)) * 1.q10
      case (2) ! shift optical range + momentum sum
        do ik=ikstr,ikend
          a_partial_sum(:,:,edisp%nbopt_min:edisp%nbopt_max,:,1) = &
                a_partial_sum(:,:,edisp%nbopt_min:edisp%nbopt_max,:,1) + &
                resp%a_full(:,:,edisp%nbopt_min:edisp%nbopt_max,:,ik) * kmesh%weightQ(ik)
        enddo
        a_partial_sum = a_partial_sum * piQ * ( echarge / (kmesh%vol * hbarevs)) * 1.q10
      case (3) ! band sum
        do iband=edisp%nbopt_min,edisp%nbopt_max
          a_partial_sum(:,:,1,:,:) = &
                a_partial_sum(:,:,1,:,:) + resp%a_full(:,:,iband,:,:)
        enddo
        a_partial_sum_dp = a_partial_sum * piQ * (echarge / (kmesh%vol * hbarevs)) * 1.q10
        deallocate(a_partial_sum)
    end select
    ! note that I put the unit factors to get the most out of the accuracy

    if (myid .eq. master) then
      select case (algo%fullOutput)
        case (1)
          allocate(a_gather(3,3,edisp%nband_max,edisp%ispin,kmesh%nkp))
        case (2)
          allocate(a_gather(3,3,edisp%nband_max,edisp%ispin,1))
        case (3)
          allocate(a_gather(3,3,1,edisp%ispin,kmesh%nkp))
      end select
    else
      allocate(a_gather(1,1,1,1,1))
    endif
#ifdef MPI
    select case (algo%fullOutput)
      case (1)
        call MPI_gatherv(a_partial_sum_dp,(ikend-ikstr+1)*9*edisp%nband_max*edisp%ispin, &
                         MPI_DOUBLE_COMPLEX, a_gather, rcounts*9*edisp%nband_max*edisp%ispin, &
                         displs*9*edisp%nband_max*edisp%ispin, MPI_DOUBLE_COMPLEX, master, &
                         MPI_COMM_WORLD, mpierr)
      case (2)
        allocate(qrksumarr(3,3,edisp%nband_max,edisp%ispin,1))
        allocate(qiksumarr(3,3,edisp%nband_max,edisp%ispin,1))
        qrksumarr = 0.q0
        qiksumarr = 0.q0
        ! we want to keep the k-sum accuracy over all cores
        do ii=1,3
          do ij=1,3
            do iband=1,edisp%nband_max
              do is=1,edisp%ispin
                call mpi_reduce_quad(real(a_partial_sum(ii,ij,iband,is,1)),qrksumarr(ii,ij,iband,is,1))
                call mpi_reduce_quad(aimag(a_partial_sum(ii,ij,iband,is,1)),qiksumarr(ii,ij,iband,is,1))
              enddo
            enddo
          enddo
        enddo
        a_gather = cmplx(real(qrksumarr,8),real(qiksumarr,8))
        deallocate(qrksumarr)
        deallocate(qiksumarr)
      case (3)
        call MPI_gatherv(a_partial_sum_dp,(ikend-ikstr+1)*9*edisp%ispin, &
                         MPI_DOUBLE_COMPLEX, a_gather, rcounts*9*edisp%ispin, &
                         displs*9*edisp%ispin, MPI_DOUBLE_COMPLEX, master, &
                         MPI_COMM_WORLD, mpierr)
    end select
#else
    if (algo%fullOutput == 2) then
      a_gather = a_partial_sum
    else
      a_gather = a_partial_sum_dp
    endif
#endif
    if (allocated(a_partial_sum)) deallocate(a_partial_sum)
    if (allocated(a_partial_sum_dp)) deallocate(a_partial_sum_dp)


    ! L12
    select case (algo%fullOutput)
      case (1) ! full
        allocate(x_partial_sum_dp(3,3,edisp%nband_max,edisp%ispin,ikstr:ikend))
      case (2) ! ksum
        allocate(x_partial_sum(3,3,edisp%nband_max,edisp%ispin,1))
        allocate(x_partial_sum_dp(3,3,edisp%nband_max,edisp%ispin,1))
        x_partial_sum = 0.q0
      case (3) ! bsum
        allocate(x_partial_sum(3,3,1,edisp%ispin,ikstr:ikend))
        allocate(x_partial_sum_dp(3,3,1,edisp%ispin,ikstr:ikend))
        x_partial_sum = 0.q0
    end select
    x_partial_sum_dp = 0.d0
    select case (algo%fullOutput)
      case (1) ! shift the optical range into the full energy band range
        x_partial_sum_dp(:,:,edisp%nbopt_min:edisp%nbopt_max,:,:) = &
            resp%x_full * piQ * (echarge / (kmesh%vol * hbarevs)) * 1.q10
      case (2) ! shift optical range + momentum sum
        do ik=ikstr,ikend
          x_partial_sum(:,:,edisp%nbopt_min:edisp%nbopt_max,:,1) = &
                x_partial_sum(:,:,edisp%nbopt_min:edisp%nbopt_max,:,1) + &
                resp%x_full(:,:,edisp%nbopt_min:edisp%nbopt_max,:,ik) * kmesh%weightQ(ik)
        enddo
        x_partial_sum = x_partial_sum * piQ * ( echarge / (kmesh%vol * hbarevs)) * 1.q10
      case (3) ! band sum
        do iband=edisp%nbopt_min,edisp%nbopt_max
          x_partial_sum(:,:,1,:,:) = &
                x_partial_sum(:,:,1,:,:) + resp%x_full(:,:,iband,:,:)
        enddo
        x_partial_sum_dp = x_partial_sum * piQ * (echarge / (kmesh%vol * hbarevs)) * 1.q10
        deallocate(x_partial_sum)
    end select
    ! note that I put the unit factors to get the most out of the accuracy

    if (myid .eq. master) then
      select case (algo%fullOutput)
        case (1)
          allocate(x_gather(3,3,edisp%nband_max,edisp%ispin,kmesh%nkp))
        case (2)
          allocate(x_gather(3,3,edisp%nband_max,edisp%ispin,1))
        case (3)
          allocate(x_gather(3,3,1,edisp%ispin,kmesh%nkp))
      end select
    else
      allocate(x_gather(1,1,1,1,1))
    endif
#ifdef MPI
    select case (algo%fullOutput)
      case (1)
        call MPI_gatherv(x_partial_sum_dp,(ikend-ikstr+1)*9*edisp%nband_max*edisp%ispin, &
                         MPI_DOUBLE_COMPLEX, x_gather, rcounts*9*edisp%nband_max*edisp%ispin, &
                         displs*9*edisp%nband_max*edisp%ispin, MPI_DOUBLE_COMPLEX, master, &
                         MPI_COMM_WORLD, mpierr)
      case (2)
        allocate(qrksumarr(3,3,edisp%nband_max,edisp%ispin,1))
        allocate(qiksumarr(3,3,edisp%nband_max,edisp%ispin,1))
        qrksumarr = 0.q0
        qiksumarr = 0.q0
        ! we want to keep the k-sum accuracy over all cores
        do ii=1,3
          do ij=1,3
            do iband=1,edisp%nband_max
              do is=1,edisp%ispin
                call mpi_reduce_quad(real(x_partial_sum(ii,ij,iband,is,1)),qrksumarr(ii,ij,iband,is,1))
                call mpi_reduce_quad(aimag(x_partial_sum(ii,ij,iband,is,1)),qiksumarr(ii,ij,iband,is,1))
              enddo
            enddo
          enddo
        enddo
        x_gather = cmplx(real(qrksumarr,8),real(qiksumarr,8))
        deallocate(qrksumarr)
        deallocate(qiksumarr)
      case (3)
        call MPI_gatherv(x_partial_sum_dp,(ikend-ikstr+1)*9*edisp%ispin, &
                         MPI_DOUBLE_COMPLEX, x_gather, rcounts*9*edisp%ispin, &
                         displs*9*edisp%ispin, MPI_DOUBLE_COMPLEX, master, &
                         MPI_COMM_WORLD, mpierr)
    end select
#else
    if (algo%fullOutput == 2) then
      x_gather = x_partial_sum
    else
      x_gather = x_partial_sum_dp
    endif
#endif
    if (allocated(x_partial_sum)) deallocate(x_partial_sum)
    if (allocated(x_partial_sum_dp)) deallocate(x_partial_sum_dp)


    if (myid .eq. master) then
      write(string,'(I6.6)') info%iStep
      string = 'step/' // trim(string) // "/L11/" // trim(adjustl(gname)) // sumstyle
      call hdf5_write_data(ifile, string, s_gather)

      write(string,'(I6.6)') info%iStep
      string = 'step/' // trim(string) // "/L12/" // trim(adjustl(gname)) // sumstyle
      call hdf5_write_data(ifile, string, a_gather)

      write(string,'(I6.6)') info%iStep
      string = 'step/' // trim(string) // "/L22/" // trim(adjustl(gname)) // sumstyle
      call hdf5_write_data(ifile, string, x_gather)
    endif

    deallocate(s_gather)
    deallocate(a_gather)
    deallocate(x_gather)
  ! full output end ... this is always for each T-point
  endif

  ! perform a local summation
  do ik = ikstr,ikend
    do iband = edisp%nbopt_min,edisp%nbopt_max
      resp%s_sum(:,:,:) = resp%s_sum(:,:,:) + resp%s_full(:,:,iband,:,ik) * kmesh%weightQ(ik)
      resp%a_sum(:,:,:) = resp%a_sum(:,:,:) + resp%a_full(:,:,iband,:,ik) * kmesh%weightQ(ik)
      resp%x_sum(:,:,:) = resp%x_sum(:,:,:) + resp%x_full(:,:,iband,:,ik) * kmesh%weightQ(ik)
    enddo
  enddo


  ! perform MPI summation
  qrsarr = 0.q0
  qraarr = 0.q0
  qrxarr = 0.q0
  qisarr = 0.q0
  qiaarr = 0.q0
  qixarr = 0.q0
#ifdef MPI
  do ii=1,3
    do ij=1,3
      do is=1,edisp%ispin
        call mpi_reduce_quad(real(resp%s_sum(ii,ij,is)),qrsarr(ii,ij,is))
        call mpi_reduce_quad(aimag(resp%s_sum(ii,ij,is)),qisarr(ii,ij,is))

        call mpi_reduce_quad(real(resp%a_sum(ii,ij,is)),qraarr(ii,ij,is))
        call mpi_reduce_quad(aimag(resp%a_sum(ii,ij,is)),qiaarr(ii,ij,is))

        call mpi_reduce_quad(real(resp%x_sum(ii,ij,is)),qrxarr(ii,ij,is))
        call mpi_reduce_quad(aimag(resp%x_sum(ii,ij,is)),qixarr(ii,ij,is))
      enddo
    enddo
  enddo
#else
  qrsarr = real(resp%s_sum)
  qisarr = aimag(resp%s_sum)

  qraarr = real(resp%a_sum)
  qiaarr = aimag(resp%a_sum)

  qrxarr = real(resp%x_sum)
  qixarr = aimag(resp%x_sum)
#endif

  qrsarr = qrsarr * piQ * ( echarge / (kmesh%vol*hbarevs)) * 1.q10
  qisarr = qisarr * piQ * ( echarge / (kmesh%vol*hbarevs)) * 1.q10
  qraarr = qraarr * piQ * ( echarge / (kmesh%vol*hbarevs)) * 1.q10
  qiaarr = qiaarr * piQ * ( echarge / (kmesh%vol*hbarevs)) * 1.q10
  qrxarr = qrxarr * piQ * ( echarge / (kmesh%vol*hbarevs)) * 1.q10
  qixarr = qixarr * piQ * ( echarge / (kmesh%vol*hbarevs)) * 1.q10

  if (myid .eq. master) then
    ! gather the data in the arrays
    resp%s_sum_range(:,:,:,info%iStep) = cmplx(real(qrsarr,8),real(qisarr,8))
    resp%a_sum_range(:,:,:,info%iStep) = cmplx(real(qraarr,8),real(qiaarr,8))
    resp%x_sum_range(:,:,:,info%iStep) = cmplx(real(qrxarr,8),real(qixarr,8))

    ! output at the last temperature step
    if ((algo%step_dir==1 .and. info%iStep==algo%steps) .or. (algo%step_dir==-1 .and. info%iStep==1)) then
      string = "/L11/" // trim(adjustl(gname)) // "/sum"
      call hdf5_write_data(ifile, string, resp%s_sum_range)
      string = "/L12/" // trim(adjustl(gname)) // "/sum"
      call hdf5_write_data(ifile, string, resp%a_sum_range)
      string = "/L22/" // trim(adjustl(gname)) // "/sum"
      call hdf5_write_data(ifile, string, resp%x_sum_range)
    endif
  endif

  if (lBoutput) then
    if (algo%fullOutput > 0) then

      ! L11B
      select case (algo%fullOutput)
        case (1) ! full
          sumstyle = '/full'
          allocate(sB_partial_sum_dp(3,3,3,edisp%nband_max,edisp%ispin,ikstr:ikend))
        case (2) ! ksum
          sumstyle = '/ksum'
          allocate(sB_partial_sum(3,3,3,edisp%nband_max,edisp%ispin,1))
          allocate(sB_partial_sum_dp(3,3,3,edisp%nband_max,edisp%ispin,1))
          sB_partial_sum = 0.q0
        case (3) ! bsum
          sumstyle = '/bsum'
          allocate(sB_partial_sum(3,3,3,1,edisp%ispin,ikstr:ikend))
          allocate(sB_partial_sum_dp(3,3,3,1,edisp%ispin,ikstr:ikend))
          sB_partial_sum = 0.q0
      end select
      sB_partial_sum_dp = 0.d0
      select case (algo%fullOutput)
        case (1) ! shift the optical range into the full energy band range
          sB_partial_sum_dp(:,:,:,edisp%nbopt_min:edisp%nbopt_max,:,:) = &
              resp%sB_full * 4.q0 / 3.q0 * piQ**2 * ( echarge / (kmesh%vol*hbarevs)) * (1.q-10 / hbarevs)
        case (2) ! shift optical range + momentum sum
          do ik=ikstr,ikend
            sB_partial_sum(:,:,:,edisp%nbopt_min:edisp%nbopt_max,:,1) = &
                  sB_partial_sum(:,:,:,edisp%nbopt_min:edisp%nbopt_max,:,1) + &
                  resp%sB_full(:,:,:,edisp%nbopt_min:edisp%nbopt_max,:,ik) * kmesh%weightQ(ik)
          enddo
          sB_partial_sum = sB_partial_sum &
              * 4.q0 / 3.q0 * piQ**2 * ( echarge / (kmesh%vol*hbarevs)) * (1.q-10 / hbarevs)
        case (3) ! band sum
          do iband=edisp%nbopt_min,edisp%nbopt_max
            sB_partial_sum(:,:,:,1,:,:) = &
                  sB_partial_sum(:,:,:,1,:,:) + resp%sB_full(:,:,:,iband,:,:)
          enddo
          sB_partial_sum_dp = sB_partial_sum &
              * 4.q0 / 3.q0 * piQ**2 * ( echarge / (kmesh%vol*hbarevs)) * (1.q-10 / hbarevs)
          deallocate(sB_partial_sum)
      end select
      ! note that I put the unit factors to get the most out of the accuracy

      if (myid .eq. master) then
        select case (algo%fullOutput)
          case (1)
            allocate(sB_gather(3,3,3,edisp%nband_max,edisp%ispin,kmesh%nkp))
          case (2)
            allocate(sB_gather(3,3,3,edisp%nband_max,edisp%ispin,1))
          case (3)
            allocate(sB_gather(3,3,3,1,edisp%ispin,kmesh%nkp))
        end select
      else
        allocate(sB_gather(1,1,1,1,1,1))
      endif
#ifdef MPI
      select case (algo%fullOutput)
        case (1)
          call MPI_gatherv(sB_partial_sum_dp,(ikend-ikstr+1)*27*edisp%nband_max*edisp%ispin, &
                           MPI_DOUBLE_COMPLEX, sB_gather, rcounts*27*edisp%nband_max*edisp%ispin, &
                           displs*27*edisp%nband_max*edisp%ispin, MPI_DOUBLE_COMPLEX, master, &
                           MPI_COMM_WORLD, mpierr)
        case (2)
          allocate(qrksumarrB(3,3,3,edisp%nband_max,edisp%ispin,1))
          allocate(qiksumarrB(3,3,3,edisp%nband_max,edisp%ispin,1))
          qrksumarrB = 0.q0
          qiksumarrB = 0.q0
          ! we want to keep the k-sum accuracy over all cores
          do ii=1,3
            do ij=1,3
              do ik=1,3
                do iband=1,edisp%nband_max
                  do is=1,edisp%ispin
                    call mpi_reduce_quad(real(sB_partial_sum(ii,ij,ik,iband,is,1)),qrksumarrB(ii,ij,ik,iband,is,1))
                    call mpi_reduce_quad(aimag(sB_partial_sum(ii,ij,ik,iband,is,1)),qiksumarrB(ii,ij,ik,iband,is,1))
                  enddo
                enddo
              enddo
            enddo
          enddo
          sB_gather = cmplx(real(qrksumarrB,8),real(qiksumarrB,8))
          deallocate(qrksumarrB)
          deallocate(qiksumarrB)
        case (3)
          call MPI_gatherv(sB_partial_sum_dp,(ikend-ikstr+1)*27*edisp%ispin, &
                           MPI_DOUBLE_COMPLEX, sB_gather, rcounts*27*edisp%ispin, &
                           displs*27*edisp%ispin, MPI_DOUBLE_COMPLEX, master, &
                           MPI_COMM_WORLD, mpierr)
      end select
#else
      if (algo%fullOutput == 2) then
        sB_gather = sB_partial_sum
      else
        sB_gather = sB_partial_sum_dp
      endif
#endif
      if (allocated(sB_partial_sum)) deallocate(sB_partial_sum)
      if (allocated(sB_partial_sum_dp)) deallocate(sB_partial_sum_dp)


      ! L12B
      select case (algo%fullOutput)
        case (1) ! full
          sumstyle = '/full'
          allocate(aB_partial_sum_dp(3,3,3,edisp%nband_max,edisp%ispin,ikstr:ikend))
        case (2) ! ksum
          sumstyle = '/ksum'
          allocate(aB_partial_sum(3,3,3,edisp%nband_max,edisp%ispin,1))
          allocate(aB_partial_sum_dp(3,3,3,edisp%nband_max,edisp%ispin,1))
          aB_partial_sum = 0.q0
        case (3) ! bsum
          sumstyle = '/bsum'
          allocate(aB_partial_sum(3,3,3,1,edisp%ispin,ikstr:ikend))
          allocate(aB_partial_sum_dp(3,3,3,1,edisp%ispin,ikstr:ikend))
          aB_partial_sum = 0.q0
      end select
      aB_partial_sum_dp = 0.d0
      select case (algo%fullOutput)
        case (1) ! shift the optical range into the full energy band range
          aB_partial_sum_dp(:,:,:,edisp%nbopt_min:edisp%nbopt_max,:,:) = &
              resp%aB_full * 4.q0 / 3.q0 * piQ**2 * ( echarge / (kmesh%vol*hbarevs)) * (1.q-10 / hbarevs)
        case (2) ! shift optical range + momentum sum
          do ik=ikstr,ikend
            aB_partial_sum(:,:,:,edisp%nbopt_min:edisp%nbopt_max,:,1) = &
                  aB_partial_sum(:,:,:,edisp%nbopt_min:edisp%nbopt_max,:,1) + &
                  resp%aB_full(:,:,:,edisp%nbopt_min:edisp%nbopt_max,:,ik) * kmesh%weightQ(ik)
          enddo
          aB_partial_sum = aB_partial_sum &
              * 4.q0 / 3.q0 * piQ**2 * ( echarge / (kmesh%vol*hbarevs)) * (1.q-10 / hbarevs)
        case (3) ! band sum
          do iband=edisp%nbopt_min,edisp%nbopt_max
            aB_partial_sum(:,:,:,1,:,:) = &
                  aB_partial_sum(:,:,:,1,:,:) + resp%aB_full(:,:,:,iband,:,:)
          enddo
          aB_partial_sum_dp = aB_partial_sum &
              * 4.q0 / 3.q0 * piQ**2 * ( echarge / (kmesh%vol*hbarevs)) * (1.q-10 / hbarevs)
          deallocate(aB_partial_sum)
      end select
      ! note that I put the unit factors to get the most out of the accuracy

      if (myid .eq. master) then
        select case (algo%fullOutput)
          case (1)
            allocate(aB_gather(3,3,3,edisp%nband_max,edisp%ispin,kmesh%nkp))
          case (2)
            allocate(aB_gather(3,3,3,edisp%nband_max,edisp%ispin,1))
          case (3)
            allocate(aB_gather(3,3,3,1,edisp%ispin,kmesh%nkp))
        end select
      else
        allocate(aB_gather(1,1,1,1,1,1))
      endif
#ifdef MPI
      select case (algo%fullOutput)
        case (1)
          call MPI_gatherv(aB_partial_sum_dp,(ikend-ikstr+1)*27*edisp%nband_max*edisp%ispin, &
                           MPI_DOUBLE_COMPLEX, aB_gather, rcounts*27*edisp%nband_max*edisp%ispin, &
                           displs*27*edisp%nband_max*edisp%ispin, MPI_DOUBLE_COMPLEX, master, &
                           MPI_COMM_WORLD, mpierr)
        case (2)
          allocate(qrksumarrB(3,3,3,edisp%nband_max,edisp%ispin,1))
          allocate(qiksumarrB(3,3,3,edisp%nband_max,edisp%ispin,1))
          qrksumarrB = 0.q0
          qiksumarrB = 0.q0
          ! we want to keep the k-sum accuracy over all cores
          do ii=1,3
            do ij=1,3
              do ik=1,3
                do iband=1,edisp%nband_max
                  do is=1,edisp%ispin
                    call mpi_reduce_quad(real(aB_partial_sum(ii,ij,ik,iband,is,1)),qrksumarrB(ii,ij,ik,iband,is,1))
                    call mpi_reduce_quad(aimag(aB_partial_sum(ii,ij,ik,iband,is,1)),qiksumarrB(ii,ij,ik,iband,is,1))
                  enddo
                enddo
              enddo
            enddo
          enddo
          aB_gather = cmplx(real(qrksumarrB,8),real(qiksumarrB,8))
          deallocate(qrksumarrB)
          deallocate(qiksumarrB)
        case (3)
          call MPI_gatherv(aB_partial_sum_dp,(ikend-ikstr+1)*27*edisp%ispin, &
                           MPI_DOUBLE_COMPLEX, aB_gather, rcounts*27*edisp%ispin, &
                           displs*27*edisp%ispin, MPI_DOUBLE_COMPLEX, master, &
                           MPI_COMM_WORLD, mpierr)
      end select
#else
      if (algo%fullOutput == 2) then
        aB_gather = aB_partial_sum
      else
        aB_gather = aB_partial_sum_dp
      endif
#endif
      if (allocated(aB_partial_sum)) deallocate(aB_partial_sum)
      if (allocated(aB_partial_sum_dp)) deallocate(aB_partial_sum_dp)

      ! L22B
      select case (algo%fullOutput)
        case (1) ! full
          sumstyle = '/full'
          allocate(xB_partial_sum_dp(3,3,3,edisp%nband_max,edisp%ispin,ikstr:ikend))
        case (2) ! ksum
          sumstyle = '/ksum'
          allocate(xB_partial_sum(3,3,3,edisp%nband_max,edisp%ispin,1))
          allocate(xB_partial_sum_dp(3,3,3,edisp%nband_max,edisp%ispin,1))
          xB_partial_sum = 0.q0
        case (3) ! bsum
          sumstyle = '/bsum'
          allocate(xB_partial_sum(3,3,3,1,edisp%ispin,ikstr:ikend))
          allocate(xB_partial_sum_dp(3,3,3,1,edisp%ispin,ikstr:ikend))
          xB_partial_sum = 0.q0
      end select
      xB_partial_sum_dp = 0.d0
      select case (algo%fullOutput)
        case (1) ! shift the optical range into the full energy band range
          xB_partial_sum_dp(:,:,:,edisp%nbopt_min:edisp%nbopt_max,:,:) = &
              resp%xB_full * 4.q0 / 3.q0 * piQ**2 * ( echarge / (kmesh%vol*hbarevs)) * (1.q-10 / hbarevs)
        case (2) ! shift optical range + momentum sum
          do ik=ikstr,ikend
            xB_partial_sum(:,:,:,edisp%nbopt_min:edisp%nbopt_max,:,1) = &
                  xB_partial_sum(:,:,:,edisp%nbopt_min:edisp%nbopt_max,:,1) + &
                  resp%xB_full(:,:,:,edisp%nbopt_min:edisp%nbopt_max,:,ik) * kmesh%weightQ(ik)
          enddo
          xB_partial_sum = xB_partial_sum &
              * 4.q0 / 3.q0 * piQ**2 * ( echarge / (kmesh%vol*hbarevs)) * (1.q-10 / hbarevs)
        case (3) ! band sum
          do iband=edisp%nbopt_min,edisp%nbopt_max
            xB_partial_sum(:,:,:,1,:,:) = &
                  xB_partial_sum(:,:,:,1,:,:) + resp%xB_full(:,:,:,iband,:,:)
          enddo
          xB_partial_sum_dp = xB_partial_sum &
              * 4.q0 / 3.q0 * piQ**2 * ( echarge / (kmesh%vol*hbarevs)) * (1.q-10 / hbarevs)
          deallocate(xB_partial_sum)
      end select
      ! note that I put the unit factors to get the most out of the accuracy

      if (myid .eq. master) then
        select case (algo%fullOutput)
          case (1)
            allocate(xB_gather(3,3,3,edisp%nband_max,edisp%ispin,kmesh%nkp))
          case (2)
            allocate(xB_gather(3,3,3,edisp%nband_max,edisp%ispin,1))
          case (3)
            allocate(xB_gather(3,3,3,1,edisp%ispin,kmesh%nkp))
        end select
      else
        allocate(xB_gather(1,1,1,1,1,1))
      endif
#ifdef MPI
      select case (algo%fullOutput)
        case (1)
          call MPI_gatherv(xB_partial_sum_dp,(ikend-ikstr+1)*27*edisp%nband_max*edisp%ispin, &
                           MPI_DOUBLE_COMPLEX, xB_gather, rcounts*27*edisp%nband_max*edisp%ispin, &
                           displs*27*edisp%nband_max*edisp%ispin, MPI_DOUBLE_COMPLEX, master, &
                           MPI_COMM_WORLD, mpierr)
        case (2)
          allocate(qrksumarrB(3,3,3,edisp%nband_max,edisp%ispin,1))
          allocate(qiksumarrB(3,3,3,edisp%nband_max,edisp%ispin,1))
          qrksumarrB = 0.q0
          qiksumarrB = 0.q0
          ! we want to keep the k-sum accuracy over all cores
          do ii=1,3
            do ij=1,3
              do ik=1,3
                do iband=1,edisp%nband_max
                  do is=1,edisp%ispin
                    call mpi_reduce_quad(real(xB_partial_sum(ii,ij,ik,iband,is,1)),qrksumarrB(ii,ij,ik,iband,is,1))
                    call mpi_reduce_quad(aimag(xB_partial_sum(ii,ij,ik,iband,is,1)),qiksumarrB(ii,ij,ik,iband,is,1))
                  enddo
                enddo
              enddo
            enddo
          enddo
          xB_gather = cmplx(real(qrksumarrB,8),real(qiksumarrB,8))
          deallocate(qrksumarrB)
          deallocate(qiksumarrB)
        case (3)
          call MPI_gatherv(xB_partial_sum_dp,(ikend-ikstr+1)*27*edisp%ispin, &
                           MPI_DOUBLE_COMPLEX, xB_gather, rcounts*27*edisp%ispin, &
                           displs*27*edisp%ispin, MPI_DOUBLE_COMPLEX, master, &
                           MPI_COMM_WORLD, mpierr)
      end select
#else
      if (algo%fullOutput == 2) then
        xB_gather = xB_partial_sum
      else
        xB_gather = xB_partial_sum_dp
      endif
#endif
      if (allocated(xB_partial_sum)) deallocate(xB_partial_sum)
      if (allocated(xB_partial_sum_dp)) deallocate(xB_partial_sum_dp)


      if (myid .eq. master) then
        write(string,'(I6.6)') info%iStep
        string = 'step/' // trim(string) // "/L11B/" // trim(adjustl(gname)) // sumstyle
        call hdf5_write_data(ifile, string, sB_gather)

        write(string,'(I6.6)') info%iStep
        string = 'step/' // trim(string) // "/L12B/" // trim(adjustl(gname)) // sumstyle
        call hdf5_write_data(ifile, string, aB_gather)

        write(string,'(I6.6)') info%iStep
        string = 'step/' // trim(string) // "/L22B/" // trim(adjustl(gname)) // sumstyle
        call hdf5_write_data(ifile, string, xB_gather)

      endif

      deallocate(sB_gather)
      deallocate(aB_gather)
      deallocate(xB_gather)
    endif ! full output

    ! perform a local summation
    do ik = ikstr,ikend
      do iband = edisp%nbopt_min,edisp%nbopt_max
        resp%sB_sum(:,:,:,:) = resp%sB_sum(:,:,:,:) + resp%sB_full(:,:,:,iband,:,ik) * kmesh%weightQ(ik)
        resp%aB_sum(:,:,:,:) = resp%aB_sum(:,:,:,:) + resp%aB_full(:,:,:,iband,:,ik) * kmesh%weightQ(ik)
        resp%xB_sum(:,:,:,:) = resp%xB_sum(:,:,:,:) + resp%xB_full(:,:,:,iband,:,ik) * kmesh%weightQ(ik)
      enddo
    enddo

    ! perform MPI summation
    qrsarrB = 0.q0
    qraarrB = 0.q0
    qrxarrB = 0.q0
    qisarrB = 0.q0
    qiaarrB = 0.q0
    qixarrB = 0.q0
#ifdef MPI
    do ii=1,3
      do ij=1,3
        do ik=1,3
          do is=1,edisp%ispin
            call mpi_reduce_quad(real(resp%sB_sum(ii,ij,ik,is)),qrsarrB(ii,ij,ik,is))
            call mpi_reduce_quad(aimag(resp%sB_sum(ii,ij,ik,is)),qisarrB(ii,ij,ik,is))

            call mpi_reduce_quad(real(resp%aB_sum(ii,ij,ik,is)),qraarrB(ii,ij,ik,is))
            call mpi_reduce_quad(aimag(resp%aB_sum(ii,ij,ik,is)),qiaarrB(ii,ij,ik,is))

            call mpi_reduce_quad(real(resp%xB_sum(ii,ij,ik,is)),qrxarrB(ii,ij,ik,is))
            call mpi_reduce_quad(aimag(resp%xB_sum(ii,ij,ik,is)),qixarrB(ii,ij,ik,is))
          enddo
        enddo
      enddo
    enddo
#else
    qrsarrB = real(resp%sB_sum)
    qisarrB = aimag(resp%sB_sum)

    qraarrB = real(resp%aB_sum)
    qiaarrB = aimag(resp%aB_sum)

    qrxarrB = real(resp%xB_sum)
    qixarrB = aimag(resp%xB_sum)
#endif

    qrsarrB = qrsarrB * 4.q0 / 3.q0 * piQ**2 * ( echarge / (kmesh%vol*hbarevs)) * (1.q-10 / hbarevs)
    qisarrB = qisarrB * 4.q0 / 3.q0 * piQ**2 * ( echarge / (kmesh%vol*hbarevs)) * (1.q-10 / hbarevs)
    qraarrB = qraarrB * 4.q0 / 3.q0 * piQ**2 * ( echarge / (kmesh%vol*hbarevs)) * (1.q-10 / hbarevs)
    qiaarrB = qiaarrB * 4.q0 / 3.q0 * piQ**2 * ( echarge / (kmesh%vol*hbarevs)) * (1.q-10 / hbarevs)
    qrxarrB = qrxarrB * 4.q0 / 3.q0 * piQ**2 * ( echarge / (kmesh%vol*hbarevs)) * (1.q-10 / hbarevs)
    qixarrB = qixarrB * 4.q0 / 3.q0 * piQ**2 * ( echarge / (kmesh%vol*hbarevs)) * (1.q-10 / hbarevs)

    ! should work=?
    if (myid .eq. master) then
      ! gather the data in the arrays
      resp%sB_sum_range(:,:,:,:,info%iStep) = cmplx(real(qrsarrB,8),real(qisarrB,8))
      resp%aB_sum_range(:,:,:,:,info%iStep) = cmplx(real(qraarrB,8),real(qiaarrB,8))
      resp%xB_sum_range(:,:,:,:,info%iStep) = cmplx(real(qrxarrB,8),real(qixarrB,8))

      ! output at the last temperature step
      if ((algo%step_dir==1 .and. info%iStep==algo%steps) .or. (algo%step_dir==-1 .and. info%iStep==1)) then
        string = "/L11B/" // trim(adjustl(gname)) // "/sum"
        call hdf5_write_data(ifile, string, resp%sB_sum_range)
        string = "/L12B/" // trim(adjustl(gname)) // "/sum"
        call hdf5_write_data(ifile, string, resp%aB_sum_range)
        string = "/L22B/" // trim(adjustl(gname)) // "/sum"
        call hdf5_write_data(ifile, string, resp%xB_sum_range)
      endif
    endif

    deallocate(qrsarrB)
    deallocate(qraarrB)
    deallocate(qrxarrB)
    deallocate(qisarrB)
    deallocate(qiaarrB)
    deallocate(qixarrB)
  endif ! Boutput

  if (myid.eq.master) then
    call hdf5_close_file(ifile)
  endif

  deallocate(qrsarr)
  deallocate(qraarr)
  deallocate(qrxarr)
  deallocate(qisarr)
  deallocate(qiaarr)
  deallocate(qixarr)

end subroutine

end module Moutput
