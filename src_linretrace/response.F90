module Mresponse
  use Mmpi_org
  use Mtypes
  use Mparams
  implicit none

  interface allocate_response
    module procedure dpresp_alloc, qpresp_alloc
  end interface

  interface calc_polygamma
    module procedure calc_polygamma_D, calc_polygamma_Q
  end interface calc_polygamma

contains

!subroutine calc_response(mu, iT, drhodT, mesh, ek, thdr, sct, dresp, dderesp, dinter, respBl, qresp)
!  implicit none
!  type(kpointmesh)        :: mesh
!  type(energydisp)        :: ek
!  type(tetramesh)         :: thdr
!  type(scattering)          :: sct
!  type(dp_resp)           :: dresp   !intraband response
!  type(dp_respinter)      :: dderesp !derivatives of intraband conductivity
!  type(dp_respinter)      :: dinter  !interband response
!  type(dp_resp)           :: respBl  !Boltzmann response
!  type(qp_resp), optional :: qresp
!  real(8) :: mu
!  real(8) :: drhodT(sct%nT)
!  integer :: iT,ib
!  integer :: itet, ik, ikk
!  integer :: ix,iy

!  complex(8),external  :: wpsipg
!  complex(16),external :: wpsipghp

!  !initialise the datatype variables
!  ! eM: incredibly the qp_resp type seems to be there also when
!  ! it is not passed by the summoning routine. I think this is a
!  ! gfortran compiler glitch in the generation of the
!  ! implicit interface, using ldebug instead of the
!  ! intrinsic fortran present(qresp) solves the issue

!  ! initialize the already allocated arrays to 0
!  call initresp (algo%lBfield, dresp)
!  call initresp (algo%lBfield, respBl)
!  call initresp (.false., dderesp)
!  call initresp (.false., dinter)
!  if (.not. algo%ldebug) then
!     call initresp_qp (algo%lBfield, qresp)
!  endif

!  ! outer k-loop
!  ! (quantities inside the tetrahedra will be interpolated over)
!  ! eM: I decided to treat the tetrahedron and the regular k-mesh cases
!  ! independently because later different parallelisation strategies can be
!  ! devised for the two approaches

!  if (algo%ltetra ) then
!     ! here the loop over bands is hidden inside the tetrahedron response routines
!     do itet=iqstr,iqend
!        !intraband transitions
!        call respintet (mu, iT, itet, thdr, ek, sct, dresp)
!        call respintet_Bl (mu, iT, itet, thdr, ek, sct, respBl)
!        if (.not.algo%ldebug) then
!           call respintet_qp (mu, iT, itet, thdr, ek, sct, qresp)
!        endif

!        !evaluate the derivatives of the response functions
!        !if the chemical potential is fixed use a semplified kernel for the
!        !derivatives (assuming also gamma to be not T dependent)
!        !if ((sct%Tstar == 0.0d0) .or. (sct%Tflat == 0.0d0)) then
!           if (algo%ltbind .and. (algo%imurestart==2)) then
!              !since one has to evaluate the derivative of mu then the first point must be skipped
!              if (iT < sct%nT) call resdertet_symm(mu, iT, itet, thdr, ek, sct, dderesp)
!           else
!              !since one has to evaluate the derivative of mu then the first point must be skipped
!              if (iT < sct%nT) call resdertet(iT, itet, thdr, ek, sct, dderesp)
!           endif
!        !endif

!        !interband transitions
!        if (algo%ldebug) then
!           !!!!!!!!!!!!! TEST
!           !if((iT==sct%nT) .and. (itet==1)) write(*,*)'test for 2-band symmetrical SC, check input parameters!!'
!           !call respintert_symm(mu, iT, itet, thdr, algo, ek, sct, dinter)
!           !!!!!!!!!!!!! TEST END
!           call respintert(mu, iT, itet, thdr, ek, sct, dinter)
!        else
!           call respintert(mu, iT, itet, thdr, ek, sct, dinter)
!        endif


!        !interpolate within tetrahedra intraband response
!        call interptra_re (iT, itet, mu, .false., mesh, ek, thdr, sct, dresp)
!        call interptra_re (iT, itet, mu, .true. , mesh, ek, thdr, sct, respBl)
!        if (.not.algo%ldebug) then
!           call interptra_re (iT, itet, mu, .false., mesh, ek, thdr, sct, dresp, qresp)
!        endif

!        !interpolate within tetrahedra intraband response derivatives
!        call interptra_re (iT, itet, mu, .false., mesh, ek, thdr, sct, dderesp)

!        !interpolate within tetrahedra interband response
!        call interptra_re (iT, itet, mu, .false., mesh, ek, thdr, sct, dinter)

!        ! add to tetrahedra-summed response functions (bands have been traced over in interptra_re )
!        !
!        ! functions had to be multiplied by 2 for spin multiplicity

!        do ix=1,lat%nalpha
!           do iy=ix,lat%nalpha
!              dresp%s_tot(ix,iy)=dresp%s_tot(ix,iy)+dresp%s_tet(ix,iy)*2.0d0
!              dresp%a_tot(ix,iy)=dresp%a_tot(ix,iy)+dresp%a_tet(ix,iy)*2.0d0

!              dderesp%s_tot(ix,iy)=dderesp%s_tot(ix,iy)+dderesp%s_tet(ix,iy)*2.0d0
!              dderesp%a_tot(ix,iy)=dderesp%a_tot(ix,iy)+dderesp%a_tet(ix,iy)*2.0d0

!              dinter%s_tot(ix,iy)=dinter%s_tot(ix,iy)+dinter%s_tet(ix,iy)*2.0d0
!              dinter%a_tot(ix,iy)=dinter%a_tot(ix,iy)+dinter%a_tet(ix,iy)*2.0d0

!              if (algo%lBfield .and. algo%ltbind ) then
!                 dresp%sB_tot(ix,iy)=dresp%sB_tot(ix,iy)+dresp%sB_tet(ix,iy)*2.0d0
!                 dresp%aB_tot(ix,iy)=dresp%aB_tot(ix,iy)+dresp%aB_tet(ix,iy)*2.0d0
!              endif

!              respBl%s_tot(ix,iy)=respBl%s_tot(ix,iy)+respBl%s_tet(ix,iy)*2.0d0
!              respBl%a_tot(ix,iy)=respBl%a_tot(ix,iy)+respBl%a_tet(ix,iy)*2.0d0
!              if (algo%lBfield .and. algo%ltbind ) then
!                 respBl%sB_tot(ix,iy)=respBl%sB_tot(ix,iy)+respBl%sB_tet(ix,iy)*2.0d0
!                 respBl%aB_tot(ix,iy)=respBl%aB_tot(ix,iy)+respBl%aB_tet(ix,iy)*2.0d0
!              endif

!              if (.not.algo%ldebug) then
!                 qresp%s_tot(ix,iy)=qresp%s_tot(ix,iy)+qresp%s_tet(ix,iy)*2.0q0
!                 qresp%a_tot(ix,iy)=qresp%a_tot(ix,iy)+qresp%a_tet(ix,iy)*2.0q0
!                 if (algo%lBfield .and. algo%ltbind) then
!                    qresp%sB_tot(ix,iy)=qresp%sB_tot(ix,iy)+qresp%sB_tet(ix,iy)*2.0q0
!                    qresp%aB_tot(ix,iy)=qresp%aB_tot(ix,iy)+qresp%aB_tet(ix,iy)*2.0q0
!                 endif
!              endif
!           enddo !iy
!        enddo    !ix
!     enddo ! loop over tetrahedra

!#ifdef MPI
!     call MPI_REDUCE(MPI_IN_PLACE, dresp%s_tot, 9, MPI_REAL8, MPI_SUM, master, MPI_COMM_WORLD, mpierr)
!     call MPI_REDUCE(MPI_IN_PLACE, dresp%a_tot, 9, MPI_REAL8, MPI_SUM, master, MPI_COMM_WORLD, mpierr)
!     if (algo%lBfield .and. algo%ltbind ) then
!        call MPI_REDUCE(MPI_IN_PLACE, dresp%sB_tot, 9, MPI_REAL8, MPI_SUM, master, MPI_COMM_WORLD, mpierr)
!        call MPI_REDUCE(MPI_IN_PLACE, dresp%aB_tot, 9, MPI_REAL8, MPI_SUM, master, MPI_COMM_WORLD, mpierr)
!     endif
!     !derivative
!     call MPI_REDUCE(MPI_IN_PLACE, dderesp%s_tot, 9, MPI_REAL8, MPI_SUM, master, MPI_COMM_WORLD, mpierr)
!     call MPI_REDUCE(MPI_IN_PLACE, dderesp%a_tot, 9, MPI_REAL8, MPI_SUM, master, MPI_COMM_WORLD, mpierr)
!     !interband
!     call MPI_REDUCE(MPI_IN_PLACE, dinter%s_tot, 9, MPI_REAL8, MPI_SUM, master, MPI_COMM_WORLD, mpierr)
!     call MPI_REDUCE(MPI_IN_PLACE, dinter%a_tot, 9, MPI_REAL8, MPI_SUM, master, MPI_COMM_WORLD, mpierr)
!     !Boltzmann
!     call MPI_REDUCE(MPI_IN_PLACE, respBl%s_tot, 9, MPI_REAL8, MPI_SUM, master, MPI_COMM_WORLD, mpierr)
!     call MPI_REDUCE(MPI_IN_PLACE, respBl%a_tot, 9, MPI_REAL8, MPI_SUM, master, MPI_COMM_WORLD, mpierr)
!     if (algo%lBfield .and. algo%ltbind ) then
!        call MPI_REDUCE(MPI_IN_PLACE, respBl%sB_tot, 9, MPI_REAL8, MPI_SUM, master, MPI_COMM_WORLD, mpierr)
!        call MPI_REDUCE(MPI_IN_PLACE, respBl%aB_tot, 9, MPI_REAL8, MPI_SUM, master, MPI_COMM_WORLD, mpierr)
!     endif
!     !intraband QP
!     if (.not.algo%ldebug) then
!        do ix=1,3
!           do iy=ix,3
!              call mpi_reduce_quad(qresp%s_tot(ix,iy),qresp%s_tot(ix,iy))
!              call mpi_reduce_quad(qresp%a_tot(ix,iy),qresp%a_tot(ix,iy))
!           enddo
!        enddo
!        if (algo%lBfield .and. algo%ltbind ) then
!           do ix=1,3
!              do iy=ix,3
!                 call mpi_reduce_quad(qresp%sB_tot(ix,iy),qresp%sB_tot(ix,iy))
!                 call mpi_reduce_quad(qresp%aB_tot(ix,iy),qresp%aB_tot(ix,iy))
!              enddo
!           enddo
!        endif
!     endif
!#endif

!  else ! no tetrahedron method

!     ! do the parallelized loop
!     do ik = iqstr,iqend
!     ! evaluate the trace over bands at each specific k-point of the mesh
!        call respinkm (mu, iT, ik, ek, sct, dresp)
!        call respinkm_Bl (mu, iT, ik, ek, sct, respBl)
!        if (.not.algo%ldebug) then
!           call respinkm_qp (mu, iT, ik, ek, sct, qresp)
!        endif
!        call respinterkm(mu, iT, ik, ek, sct, dinter)
!     enddo

!        !evaluate the derivatives of the response functions
!        !if the chemical potential is fixed use a semplified kernel for the
!        !derivatives (assuming also gamma to be not T dependent)
!        !if ((sct%Tstar == 0.0d0) .or. (sct%Tflat == 0.0d0)) then
!        !   if (algo%ltbind .and. (algo%imurestart==2)) then
!        !      if (iT < sct%nT) call resderkm_symm(mu, iT, ik, ek, sct, dderesp)
!        !   else
!        !      !since one has to evaluate the derivative of mu then the first point must be skipped
!        !      if (iT < sct%nT) call resderkm(iT, ik, ek, sct, dderesp)
!        !   endif
!        !endif


!     ! mP note:
!     ! now we have all the data for each optical band and k-point
!     ! we do total summation in two steps, since we want to have the
!     ! band-resolved quantities as well

!     ! perform partial k-sum on each core
!     do ix=1,lat%nalpha
!        do iy=ix,lat%nalpha
!           do ib=1,ek%nband_max
!              if(ib<ek%nbopt_min) cycle
!              if(ib>ek%nbopt_max) cycle

!              ! partial k-sum (in the single-core application this is the whole BZ)
!              do ik=iqstr,iqend
!                 ! ikk = symm%symop_id(1,ik) ! k-point of the corresponding saved data
!                 ikk = ik
!                 !multiply by a factor that includes spin multiplicity and the term
!                 !beta/gamma (beta^2/gamma) for s (a) in presence of magnetic field gamma --> gamma^2
!                 !for the Boltzmann response beta --> 1
!                 if (allocated(sct%ykb)) then
!                    dresp%gamma=real(ek%z(ikk,ib)*(sct%gam(iT)+sct%ykb(iT,ikk,ib)),8)
!                    if (.not.algo%ldebug) qresp%gamma=real(ek%z(ikk,ib)*(sct%gam(iT)+sct%ykb(iT,ikk,ib)),16)
!                 else
!                    dresp%gamma=real(ek%z(ikk,ib)*sct%gam(iT),8)
!                    if (.not.algo%ldebug) qresp%gamma=real(ek%z(ikk,ib)*sct%gam(iT),16)
!                 endif
!                 dresp%s_local(ib,ix,iy) = dresp%s_local(ib,ix,iy) &
!                                         + dresp%s_tmp(ik,ib,ix,iy)*2.0d0*beta/dresp%gamma * mesh%weight(ikk)
!                 dresp%a_local(ib,ix,iy) = dresp%a_local(ib,ix,iy) &
!                                         + dresp%a_tmp(ik,ib,ix,iy)*2.0d0*(beta**2)/dresp%gamma * mesh%weight(ikk)
!                 dderesp%s_local(ib,ix,iy) = dderesp%s_local(ib,ix,iy) &
!                                           + dderesp%s_tmp(ik,ib,ix,iy)*2.0d0 * mesh%weight(ikk)
!                 dderesp%a_local(ib,ix,iy) = dderesp%a_local(ib,ix,iy) &
!                                           + dderesp%a_tmp(ik,ib,ix,iy)*2.0d0 * mesh%weight(ikk)
!                 dinter%s_local(ib,ix,iy) = dinter%s_local(ib,ix,iy) &
!                                          + dinter%s_tmp(ik,ib,ix,iy)*2.0d0 * mesh%weight(ikk)
!                 dinter%a_local(ib,ix,iy) = dinter%a_local(ib,ix,iy) &
!                                          + dinter%a_tmp(ik,ib,ix,iy)*2.0d0 * mesh%weight(ikk)
!                 if (algo%lBfield .and. algo%ltbind ) then
!                    dresp%sB_local(ib,ix,iy) = dresp%sB_local(ib,ix,iy) &
!                                             + dresp%sB_tmp(ik,ib,ix,iy)*2.0d0*beta/(dresp%gamma**2) * mesh%weight(ikk)
!                    dresp%aB_local(ib,ix,iy) = dresp%aB_local(ib,ix,iy) &
!                                             + dresp%aB_tmp(ik,ib,ix,iy)*2.0d0*(beta/dresp%gamma)**2 * mesh%weight(ikk)
!                 endif

!                 respBl%s_local(ib,ix,iy) = respBl%s_local(ib,ix,iy) &
!                                          + respBl%s_tmp(ik,ib,ix,iy)*2.0d0/dresp%gamma * mesh%weight(ikk)
!                 respBl%a_local(ib,ix,iy) = respBl%a_local(ib,ix,iy) &
!                                          + respBl%a_tmp(ik,ib,ix,iy)*2.0d0/dresp%gamma * mesh%weight(ikk)
!                 if (algo%lBfield .and. algo%ltbind ) then
!                    respBl%sB_local(ib,ix,iy) = respBl%sB_local(ib,ix,iy) &
!                                              + respBl%sB_tmp(ik,ib,ix,iy)*2.0d0/(dresp%gamma**2) * mesh%weight(ikk)
!                    respBl%aB_local(ib,ix,iy) = respBl%aB_local(ib,ix,iy) &
!                                              + respBl%aB_tmp(ik,ib,ix,iy)*2.0d0/(dresp%gamma**2) * mesh%weight(ikk)
!                 endif

!                 if (.not.algo%ldebug) then
!                    qresp%s_local(ib,ix,iy) = qresp%s_local(ib,ix,iy) &
!                                            + qresp%s_tmp(ik,ib,ix,iy)*2.0q0*betaQ/qresp%gamma * mesh%weight(ikk)
!                    qresp%a_local(ib,ix,iy) = qresp%a_local(ib,ix,iy) &
!                                            + qresp%a_tmp(ik,ib,ix,iy)*2.0q0*(betaQ**2)/qresp%gamma * mesh%weight(ikk)
!                    if (algo%lBfield .and. algo%ltbind ) then
!                       qresp%sB_local(ib,ix,iy) = qresp%sB_local(ib,ix,iy) &
!                                                + qresp%sB_tmp(ik,ib,ix,iy)*2.0q0*betaQ/(qresp%gamma**2) * mesh%weight(ikk)
!                       qresp%aB_local(ib,ix,iy) = qresp%aB_local(ib,ix,iy) &
!                                                + qresp%aB_tmp(ik,ib,ix,iy)*2.0q0*(betaQ/qresp%gamma)**2 * mesh%weight(ikk)
!                    endif
!                 endif
!              enddo
!           enddo !nbands
!        enddo !iy
!     enddo !ix

!     ! perform the total k-summation over all cores
!     ! this is gathered only on the master node.
!#ifdef MPI
!     call MPI_REDUCE(dresp%s_local, dresp%s, ek%nband_max*9, MPI_REAL8, MPI_SUM, master, MPI_COMM_WORLD, mpierr)
!     call MPI_REDUCE(dresp%a_local, dresp%a, ek%nband_max*9, MPI_REAL8, MPI_SUM, master, MPI_COMM_WORLD, mpierr)
!     if (algo%lBfield .and. algo%ltbind ) then
!        call MPI_REDUCE(dresp%sB_local, dresp%sB, ek%nband_max*9, MPI_REAL8, MPI_SUM, master, MPI_COMM_WORLD, mpierr)
!        call MPI_REDUCE(dresp%aB_local, dresp%aB, ek%nband_max*9, MPI_REAL8, MPI_SUM, master, MPI_COMM_WORLD, mpierr)
!     endif
!     !derivative
!     call MPI_REDUCE(dderesp%s_local, dderesp%s, ek%nband_max*9, MPI_REAL8, MPI_SUM, master, MPI_COMM_WORLD, mpierr)
!     call MPI_REDUCE(dderesp%a_local, dderesp%a, ek%nband_max*9, MPI_REAL8, MPI_SUM, master, MPI_COMM_WORLD, mpierr)
!     !interband
!     call MPI_REDUCE(dinter%s_local, dinter%s, ek%nband_max*9, MPI_REAL8, MPI_SUM, master, MPI_COMM_WORLD, mpierr)
!     call MPI_REDUCE(dinter%a_local, dinter%a, ek%nband_max*9, MPI_REAL8, MPI_SUM, master, MPI_COMM_WORLD, mpierr)
!     !Boltzmann
!     call MPI_REDUCE(respBl%s_local, respBl%s, ek%nband_max*9, MPI_REAL8, MPI_SUM, master, MPI_COMM_WORLD, mpierr)
!     call MPI_REDUCE(respBl%a_local, respBl%a, ek%nband_max*9, MPI_REAL8, MPI_SUM, master, MPI_COMM_WORLD, mpierr)
!     if (algo%lBfield .and. algo%ltbind ) then
!        call MPI_REDUCE(respBl%sB_local, respBl%sB, ek%nband_max*9, MPI_REAL8, MPI_SUM, master, MPI_COMM_WORLD, mpierr)
!        call MPI_REDUCE(respBl%aB_local, respBl%aB, ek%nband_max*9, MPI_REAL8, MPI_SUM, master, MPI_COMM_WORLD, mpierr)
!     endif

!     ! intraband QP
!     ! here we unfortunately have to work with the slow calls
!     ! because im not interested in writing a generic quad version of mpi reduce
!     if (.not.algo%ldebug) then
!        do ib=ek%nbopt_min,ek%nbopt_max
!           do ix=1,lat%nalpha
!              do iy=ix,lat%nalpha
!                 call mpi_reduce_quad(qresp%s_local(ib,ix,iy),qresp%s(ib,ix,iy))
!                 call mpi_reduce_quad(qresp%a_local(ib,ix,iy),qresp%a(ib,ix,iy))
!                 if (algo%lBfield .and. algo%ltbind ) then
!                    call mpi_reduce_quad(qresp%sB_local(ib,ix,iy),qresp%sB(ib,ix,iy))
!                    call mpi_reduce_quad(qresp%aB_local(ib,ix,iy),qresp%aB(ib,ix,iy))
!                 endif
!              enddo
!           enddo
!        enddo
!     endif
!#else
!     ! in the single core application the local is simply the total one
!     drsep%s = dresp%s_local
!     dresp%a = dresp%a_local
!     if (algo%lBfield .and. algo%ltbind) then
!        dresp%sB = dresp%sB_local
!        dresp%aB = dresp%aB_local
!        respBl%sB = respBl%sB_local
!        respBl%aB = respBl%aB_local
!     endif
!     ddersp%s = ddersp%s_local
!     ddresp%a = ddersp%a_local
!     dinter%s = dinter%s_local
!     dinter%a = dinter%a_local
!     respBl%s = respBl%s_local
!     respBl%a = respBl%a_local
!#endif

!     ! perform the band summation
!     if (myid.eq.master) then
!        do ib=1,ek%nband_max
!           if (ib < ek%nbopt_min) cycle
!           if (ib > ek%nbopt_max) cycle

!           dresp%s_tot(:,:)  = dresp%s_tot(:,:) + dresp%s(ib,:,:)
!           dresp%a_tot(:,:)  = dresp%a_tot(:,:) + dresp%a(ib,:,:)
!           !same treatment for the derivatives...
!           dderesp%s_tot(:,:)= dderesp%s_tot(:,:) + dderesp%s(ib,:,:)
!           dderesp%a_tot(:,:)= dderesp%a_tot(:,:) + dderesp%a(ib,:,:)
!           !and for the interband response
!           dinter%s_tot(:,:)= dinter%s_tot(:,:) + dinter%s(ib,:,:)
!           dinter%a_tot(:,:)= dinter%a_tot(:,:) + dinter%a(ib,:,:)
!           if (algo%lBfield .and. algo%ltbind ) then
!              dresp%sB_tot(:,:)= dresp%sB_tot(:,:) + dresp%sB(ib,:,:)
!              dresp%aB_tot(:,:)= dresp%aB_tot(:,:) + dresp%aB(ib,:,:)
!           endif

!           respBl%s_tot(:,:)= respBl%s_tot(:,:) + respBl%s(ib,:,:)
!           respBl%a_tot(:,:)= respBl%a_tot(:,:) + respBl%a(ib,:,:)
!           if (algo%lBfield .and. algo%ltbind ) then
!              respBl%sB_tot(:,:)= respBl%sB_tot(:,:) + dresp%sB(ib,:,:)
!              respBl%aB_tot(:,:)= respBl%aB_tot(:,:) + dresp%aB(ib,:,:)
!           endif

!           if (.not.algo%ldebug) then
!              qresp%s_tot(:,:)= qresp%s_tot(:,:) + qresp%s(ib,:,:)
!              qresp%a_tot(:,:)= qresp%a_tot(:,:) + qresp%a(ib,:,:)
!              if (algo%lBfield .and. algo%ltbind ) then
!                 qresp%sB_tot(:,:)= qresp%sB_tot(:,:) + qresp%sB(ib,:,:)
!                 qresp%aB_tot(:,:)= qresp%aB_tot(:,:) + qresp%aB(ib,:,:)
!              endif
!           endif

!        enddo !over bands
!     endif

!     ! band-resolved quantities in resp%(s/a/sB/aB)
!     ! band-summed quantities   in resp%(s_tot/a_tot/sB_tot/aB_tot)

!  endif !ltetra

!  if (myid.eq.master) then
!  ! At this stage solve the equation (d^2 rho)/(d beta^2) =0
!  ! since the following multiplications in globfac affect equally the
!  ! conductivity and its derivatives there is  no impact on the
!  ! resulting T* temperature (T at which the resistivity has an inflection)
!     !if (sct%Tstar == 0.0d0) then
!        if (iT < sct%nT) call findrhoflex (iT, dresp, dderesp, sct)
!     !endif
!     if (sct%Tflat == 0.0d0) then
!        if (iT < sct%nT) call findrhoflat (iT, dresp, dderesp, sct)
!     endif
!     !if ((iT<sct%nT) .and. (iT>1)) write(*,*)'calling finddrhomax'
!     if ((iT<sct%nT) .and. (iT>1)) call finddrhomax (iT, dresp, dderesp, sct, drhodT)

!  ! At this point the values in the response datatypes should be consistent (in terms of prefactors/dimensionality)
!  ! there are some global factors missing that are taken care of in the following routine
!     if (algo%ldebug) then
!        call globfac(mesh, dresp)
!     else
!        call globfac(mesh, dresp, qresp)
!     endif
!     call globfac(mesh, respBl)
!     call globfac(mesh, dinter)

!  ! The derived variables (Seebeck, Nernst, R_H) are computed in the following routine
!     if (algo%ldebug) then
!        call derresp(dresp)
!     else
!        call derresp(dresp, qresp)
!     endif
!     call derresp(respBl)
!     call derresp(dinter)

!     if (algo%ldebug) then
!        call wrtresp(iT, sct, dresp, dinter, respBl)
!     else
!        call wrtresp(iT, sct, dresp, dinter, respBl, qresp)
!     endif
!  endif

!end subroutine calc_response

!  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  ! INITRESP
!  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  ! This subroutine initialises the datatypes
!  !
!  subroutine initresp (lBfield, dresp)
!    implicit none
!    logical :: lBfield
!    class(dp_resp) :: dresp

!    dresp%s_tot = 0.d0; dresp%s = 0.d0; dresp%s_tmp = 0.d0; dresp%s_local = 0.d0
!    dresp%a_tot = 0.d0; dresp%a = 0.d0; dresp%a_tmp = 0.d0; dresp%a_local = 0.d0
!    if (lBfield) then
!       dresp%sB_tot = 0.0d0; dresp%sB = 0.d0; dresp%sB_tmp = 0.d0; dresp%sB_local = 0.d0
!       dresp%aB_tot = 0.0d0; dresp%aB = 0.d0; dresp%aB_tmp = 0.d0; dresp%aB_local = 0.d0
!    endif
!  end subroutine initresp

!  subroutine initresp_qp (lBfield, dresp)
!    implicit none
!    logical :: lBfield
!    class(qp_resp) :: dresp

!    dresp%s_tot = 0.q0; dresp%s = 0.q0; dresp%s_tmp = 0.q0; dresp%s_local = 0.q0
!    dresp%a_tot = 0.q0; dresp%a = 0.q0; dresp%a_tmp = 0.q0; dresp%a_local = 0.q0
!    if (lBfield) then
!       dresp%sB_tot = 0.0q0; dresp%sB = 0.q0; dresp%sB_tmp = 0.q0; dresp%sB_local = 0.q0
!       dresp%aB_tot = 0.0q0; dresp%aB = 0.q0; dresp%aB_tmp = 0.q0; dresp%aB_local = 0.q0
!    endif
!  end subroutine initresp_qp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! RESPINTET
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! This routine (and the following one in QP)
!! evaluate the conduction kernel, and the
!! full response functions on the vertices of a
!! tetrahedron in the KUBO formalism
!!
!subroutine respintet(mu, iT, itet, thdr, ek, sct, resp)
!  implicit none
!  type (dp_resp)    :: resp ! dynamical datatype allow only for an inclusion through extension of the parent type,
!                         ! I could have declared a dummy class pointer that would be assigned to either dp or qp response type.
!                         ! To do so, the varaibles defined in the individual types must have had different names and since I
!                         ! REALLY dislike having that extra Q for each varible I decided to stick to static datatypes.
!  type (tetramesh)  :: thdr
!  type (energydisp) :: ek
!  type (scattering)   :: sct
!  real(8), intent(in) :: mu
!  integer, intent(in) :: iT
!  integer, intent(in) :: itet
!  integer :: iband, ipg
!  integer :: ix,iy
!  integer :: ik, ikk, iktet
!  complex(8),external  :: wpsipg
!  complex(16),external :: wpsipghp
!!local variables
!  real(8), allocatable :: s_tmp_tetra(:,:,:,:),  a_tmp_tetra(:,:,:,:)
!  real(8), allocatable :: sB_tmp_tetra(:,:,:,:), aB_tmp_tetra(:,:,:,:)
!  real(8) :: Mopt(6, ek%nbopt_min:ek%nbopt_max, ek%nbopt_min:ek%nbopt_max)

!  !allocation
!  if(.not. allocated(s_tmp_tetra)) allocate(s_tmp_tetra(4,ek%nband_max,3,3))
!  if(.not. allocated(a_tmp_tetra)) allocate(a_tmp_tetra(4,ek%nband_max,3,3))
!  if (algo%lBfield .and. algo%ltbind) then
!     if(.not. allocated(sB_tmp_tetra)) allocate(sB_tmp_tetra(4,ek%nband_max,3,3))
!     if(.not. allocated(aB_tmp_tetra)) allocate(aB_tmp_tetra(4,ek%nband_max,3,3))
!  endif
!  !initialisation
!  s_tmp_tetra=0.0d0 ; a_tmp_tetra=0.0d0
!  if (algo%lBfield .and. algo%ltbind) then
!     sB_tmp_tetra=0.0d0 ; aB_tmp_tetra=0.0d0
!  endif

!   do iktet=1,4  !loop over corners of the tetrahedron

!      ik = thdr%idtet(iktet,itet)
!      ! ikk = symm%symop_id(1,ik)
!      ! call getmopt(ek, ik, Mopt, .false.) ! intra

!      do iband=1,ek%nband_max !loop over bands (these will be traced over)

!         if (iband < ek%nbopt_min) cycle
!         if (iband > ek%nbopt_max) cycle
!         resp%z=real(ek%z(ikk,iband),8)
!         if (allocated(sct%ykb)) then
!            resp%gamma=resp%z*real(sct%gam(iT)*sct%ykb(iT,ikk,iband),8)
!         else
!            resp%gamma=resp%z*real(sct%gam(iT),8)
!         endif
!         resp%aqp=resp%z*real(ek%band(ikk,iband)-mu,8)
!         ! pre-compute all needed digamma functions
!         resp%zarg=0.5d0+beta2p*(ci*resp%aqp+resp%gamma)
!         do ipg=1,3 ! XXX need 0 for alphaxy ????
!            resp%ctmp=wpsipg(resp%zarg,ipg)
!            resp%RePolyGamma(ipg)=real(resp%ctmp,8)
!            resp%ImPolyGamma(ipg)=imag(resp%ctmp)
!         enddo

!         ! compute transport kernels (omega-part)
!         !
!         resp%tmp=resp%z**2 / (4.d0*pi**3) ! missing: beta/gamma (multiplied later to keep numbers reasonable here)
!         resp%s_ker = resp%tmp * (resp%RePolyGamma(1) - resp%gamma*beta2p * resp%RePolyGamma(2) )
!         resp%a_ker = resp%tmp * ( resp%aqp * resp%RePolyGamma(1) - resp%aqp*resp%gamma*beta2p*resp%RePolyGamma(2) &
!                    - resp%gamma**2.d0 * beta2p * resp%ImPolyGamma(2) )


!         resp%tmp=resp%tmp*3.d0*resp%z/(4.d0*pi) ! additionally missing: 1/gamma (multiplied later) XXX

!         if(algo%lBfield .and. algo%ltbind) then
!            resp%sB_ker = resp%tmp*(-resp%RePolyGamma(1)-resp%gamma*beta2p*resp%RePolyGamma(2)-(beta2p*resp%gamma)**2.d0/3.d0 &
!                          * resp%RePolyGamma(3))

!            resp%aB_ker=resp%tmp*(resp%aqp*resp%RePolyGamma(1)-resp%aqp*beta2p*resp%gamma*resp%RePolyGamma(2)+resp%gamma**2/3.d0 &
!                 * beta2p * resp%ImPolyGamma(2) &
!                 - resp%aqp*resp%gamma**2.d0 / 3.d0 * beta2p**2.d0 * resp%ImPolyGamma(3) &
!                 + resp%gamma**3.d0 / 3.d0 * beta2p**2.d0 * resp%RePolyGamma(3) )
!          endif

!         ! B = 0
!         !tmp=vka(ik,ib,ix)*vka(ik,ib,iy)
!         do ix=1,lat%nalpha

!            if (algo%ltbind) then
!               resp%tmp=ek%Mopt(ix,thdr%idtet(ik,itet), iband, iband)*ek%Mopt(ix,thdr%idtet(ik,itet), iband, iband)
!            else
!               ! resp%tmp=ek%Mopt(ix,thdr%idtet(ik,itet), iband, iband) !the optical matrix elements given by Wien2k are squared already
!               resp%tmp=Mopt(ix,iband,iband)
!            endif

!            s_tmp_tetra(iktet,iband,ix,ix)=resp%s_ker * resp%tmp
!            a_tmp_tetra(iktet,iband,ix,ix)=resp%a_ker * resp%tmp

!            do iy=ix+1,lat%nalpha
!               ! resp%tmp=ek%Mopt(ix+iy+1,thdr%idtet(ik,itet), iband, iband) !the optical matrix elements given by Wien2k are squared already
!               resp%tmp=Mopt(ix+iy+1,iband,iband)
!               s_tmp_tetra(iktet,iband,ix,iy)=resp%s_ker * resp%tmp
!               a_tmp_tetra(iktet,iband,ix,iy)=resp%a_ker * resp%tmp
!            enddo !iy
!         enddo ! ix

!         ! B .ne. 0
!         !tmp=vka(ik,ib,ix)*( vkab(ik,ib,iy,ix)*vka(ik,ib,iy) - vkab(ik,ib,iy,iy)*vka(ik,ib,ix)    )
!         if (algo%lBfield .and. algo%ltbind ) then
!            do ix=1,lat%nalpha
!               do iy=ix+1,3

!                  resp%tmp =ek%Mopt(ix,thdr%idtet(ik,itet), iband, iband)* &
!                    (ek%M2(iy, ix, thdr%idtet(ik,itet), iband)*ek%Mopt(ix,thdr%idtet(ik,itet), iband, iband) - &
!                    ek%M2(iy, iy, thdr%idtet(ik,itet), iband)* ek%Mopt(ix,thdr%idtet(ik,itet), iband, iband) )
!                  sB_tmp_tetra(iktet,iband,ix,iy)=resp%sB_ker * resp%tmp
!                  aB_tmp_tetra(iktet,iband,ix,iy)=resp%aB_ker * resp%tmp

!               enddo !iy
!            enddo ! ix
!         endif !lBfield

!         ! Now copy the local variable into the datastructure that will be passed to the interptra_re
!         resp%s_tmp(iktet,iband,:,:) = s_tmp_tetra(iktet,iband,:,:)
!         resp%a_tmp(iktet,iband,:,:) = a_tmp_tetra(iktet,iband,:,:)
!         if(algo%lBfield .and. algo%ltbind ) then
!            resp%sB_tmp(iktet,iband,:,:) = sB_tmp_tetra(iktet,iband,:,:)
!            resp%aB_tmp(iktet,iband,:,:) = aB_tmp_tetra(iktet,iband,:,:)
!         endif

!      enddo ! iband

!   enddo ! ik

!  !!!!!!!!!!!!! TEST
!  ! if ((mod(itet,100) ==0) .and. (sct%TT(iT)==90.00)) then
!  !   write(itet+1,*) s_tmp_tetra(1,:,1,1), a_tmp_tetra(1,:,1,1)
!  !   write(itet+2,*) s_tmp_tetra(2,:,1,1), a_tmp_tetra(2,:,1,1)
!  !   write(itet+3,*) s_tmp_tetra(3,:,1,1), a_tmp_tetra(3,:,1,1)
!  !   write(itet+4,*) s_tmp_tetra(4,:,1,1), a_tmp_tetra(4,:,1,1)
!  !   write(itet+10,*) resp%s_tmp(1,:,1,1), resp%a_tmp(1,:,1,1)
!  !   write(itet+20,*) resp%s_tmp(2,:,1,1), resp%a_tmp(2,:,1,1)
!  !   write(itet+30,*) resp%s_tmp(3,:,1,1), resp%a_tmp(3,:,1,1)
!  !   write(itet+40,*) resp%s_tmp(4,:,1,1), resp%a_tmp(4,:,1,1)
!  ! endif
!  !!!!!!!!!!!!! TEST END

!  deallocate(s_tmp_tetra)
!  deallocate(a_tmp_tetra)
!  if (allocated(sB_tmp_tetra)) then
!     deallocate(sB_tmp_tetra)
!     deallocate(aB_tmp_tetra)
!  endif

!end subroutine respintet

!subroutine respintet_qp(mu, iT, itet, thdr, ek, sct, resp)
!  implicit none
!  type (qp_resp) :: resp ! dynamical datatype allow only for an inclusion through extension of the parenttype
!                         ! I could have declared a dummy class pointer that would be assigned to either dp or qp response type
!                         ! to do so the varaibles defined in the individual types must have had different names and since I
!                         ! REALLY dislike having that extra Q for each varible I decided to stick to static datatypes
!  type (tetramesh) :: thdr
!  type (energydisp) :: ek
!  type (scattering) :: sct
!  real(8), intent(in) :: mu
!  integer, intent(in) :: iT
!  integer, intent(in) :: itet
!  integer :: iband, ipg
!  integer :: ix,iy
!  integer :: ik, ikk, iktet
!  complex(8),external  :: wpsipg
!  complex(16),external :: wpsipghp
!!local variables
!  real(16), allocatable :: s_tmp_tetra(:,:,:,:),  a_tmp_tetra(:,:,:,:)
!  real(16), allocatable :: sB_tmp_tetra(:,:,:,:), aB_tmp_tetra(:,:,:,:)
!  real(8) :: Mopt(6,ek%nbopt_min:ek%nbopt_max, ek%nbopt_min:ek%nbopt_max)

!  !allocation
!  if(.not. allocated(s_tmp_tetra)) allocate(s_tmp_tetra(4,ek%nband_max,3,3))
!  if(.not. allocated(a_tmp_tetra)) allocate(a_tmp_tetra(4,ek%nband_max,3,3))
!  if (algo%lBfield .and. algo%ltbind) then
!     if(.not. allocated(sB_tmp_tetra)) allocate(sB_tmp_tetra(4,ek%nband_max,3,3))
!     if(.not. allocated(aB_tmp_tetra)) allocate(aB_tmp_tetra(4,ek%nband_max,3,3))
!  endif
!  !initialisation
!  s_tmp_tetra=0.0q0 ; a_tmp_tetra=0.0q0
!  if (algo%lBfield .and. algo%ltbind) then
!     sB_tmp_tetra=0.0q0 ; aB_tmp_tetra=0.0q0
!  endif

!  do iktet=1,4  !loop over corners of the tetrahedron

!     ik = thdr%idtet(iktet,itet)
!     ! ikk = symm%symop_id(1,ik)
!     ! call getmopt(ek, ik, Mopt, .false.) ! intra

!     do iband=1,ek%nband_max !loop over bands (these will be traced over)

!        ! if the band is not contained in the optical matrices just do nothing
!        if (iband < ek%nbopt_min) cycle
!        if (iband > ek%nbopt_max) cycle
!        resp%z=real(ek%z(ikk,iband),16)
!        if (allocated(sct%ykb)) then
!           resp%gamma=resp%z*real(sct%gam(iT)+sct%ykb(iT,ikk,iband),16)
!        else
!           resp%gamma=resp%z*real(sct%gam(iT),16)
!        endif
!        ! pre-compute all needed digamma functions
!        resp%aqp=resp%z*real(ek%band(ikk,iband)-mu,16)
!        resp%zarg=0.5q0+beta2pQ*(ciQ*resp%aqp+resp%gamma)
!        do ipg=1,3 ! XXX need 0 for alphaxy ????
!           resp%ctmp=wpsipghp(resp%zarg,ipg)
!           resp%RePolyGamma(ipg)=real(resp%ctmp,16)
!           resp%ImPolyGamma(ipg)=imag(resp%ctmp)
!        enddo

!        ! compute transport kernels (omega-part)
!        !
!        resp%tmp=resp%z**2 / (4.q0*piQ**3) ! missing: beta/gamma (multiplied later to keep number reasonable here)
!        resp%s_ker = resp%tmp * ( resp%RePolyGamma(1) - resp%gamma*beta2pQ * resp%RePolyGamma(2) )
!        resp%a_ker = resp%tmp * ( resp%aqp * resp%RePolyGamma(1) - resp%aqp*resp%gamma*beta2pQ*resp%RePolyGamma(2) &
!                   - resp%gamma**2.q0 * beta2pQ * resp%ImPolyGamma(2) )


!        resp%tmp=resp%tmp*3.q0*resp%z/(4.q0*piQ) ! additionally missing: 1/gamma (multiplied later) XXX

!        if(algo%lBfield) then
!           resp%sB_ker=resp%tmp*(-resp%RePolyGamma(1)-resp%gamma*beta2pQ*resp%RePolyGamma(2)-(beta2pQ*resp%gamma)**2/3.q0 &
!                      * resp%RePolyGamma(3))


!           resp%aB_ker=resp%tmp*(resp%aqp*resp%RePolyGamma(1)-resp%aqp*beta2pQ*resp%gamma*resp%RePolyGamma(2)+resp%gamma**2/3.q0 &
!                * beta2pQ * resp%ImPolyGamma(2) &
!                - resp%aqp*resp%gamma**2.q0 / 3.q0 * beta2pQ**2.q0 * resp%ImPolyGamma(3) &
!                + resp%gamma**3.q0 / 3.q0 * beta2pQ**2.q0 * resp%RePolyGamma(3) )
!         endif


!        ! B = 0
!        !tmp=vka(ik,ib,ix)*vka(ik,ib,iy)
!        do ix=1,lat%nalpha

!           if (algo%ltbind) then
!              resp%tmp=real(ek%Mopt(ix,thdr%idtet(ik,itet),iband,iband)*ek%Mopt(ix,thdr%idtet(ik,itet),iband,iband),16)
!           else
!              ! resp%tmp=real(ek%Mopt(ix,thdr%idtet(ik,itet), iband, iband),16) !the optical matrix elements given by Wien2k are squared already
!              resp%tmp=real(Mopt(ix,iband,iband),16) !the optical matrix elements given by Wien2k are squared already
!           endif

!           s_tmp_tetra(iktet,iband,ix,ix)=resp%s_ker * resp%tmp
!           a_tmp_tetra(iktet,iband,ix,ix)=resp%a_ker * resp%tmp

!           do iy=ix+1,lat%nalpha
!              ! resp%tmp=real(ek%Mopt(ix+iy+1,thdr%idtet(ik,itet), iband, iband),16)
!              resp%tmp=real(Mopt(ix+iy+1,iband,iband),16)
!              s_tmp_tetra(iktet,iband,ix,iy)=resp%s_ker * resp%tmp
!              a_tmp_tetra(iktet,iband,ix,iy)=resp%a_ker * resp%tmp
!           enddo !iy
!        enddo ! ix


!        ! B .ne. 0
!        if (algo%lBfield .and. algo%ltbind ) then
!           do ix=1,lat%nalpha
!              do iy=ix+1,lat%nalpha

!                 !tmp=vka(ik,ib,ix)*( vkab(ik,ib,iy,ix)*vka(ik,ib,iy) - vkab(ik,ib,iy,iy)*vka(ik,ib,ix)    )
!                 resp%tmp = real(ek%Mopt(ix,thdr%idtet(ik,itet), iband, iband)* &
!                   (ek%M2(iy, ix, thdr%idtet(ik,itet), iband)*ek%Mopt(ix,thdr%idtet(ik,itet), iband, iband) - &
!                   ek%M2(iy, iy, thdr%idtet(ik,itet), iband)* ek%Mopt(ix,thdr%idtet(ik,itet), iband, iband) ),16)
!                 sB_tmp_tetra(iktet,iband,ix,iy)=resp%sB_ker * resp%tmp
!                 aB_tmp_tetra(iktet,iband,ix,iy)=resp%aB_ker * resp%tmp

!              enddo !iy
!           enddo ! ix
!        endif !lBfield

!        ! Now copy the local variable into the datastructure that will be passed to the interptra_re
!        resp%s_tmp(iktet,iband,:,:) = s_tmp_tetra(iktet,iband,:,:)
!        resp%a_tmp(iktet,iband,:,:) = a_tmp_tetra(iktet,iband,:,:)
!        if(algo%lBfield .and. algo%ltbind ) then
!           resp%sB_tmp(iktet,iband,:,:) = sB_tmp_tetra(iktet,iband,:,:)
!           resp%aB_tmp(iktet,iband,:,:) = aB_tmp_tetra(iktet,iband,:,:)
!        endif

!     enddo ! iband
!  enddo ! iktet

!end subroutine respintet_qp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! RESPINTET_BL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! This routine (and the following one in QP)
!! evaluate the conduction kernel, and the
!! full response functions on the vertices of a
!! tetrahedron in the BOLTZMANN formalism
!!
!subroutine respintet_Bl(mu, iT, itet, thdr, ek, sct, resp)
!  implicit none
!  type (dp_resp) :: resp ! dynamical datatype allow only for an inclusion through extension of the parenttype
!                         ! I could have declared a dummy class pointer that would be assigned to either dp or qp response type
!                         ! to do so the varaibles defined in the individual types must have had different names and since I
!                         ! REALLY dislike having that extra Q for each varible I decided to stick to static datatypes
!  type (tetramesh) :: thdr
!  type (energydisp) :: ek
!  type (scattering) :: sct
!  real(8), intent(in) :: mu
!  integer, intent(in) :: iT
!  integer, intent(in) :: itet
!  integer :: iband, ipg
!  integer :: ix,iy
!  integer :: ik, ikk, iktet
!!local variables
!  real(8), allocatable :: s_tmp_tetra(:,:,:,:),  a_tmp_tetra(:,:,:,:)
!  real(8), allocatable :: sB_tmp_tetra(:,:,:,:), aB_tmp_tetra(:,:,:,:)
!!external variables
!  real(8) :: Mopt(6, ek%nbopt_min:ek%nbopt_max, ek%nbopt_min:ek%nbopt_max)

!  !allocation
!  if(.not. allocated(s_tmp_tetra)) allocate(s_tmp_tetra(4,ek%nband_max,3,3))
!  if(.not. allocated(a_tmp_tetra)) allocate(a_tmp_tetra(4,ek%nband_max,3,3))
!  if (algo%lBfield .and. algo%ltbind) then
!     if(.not. allocated(sB_tmp_tetra)) allocate(sB_tmp_tetra(4,ek%nband_max,3,3))
!     if(.not. allocated(aB_tmp_tetra)) allocate(aB_tmp_tetra(4,ek%nband_max,3,3))
!  endif
!  !initialisation
!  s_tmp_tetra=0.0d0 ; a_tmp_tetra=0.0d0
!  if (algo%lBfield .and. algo%ltbind) then
!     sB_tmp_tetra=0.0d0 ; aB_tmp_tetra=0.0d0
!  endif


!   do iktet=1,4  !loop over corners of the tetrahedron

!      ik = thdr%idtet(iktet,itet) ! itet -> tetrahedron number which is parallelized, iktet the corresponding 4 corners
!      ! ikk = symm%symop_id(1,ik)
!      ! call getmopt(ek, ik, Mopt, .false.) ! intra

!      do iband=1,ek%nband_max !loop over bands (these will be traced over)

!         if (iband < ek%nbopt_min) cycle
!         if (iband > ek%nbopt_max) cycle

!         ! if the band is not contained in the optical matrices just do nothing
!         resp%z=real(ek%z(ikk,iband),8)
!         if (allocated(sct%ykb)) then
!            resp%gamma=resp%z*real(sct%gam(iT)+sct%ykb(iT,ik,iband),8)
!         else
!            resp%gamma=resp%z*real(sct%gam(iT),8)
!         endif
!         resp%aqp=resp%z*real(ek%band(ikk,iband)-mu,8)

!         ! compute transport kernels (omega-part)
!         !
!         resp%tmp= (-1.d0) * resp%z**2 * beta2p ! missing: 1/gamma (multiplied later to keep number reasonable here)
!         resp%s_ker = resp%tmp * dfermi(resp%aqp, beta)
!         resp%a_ker = resp%tmp * resp%aqp * dfermi(resp%aqp, beta)


!         resp%tmp=resp%tmp*3.d0*resp%z/(4.d0*pi) ! additionally missing: 1/gamma**2 (multiplied later) XXX

!         if(algo%lBfield .and. algo%ltbind) then
!            resp%sB_ker = resp%tmp * dfermi(resp%aqp, beta)
!            resp%aB_ker = resp%tmp * resp%aqp * dfermi(resp%aqp, beta)
!         endif


!         ! B = 0
!         !tmp=vka(ik,ib,ix)*vka(ik,ib,iy)
!         do ix=1,lat%nalpha

!            if (algo%ltbind) then
!               resp%tmp=ek%Mopt(ix,thdr%idtet(ik,itet), iband, iband)*ek%Mopt(ix,thdr%idtet(ik,itet), iband, iband)
!            else
!               ! resp%tmp=ek%Mopt(ix,thdr%idtet(ik,itet), iband, iband) !the optical matrix elements given by Wien2k are squared already
!               resp%tmp=Mopt(ix,iband,iband) !the optical matrix elements given by Wien2k are squared already
!            endif
!            s_tmp_tetra(iktet,iband,ix,ix)=resp%s_ker * resp%tmp
!            a_tmp_tetra(iktet,iband,ix,ix)=resp%a_ker * resp%tmp

!            do iy=ix+1,lat%nalpha
!               ! resp%tmp=ek%Mopt(ix+iy+1,thdr%idtet(ik,itet), iband, iband) !the optical matrix elements given by Wien2k are squared already
!               resp%tmp=Mopt(ix+iy+1,iband,iband)
!               s_tmp_tetra(iktet,iband,ix,iy)=resp%s_ker * resp%tmp
!               a_tmp_tetra(iktet,iband,ix,iy)=resp%a_ker * resp%tmp
!            enddo !iy
!         enddo ! ix

!         ! B .ne. 0
!         if (algo%lBfield .and. algo%ltbind ) then
!            do ix=1,lat%nalpha
!               do iy=ix+1,lat%nalpha

!                  !tmp=vka(ik,ib,ix)*( vkab(ik,ib,iy,ix)*vka(ik,ib,iy) - vkab(ik,ib,iy,iy)*vka(ik,ib,ix)    )
!                  resp%tmp =ek%Mopt(ix,thdr%idtet(ik,itet), iband, iband)* &
!                    (ek%M2(iy, ix, thdr%idtet(ik,itet), iband)*ek%Mopt(ix,thdr%idtet(ik,itet), iband, iband) - &
!                    ek%M2(iy, iy, thdr%idtet(ik,itet), iband)* ek%Mopt(ix,thdr%idtet(ik,itet), iband, iband) )

!                  sB_tmp_tetra(iktet,iband,ix,iy)=resp%sB_ker * resp%tmp
!                  aB_tmp_tetra(iktet,iband,ix,iy)=resp%aB_ker * resp%tmp

!               enddo !iy
!            enddo ! ix
!         endif !lBfield

!         ! Now copy the local variable into the datastructure that will be passed to the interptra_re
!         resp%s_tmp(iktet,iband,:,:) = s_tmp_tetra(iktet,iband,:,:)
!         resp%a_tmp(iktet,iband,:,:) = a_tmp_tetra(iktet,iband,:,:)
!         if(algo%lBfield .and. algo%ltbind ) then
!            resp%sB_tmp(iktet,iband,:,:) = sB_tmp_tetra(iktet,iband,:,:)
!            resp%aB_tmp(iktet,iband,:,:) = aB_tmp_tetra(iktet,iband,:,:)
!         endif

!      enddo ! iband
!   enddo ! iktet - 4 corners

!end subroutine respintet_Bl

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! RESPINTERT
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! This routine (QP counterpart missing)
!! evaluates the conduction kernel, and the
!! full response functions on the vertices of a
!! tetrahedron in the KUBO formalism
!! for interband transitions
!! singularities might arise at conical intersections
!! between the two bands
!!
!subroutine respintert(mu, iT, itet, thdr, ek, sct, resp)
!  implicit none
!  type (dp_respinter) :: resp
!  type (tetramesh) :: thdr
!  type (energydisp) :: ek
!  type (scattering) :: sct
!  real(8), intent(in) :: mu
!  integer, intent(in) :: iT
!  integer, intent(in) :: itet
!  integer :: ib1, ib2, ipg !band1, band2, k-point,
!  integer :: ik, ikk, iktet
!  integer :: ix,iy
!  complex(8),external  :: wpsipg
!  complex(16),external :: wpsipghp
!!local variables
!  real(8), allocatable :: s_tmp_tetra(:,:,:,:),  a_tmp_tetra(:,:,:,:)
!  real(8) :: Dqp, Dgamma, Ggamma !qp energy difference, scattering rate difference and sum
!  real(8) :: DD1, DD2   !denominators
!  real(8) :: ReK, ImK, tmp_s, tmp_a
!  real(8) :: Mopt(6, ek%nbopt_min:ek%nbopt_max, ek%nbopt_min:ek%nbopt_max)

!  !allocation
!  if(.not. allocated(s_tmp_tetra)) allocate(s_tmp_tetra(4,ek%nband_max,3,3))
!  if(.not. allocated(a_tmp_tetra)) allocate(a_tmp_tetra(4,ek%nband_max,3,3))
!  !initialisation
!  s_tmp_tetra=0.0d0 ; a_tmp_tetra=0.0d0

!   do iktet=1,4  !loop over corners of the tetrahedron

!      ik = thdr%idtet(iktet,itet)
!      ! ikk = symm%symop_id(1,ik)
!      ! call getmopt(ek, ik, Mopt, .true.) ! inter

!      do ib1=1,ek%nband_max !loop over bands (these will be traced over)

!         ! if the band is not contained in the optical matrices just do nothing
!         if (ib1 < ek%nbopt_min) cycle
!         if (ib1 > ek%nbopt_max) cycle

!         resp%z1=real(ek%z(ikk,ib1),8)
!         if (allocated(sct%ykb)) then
!            resp%gamma1=resp%z1*real(sct%gam(iT)+sct%ykb(iT,ikk,ib1),8)
!         else
!            resp%gamma1=resp%z1*real(sct%gam(iT),8)
!         endif
!         resp%aqp1=resp%z1*real(ek%band(ikk,ib1)-mu,8)
!         !the first state has to belong to the occupied manifold (is this so? I'm not certain anymore...)
!         if (resp%aqp1 > 0.0d0) cycle
!         resp%zarg=0.5d0+beta2p*(ci*resp%aqp1+resp%gamma1)
!         do ipg=1,1
!            resp%ctmp=wpsipg(resp%zarg,ipg)
!            resp%RePolyGamma1(ipg)=real(resp%ctmp,8)
!            resp%ImPolyGamma1(ipg)=imag(resp%ctmp)
!         enddo

!         ! compute transport kernels (omega-part)
!         do ib2=1,ek%nband_max
!            if (ib2 < ek%nbopt_min) cycle
!            if (ib2 > ek%nbopt_max) cycle
!            if (ib2 == ib1 ) cycle

!            !singularities might arise if ek%band1 = ek%band2

!            !second band variables and derived quantities
!            resp%z2=real(ek%z(ikk,ib2),8)
!            if (allocated(sct%ykb)) then
!               resp%gamma2=resp%z2*real(sct%gam(iT)+sct%ykb(iT,ikk,ib2),8)
!            else
!               resp%gamma2=resp%z2*real(sct%gam(iT),8)
!            endif
!            resp%aqp2=resp%z2*real(ek%band(ikk,ib2)-mu,8)
!            !the second state has to belong to the unoccupied manifold
!            if (resp%aqp2 < 0.0d0) cycle
!            resp%zarg=0.5d0+beta2p*(ci*resp%aqp2+resp%gamma2)
!            do ipg=1,1
!               resp%ctmp=wpsipg(resp%zarg,ipg)
!               resp%RePolyGamma2(ipg)=real(resp%ctmp,8)
!               resp%ImPolyGamma2(ipg)=imag(resp%ctmp)
!            enddo

!            Dqp    = resp%aqp1 - resp%aqp2     !Delta csi in eq
!            Dgamma = resp%gamma1 - resp%gamma2 !Delta in eq
!            Ggamma = resp%gamma1 + resp%gamma2 !Gamma in eq
!            DD1 = 1.0d0/(Dqp**2 + Ggamma**2)
!            DD2 = 1.0d0/(Dqp**2 + Dgamma**2)

!            ReK = 2.0d0*resp%gamma1*resp%gamma2*DD2*( (resp%gamma2*resp%aqp1) - (resp%gamma1*resp%aqp2) )
!            ImK = resp%gamma1*resp%gamma2*DD2*( (Ggamma*Dgamma) + (resp%aqp1+resp%aqp2)*Dqp )

!            tmp_s = (resp%z1*resp%z2 * resp%gamma1*resp%gamma2)*DD1*beta/(pi**3)
!            tmp_a = 0.5d0*resp%z1*resp%z2*DD1*(beta**2)/(pi**3)


!            resp%s_ker = tmp_s * ( ((DD2*Dgamma + 0.5d0/resp%gamma2)*resp%RePolyGamma2(1)) &
!                       - ((DD2*Dgamma - 0.5d0/resp%gamma1)*resp%RePolyGamma1(1)) &
!                       + (DD2*Dqp*(resp%ImPolyGamma2(1) - resp%ImPolyGamma1(1))) )

!            resp%a_ker = tmp_a * ( (resp%aqp1*resp%gamma2*resp%RePolyGamma1(1)) + (resp%aqp2*resp%gamma1*resp%RePolyGamma2(1)) &
!                       + (ReK*(resp%RePolyGamma1(1) - resp%RePolyGamma2(1))) + (ImK*(resp%ImPolyGamma2(1) - resp%ImPolyGamma1(1))) )

!            ! B = 0
!            ! tmp=vka(ik,ib,ix)*vka(ik,ib,iy)
!            do ix=1,lat%nalpha
!               if (algo%ltbind) then
!                  resp%tmp=ek%Mopt(ix,thdr%idtet(ik,itet), ib1, ib2)*ek%Mopt(ix,thdr%idtet(ik,itet), ib1, ib2)
!               else
!                  ! resp%tmp=ek%Mopt(ix,thdr%idtet(ik,itet), ib1, ib2) !the optical matrix elements given by Wien2k are squared already
!                  resp%tmp=Mopt(ix,ib1,ib2)
!               endif

!               s_tmp_tetra(iktet,ib1,ix,ix)=s_tmp_tetra(iktet,ib1,ix,ix) + (resp%s_ker * resp%tmp)
!               a_tmp_tetra(iktet,ib1,ix,ix)=a_tmp_tetra(iktet,ib1,ix,ix) + (resp%a_ker * resp%tmp)

!               do iy=ix+1,lat%nalpha
!                  ! resp%tmp=ek%Mopt(ix+iy+1,thdr%idtet(ik,itet), ib1, ib2)
!                  resp%tmp=Mopt(ix+iy+1,ib1,ib2)
!                  s_tmp_tetra(iktet,ib1,ix,iy)=s_tmp_tetra(iktet,ib1,ix,iy) + (resp%s_ker * resp%tmp)
!                  a_tmp_tetra(iktet,ib1,ix,iy)=a_tmp_tetra(iktet,ib1,ix,iy) + (resp%a_ker * resp%tmp)
!               enddo !iy
!            enddo ! ix
!         enddo !ib2

!         ! Now copy the local variable into the datastructure that will be passed to the interptra_re
!         resp%s_tmp(iktet,ib1,:,:) = s_tmp_tetra(iktet,ib1,:,:)
!         resp%a_tmp(iktet,ib1,:,:) = a_tmp_tetra(iktet,ib1,:,:)

!      enddo ! ib1
!   enddo ! ik

!end subroutine respintert

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! RESPINTERT_SYMM
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! This routine (QP counterpart missing)
!! evaluates the conduction kernel, and the
!! full response functions on the vertices of a
!! tetrahedron in the KUBO formalism
!! for interband transitions.
!! It uses a semplified kernel obtained for a
!! 2 band symmetrical semiconductor with
!! a given band gap !
!! Expressions valid only at the Gamma point
!! have been commented out
!!
!subroutine respintert_symm(mu, iT, itet, thdr, ek, sct, resp)
!  implicit none
!  type (dp_respinter) :: resp
!  type (tetramesh) :: thdr
!  type (energydisp) :: ek
!  type (scattering) :: sct
!  real(8), intent(in) :: mu
!  !real(8), intent(in) :: gap
!  integer, intent(in) :: iT
!  integer, intent(in) :: itet
!  integer :: ib1, ib2, ipg !band1, band2, k-point,
!  integer :: ik, ikk, iktet
!  integer :: ix,iy
!  complex(8),external  :: wpsipg
!  complex(16),external :: wpsipghp
!!local variables
!  real(8), allocatable :: s_tmp_tetra(:,:,:,:),  a_tmp_tetra(:,:,:,:)
!  real(8) :: Dqp, Dgamma, Ggamma !qp energy difference, scattering rate difference and sum
!  real(8) :: DD1, DD2   !denominators
!  real(8) :: ReK, ImK, tmp_s, tmp_a
!  real(8) :: Mopt(6, ek%nbopt_min:ek%nbopt_max, ek%nbopt_min:ek%nbopt_max)

!  !allocation
!  if(.not. allocated(s_tmp_tetra)) allocate(s_tmp_tetra(4,ek%nband_max,3,3))
!  if(.not. allocated(a_tmp_tetra)) allocate(a_tmp_tetra(4,ek%nband_max,3,3))
!  !initialisation
!  s_tmp_tetra=0.0d0 ; a_tmp_tetra=0.0d0

!   do iktet=1,4  !loop over corners of the tetrahedron

!      ik = thdr%idtet(iktet,itet)
!      ! ikk = symm%symop_id(1,ik)
!      ! call getmopt(ek, ik, Mopt, .true.) ! inter

!      do ib1=1,ek%nband_max !loop over bands (these will be traced over)

!         ! if the band is not contained in the optical matrices just do nothing
!         if (ib1 < ek%nbopt_min) cycle
!         if (ib1 > ek%nbopt_max) cycle
!         resp%z1=real(ek%z(ikk,ib1),8)
!         if (allocated(sct%ykb)) then
!            resp%gamma1=resp%z1*real(sct%gam(iT)+sct%ykb(iT,ikk,ib1),8)
!         else
!            resp%gamma1=resp%z1*real(sct%gam(iT),8)
!         endif
!         resp%aqp1=resp%z1*real(ek%band(ikk,ib1),8) !in a symmetric SC mu=0
!         ! if the band is unoccupied cycle
!         if(resp%aqp1 > mu) cycle
!         resp%zarg=0.5d0+beta2p*((ci*resp%aqp1)+resp%gamma1)
!         do ipg=1,1
!            resp%ctmp=wpsipg(resp%zarg,ipg)
!            resp%RePolyGamma1(ipg)=real(resp%ctmp,8)
!            resp%ImPolyGamma1(ipg)=imag(resp%ctmp)
!         enddo

!         ! compute transport kernels (omega-part)
!         !
!         do ib2=1,ek%nband_max
!            if (ib2 < ek%nbopt_min) cycle
!            if (ib2 > ek%nbopt_max) cycle
!            if (ib2 == ib1 ) cycle

!            !second band variables and derived quantities
!            resp%z2=real(ek%z(ikk,ib2),8)
!            resp%gamma2=resp%gamma1   !only one gamma required !real(sct%gam(iT,ib2),8)
!            resp%aqp2=resp%z2*real(ek%band(ikk,ib2),8) !in a symmetric SC mu=0
!            ! if the second state is occupied cycle (interband contribution)
!            if(resp%aqp2 < mu) cycle
!            resp%zarg=0.5d0+beta2p*(ci*resp%aqp2+resp%gamma2)
!            do ipg=1,1
!               resp%ctmp=wpsipg(resp%zarg,ipg)
!               resp%RePolyGamma2(ipg)=real(resp%ctmp,8)
!               resp%ImPolyGamma2(ipg)=imag(resp%ctmp)
!            enddo

!            Dqp    = resp%aqp1 - resp%aqp2     !Delta csi in eq
!            !DD1 = 1.0d0/(gap**2 + 4.0d0*(resp%gamma1**2) )
!            DD1 = 1.0d0/(Dqp**2 + 4.0d0*(resp%gamma1**2) )

!            tmp_s = DD1*((resp%z1 * resp%gamma1)**2)*beta/(pi**3)
!            tmp_a = DD1*((resp%z1 * beta)**2)/(2.0d0*(pi**3))


!            resp%s_ker = tmp_s * ( (resp%RePolyGamma2(1) + resp%RePolyGamma1(1))/(2.0d0*resp%gamma1) &
!                       + (resp%ImPolyGamma2(1) - resp%ImPolyGamma1(1))/Dqp )
!                       !- (resp%ImPolyGamma2(1) - resp%ImPolyGamma1(1))/gap ) !only at the Gamma point!!
!            resp%a_ker = tmp_a * ( resp%gamma1*(resp%aqp1*resp%RePolyGamma1(1) + resp%aqp2*resp%RePolyGamma2(1)) &
!                       + (resp%gamma1**2)*(resp%aqp1+resp%aqp2)*(resp%ImPolyGamma2(1)-resp%ImPolyGamma1(1))/Dqp  &
!                       + (resp%gamma1**3)*2.0d0*(resp%RePolyGamma1(1) - resp%RePolyGamma2(1))/Dqp )

!            !only at the Gamma point!!
!            !resp%a_ker = tmp_a * ( resp%gamma1*abs(resp%aqp1)*(resp%RePolyGamma2(1)-resp%RePolyGamma1(1)) &
!            !           + abs(resp%aqp1)*(resp%gamma1**3)*(resp%RePolyGamma2(1)-resp%RePolyGamma1(1))/(resp%aqp1**2) )

!            ! B = 0
!            !tmp=vka(ik,ib,ix)*vka(ik,ib,iy)
!            do ix=1,lat%nalpha
!               if (algo%ltbind) then
!                  resp%tmp=ek%Mopt(ix,thdr%idtet(ik,itet), ib1, ib2)*ek%Mopt(ix,thdr%idtet(ik,itet), ib1, ib2)
!               else
!                  ! resp%tmp=ek%Mopt(ix,thdr%idtet(ik,itet), ib1, ib2) !the optical matrix elements given by Wien2k are squared already
!                  resp%tmp=Mopt(ix,ib1,ib2)
!               endif

!               s_tmp_tetra(iktet,ib1,ix,ix)=s_tmp_tetra(iktet,ib1,ix,ix) + (resp%s_ker * resp%tmp)
!               a_tmp_tetra(iktet,ib1,ix,ix)=a_tmp_tetra(iktet,ib1,ix,ix) + (resp%a_ker * resp%tmp)

!               do iy=ix+1,lat%nalpha
!                  ! resp%tmp=ek%Mopt(ix+iy+1,thdr%idtet(ik,itet), ib1, ib2)
!                  resp%tmp=Mopt(ix+iy+1,ib1,ib2)
!                  s_tmp_tetra(iktet,ib1,ix,iy)=s_tmp_tetra(iktet,ib1,ix,iy) + (resp%s_ker * resp%tmp)
!                  a_tmp_tetra(iktet,ib1,ix,iy)=a_tmp_tetra(iktet,ib1,ix,iy) + (resp%a_ker * resp%tmp)
!               enddo !iy
!            enddo ! ix
!         enddo !ib2

!         ! Now copy the local variable into the datastructure that will be passed to the interptra_re
!         resp%s_tmp(iktet,ib1,:,:) = s_tmp_tetra(iktet,ib1,:,:)
!         resp%a_tmp(iktet,ib1,:,:) = a_tmp_tetra(iktet,ib1,:,:)

!      enddo ! ib1
!   enddo ! iktet

!end subroutine respintert_symm

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! RESPINKM
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! This routine evaluates the conduction kernel
!! for a given k-point in the KUBO formalism
!! for intraband transitions.
!!
!subroutine respinkm(mu, iT, ik, ek, sct, resp)
!  implicit none
!  type (dp_resp) :: resp ! dynamical datatype allow only for an inclusion through extension of the parenttype
!                         ! I could have declared a dummy class pointer that would be assigned to either dp or qp response type
!                         ! to do so the varaibles defined in the individual types must have had different names and since I
!                         ! REALLY dislike having that extra Q for each varible I decided to stick to static datatypes
!  type (energydisp) :: ek
!  type (scattering) :: sct
!  real(8), intent(in) :: mu
!  integer, intent(in) :: iT
!  integer, intent(in) :: ik
!  integer :: iband, ipg
!  integer :: ix,iy
!  complex(8),external  :: wpsipg
!  complex(16),external :: wpsipghp

!  integer :: ikk
!  real(8) :: Mopt(6, ek%nbopt_min:ek%nbopt_max, ek%nbopt_min:ek%nbopt_max)

!  ! ikk  = symm%symop_id(1,ik)
!  ikk  = ik
!  ! get the optical element and save it into Mopt
!  ! if we have a reducible element we just copy it from ek
!  ! otherwise we have to rotate the element of the corresponding irreducible k-point
!  ! call getmopt(ek, ik, Mopt, .false.) ! intra


!  ! TODO optimize this ... when needed -> read in ??
!  if (lat%lortho) then
!     Mopt = 0.d0
!     Mopt(:3, :, :) = ek%Mopt(:,ikk,:,:)
!  else
!     Mopt(:, :, :) = ek%Mopt(:,ikk,:,:)
!  endif


!  !loop over k-points is external
!  do iband=1,ek%nband_max !loop over bands (these will be traced over)
!    ! if the band is not contained in the optical matrices just do nothing
!    if (iband < ek%nbopt_min) cycle
!    if (iband > ek%nbopt_max) cycle


!    resp%z=real(ek%z(ikk,iband),8)
!    if (allocated(sct%ykb)) then
!       resp%gamma=resp%z*real(sct%gam(iT)+sct%ykb(iT,ikk,iband),8)
!    else
!       resp%gamma=resp%z*real(sct%gam(iT),8)
!    endif
!    !resp%aqp=resp%z*real(ek%band(ik,iband)+selfnrg%Re(ik,iband)-mu,8)
!    resp%aqp=resp%z*real(ek%band(ikk,iband)-mu,8)
!    resp%zarg=0.5d0+beta2p*(ci*resp%aqp+resp%gamma)
!    ! pre-compute all needed digamma functions
!    do ipg=1,3 ! XXX need 0 for alphaxy ????
!       resp%ctmp=wpsipg(resp%zarg,ipg)
!       resp%RePolyGamma(ipg)=real(resp%ctmp)
!       resp%ImPolyGamma(ipg)=aimag(resp%ctmp)
!    enddo

!    ! compute transport kernels (omega-part)
!    !
!    resp%tmp=resp%z**2 / (4.d0*pi**3) ! missing: beta/gamma (multiplied later to keep number reasonable here)
!    resp%s_ker = resp%tmp * (resp%RePolyGamma(1) - resp%gamma*beta2p * resp%RePolyGamma(2) )
!    resp%a_ker = resp%tmp * (resp%aqp * resp%RePolyGamma(1) - resp%aqp*resp%gamma*beta2p*resp%RePolyGamma(2) &
!               - resp%gamma**2.d0 * beta2p * resp%ImPolyGamma(2) )


!    resp%tmp=resp%tmp*3.d0*resp%z/(4.d0*pi) ! additionally missing: 1/gamma (multiplied later) XXX

!    if(algo%lBfield .and. algo%ltbind) then
!       resp%sB_ker=resp%tmp*(-resp%RePolyGamma(1)-resp%gamma*beta2p*resp%RePolyGamma(2)-(beta2p*resp%gamma)**2.d0/3.d0 &
!                  * resp%RePolyGamma(3))

!       resp%aB_ker=resp%tmp*(resp%aqp*resp%RePolyGamma(1)-resp%aqp*beta2p*resp%gamma*resp%RePolyGamma(2)+resp%gamma**2.d0/3.d0  &
!            * beta2p * resp%ImPolyGamma(2) &
!            - resp%aqp*resp%gamma**2.d0 / 3.d0 * beta2p**2.d0 * resp%ImPolyGamma(3) &
!            + resp%gamma**3.d0 / 3.d0 * beta2p**2.d0 * resp%RePolyGamma(3) )
!     endif


!    ! B = 0
!    !tmp=vka(ik,ib,ix)*vka(ik,ib,iy)
!    do ix=1,lat%nalpha
!       if (algo%ltbind) then
!          resp%tmp=ek%Mopt(ix,ik, iband, iband)*ek%Mopt(ix,ik, iband, iband)
!       else
!          !the expression requires only the diagonal of the optical matrix elements because a trace is evaluated
!          ! resp%tmp=ek%Mopt(ix,ik, iband, iband) !the optical matrix elements given by Wien2k are squared already
!          resp%tmp=Mopt(ix, iband, iband)
!       endif

!       resp%s_tmp(ik,iband,ix,ix)=resp%s_ker * resp%tmp
!       resp%a_tmp(ik,iband,ix,ix)=resp%a_ker * resp%tmp

!    ! TESTED 25.05.2018
!       do iy=ix+1,lat%nalpha
!          ! resp%tmp=ek%Mopt(ix+iy+1,ik, iband, iband)
!          resp%tmp=Mopt(ix+iy+1, iband, iband)
!          resp%s_tmp(ik,iband,ix,iy)=resp%s_ker * resp%tmp
!          resp%a_tmp(ik,iband,ix,iy)=resp%a_ker * resp%tmp
!       enddo !iy
!    enddo ! ix


!    ! B .ne. 0
!    if (algo%lBfield .and. algo%ltbind ) then
!       do ix=1,lat%nalpha
!          do iy=ix+1,lat%nalpha

!             !tmp=vka(ik,ib,ix)*( vkab(ik,ib,iy,ix)*vka(ik,ib,iy) - vkab(ik,ib,iy,iy)*vka(ik,ib,ix)    )
!             resp%tmp =ek%Mopt(ix, ik, iband, iband)*( ek%M2(iy, ix, ik, iband)*ek%Mopt(ix, ik, iband, iband) - &
!                ek%M2(iy, iy, ik, iband)* ek%Mopt(ix, ik, iband, iband) )
!                 resp%sB_tmp(ik,iband,ix,iy)=resp%sB_ker * resp%tmp
!                 resp%aB_tmp(ik,iband,ix,iy)=resp%aB_ker * resp%tmp

!          enddo !iy
!       enddo ! ix
!    endif !lBfield

!  enddo ! iband

!end subroutine respinkm

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! RESPINKM_QP
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! This routine evaluates the conduction kernel
!! for a given k-point  in the KUBO formalism
!! for intraband transitions.
!! Quadruple precision conterpart of the
!! preceeding routine
!!
!subroutine respinkm_qp(mu, iT, ik, ek, sct, resp)
!  implicit none
!  type (qp_resp) :: resp ! dynamical datatype allow only for an inclusion through extension of the parenttype
!                         ! I could have declared a dummy class pointer that would be assigned to either dp or qp response type
!                         ! to do so the varaibles defined in the individual types must have had different names and since I
!                         ! REALLY dislike having that extra Q for each varible I decided to stick to static datatypes
!  type (energydisp) :: ek
!  type (scattering) :: sct
!  real(8), intent(in) :: mu
!  integer, intent(in) :: iT
!  integer, intent(in) :: ik
!  integer :: iband, ipg
!  integer :: ix,iy
!  complex(8),external  :: wpsipg
!  complex(16),external :: wpsipghp

!  integer :: ikk
!  real(8) :: Mopt(6,ek%nbopt_min:ek%nbopt_max, ek%nbopt_min:ek%nbopt_max)

!  ikk = ik
!  if (lat%lortho) then
!     Mopt = 0.d0
!     Mopt(:3, :, :) = ek%Mopt(:,ikk,:,:)
!  else
!     Mopt(:, :, :) = ek%Mopt(:,ikk,:,:)
!  endif
!  ! ikk  = symm%symop_id(1,ik)
!  ! call getmopt(ek, ik, Mopt, .false.) ! intra

!  do iband=1,ek%nband_max !loop over bands (these will be traced over)


!     if (iband < ek%nbopt_min) cycle
!     if (iband > ek%nbopt_max) cycle


!     ! if the band is not contained in the optical matrices just do nothing
!     resp%z=real(ek%z(ikk,iband),16)
!     if (allocated(sct%ykb)) then
!        resp%gamma=resp%z*real(sct%gam(iT)+sct%ykb(iT,ikk,iband),16)
!     else
!        resp%gamma=resp%z*real(sct%gam(iT),16)
!     endif
!     ! pre-compute all needed digamma functions
!     resp%aqp=resp%z*real(ek%band(ikk,iband)-mu,16)
!     resp%zarg=0.5q0+beta2pQ*(ciQ*resp%aqp+resp%gamma)
!     do ipg=1,3 ! XXX need 0 for alphaxy ????
!        resp%ctmp=wpsipghp(resp%zarg,ipg)
!        resp%RePolyGamma(ipg)=real(resp%ctmp,16)
!        resp%ImPolyGamma(ipg)=imag(resp%ctmp)
!     enddo

!     ! compute transport kernels (omega-part)
!     !
!     resp%tmp=resp%z**2 / (4.q0*piQ**3) ! missing: beta/gamma (multiplied later to keep number reasonable here)
!     resp%s_ker=resp%tmp*(resp%RePolyGamma(1)-resp%gamma*beta2pQ*resp%RePolyGamma(2) )
!     resp%a_ker=resp%tmp*(resp%aqp*resp%RePolyGamma(1)-resp%aqp*resp%gamma*beta2pQ*resp%RePolyGamma(2) &
!               -resp%gamma**2*beta2pQ*resp%ImPolyGamma(2))


!     resp%tmp=resp%tmp*3.q0*resp%z/(4.q0*piQ) ! additionally missing: 1/gamma (multiplied later) XXX

!     if(algo%lBfield) then
!        resp%sB_ker=resp%tmp*(-resp%RePolyGamma(1)-resp%gamma*beta2pQ*resp%RePolyGamma(2)-(beta2pQ*resp%gamma)**2/3.q0 &
!                   *resp%RePolyGamma(3))


!        resp%aB_ker=resp%tmp*(resp%aqp*resp%RePolyGamma(1)-resp%aqp*beta2pQ*resp%gamma*resp%RePolyGamma(2)+resp%gamma**2.q0/3.q0  &
!             * beta2pQ * resp%ImPolyGamma(2) &
!             - resp%aqp*resp%gamma**2.q0 / 3.q0 * beta2pQ**2.q0 * resp%ImPolyGamma(3) &
!             + resp%gamma**3.q0 / 3.q0 * beta2pQ**2.q0 * resp%RePolyGamma(3) )
!      endif


!     ! B = 0
!     !tmp=vka(ik,ib,ix)*vka(ik,ib,iy)
!     do ix=1,lat%nalpha

!        if (algo%ltbind) then
!           resp%tmp=real(ek%Mopt(ix,ik, iband, iband)*ek%Mopt(ix,ik, iband, iband),16)
!        else
!           ! resp%tmp=real(ek%Mopt(ix,ik, iband, iband),16) !the optical matrix elements given by Wien2k are squared already
!           resp%tmp=real(Mopt(ix,iband,iband),16) !the optical matrix elements given by Wien2k are squared already
!        endif

!        resp%s_tmp(ik,iband,ix,ix)=resp%s_ker * resp%tmp
!        resp%a_tmp(ik,iband,ix,ix)=resp%a_ker * resp%tmp

!        do iy=ix+1,lat%nalpha
!           ! resp%tmp=real(ek%Mopt(ix+iy+1,ik, iband, iband),16)
!           resp%tmp=real(Mopt(ix+iy+1,iband,iband),16)
!           resp%s_tmp(ik,iband,ix,iy)=resp%s_ker * resp%tmp
!           resp%a_tmp(ik,iband,ix,iy)=resp%a_ker * resp%tmp
!        enddo !iy
!     enddo ! ix


!     ! B .ne. 0
!     if (algo%lBfield .and. algo%ltbind ) then
!        do ix=1,lat%nalpha
!           do iy=ix+1,lat%nalpha

!              !tmp=vka(ik,ib,ix)*( vkab(ik,ib,iy,ix)*vka(ik,ib,iy) - vkab(ik,ib,iy,iy)*vka(ik,ib,ix)    )
!              resp%tmp = real(ek%Mopt(ix,ik,iband,iband)*(ek%M2(iy,ix,ik,iband)*ek%Mopt(ix,ik,iband,iband) &
!                       - ek%M2(iy, iy, ik, iband)* ek%Mopt(ix, ik, iband, iband) ),16)
!              resp%sB_tmp(ik,iband,ix,iy)=resp%sB_ker * resp%tmp
!              resp%aB_tmp(ik,iband,ix,iy)=resp%aB_ker * resp%tmp

!           enddo !iy
!        enddo ! ix
!     endif !lBfield

!  enddo ! iband

!end subroutine respinkm_qp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! RESPINKM_BL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! This routine evaluates the conduction kernel
!! for a given k-point in the BOLTZMANN formalism
!! for intraband transitions.
!!
!subroutine respinkm_Bl(mu, iT, ik, ek, sct, resp)
!  implicit none
!  type (dp_resp) :: resp ! dynamical datatype allow only for an inclusion through extension of the parenttype
!                         ! I could have declared a dummy class pointer that would be assigned to either dp or qp response type
!                         ! to do so the varaibles defined in the individual types must have had different names and since I
!                         ! REALLY dislike having that extra Q for each varible I decided to stick to static datatypes
!  type (energydisp) :: ek
!  type (scattering) :: sct
!  real(8), intent(in) :: mu
!  integer, intent(in) :: iT
!  integer, intent(in) :: ik
!  integer :: iband, ipg
!  integer :: ix,iy

!  integer :: ikk
!  real(8) :: Mopt(6,ek%nbopt_min:ek%nbopt_max,ek%nbopt_min:ek%nbopt_max)

!  ikk = ik
!  if (lat%lortho) then
!     Mopt = 0.d0
!     Mopt(:3, :, :) = ek%Mopt(:,ikk,:,:)
!  else
!     Mopt(:, :, :) = ek%Mopt(:,ikk,:,:)
!  endif
!  ! ikk = symm%symop_id(1,ik)
!  ! call getmopt(ek, ik, Mopt, .false.) ! intra

!  do iband=1,ek%nband_max !loop over bands (these will be traced over)

!    ! if the band is not contained in the optical matrices just do nothing
!    if (iband < ek%nbopt_min) cycle
!    if (iband > ek%nbopt_max) cycle


!    resp%z=real(ek%z(ikk,iband),8)
!    if (allocated(sct%ykb)) then
!       resp%gamma=resp%z*real(sct%gam(iT)+sct%ykb(iT,ikk,iband),8)
!    else
!       resp%gamma=resp%z*real(sct%gam(iT),8)
!    endif
!    ! pre-compute all needed digamma functions
!    resp%aqp=resp%z*real(ek%band(ikk,iband)-mu,8)

!    ! compute transport kernels (omega-part)
!    !
!    resp%tmp=resp%z**2 / (2.d0*pi) ! missing: 1/gamma (multiplied later to keep number reasonable here)
!    resp%s_ker = resp%tmp * dfermi(resp%aqp, beta)
!    resp%a_ker = resp%tmp * resp%aqp * dfermi(resp%aqp, beta)


!    resp%tmp=resp%tmp*3.d0*resp%z/(4.d0*pi) ! additionally missing: 1/gamma (multiplied later) XXX

!    if(algo%lBfield .and. algo%ltbind) then
!       resp%sB_ker = resp%tmp * ( -dfermi(resp%aqp, beta)  )
!       resp%aB_ker = resp%tmp * ( resp%aqp * dfermi(resp%aqp, beta) )
!     endif


!    ! B = 0
!    !tmp=vka(ik,ib,ix)*vka(ik,ib,iy)
!    do ix=1,lat%nalpha
!       if (algo%ltbind) then
!          resp%tmp=ek%Mopt(ix,ik, iband, iband)*ek%Mopt(ix,ik, iband, iband)
!       else
!          ! resp%tmp=ek%Mopt(ix,ik, iband, iband) !the optical matrix elements given by Wien2k are squared already
!          resp%tmp=Mopt(ix,iband,iband)
!       endif

!       resp%s_tmp(ik,iband,ix,ix)=resp%s_ker * resp%tmp
!       resp%a_tmp(ik,iband,ix,ix)=resp%a_ker * resp%tmp

!       do iy=ix+1,lat%nalpha
!          ! resp%tmp=ek%Mopt(ix+iy+1,ik, iband, iband) !the optical matrix elements given by Wien2k are squared already
!          resp%tmp=Mopt(ix+iy+1,iband,iband)
!          resp%s_tmp(ik,iband,ix,iy)=resp%s_ker * resp%tmp
!          resp%a_tmp(ik,iband,ix,iy)=resp%a_ker * resp%tmp
!       enddo !iy
!    enddo ! ix


!    ! B .ne. 0
!    if (algo%lBfield .and. algo%ltbind ) then
!       do ix=1,lat%nalpha
!          do iy=ix+1,lat%nalpha

!             !tmp=vka(ik,ib,ix)*( vkab(ik,ib,iy,ix)*vka(ik,ib,iy) - vkab(ik,ib,iy,iy)*vka(ik,ib,ix)    )
!             resp%tmp =ek%Mopt(ix, ik, iband, iband)*( ek%M2(iy, ix, ik, iband)*ek%Mopt(ix, ik, iband, iband) - &
!                ek%M2(iy, iy, ik, iband)* ek%Mopt(ix, ik, iband, iband) )
!                 resp%sB_tmp(ik,iband,ix,iy)=resp%sB_ker * resp%tmp
!                 resp%aB_tmp(ik,iband,ix,iy)=resp%aB_ker * resp%tmp

!          enddo !iy
!       enddo ! ix
!    endif !lBfield

!  enddo ! iband

!end subroutine respinkm_Bl

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! RESPINTERKM
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! This routine (QP counterpart missing)
!! evaluates the conduction kernel, and the
!! full response functions for a given k-point
!! in the KUBO formalism for interband transitions
!! singularities might arise at conical intersections
!! between the two bands
!!
!subroutine respinterkm(mu, iT, ik, ek, sct, resp)
!  implicit none
!  type (dp_respinter) :: resp
!  type (energydisp) :: ek
!  type (scattering) :: sct
!  real(8), intent(in) :: mu
!  integer, intent(in) :: iT
!  integer, intent(in) :: ik
!  integer :: ib1, ib2, ipg !band1, band2, degree of Polygamma f'ns
!  integer :: ix,iy
!  complex(8),external  :: wpsipg
!  complex(16),external :: wpsipghp
!!local variables
!  real(8) :: Dqp, Dgamma, Ggamma !qp energy difference, scattering rate difference and sum
!  real(8) :: DD1, DD2   !denominators
!  real(8) :: ReK, ImK, tmp_s, tmp_a
!  integer :: ikk
!  real(8) :: Mopt(6, ek%nbopt_min:ek%nbopt_max, ek%nbopt_min:ek%nbopt_max)


!  ikk = ik
!  if (lat%lortho) then
!     Mopt = 0.d0
!     Mopt(:3, :, :) = ek%Mopt(:,ikk,:,:)
!  else
!     Mopt(:, :, :) = ek%Mopt(:,ikk,:,:)
!  endif
!  ! ikk = symm%symop_id(1,ik)
!  ! call getmopt(ek, ik, Mopt, .true.) ! inter

!  do ib1=1,ek%nband_max !loop over bands (these will be traced over)
!     ! if the band is not contained in the optical matrices just do nothing
!     if (ib1 < ek%nbopt_min) cycle
!     if (ib1 > ek%nbopt_max) cycle


!     resp%z1=real(ek%z(ikk,ib1),8)
!     if (allocated(sct%ykb)) then
!        resp%gamma1=resp%z1*real(sct%gam(iT)+sct%ykb(iT,ikk,ib1),8)
!     else
!        resp%gamma1=resp%z1*real(sct%gam(iT),8)
!     endif
!     resp%aqp1=resp%z1*real(ek%band(ikk,ib1)-mu,8)
!     !the first state has to belong to the occupied manifold
!     if (resp%aqp1 > 0.0d0) cycle
!     resp%zarg=0.5d0+beta2p*(ci*resp%aqp1+resp%gamma1)
!     do ipg=1,1
!        resp%ctmp=wpsipg(resp%zarg,ipg)
!        resp%RePolyGamma1(ipg)=real(resp%ctmp,8)
!        resp%ImPolyGamma1(ipg)=imag(resp%ctmp)
!     enddo

!     ! compute transport kernels (omega-part)
!     !
!     !do ib2=ib1+1,ek%nband_max
!     do ib2=1,ek%nband_max
!        if (ib2 < ek%nbopt_min) cycle
!        if (ib2 > ek%nbopt_max) cycle
!        if (ib2 == ib1 ) cycle
!        !singularities might arise if ek%band1 = ek%band2


!        !second band variables and derived quantities
!        resp%z2=real(ek%z(ikk,ib2),8)
!        if (allocated(sct%ykb)) then
!           resp%gamma2=resp%z2*real(sct%gam(iT)+sct%ykb(iT,ikk,ib2),8)
!        else
!           resp%gamma2=resp%z2*real(sct%gam(iT),8)
!        endif
!        resp%aqp2=resp%z2*real(ek%band(ikk,ib2)-mu,8)
!        !the second state has to belong to the unoccupied manifold
!        if (resp%aqp2 < 0.0d0) cycle
!        resp%zarg=0.5d0+beta2p*(ci*resp%aqp2+resp%gamma2)
!        do ipg=1,1
!           resp%ctmp=wpsipg(resp%zarg,ipg)
!           resp%RePolyGamma2(ipg)=real(resp%ctmp,8)
!           resp%ImPolyGamma2(ipg)=imag(resp%ctmp)
!        enddo

!        Dqp    = resp%aqp1 - resp%aqp2     !Delta csi in eq
!        Dgamma = resp%gamma1 - resp%gamma2 !Delta in eq
!        Ggamma = resp%gamma1 + resp%gamma2 !Gamma in eq
!        DD1 = 1.0d0/(Dqp**2 + Ggamma**2)
!        DD2 = 1.0d0/(Dqp**2 + Dgamma**2)

!        ReK = 2.0d0*resp%gamma1*resp%gamma2*DD2*( (resp%gamma2*resp%aqp1) - (resp%gamma1*resp%aqp2) )
!        ImK = resp%gamma1*resp%gamma2*DD2*( (Ggamma*Dgamma) + (resp%aqp1+resp%aqp2)*Dqp )

!        tmp_s = (resp%z1*resp%z2 * resp%gamma1*resp%gamma2)*DD1*beta/(pi**3)
!        tmp_a = 0.5d0*resp%z1*resp%z2*DD1*(beta**2)/(pi**3)


!        resp%s_ker = tmp_s * ( ((DD2*Dgamma + 0.5d0/resp%gamma2)*resp%RePolyGamma2(1)) &
!                   - ((DD2*Dgamma - 0.5d0/resp%gamma1)*resp%RePolyGamma1(1)) &
!                   + (DD2*Dqp*(resp%ImPolyGamma2(1) - resp%ImPolyGamma1(1))) )

!        resp%a_ker = tmp_a * ( (resp%aqp1*resp%gamma2*resp%RePolyGamma1(1)) + (resp%aqp2*resp%gamma1*resp%RePolyGamma2(1)) &
!                   + (ReK*(resp%RePolyGamma1(1) - resp%RePolyGamma2(1))) + (ImK*(resp%ImPolyGamma2(1) - resp%ImPolyGamma1(1))) )

!        ! B = 0
!        !tmp=vka(ik,ib,ix)*vka(ik,ib,iy)
!        do ix=1,lat%nalpha
!           if (algo%ltbind) then
!              resp%tmp=ek%Mopt(ix,ik, ib1, ib2)*ek%Mopt(ix,ik, ib1, ib2)
!           else
!              ! resp%tmp=ek%Mopt(ix,ik, ib1, ib2) !the optical matrix elements given by Wien2k are squared already
!              resp%tmp=Mopt(ix,ib1,ib2)
!           endif

!           resp%s_tmp(ik,ib1,ix,ix)=resp%s_tmp(ik,ib1,ix,ix) + (resp%s_ker * resp%tmp)
!           resp%a_tmp(ik,ib1,ix,ix)=resp%a_tmp(ik,ib1,ix,ix) + (resp%a_ker * resp%tmp)

!           do iy=ix+1,lat%nalpha
!              ! resp%tmp=ek%Mopt(ix+iy+1,ik, ib1, ib2) !the optical matrix elements given by Wien2k are squared already
!              resp%tmp=Mopt(ix+iy+1,ib1,ib2)
!              resp%s_tmp(ik,ib1,ix,iy)=resp%s_tmp(ik,ib1,ix,iy) + (resp%s_ker * resp%tmp)
!              resp%a_tmp(ik,ib1,ix,iy)=resp%a_tmp(ik,ib1,ix,iy) + (resp%a_ker * resp%tmp)
!           enddo !iy
!        enddo ! ix

!     enddo !ib2
!  enddo ! ib1

!end subroutine respinterkm

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! RESPINTERKM_SYMM
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! This routine (QP counterpart missing)
!! evaluates the conduction kernel, and the
!! full response functions at a given k-point
!! in the KUBO formalism for interband transitions.
!! It uses a semplified kernel obtained for a
!! 2 band symmetrical semiconductor with
!! a given band gap !
!! Expressions valid only at the Gamma point
!! have been commented out
!!
!subroutine respinterkm_symm(mu, iT, ik, ek, sct, resp)
!  implicit none
!  type (dp_respinter) :: resp
!  type (energydisp) :: ek
!  type (scattering) :: sct
!  real(8), intent(in) :: mu
!  !real(8), intent(in) :: gap
!  integer, intent(in) :: iT
!  integer, intent(in) :: ik
!  integer :: ib1, ib2, ipg !band1, band2, degree of Polygamma f'ns
!  integer :: ix,iy
!  complex(8),external  :: wpsipg
!  complex(16),external :: wpsipghp
!!local variables
!  real(8) :: Dqp, Dgamma, Ggamma !qp energy difference, scattering rate difference and sum
!  real(8) :: DD1, DD2   !denominators
!  real(8) :: ReK, ImK, tmp_s, tmp_a
!  integer :: ikk
!  real(8) :: Mopt(6, ek%nbopt_min:ek%nbopt_max, ek%nbopt_min:ek%nbopt_max)

!  ikk = ik
!  if (lat%lortho) then
!     Mopt = 0.d0
!     Mopt(:3, :, :) = ek%Mopt(:,ikk,:,:)
!  else
!     Mopt(:, :, :) = ek%Mopt(:,ikk,:,:)
!  endif
!  ! ikk = symm%symop_id(1,ik)
!  ! call getmopt(ek, ik, Mopt, .true.) ! inter

!  do ib1=1,ek%nband_max !loop over bands (these will be traced over)
!     ! if the band is not contained in the optical matrices just do nothing
!     if (ib1 < ek%nbopt_min) cycle
!     if (ib1 > ek%nbopt_max) cycle

!     resp%z1=real(ek%z(ikk,ib1),8)
!     if (allocated(sct%ykb)) then
!        resp%gamma1=resp%z1*real(sct%gam(iT)+sct%ykb(iT,ikk,ib1),8)
!     else
!        resp%gamma1=resp%z1*real(sct%gam(iT),8)
!     endif
!     resp%aqp1=resp%z1*real(ek%band(ikk,ib1),8) !in a symmetric SC mu=0
!     ! if the band is unoccupied cycle
!     if(resp%aqp1 > mu) cycle
!     resp%zarg=0.5d0+beta2p*((ci*resp%aqp1)+resp%gamma1)
!     do ipg=1,1
!        resp%ctmp=wpsipg(resp%zarg,ipg)
!        resp%RePolyGamma1(ipg)=real(resp%ctmp,8)
!        resp%ImPolyGamma1(ipg)=imag(resp%ctmp)
!     enddo

!     ! compute transport kernels (omega-part)
!     !
!     do ib2=1,ek%nband_max
!        if (ib2 < ek%nbopt_min) cycle
!        if (ib2 > ek%nbopt_max) cycle
!        if (ib2 == ib1 ) cycle


!        !second band variables and derived quantities
!        resp%z2=real(ek%z(ikk,ib2),8)
!        resp%gamma2=resp%gamma1   !only one gamma required !real(sct%gam(iT,ib2),8)
!        resp%aqp2=resp%z2*real(ek%band(ikk,ib2),8) !in a symmetric SC mu=0
!        ! if the second state is occupied cycle (interband contribution)
!        if(resp%aqp2 < mu) cycle
!        resp%zarg=0.5d0+beta2p*(ci*resp%aqp2+resp%gamma2)
!        do ipg=1,1
!           resp%ctmp=wpsipg(resp%zarg,ipg)
!           resp%RePolyGamma2(ipg)=real(resp%ctmp,8)
!           resp%ImPolyGamma2(ipg)=imag(resp%ctmp)
!        enddo

!        Dqp    = resp%aqp1 - resp%aqp2     !Delta csi in eq
!        !DD1 = 1.0d0/(gap**2 + 4.0d0*(resp%gamma1**2) )
!        DD1 = 1.0d0/(Dqp**2 + 4.0d0*(resp%gamma1**2) )

!        tmp_s = DD1*((resp%z1 * resp%gamma1)**2)*beta/(pi**3)
!        tmp_a = DD1*((resp%z1 * beta)**2)/(2.0d0*(pi**3))


!        resp%s_ker = tmp_s * ( (resp%RePolyGamma2(1) + resp%RePolyGamma1(1))/(2.0d0*resp%gamma1) &
!                   + (resp%ImPolyGamma2(1) - resp%ImPolyGamma1(1))/Dqp )
!                   !- (resp%ImPolyGamma2(1) - resp%ImPolyGamma1(1))/gap ) !only at the Gamma point!!
!        resp%a_ker = tmp_a * ( resp%gamma1*(resp%aqp1*resp%RePolyGamma1(1) + resp%aqp2*resp%RePolyGamma2(1)) &
!                   + (resp%gamma1**2)*(resp%aqp1+resp%aqp2)*(resp%ImPolyGamma2(1)-resp%ImPolyGamma1(1))/Dqp  &
!                   + (resp%gamma1**3)*2.0d0*(resp%RePolyGamma1(1) - resp%RePolyGamma2(1))/Dqp )

!        !only at the Gamma point!!
!        !resp%a_ker = tmp_a * ( resp%gamma1*abs(resp%aqp1)*(resp%RePolyGamma2(1)-resp%RePolyGamma1(1)) &
!        !           + abs(resp%aqp1)*(resp%gamma1**3)*(resp%RePolyGamma2(1)-resp%RePolyGamma1(1))/(resp%aqp1**2) )

!        ! B = 0
!        !tmp=vka(ik,ib,ix)*vka(ik,ib,iy)
!        do ix=1,lat%nalpha
!           if (algo%ltbind) then
!              resp%tmp=ek%Mopt(ix,ik, ib1, ib2)*ek%Mopt(ix,ik, ib1, ib2)
!           else
!              ! resp%tmp=ek%Mopt(ix,ik, ib1, ib2) !the optical matrix elements given by Wien2k are squared already
!              resp%tmp=Mopt(ix,ib1,ib2)
!           endif

!           resp%s_tmp(ik,ib1,ix,ix)=resp%s_tmp(ik,ib1,ix,ix) + (resp%s_ker * resp%tmp)
!           resp%a_tmp(ik,ib1,ix,ix)=resp%a_tmp(ik,ib1,ix,ix) + (resp%a_ker * resp%tmp)

!           do iy=ix+1,lat%nalpha
!              ! resp%tmp=ek%Mopt(ix+iy+1,ik, ib1, ib2)
!              resp%tmp=Mopt(ix+iy+1,ib1,ib2)
!              resp%s_tmp(ik,ib1,ix,iy)=resp%s_tmp(ik,ib1,ix,iy) + (resp%s_ker * resp%tmp)
!              resp%a_tmp(ik,ib1,ix,iy)=resp%a_tmp(ik,ib1,ix,iy) + (resp%a_ker * resp%tmp)
!           enddo !iy
!        enddo ! ix

!     enddo !ib2
!  enddo ! ib1

!end subroutine respinterkm_symm

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! RESDERTET_SYMM
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! This routine
!! evaluates the conduction kernel derivatives
!! with respect to beta=1/(kB*T) under the assumption
!! that both the chemical potential and the
!! scattering rate are temperature independent
!! (hence it can give meaningful results only for
!! a symmetric semiconductor).
!! the first derivative is saved in resp%s
!! the second derivative is saved in resp%a
!!
!subroutine resdertet_symm(mu, iT, itet, thdr, ek, sct, resp)
!  implicit none
!  type (dp_respinter) :: resp
!  type (tetramesh) :: thdr
!  type (energydisp) :: ek
!  type (scattering) :: sct
!  real(8), intent(in) :: mu
!  integer, intent(in) :: iT
!  integer, intent(in) :: itet
!  integer :: iband, ipg
!  integer :: ik, ikk, iktet
!  integer :: ix,iy
!  complex(8),external  :: wpsipg
!!local variables
!  real(8), allocatable :: s_tmp_tetra(:,:,:,:),  a_tmp_tetra(:,:,:,:)
!  real(8) :: Mopt(6, ek%nbopt_min:ek%nbopt_max, ek%nbopt_min:ek%nbopt_max)

!  !allocation
!  if(.not. allocated(s_tmp_tetra)) allocate(s_tmp_tetra(4,ek%nband_max,3,3))
!  if(.not. allocated(a_tmp_tetra)) allocate(a_tmp_tetra(4,ek%nband_max,3,3))
!  !initialisation
!  s_tmp_tetra=0.0d0 ; a_tmp_tetra=0.0d0

!   do iktet=1,4  !loop over corners of the tetrahedron

!      ik = thdr%idtet(iktet,itet)
!      ! ikk = symm%symop_id(1,ik)
!      ! call getmopt(ek, ik, Mopt, .false.) ! intra

!      do iband=1,ek%nband_max !loop over bands (these will be traced over)

!         ! if the band is not contained in the optical matrices just do nothing
!         if (iband < ek%nbopt_min) cycle
!         if (iband > ek%nbopt_max) cycle
!         resp%z=real(ek%z(ikk,iband),8)
!         if (allocated(sct%ykb)) then
!            resp%gamma=resp%z*real(sct%gam(iT)+sct%ykb(iT,ikk,iband),8)
!         else
!            resp%gamma=resp%z*real(sct%gam(iT),8)
!         endif
!         ! pre-compute all needed digamma functions
!         resp%aqp=resp%z*real(ek%band(ikk,iband)-mu,8)
!         resp%zarg=0.5d0+beta2p*(ci*resp%aqp+resp%gamma)
!         do ipg=1,4
!            resp%ctmp=wpsipg(resp%zarg,ipg)
!            resp%RePolyGamma(ipg)=real(resp%ctmp,8)
!            resp%ImPolyGamma(ipg)=imag(resp%ctmp)
!         enddo

!         ! compute transport kernel derivatives (omega-part)
!         !
!         ! 1st derivative w.r.t. beta
!         resp%tmp=resp%z**2 / (4.d0*pi**3) ! for the 2nd derivative there is a factor 1/pi missing
!         resp%s_ker = resp%tmp * ((1.0d0/resp%gamma)*resp%RePolyGamma(1) - beta2p*resp%RePolyGamma(2) &
!                    - beta2p*(resp%aqp/resp%gamma)*resp%ImPolyGamma(2) + resp%aqp*(beta2p**2)*resp%ImPolyGamma(3) &
!                    - resp%gamma*(beta2p**2)*resp%RePolyGamma(3) )

!         !
!         ! 2nd derivative w.r.t. beta
!         resp%tmp=resp%z**2 / (4.d0*pi**4)
!         resp%a_ker = resp%tmp * (-(resp%aqp/resp%gamma)*resp%ImPolyGamma(2) &
!                    - 0.5d0*resp%gamma*beta2p*(3.0d0+(resp%aqp/resp%gamma)**2 )*resp%RePolyGamma(3) &
!                    + beta2p*resp%aqp*resp%ImPolyGamma(3) + (beta2p**2)*resp%aqp*resp%gamma*resp%ImPolyGamma(4) &
!                    + 0.5d0*(beta2p**2)*(resp%aqp**2 - resp%gamma**2)*resp%RePolyGamma(4) )

!         !only the xx component has been evaluated
!         do ix=1,1
!            do iy=1,1

!               !tmp=vka(ik,ib,ix)*vka(ik,ib,iy)
!               if (algo%ltbind) then
!                  resp%tmp=ek%Mopt(ix,thdr%idtet(ik,itet), iband, iband)*ek%Mopt(ix,thdr%idtet(ik,itet), iband, iband)
!               else
!                  write(*,*) 'resdertet_symm: the expression for the derivatives is only valid for a symmetric SC'
!                  STOP
!               endif

!               s_tmp_tetra(ik,iband,ix,iy)=resp%s_ker * resp%tmp
!               a_tmp_tetra(ik,iband,ix,iy)=resp%a_ker * resp%tmp

!            enddo !iy
!         enddo ! ix

!         ! Now copy the local variable into the datastructure that will be passed to the interptra_re
!         resp%s_tmp(ik,iband,1,1) = s_tmp_tetra(ik,iband,1,1)
!         resp%a_tmp(ik,iband,1,1) = a_tmp_tetra(ik,iband,1,1)

!      enddo ! iband
!   enddo ! ik

!end subroutine resdertet_symm

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! RESDERTET
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! This routine
!! evaluates the conduction kernel derivatives
!! with respect to beta=1/(kB*T) with a temperature
!! dependent chemical potential and
!! scattering rate. The 2nd derivative of
!! the chemical potential is neglected
!! for the scattering rate it is evaluated assuming
!! gamma(T) = gc0 + gc2*T^2, so
!! gam2dot = 6*gc2/(kB^2 * beta^4)
!! the first derivative of the conductivity is saved in resp%s
!! the second derivative of the conductivity is saved in resp%a
!!
!subroutine resdertet(iT, itet, thdr, ek, sct, resp)
!  implicit none
!  type (dp_respinter) :: resp
!  type (tetramesh) :: thdr
!  type (energydisp) :: ek
!  type (scattering) :: sct
!  integer, intent(in) :: iT
!  integer, intent(in) :: itet
!  complex(8),external  :: wpsipg
!!local variables
!  real(8), allocatable :: s_tmp_tetra(:,:,:,:),  a_tmp_tetra(:,:,:,:)
!  real(8) :: muder, mudot   ! derivative of the chemical potential w.r.t. T, beta
!  real(8) :: gamder, gamdot ! derivative of the scattering rate w.r.t. T, beta
!  real(8) :: dlogg          ! gammadot/gamma
!  real(8) :: gam2dot        ! 2nd derivative of gamma w.r.t. beta, assuming gamma(T) = gc0 + gc2*T^2
!  real(8) :: csim           ! aqp/beta - mudot
!  integer :: iband, ipg
!  integer :: ik, ikk, iktet
!  integer :: ix,iy
!  real(8) :: Mopt(6, ek%nbopt_min:ek%nbopt_max, ek%nbopt_min:ek%nbopt_max)

!  !allocation
!  if(.not. allocated(s_tmp_tetra)) allocate(s_tmp_tetra(4,ek%nband_max,3,3))
!  if(.not. allocated(a_tmp_tetra)) allocate(a_tmp_tetra(4,ek%nband_max,3,3))
!  !initialisation
!  s_tmp_tetra=0.0d0 ; a_tmp_tetra=0.0d0

!  !need to call this subroutine for iT<nT
!  if(iT == sct%nT) STOP !hopefully this case has been taken care of before calling the routine
!  muder = (sct%mu(iT+1)-sct%mu(iT))/sct%dT
!  mudot = -kB*muder*((sct%TT(iT))**2)

!   do iktet=1,4  !loop over corners of the tetrahedron

!      ik = thdr%idtet(iktet,itet)
!      ! ikk = symm%symop_id(1,ik)
!      ! call getmopt(ek, ik, Mopt, .false.) ! intra

!      do iband=1,ek%nband_max !loop over bands (these will be traced over)

!         ! if the band is not contained in the optical matrices just do nothing
!         if (iband < ek%nbopt_min) cycle
!         if (iband > ek%nbopt_max) cycle
!         resp%z=real(ek%z(ikk,iband),8)
!         if (allocated(sct%ykb)) then
!            resp%gamma=resp%z*real(sct%gam(iT)+sct%ykb(iT,ikk,iband),8)
!         else
!            resp%gamma=resp%z*real(sct%gam(iT),8)
!         endif
!         !TODO: replace these values with calls to derrich subroutine
!         gamder = (sct%gam(iT+1)-sct%gam(iT))/sct%dT !safeguard condition set above
!         gamdot = -kB*gamder*((sct%TT(iT))**2)
!         dlogg  = gamdot/resp%gamma
!         if (sct%ng ==2) then
!            gam2dot = 6.0d0*sct%gc(2)/((kB**2)*(beta**4))
!         else
!            gam2dot = 0.0d0
!         endif

!         !resp%aqp=resp%z*real(ek%band(thdr%idtet(ik,itet),iband)+selfnrg%Re(thdr%idtet(ik,itet),iband)-sct%mu(iT),8)
!         resp%aqp=resp%z*real(ek%band(ikk,iband)-sct%mu(iT),8)
!         csim = (resp%aqp/beta)-mudot
!         ! pre-compute all needed digamma functions
!         resp%zarg=0.5d0+beta2p*(ci*resp%aqp+resp%gamma)
!         do ipg=1,4
!            resp%ctmp=wpsipg(resp%zarg,ipg)
!            resp%RePolyGamma(ipg)=real(resp%ctmp,8)
!            resp%ImPolyGamma(ipg)=imag(resp%ctmp)
!         enddo

!         ! compute transport kernel derivatives

!         ! 1st derivative w.r.t. beta (there is a term 2/beta * sigma that needs to be added up)
!         resp%tmp=((resp%z*beta)**2) / (8.d0*pi**4)
!         resp%s_ker = resp%tmp*((dlogg+1.0d0/beta)*(resp%RePolyGamma(2)-(1.0d0/(beta2p*resp%gamma))*resp%RePolyGamma(1)) &
!                    - (csim/resp%gamma)*resp%ImPolyGamma(2) + beta2p*csim*resp%ImPolyGamma(3) &
!                    - (dlogg + 1.0d0/beta)*beta2p*resp%gamma*resp%RePolyGamma(3) )


!         ! 2nd derivative w.r.t. beta (there is a term 2/beta^2 * sigma that needs to be added up)
!         resp%tmp=((resp%z*beta)**2) / (8.d0*pi**4)
!         resp%a_ker = resp%tmp*( ((2.0d0*dlogg/beta) + (2.0d0/(beta**2)) + (gam2dot/resp%gamma) ) &
!                    * (resp%RePolyGamma(2) - (1.0d0/(beta2p*resp%gamma))*resp%RePolyGamma(1)) &
!                    + ((csim*(2.0d0*dlogg - 3.0d0/beta)/resp%gamma) + (csim/(beta*resp%gamma)))*resp%ImPolyGamma(2) &
!                    - beta2p*(resp%gamma*(3.0d0/(beta**2) + 4.0d0*dlogg/beta + gam2dot/resp%gamma - dlogg**2 ) &
!                              + (csim**2)/resp%gamma )*resp%RePolyGamma(3) &
!                    - (2.0d0*beta2p*csim*dlogg + (4.0d0*mudot -2.0d0*resp%aqp/beta)/(2.0d0*pi) )*resp%ImPolyGamma(3) &
!                    + (beta2p**2)*(csim**2 - (gamdot + resp%gamma/beta)**2)*resp%RePolyGamma(4) &
!                    + 2.0d0*(beta2p**2)*(csim*(gamdot + resp%gamma/beta))*resp%ImPolyGamma(4) )

!         ! Now add the missing terms:
!         resp%tmp=(((resp%z*beta)**2)/(8.d0*pi**4)) * ((1.0d0/(beta2p*resp%gamma))*resp%RePolyGamma(1) - resp%RePolyGamma(2))
!         resp%s_ker = resp%s_ker + (2.0d0*resp%tmp/beta)
!         resp%a_ker = resp%a_ker + (2.0d0*resp%tmp/(beta**2))

!         !tmp=vka(ik,ib,ix)*vka(ik,ib,iy)
!         do ix=1,lat%nalpha

!            if (algo%ltbind) then
!               resp%tmp=ek%Mopt(ix,thdr%idtet(ik,itet), iband, iband)*ek%Mopt(ix,thdr%idtet(ik,itet), iband, iband)
!            else
!               ! resp%tmp=ek%Mopt(ix,thdr%idtet(ik,itet), iband, iband) !the optical matrix elements given by Wien2k are squared already
!               resp%tmp=Mopt(ix,iband,iband)
!            endif

!            s_tmp_tetra(iktet,iband,ix,ix)=resp%s_ker * resp%tmp
!            a_tmp_tetra(iktet,iband,ix,ix)=resp%a_ker * resp%tmp

!            do iy=ix+1,lat%nalpha
!               ! resp%tmp=ek%Mopt(ix+iy+1,thdr%idtet(ik,itet), iband, iband)
!               resp%tmp=Mopt(ix+iy+1,iband,iband)
!               s_tmp_tetra(iktet,iband,ix,iy)=resp%s_ker * resp%tmp
!               a_tmp_tetra(iktet,iband,ix,iy)=resp%a_ker * resp%tmp
!            enddo !iy
!         enddo ! ix

!         ! Now copy the local variable into the datastructure that will be passed to the interptra_re
!         resp%s_tmp(iktet,iband,:,:) = s_tmp_tetra(iktet,iband,:,:)
!         resp%a_tmp(iktet,iband,:,:) = a_tmp_tetra(iktet,iband,:,:)

!      enddo ! iband
!   enddo ! iktet

!end subroutine resdertet

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! RESDERKM_SYMM
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! This routine
!! evaluates the conduction kernel derivatives
!! with respect to beta=1/(kB*T) under the assumption
!! that both the chemical potential and the
!! scattering rate are temperature independent
!! (hence it can give meaningful results only for
!! a symmetric semiconductor).
!! the first derivative is saved in resp%s
!! the second derivative is saved in resp%a
!!
!subroutine resderkm_symm(mu, iT, ik, ek, sct, resp)
!  implicit none
!  type (dp_respinter) :: resp
!  type (energydisp) :: ek
!  type (scattering) :: sct
!  real(8), intent(in) :: mu
!  integer, intent(in) :: iT
!  integer, intent(in) :: ik
!  integer :: iband, ipg
!  integer :: ix,iy,ikk
!  real(8) :: Mopt(6,ek%nbopt_min:ek%nbopt_max, ek%nbopt_min:ek%nbopt_max)
!  complex(8),external  :: wpsipg


!  ikk = ik
!  if (lat%lortho) then
!     Mopt = 0.d0
!     Mopt(:3, :, :) = ek%Mopt(:,ikk,:,:)
!  else
!     Mopt(:, :, :) = ek%Mopt(:,ikk,:,:)
!  endif
!  ! ikk = symm%symop_id(1,ik)
!  ! call getmopt(ek, ik, Mopt, .false.) ! intra

!  do iband=1,ek%nband_max !loop over bands (these will be traced over)

!    ! if the band is not contained in the optical matrices just do nothing
!     if (iband < ek%nbopt_min) cycle
!     if (iband > ek%nbopt_max) cycle

!     resp%z=real(ek%z(ikk,iband),8)
!     if (allocated(sct%ykb))then
!        resp%gamma=resp%z*real(sct%gam(iT)+sct%ykb(iT,ikk,iband),8)
!     else
!        resp%gamma=resp%z*real(sct%gam(iT),8)
!     endif
!     ! pre-compute all needed digamma functions
!     resp%aqp=resp%z*real(ek%band(ikk,iband)-mu,8)
!     resp%zarg=0.5d0+beta2p*(ci*resp%aqp+resp%gamma)
!     do ipg=1,4
!        resp%ctmp=wpsipg(resp%zarg,ipg)
!        resp%RePolyGamma(ipg)=real(resp%ctmp,8)
!        resp%ImPolyGamma(ipg)=imag(resp%ctmp)
!     enddo

!     ! compute transport kernel derivatives (omega-part)

!     ! 1st derivative w.r.t. beta
!     resp%tmp=resp%z**2 / (4.d0*pi**3) ! for the 2nd derivative there is a factor 1/pi missing
!     resp%s_ker = resp%tmp * ((1.0d0/resp%gamma)*resp%RePolyGamma(1) - beta2p*resp%RePolyGamma(2) &
!                - beta2p*(resp%aqp/resp%gamma)*resp%ImPolyGamma(2) + resp%aqp*(beta2p**2)*resp%ImPolyGamma(3) &
!                - resp%gamma*(beta2p**2)*resp%RePolyGamma(3) )

!     ! 2nd derivative w.r.t. beta
!     resp%tmp=resp%z**2 / (4.d0*pi**4)
!     resp%a_ker = resp%tmp * (-(resp%aqp/resp%gamma)*resp%ImPolyGamma(2) &
!                - 0.5d0*resp%gamma*beta2p*(3.0d0+(resp%aqp/resp%gamma)**2 )*resp%RePolyGamma(3) &
!                + beta2p*resp%aqp*resp%ImPolyGamma(3) + (beta2p**2)*resp%aqp*resp%gamma*resp%ImPolyGamma(4) &
!                + 0.5d0*(beta2p**2)*(resp%aqp**2 - resp%gamma**2)*resp%RePolyGamma(4) )

!     !only the xx component has been evaluated
!     do ix=1,1
!        do iy=1,1

!           !tmp=vka(ik,ib,ix)*vka(ik,ib,iy)
!           if (algo%ltbind) then
!              resp%tmp=ek%Mopt(ix, ik, iband, iband)*ek%Mopt(ix, ik, iband, iband)
!           else
!              write(*,*) 'resderkm_symm: the expression for the derivatives is only valid for a symmetric SC'
!              STOP
!           endif

!           resp%s_tmp(ik,iband,ix,iy)=resp%s_ker * resp%tmp
!           resp%a_tmp(ik,iband,ix,iy)=resp%a_ker * resp%tmp

!        enddo !iy
!     enddo ! ix

!  enddo ! iband

!end subroutine resderkm_symm

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! RESDERKM
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! This routine
!! evaluates the conduction kernel derivatives
!! with respect to beta=1/(kB*T) with a temperature
!! dependent chemical potential and
!! scattering rate. The 2nd derivatives of both
!! the chemical potential and gamma are neglected
!! the first derivative of the conductivity is saved in resp%s
!! the second derivative of the conductivity is saved in resp%a
!!
!subroutine resderkm(iT, ik, ek, sct, resp)
!  implicit none
!  type (dp_respinter) :: resp
!  type (energydisp) :: ek
!  type (scattering) :: sct
!  integer, intent(in) :: iT
!  integer, intent(in) :: ik
!  complex(8),external  :: wpsipg
!  !complex(16),external :: wpsipghp
!!local variables
!  real(8) :: muder, mudot   ! derivative of the chemical potential w.r.t. T, beta
!  real(8) :: gamder, gamdot ! derivative of the scattering rate w.r.t. T, beta
!  real(8) :: dlogg          ! gammadot/gamma
!  real(8) :: csim           ! aqp/beta - mudot
!  integer :: iband, ipg
!  integer :: ix,iy
!  integer :: ikk
!  real(8) :: Mopt(6, ek%nbopt_min:ek%nbopt_max, ek%nbopt_min:ek%nbopt_max)


!  !need to call this subroutine for iT<nT
!  if(iT == sct%nT) stop !hopefully this case has been taken care of before calling the routine
!  muder = (sct%mu(iT+1)-sct%mu(iT))/sct%dT
!  mudot = -kB*muder*((sct%TT(iT))**2)


!  ikk = ik
!  ! ikk = symm%symop_id(1,ik)
!  ! call getmopt(ek, ik, Mopt, .false.) !intra

!  do iband=1,ek%nband_max !loop over bands (these will be traced over)

!     ! if the band is not contained in the optical matrices just do nothing
!     if (iband < ek%nbopt_min) cycle
!     if (iband > ek%nbopt_max) cycle


!     resp%z=real(ek%z(ikk,iband),8)
!     if (allocated(sct%ykb))then
!        resp%gamma=resp%z*real(sct%gam(iT)+sct%ykb(iT,ikk,iband),8)
!     else
!        resp%gamma=resp%z*real(sct%gam(iT),8)
!     endif
!     !TODO: replace these values with calls to derrich subroutine
!     gamder = (sct%gam(iT+1)-sct%gam(iT))/sct%dT !safeguard condition set above
!     gamdot = -kB*gamder*((sct%TT(iT))**2)
!     dlogg  = gamdot/resp%gamma

!     !resp%aqp=resp%z*real(ek%band(ik,iband)+selfnrg%Re(ik,iband)-sct%mu(iT),8)
!     resp%aqp=resp%z*real(ek%band(ikk,iband)-sct%mu(iT),8)
!     csim = (resp%aqp/beta)-mudot
!     ! pre-compute all needed digamma functions
!     resp%zarg=0.5d0+beta2p*(ci*resp%aqp+resp%gamma)
!     do ipg=1,4
!        resp%ctmp=wpsipg(resp%zarg,ipg)
!        resp%RePolyGamma(ipg)=real(resp%ctmp,8)
!        resp%ImPolyGamma(ipg)=imag(resp%ctmp)
!     enddo

!     ! compute transport kernel derivatives

!     ! 1st derivative w.r.t. beta (there is a term 2/beta * sigma that needs to be added up)
!     resp%tmp=((resp%z*beta)**2) / (8.d0*pi**4)
!     resp%s_ker = resp%tmp*((dlogg+1.0d0/beta)*(resp%RePolyGamma(2)-(1.0d0/(beta2p*resp%gamma))*resp%RePolyGamma(1)) &
!                - (csim/resp%gamma)*resp%ImPolyGamma(2) + beta2p*csim*resp%ImPolyGamma(3) &
!                - (dlogg + 1.0d0/beta)*beta2p*resp%gamma*resp%RePolyGamma(3) )

!     ! 2nd derivative w.r.t. beta (there is a term 2/beta^2 * sigma that needs to be added up)
!     resp%tmp=((resp%z*beta)**2) / (8.d0*pi**4)
!     resp%a_ker = resp%tmp*( ((2.0d0*dlogg/beta) + (2.0d0/(beta**2))) &
!                * (resp%RePolyGamma(2) - (1.0d0/(beta2p*resp%gamma))*resp%RePolyGamma(1)) &
!                + ((csim*(2.0d0*dlogg - 3.0d0/beta)/resp%gamma) + (csim/(beta*resp%gamma)))*resp%ImPolyGamma(2) &
!                - beta2p*(resp%gamma*(3.0d0/(beta**2) + 4.0d0*dlogg/beta) + (csim**2)/resp%gamma )*resp%RePolyGamma(3) &
!                - (2.0d0*beta2p*csim*dlogg + (4.0d0*mudot -2.0d0*resp%aqp/beta)/(2.0d0*pi) )*resp%ImPolyGamma(3) &
!                + (beta2p**2)*(csim**2 - (gamdot + resp%gamma/beta)**2)*resp%RePolyGamma(4) &
!                + 2.0d0*(beta2p**2)*(csim*(gamdot + resp%gamma/beta))*resp%ImPolyGamma(4) )


!     ! Now add the missing terms:
!     resp%tmp=(((resp%z*beta)**2)/(8.d0*pi**4)) * ((1.0d0/(beta2p*resp%gamma))*resp%RePolyGamma(1) - resp%RePolyGamma(2))
!     resp%s_ker = resp%s_ker + (2.0d0*resp%tmp/beta)
!     resp%a_ker = resp%a_ker + (2.0d0*resp%tmp/(beta**2))

!     do ix=1,lat%nalpha
!           if (algo%ltbind) then
!              !the expression requires only the diagonal of the optical matrix elements because a trace is evaluated
!              resp%tmp=ek%Mopt(ix, ik, iband, iband)*ek%Mopt(ix, ik, iband, iband)
!           else
!              !the expression requires only the diagonal of the optical matrix elements because a trace is evaluated
!              ! resp%tmp=ek%Mopt(ix, ik, iband, iband) !the optical matrix elements given by Wien2k are squared already
!              resp%tmp=Mopt(ix,iband,iband)
!           endif

!           resp%s_tmp(ik,iband,ix,ix)=resp%s_ker * resp%tmp
!           resp%a_tmp(ik,iband,ix,ix)=resp%a_ker * resp%tmp

!        do iy=ix+1,lat%nalpha
!           ! resp%tmp=ek%Mopt(ix+iy+1,ik, iband, iband)
!           resp%tmp=Mopt(ix+iy+1,iband,iband)
!           resp%s_tmp(ik,iband,ix,iy)=resp%s_ker * resp%tmp
!           resp%a_tmp(ik,iband,ix,iy)=resp%a_ker * resp%tmp
!        enddo !iy
!     enddo ! ix

!  enddo ! iband

!end subroutine resderkm

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! FINDRHOFLEX
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! This routine tries to solve the equation
!! 2*(d sigma/d beta)^2 - sigma*(d^2 sigma/d beta^2) = 0
!! for a given temperature and stores away the first
!! temperature that fulfills it (in sct%Tstar)
!! since the equation is difficult to solve point by point
!! (i.e. w/o knowing the global behaviour with T)
!! I neglect the term 2*(d sigma/d beta)^2 and find the T
!! for which the 2nd derivative changes sign
!! the first derivative of the conductivity saved in resp2%s
!! the second derivative of the conductivity saved in resp2%a
!!
!subroutine findrhoflex(iT, resp1, resp2, sct)
!   implicit none
!   integer, intent(in) :: iT
!   type(dp_resp) :: resp1      !contains the intraband conductivity
!   type(dp_respinter) :: resp2 !contains its derivatives
!   type(scattering):: sct
!   !local variables
!   real(8), parameter :: tol=1.0d-12
!   real(8) :: tmp, tmp1, tmp2

!   if (.not. allocated(sct%d0)) then
!      allocate(sct%d0(1:sct%nT))
!      allocate(sct%d1(1:sct%nT))
!      allocate(sct%d2(1:sct%nT))
!   endif

!   tmp1 = 2.0d0*((resp2%s_tot(1,1))**2)
!   tmp2 = resp1%s_tot(1,1)*resp2%a_tot(1,1)
!   tmp  = tmp1 - tmp2
!   sct%d0(iT) = tmp
!   sct%d1(iT) = tmp1
!   sct%d2(iT) = tmp2

!end subroutine findrhoflex

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! FINDRHOFLAT
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! This routine finds the value of T for which
!! the derivative  d sigma/d beta changes sign
!! (onset of saturation for the resistivity)
!! the first derivative of the conductivity saved in resp2%s
!! the second derivative of the conductivity saved in resp2%a
!!
!subroutine findrhoflat(iT, resp1, resp2, sct)
!   implicit none
!   integer, intent(in) :: iT
!   type(dp_resp) :: resp1      !contains the intraband conductivity
!   type(dp_respinter) :: resp2 !contains its derivatives
!   type(scattering):: sct
!   !local variables
!   real(8), parameter :: tol=1.0d-12
!   real(8) :: tmp

!   tmp = -resp2%s_tot(1,1)
!   if (tmp < tol) then
!      sct%Tflat = sct%TT(iT)
!   else
!      sct%Tflat = 0.0d0
!   endif

!end subroutine findrhoflat

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! FINDDRHOMAX
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! This routine evaluates the 1st derivative of the
!! resistivity at a given temperature and saves it
!! in the variable sct%drhodT
!! the first derivative of the conductivity saved in resp2%s
!! the second derivative of the conductivity saved in resp2%a
!!
!subroutine finddrhomax(iT, resp1, resp2, sct, drhodT)
!   implicit none

!   integer, intent(in) :: iT
!   type(dp_resp) :: resp1      !contains the intraband conductivity
!   type(dp_respinter) :: resp2 !contains its derivatives
!   type(scattering):: sct
!   real(8), intent(out) :: drhodT(sct%nT)

!   !drhodT(iT) = resp2%s_tot(1,1)/(kB*((sct%TT(iT)*resp1%s_tot(1,1))**2))
!   drhodT(iT) = -resp2%s_tot(1,1)/(resp1%s_tot(1,1)**2)

!end subroutine finddrhomax

!subroutine globfac(kmesh, resp, hpresp)
!   implicit none
!   type(kpointmesh) :: kmesh
!   class(dp_resp) :: resp
!   type(qp_resp), optional ::hpresp
!!local variables
!   integer :: ktot
!   real(8) :: fac,facB
!   real(16):: facQ,facBQ

!   fac   = 2.d0 * pi * ( echarge / (lat%vol*hbarevs)) * 1.d10
!   facB  = 2.d0 * pi**2 * ( echarge / (lat%vol*hbarevs) ) * (1.d-10 / hbareVs)
!   facQ  = 2.q0 * piQ * ( real(echarge,16) / real(lat%vol*hbarevs,16)) * 1.q10
!   facBQ = 2.q0 * piQ**2 * ( real(echarge,16) / real(lat%vol*hbarevs,16)) *  (1.q-10 / real(hbareVs,16))

!   resp%s  = resp%s * fac ! --> sigma in 1/(Ohm m)     [vk's are in eV*Angstroem]
!   resp%a  = resp%a * fac * ( - beta * kb)  ! --> S=alpha/sigma in units V/K (below conversion to mV/K for output of S)
!   resp%s_tot  = resp%s_tot * fac
!   resp%a_tot  = resp%a_tot * fac * ( - beta * kb)
!   if(algo%lBfield) then
!      resp%sB = resp%sB * facB
!      resp%aB = resp%aB * facB * ( - beta * kb)
!      resp%sB_tot = resp%sB_tot * facB
!      resp%aB_tot = resp%aB_tot * facB * ( - beta * kb)
!   endif

!   if (present(hpresp)) then
!      hpresp%s  = hpresp%s * facQ ! --> sigma in 1/(Ohm m)     [vk's are in eV*Angstroem]
!      hpresp%a  = hpresp%a * facQ * ( - betaQ * kbQ)  ! --> S=alpha/sigma in units V/K (below conversion to mV/K for output of S)
!      hpresp%s_tot  = hpresp%s_tot * facQ
!      hpresp%a_tot  = hpresp%a_tot * facQ * ( - betaQ * kbQ)
!      if(algo%lBfield) then
!         hpresp%sB = hpresp%sB * facBQ
!         hpresp%aB = hpresp%aB * facBQ * ( - betaQ * kbQ)
!         hpresp%sB_tot = hpresp%sB_tot * facBQ
!         hpresp%aB_tot = hpresp%aB_tot * facBQ * ( - betaQ * kbQ)
!      endif
!   endif

!end subroutine globfac

!subroutine derresp(resp, hpresp)
!   implicit none
!   class(dp_resp)  :: resp
!   type(qp_resp),optional :: hpresp
!!local variables
!   integer :: ix

!! In Seebeck: *1000 so as to yield [S]=mV/K

!     do ix=1,3
!        if (.not. lat%lortho) then
!           resp%Seebeck(ix)=1000.d0*resp%a_tot(ix,ix)/resp%s_tot(ix,ix)
!           if (present(hpresp)) hpresp%Seebeck(ix)=1000.q0*hpresp%a_tot(ix,ix)/hpresp%s_tot(ix,ix)
!        else
!           resp%Seebeck(ix)=1000.d0*resp%a_tot(1,1)/resp%s_tot(1,1)
!           if (present(hpresp)) hpresp%Seebeck(ix)=1000.q0*hpresp%a_tot(1,1)/hpresp%s_tot(1,1)
!        endif
!     enddo


!!     1 = xy
!!     2 = xz
!!     3 = yz

!     if (algo%lBfield) then
!        resp%Nernst(1) = (resp%aB_tot(1,2)*resp%s_tot(1,1)-resp%a_tot(1,1)*resp%sB_tot(1,2))/(resp%s_tot(1,1)**2)
!        resp%Nernst(2) = (resp%aB_tot(1,3)*resp%s_tot(1,1)-resp%a_tot(1,1)*resp%sB_tot(1,3))/(resp%s_tot(1,1)**2)
!        resp%Nernst(3) = (resp%aB_tot(2,3)*resp%s_tot(2,2)-resp%a_tot(2,2)*resp%sB_tot(2,3))/(resp%s_tot(2,2)**2)
!        resp%Nernst = resp%Nernst * 1000.d0 ! V/K --> mV/K
!        if (present(hpresp)) then
!           hpresp%Nernst(1) = (hpresp%aB_tot(1,2)*hpresp%s_tot(1,1)-hpresp%a_tot(1,1)*hpresp%sB_tot(1,2))/(hpresp%s_tot(1,1)**2)
!           hpresp%Nernst(2) = (hpresp%aB_tot(1,3)*hpresp%s_tot(1,1)-hpresp%a_tot(1,1)*hpresp%sB_tot(1,3))/(hpresp%s_tot(1,1)**2)
!           hpresp%Nernst(3) = (hpresp%aB_tot(2,3)*hpresp%s_tot(2,2)-hpresp%a_tot(2,2)*hpresp%sB_tot(2,3))/(hpresp%s_tot(2,2)**2)
!           hpresp%Nernst = hpresp%Nernst * 1000.q0 ! V/K --> mV/K

!           hpresp%Nernstpart(1) = ( hpresp%aB_tot(1,2)*hpresp%s_tot(1,1))/(hpresp%s_tot(1,1)**2)*1000.q0
!           hpresp%Nernstpart(2) = (-hpresp%a_tot(1,1)*hpresp%sB_tot(1,2))/(hpresp%s_tot(1,1)**2)*1000.q0
!        endif
!        resp%RH(1) = -resp%sB_tot(1,2)/(resp%s_tot(1,1)*resp%s_tot(2,2))
!        resp%RH(2) = -resp%sB_tot(1,3)/(resp%s_tot(1,1)*resp%s_tot(3,3))
!        resp%RH(3) = -resp%sB_tot(2,3)/(resp%s_tot(3,3)*resp%s_tot(3,3))
!        resp%RH = resp%RH * 1.d+7 ! --> 10^-7 m^3/C
!        if (present(hpresp)) then
!           hpresp%RH(1) = -hpresp%sB_tot(1,2)/(hpresp%s_tot(1,1)*hpresp%s_tot(2,2))
!           hpresp%RH(2) = -hpresp%sB_tot(1,3)/(hpresp%s_tot(1,1)*hpresp%s_tot(3,3))
!           hpresp%RH(3) = -hpresp%sB_tot(2,3)/(hpresp%s_tot(3,3)*hpresp%s_tot(3,3))
!           hpresp%RH = hpresp%RH * 1.q+7 ! --> 10^-7 m^3/C
!        endif
!     endif
!end subroutine derresp

!subroutine wrtresp(iT, sct, resp, respinter, respBl, hpresp)
!   implicit none
!   integer, intent(in) :: iT
!   type(scattering)  :: sct
!   type(dp_resp)   :: resp
!   type(dp_respinter) :: respinter
!   type(dp_resp)   :: respBl
!   type(qp_resp),optional :: hpresp
!!local variable
!   integer :: ia, ib
!   real(8) :: T, det, tmp, Mtmp(3,3)
!   real(16):: dpdet, dptmp

!   T = sct%TT(iT)

!   write(30,'(100E20.12)') T,(resp%s_tot(ia,ia), ia=1,lat%nalpha )
!   do ib =1,size(resp%s,1)
!   !write(31,'(100E20.12)') T, (resp%s(:,ia,ia),   ia=1,lat%nalpha )
!   ! this creates a compiler bug
!   !
!   ! write(31,'(E20.12, I6, E20.12)') T,ib, (resp%s(ib,ia,ia),   ia=1,lat%nalpha )
!   !
!   enddo
!   write(40,'(100E20.12)') T,(resp%a_tot(ia,ia), ia=1,lat%nalpha )
!   write(41,'(100E20.12)') T,(resp%a(:,ia,ia),   ia=1,lat%nalpha )
!   if (present(hpresp)) then
!      write(130,'(100E20.12)') T,(hpresp%s_tot(ia,ia), ia=1,lat%nalpha )
!      write(131,'(100E20.12)') T,(hpresp%s(:,ia,ia),   ia=1,lat%nalpha )
!      write(140,'(100E20.12)') T,(hpresp%a_tot(ia,ia), ia=1,lat%nalpha )
!      write(141,'(100E20.12)') T,(hpresp%a(:,ia,ia),   ia=1,lat%nalpha )
!   endif
!   write(230,'(100E20.12)') T,(respBl%s_tot(ia,ia), ia=1,lat%nalpha )
!   write(231,'(100E20.12)') T,(respBl%s(:,ia,ia),   ia=1,lat%nalpha )
!   write(240,'(100E20.12)') T,(respBl%a_tot(ia,ia), ia=1,lat%nalpha )
!   write(241,'(100E20.12)') T,(respBl%a(:,ia,ia),   ia=1,lat%nalpha )

!   write(330,'(100E20.12)') T,(respinter%s_tot(ia,ia), ia=1,lat%nalpha )
!   write(331,'(100E20.12)') T,(respinter%s(:,ia,ia),   ia=1,lat%nalpha )
!   write(340,'(100E20.12)') T,(respinter%a_tot(ia,ia), ia=1,lat%nalpha )
!   write(341,'(100E20.12)') T,(respinter%a(:,ia,ia),   ia=1,lat%nalpha )

!! XXX ACHTUNG not writing out all combinations...
!   if (algo%lBfield) then
!      write(50,'(100E20.12)') T,resp%sB_tot(1,2)
!      write(51,'(100E20.12)') T,resp%sB(:,1,2)
!      write(60,'(100E20.12)') T,resp%aB_tot(1,2)
!      write(61,'(100E20.12)') T,resp%aB(:,1,2)
!      if (present(hpresp)) then
!         write(150,'(100E20.12)') T,hpresp%sB_tot(1,2)
!         write(151,'(100E20.12)') T,hpresp%sB(:,1,2)
!         write(160,'(100E20.12)') T,hpresp%aB_tot(1,2)
!         write(161,'(100E20.12)') T,hpresp%aB(:,1,2)
!      endif
!      write(250,'(100E20.12)') T,respBl%sB_tot(1,2)
!      write(251,'(100E20.12)') T,respBl%sB(:,1,2)
!      write(260,'(100E20.12)') T,respBl%aB_tot(1,2)
!      write(261,'(100E20.12)') T,respBl%aB(:,1,2)
!   endif

!   write(70,'(100E20.12)') T,resp%Seebeck(:) ! in mV/K
!   write(370,'(100E20.12)') T,respinter%Seebeck(:) ! in mV/K
!   if (algo%lBfield) then
!      write(71,'(100E20.12)') T, resp%Nernst(:)  ! in mV/KT
!      write(72,'(100E20.12)') T, resp%RH(:) ! in 10^-7 m^3/C
!      write(73,'(100E20.12)') T, resp%sB_tot(1,2) / resp%s_tot(1,1) ! in 1/T
!      write(74,'(100E20.12)') T, resp%aB_tot(1,2) / resp%a_tot(1,1) ! in 1/T
!   endif
!   if (lat%nalpha==1) then
!   ! diagonal conductivity tensor
!      write(75,'(100E20.12)') T, (1.d0/resp%s_tot(1,1)) ! resistivity in Ohm m
!      write(375,'(100E20.12)') T, (1.d0/respinter%s_tot(1,1)) ! resistivity in Ohm m
!   else
!      if (algo%ldebug) then
!         write(75,'(100E20.12)') T, (1.d0/resp%s_tot(ia,ia),ia=1,lat%nalpha) ! resistivity in Ohm m
!         write(375,'(100E20.12)') T, (1.d0/respinter%s_tot(ia,ia),ia=1,lat%nalpha) ! resistivity in Ohm m
!      else
!   ! evaluate the inverse of conductivity tensor
!         det = (resp%s_tot(1,1)*resp%s_tot(2,2)*resp%s_tot(3,3)) &
!               + (2.0d0*resp%s_tot(1,2)*resp%s_tot(2,3)*resp%s_tot(1,3)) &
!               - (resp%s_tot(1,1)*(resp%s_tot(2,3)**2)) &
!               - (resp%s_tot(2,2)*(resp%s_tot(1,3)**2)) &
!               - (resp%s_tot(3,3)*(resp%s_tot(1,2)**2))
!         tmp = ((resp%s_tot(3,3)*resp%s_tot(2,2)) - (resp%s_tot(2,3)**2))/det
!         write(75,'(100E20.12)') T, tmp, det, resp%s_tot(1,2),resp%s_tot(2,3),resp%s_tot(1,3) ! resistivity in Ohm m (xx component of the inverted matrix)
!         Mtmp(:,:) = resp%s_tot(:,:)+respinter%s_tot(:,:)
!         det = (Mtmp(1,1)*Mtmp(2,2)*Mtmp(3,3)) &
!               + (2.0d0*Mtmp(1,2)*Mtmp(2,3)*Mtmp(1,3)) &
!               - (Mtmp(1,1)*(Mtmp(2,3)**2)) &
!               - (Mtmp(2,2)*(Mtmp(1,3)**2)) &
!               - (Mtmp(3,3)*(Mtmp(1,2)**2))
!         tmp = ((Mtmp(3,3)*Mtmp(2,2)) - (Mtmp(2,3)**2))/det
!         write(375,'(100E20.12)') T, tmp ! resistivity in Ohm m (xx component of the inverted matrix)
!      endif
!   endif

!   if (present(hpresp)) then
!      write(170,'(100E20.12)') T, real(hpresp%Seebeck(:),8) ! in mV/K
!      if (algo%lBfield) then
!         write(171,'(100E20.12)') T, real(hpresp%Nernst(:),8)  ! in mV/KT
!         write(172,'(100E20.12)') T, real(hpresp%RH(:),8)      ! in 10^-7 m^3/C
!         write(173,'(100E20.12)') T, real(hpresp%sB_tot(1,2) / resp%s_tot(1,1) ,8)    ! in 1/T
!         write(174,'(100E20.12)') T, real(hpresp%aB_tot(1,2) / resp%a_tot(1,1) ,8)    ! in 1/T

!         write(180,'(100E20.12)') T,real(hpresp%Nernstpart(1),8)  ! in mV/K
!         write(181,'(100E20.12)') T,real(hpresp%Nernstpart(2),8)  ! in mV/K

!      endif
!      if (lat%nalpha==1) then
!      ! diagonal conductivity tensor
!         write(175,'(100E20.12)') T, (real(1.q0/hpresp%s_tot(ia,ia),16),ia=1,lat%nalpha) ! resistivity in Ohm m
!      else
!         if (algo%ldebug) then
!            write(175,'(100E20.12)') T, (real(1.q0/hpresp%s_tot(ia,ia),16),ia=1,lat%nalpha) ! resistivity in Ohm m
!         else
!      ! evaluate the inverse of conductivity tensor
!            dpdet = (hpresp%s_tot(1,1)*hpresp%s_tot(2,2)*hpresp%s_tot(3,3)) &
!                    + (2.0q0*hpresp%s_tot(1,2)*hpresp%s_tot(2,3)*hpresp%s_tot(1,3)) &
!                    - (hpresp%s_tot(1,1)*(hpresp%s_tot(2,3)**2)) &
!                    - (hpresp%s_tot(2,2)*(hpresp%s_tot(1,3)**2)) &
!                    - (hpresp%s_tot(3,3)*(hpresp%s_tot(1,2)**2))
!            dptmp = ((hpresp%s_tot(3,3)*hpresp%s_tot(2,2)) - (hpresp%s_tot(2,3)**2))/dpdet
!            write(175,'(100E20.12)') T, real(dptmp,16),dpdet,hpresp%s_tot(1,2),hpresp%s_tot(2,3),hpresp%s_tot(1,3)! resistivity in Ohm m (xx component of the inverted matrix)
!         endif
!      endif
!   endif

!   write(270,'(100E20.12)') T,respBl%Seebeck(:) ! in mV/K
!   if (algo%lBfield) then
!      write(271,'(100E20.12)') T, respBl%Nernst(:)  ! in mV/KT
!      write(272,'(100E20.12)') T, respBl%RH(:) ! in 10^-7 m^3/C
!      write(273,'(100E20.12)') T, respBl%sB_tot(1,2) / resp%s_tot(1,1) ! in 1/T
!      write(274,'(100E20.12)') T, respBl%aB_tot(1,2) / resp%a_tot(1,1) ! in 1/T
!   endif
!   if (lat%nalpha==1) then
!   ! diagonal conductivity tensor
!      write(275,'(100E20.12)') T, (1.d0/respBl%s_tot(ia,ia),ia=1,lat%nalpha) ! resistivity in Ohm m
!   else
!      if (algo%ldebug) then
!         write(275,'(100E20.12)') T, (1.d0/respBl%s_tot(ia,ia),ia=1,lat%nalpha) ! resistivity in Ohm m
!      else
!   ! evaluate the inverse of conductivity tensor
!         det = (respBl%s_tot(1,1)*respBl%s_tot(2,2)*respBl%s_tot(3,3)) &
!               + (2.0d0*respBl%s_tot(1,2)*respBl%s_tot(2,3)*respBl%s_tot(1,3)) &
!               - (respBl%s_tot(1,1)*respBl%s_tot(2,3)**2) &
!               - (respBl%s_tot(2,2)*respBl%s_tot(1,3)**2) &
!               - (respBl%s_tot(3,3)*respBl%s_tot(1,2)**2)
!         tmp = ((respBl%s_tot(3,3)*respBl%s_tot(2,2)) - (respBl%s_tot(2,3)**2))/det
!         write(275,'(100E20.12)') T, tmp ! resistivity in Ohm m (xx component of the inverted matrix)
!      endif
!   endif

!end subroutine wrtresp

!subroutine response_open_files()
!   implicit none

!   open(30,file='sigma_tot.dat',status='unknown')
!   open(31,file='sigma_band_xx.dat',status='unknown')
!   open(40,file='peltier_tot.dat',status='unknown')
!   open(41,file='peltier_band_xx.dat',status='unknown')
!! 50 for sxy
!   if (algo%lBfield) then
!   open(50,file='sigmaB_tot.dat',status='unknown')
!   open(51,file='sigmaB_band_xy.dat',status='unknown')
!! 60 for axy
!   open(60,file='peltierB_tot.dat',status='unknown')
!   open(61,file='peltierB_band_xy.dat',status='unknown')

!   open(71,file='Nernst.dat',status='unknown')
!   open(72,file='RH.dat',status='unknown')
!   open(73,file='muH.dat',status='unknown') ! Hall mobility
!   open(74,file='mut.dat',status='unknown') ! thermal counterpart of Hall mobility
!   endif

!! 70 for aux's: Seebeck, Hall, Nernst...
!   open(70,file='Seebeck.dat',status='unknown')

!   open(75,file='resistivity.dat',status='unknown')


!! QUAD PRECISION
!   open(130,file='sigma_tot.datQ',status='unknown')
!   open(131,file='sigma_band_xx.datQ',status='unknown')
!   open(140,file='peltier_tot.datQ',status='unknown')
!   open(141,file='peltier_band_xx.datQ',status='unknown')
!! 50 for sxy
!   if (algo%lBfield) then
!   open(150,file='sigmaB_tot.datQ',status='unknown')
!   open(151,file='sigmaB_band_xy.datQ',status='unknown')
!! 60 for axy
!   open(160,file='peltierB_tot.datQ',status='unknown')
!   open(161,file='peltierB_band_xy.datQ',status='unknown')
!   open(171,file='Nernst.datQ',status='unknown')
!   open(172,file='RH.datQ',status='unknown')
!   open(173,file='muH.datQ',status='unknown') ! Hall mobility
!   open(174,file='mut.datQ',status='unknown') ! thermal counterpart of Hall mobility
!   open(180,file='Nernst_part1.datQ',status='unknown')
!   open(181,file='Nernst_part2.datQ',status='unknown')
!   endif
!! 70 for aux's: Seebeck, Hall, Nernst...
!   open(170,file='Seebeck.datQ',status='unknown')

!   open(175,file='resistivity.datQ',status='unknown')

!! DOUBLE PRECISION, BOLTZMANN MODE
!   open(230,file='sigma_tot.datB',status='unknown')
!   open(231,file='sigma_band_xx.datB',status='unknown')
!   open(240,file='peltier_tot.datB',status='unknown')
!   open(241,file='peltier_band_xx.datB',status='unknown')
!! 50 for sxy
!   if (algo%lBfield) then
!   open(250,file='sigmaB_tot.datB',status='unknown')
!   open(251,file='sigmaB_band_xy.datB',status='unknown')
!! 60 for axy
!   open(260,file='peltierB_tot.datB',status='unknown')
!   open(261,file='peltierB_band_xy.datB',status='unknown')
!   open(271,file='Nernst.datB',status='unknown')
!   open(272,file='RH.datB',status='unknown')
!   open(273,file='muH.datB',status='unknown') ! Hall mobility
!   open(274,file='mut.datB',status='unknown') ! thermal counterpart of Hall mobility
!   !open(280,file='Nernst_part1.datBQ',status='unknown')
!   !open(281,file='Nernst_part2.datBQ',status='unknown')
!   endif
!! 70 for aux's: Seebeck, Hall, Nernst...
!   open(270,file='Seebeck.datB',status='unknown')
!   open(275,file='resistivity.datB',status='unknown')

!! DOUBLE PRECISION, INTERBAND RESPONSE
!   open(330,file='sigma_inter_tot.dat',status='unknown')
!   open(331,file='sigma_interband_xx.dat',status='unknown')
!   open(340,file='peltier_inter_tot.dat',status='unknown')
!   open(341,file='peltier_interband_xx.dat',status='unknown')

!! 70 for aux's: Seebeck, Hall, Nernst...
!   open(370,file='Seebeck_inter.dat',status='unknown')

!   open(375,file='resistivity_inter.dat',status='unknown')

!end subroutine response_open_files

!subroutine response_close_files()
!   implicit none

!   close(30)
!   close(31)
!   close(40)
!   close(41)
!   if (algo%lBfield) then
!     close(50)
!     close(51)
!     close(60)
!     close(61)
!     close(71)
!     close(72)
!     close(73)
!     close(74)
!   endif
!   close(70)
!   close(75)

!   close(130)
!   close(131)
!   close(140)
!   close(141)
!   if (algo%lBfield) then
!     close(150)
!     close(151)
!     close(160)
!     close(161)
!     close(171)
!     close(172)
!     close(173)
!     close(174)
!   endif

!   close(170)
!   close(175)

!   close(180)
!   close(181)

!!Boltzmann
!   close(230)
!   close(231)
!   close(240)
!   close(241)
!   if (algo%lBfield) then
!     close(250)
!     close(251)
!     close(260)
!     close(261)
!     close(271)
!     close(272)
!     close(273)
!     close(274)
!   endif
!   close(270)
!   close(275)
!   !close(280)
!   !close(281)

!!Interband
!   close(330)
!   close(331)
!   close(340)
!   close(341)
!   close(370)
!   close(375)

!end subroutine response_close_files


subroutine dpresp_alloc(algo, edisp, dpresp)
  implicit none
  type(algorithm)   :: algo
  type(energydisp)  :: edisp
  type(response_dp) :: dpresp

  ! allocate transport variables
  allocate(dpresp%s_full(3,3,edisp%nbopt_min:edisp%nbopt_max,edisp%ispin,ikstr:ikend))
  allocate(dpresp%a_full(3,3,edisp%nbopt_min:edisp%nbopt_max,edisp%ispin,ikstr:ikend))
  allocate(dpresp%s_sum(3,3,edisp%iSpin))
  allocate(dpresp%a_sum(3,3,edisp%iSpin))

  if (algo%lBfield) then
     allocate(dpresp%sB_full(3,3,edisp%nbopt_min:edisp%nbopt_max,edisp%ispin,ikstr:ikend))
     allocate(dpresp%aB_full(3,3,edisp%nbopt_min:edisp%nbopt_max,edisp%ispin,ikstr:ikend))
     allocate(dpresp%sB_sum(3,3,edisp%iSpin))
     allocate(dpresp%aB_sum(3,3,edisp%iSpin))
  endif

end subroutine dpresp_alloc

subroutine qpresp_alloc(algo, edisp, qpresp)
  implicit none
  type(algorithm)   :: algo
  type(energydisp)  :: edisp
  type(response_qp) :: qpresp

  ! allocate transport variables
  allocate(qpresp%s_full(3,3,edisp%nbopt_min:edisp%nbopt_max,edisp%ispin,ikstr:ikend))
  allocate(qpresp%a_full(3,3,edisp%nbopt_min:edisp%nbopt_max,edisp%ispin,ikstr:ikend))
  allocate(qpresp%s_sum(3,3,edisp%iSpin))
  allocate(qpresp%a_sum(3,3,edisp%iSpin))

  if (algo%lBfield) then
     allocate(qpresp%sB_full(3,3,edisp%nbopt_min:edisp%nbopt_max,edisp%ispin,ikstr:ikend))
     allocate(qpresp%aB_full(3,3,edisp%nbopt_min:edisp%nbopt_max,edisp%ispin,ikstr:ikend))
     allocate(qpresp%sB_sum(3,3,edisp%iSpin))
     allocate(qpresp%aB_sum(3,3,edisp%iSpin))
  endif
end subroutine qpresp_alloc

subroutine calc_polygamma_D(PolyGamma, mu, edisp, sct, kmesh, algo, info)
  implicit none
  type(algorithm)  :: algo
  type(energydisp) :: edisp
  type(kpointmesh) :: kmesh
  type(scattering) :: sct
  type(runinfo)    :: info

  real(8), intent(in) :: mu
  complex(8) :: PolyGamma(3,edisp%nbopt_min:edisp%nbopt_max, ikstr:ikend, edisp%ispin)

  complex(8), external :: wpsipg
  complex(8), allocatable :: to_evaluate(:,:,:)
  integer :: ipg, iband, ik, is

  allocate(to_evaluate(edisp%nbopt_min:edisp%nbopt_max, ikstr:ikend, edisp%ispin))

  if (algo%lScatteringFile) then
    to_evaluate = 0.5d0 + info%beta2p * &
                  (sct%gam(:,ikstr:ikend,:) + ci*(sct%zqp(:,ikstr:ikend,:)*edisp%band(:,ikstr:ikend,:) - mu))
  else
    to_evaluate = 0.5d0 + info%beta2p * &
                  (sct%gamscalar + ci*(sct%zqpscalar*edisp%band(:,ikstr:ikend,:) - mu))
  endif

  do is = 1,edisp%ispin
    do ik = ikstr, ikend
      do iband = edisp%nbopt_min,edisp%nbopt_max
        do ipg = 1,3
          PolyGamma(ipg,iband,ik,is) = wpsipg(to_evaluate(iband,ik,is),ipg)
        enddo
      enddo
    enddo
  enddo

  deallocate(to_evaluate)

end subroutine calc_polygamma_D

subroutine calc_polygamma_Q(PolyGamma, mu, edisp, sct, kmesh, algo, info)
  implicit none
  type(algorithm)  :: algo
  type(energydisp) :: edisp
  type(kpointmesh) :: kmesh
  type(scattering) :: sct
  type(runinfo)    :: info

  real(8), intent(in) :: mu
  complex(16) :: PolyGamma(3,edisp%nbopt_min:edisp%nbopt_max, ikstr:ikend, edisp%ispin)

  complex(16), external :: wpsipghp
  complex(16), allocatable :: to_evaluate(:,:,:)
  integer :: ipg, iband, ik, is

  allocate(to_evaluate(edisp%nbopt_min:edisp%nbopt_max, ikstr:ikend, edisp%ispin))

  if (algo%lScatteringFile) then
    to_evaluate = 0.5q0 + info%beta2pQ * &
                  (sct%gam(:,ikstr:ikend,:) + ciQ*(sct%zqp(:,ikstr:ikend,:)*edisp%band(:,ikstr:ikend,:) - mu))
  else
    to_evaluate = 0.5q0 + info%beta2pQ * &
                  (sct%gamscalar + ci*(sct%zqpscalar*edisp%band(:,ikstr:ikend,:) - mu))
  endif

  do is = 1,edisp%ispin
    do ik = ikstr, ikend
      do iband = edisp%nbopt_min,edisp%nbopt_max
        do ipg = 1,3
          PolyGamma(ipg,iband,ik,is) = wpsipghp(to_evaluate(iband,ik,is),ipg)
        enddo
      enddo
    enddo
  enddo

  deallocate(to_evaluate)

end subroutine calc_polygamma_Q





 ! INPUT: (renolmalised) bandstructure: (sct) ek
 !        chemical potential
 !        energy window (to be defined)
 ! OUTPUT: integrated DOS (summed over bands and over k-points)
 !
!subroutine intldos(iT, dos, kmesh, ek, sct)
!  !passed variables
!  integer :: iT
!  type(dosgrid) :: dos
!  type(kpointmesh) :: kmesh !reducible k-mesh
!  type(energydisp) :: ek
!  type(scattering) :: sct
!  !local variables
!  integer :: ee, ik, ibn, ikk
!  real(8) :: eps, eta
!  real(8) :: s, dee
!  real(8), allocatable :: AA(:)
!  complex :: iu, G0

!  dee= (dos%emax-dos%emin)/real(dos%nnrg)
!  eta=dee
!  allocate(AA(dos%nnrg))
!  AA = 0.0d0
!  do ee=1,dos%nnrg
!     do ibn=1, ek%nband_max
!        do ikk=1, kmesh%ktot
!           if (ek%band(ikk,ibn) > band_fill_value) cycle
!           eps = ek%z(ikk,ibn)*ek%band(ikk,ibn)-sct%mu(iT)
!           if (eps > 0.0d0) exit !occupied bands only -> we skip to the next k-point
!           G0=1.0d0/(dos%enrg(ee) - eps - ci*sct%gam(iT))
!           ! A = -1/pi  Im G0
!           AA(ee) = AA(ee) - aimag(G0)*kmesh%weight(ikk)
!        enddo
!     enddo
!  enddo !ee
!  AA = AA*2.d0*dee/(pi) ! spin, spacing, pi from spectralfunction, kmesh normalization
!  ! AA = AA*2.d0*dee/(pi*kmesh%kred) ! spin, spacing, pi from spectralfunction, kmesh normalization
!  !use trapezoidal rule to evaluate the number of electrons
!  s=0.5d0*(AA(1)+AA(dos%nnrg))
!  do ee=2,dos%nnrg-1
!     s=s+AA(ee)
!  enddo

!  write(*,120)'integrated spectral density ',s,' electron number ',ek%nelect, 'diff',ek%nelect-s
!  120 FORMAT (A,f12.6,A,f12.6)
!  do ee=1,dos%nnrg
!     write(120,'(E12.6,3x,E12.6)') dos%enrg(ee),AA(ee)
!  enddo

!end subroutine intldos


end module Mresponse
