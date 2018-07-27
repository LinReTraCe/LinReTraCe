

module response

use mpi_org
implicit none

! DOUBLE PRECISION
type dp_resp
!integer,parameter :: iq=8
!integer, kind :: iq=8
!kernels
real(8) ::  s_ker  ! for conductivity
real(8) ::  sB_ker ! for conductivity in B-field
real(8) ::  a_ker  ! for Peltier
real(8) ::  aB_ker ! for Peltier in B-field          
!response functions...
real(8), allocatable :: s_tmp(:,:,:,:)   ! nk,nband,3,3 for conductivity
real(8), allocatable :: sB_tmp(:,:,:,:)  ! nk,nband,3,3 for conductivity in B-field
real(8), allocatable :: a_tmp(:,:,:,:)   ! nk,nband,3,3 for Peltier
real(8), allocatable :: aB_tmp(:,:,:,:)  ! nk,nband,3,3 for Peltier in B-field     
real(8), allocatable :: s(:,:,:)         ! nband,3,3 for conductivity
real(8), allocatable :: sB(:,:,:)        ! nband,3,3 for conductivity in B-field
real(8), allocatable :: a(:,:,:)         ! nband,3,3 for Peltier
real(8), allocatable :: aB(:,:,:)        ! nband,3,3 for Peltier in B-field           
real(8), allocatable :: s_local(:,:,:)   ! used only in the mpi version
real(8), allocatable :: sB_local(:,:,:)  ! 
real(8), allocatable :: a_local(:,:,:)   ! 
real(8), allocatable :: aB_local(:,:,:)  ! 

real(8) :: s_tot(3,3),  s_tet(3,3)
real(8) :: sB_tot(3,3), sB_tet(3,3)
real(8) :: a_tot(3,3),  a_tet(3,3)
real(8) :: aB_tot(3,3), aB_tet(3,3)
real(8) :: Seebeck(3),Nernst(3),RH(3)

real(8) :: RePolyGamma(0:4),ImPolyGamma(0:4),gamma,aqp,z,tmp
complex(8) :: ctmp,zarg
end type

!interband transitions functions
type, extends(dp_resp) :: dp_respinter
!variables for the second band
real(8) :: RePolyGamma1(0:4),ImPolyGamma1(0:4)
real(8) :: RePolyGamma2(0:4),ImPolyGamma2(0:4)
real(8) :: gamma1,gamma2, aqp1,aqp2, z1,z2 
end type

type qp_resp
!integer,parameter :: iq=16
!integer, kind :: iq=16
!kernels
real(16) ::  s_ker  ! for conductivity
real(16) ::  sB_ker ! nband,3,3 for conductivity in B-field
real(16) ::  a_ker  ! nband,3,3 for Peltier
real(16) ::  aB_ker ! nband,3,3 for Peltier in B-field          
!response functions...
real(16), allocatable :: s_tmp(:,:,:,:)   ! nk,nband,3,3 for conductivity
real(16), allocatable :: sB_tmp(:,:,:,:)  ! nk,nband,3,3 for conductivity in B-field
real(16), allocatable :: a_tmp(:,:,:,:)   ! nk,nband,3,3 for Peltier
real(16), allocatable :: aB_tmp(:,:,:,:)  ! nk,nband,3,3 for Peltier in B-field     
real(16), allocatable :: s(:,:,:)         ! nband,3,3 for conductivity
real(16), allocatable :: sB(:,:,:)        ! nband,3,3 for conductivity in B-field
real(16), allocatable :: a(:,:,:)         ! nband,3,3 for Peltier
real(16), allocatable :: aB(:,:,:)        ! nband,3,3 for Peltier in B-field           
real(16), allocatable :: s_local(:,:,:)   ! used only in the mpi version
real(16), allocatable :: sB_local(:,:,:)  ! 
real(16), allocatable :: a_local(:,:,:)   ! 
real(16), allocatable :: aB_local(:,:,:)  ! 

real(16) :: s_tot(3,3),  s_tet(3,3)
real(16) :: sB_tot(3,3), sB_tet(3,3)
real(16) :: a_tot(3,3),  a_tet(3,3)
real(16) :: aB_tot(3,3), aB_tet(3,3)
real(16) :: Seebeck(3),Nernst(3),RH(3),Nernstpart(2) ! nernstpart is axy sxx bzw axx sxy / sxx^2 

real(16) :: RePolyGamma(0:4),ImPolyGamma(0:4),gamma,aqp,z,tmp
complex(16) :: ctmp,zarg
end type

end module response


subroutine calc_response(mu, iT, drhodT, algo, mesh, ek, thdr, sct, dresp, dderesp, dinter, respBl, qresp)
  use response
  use params 
  use types
  use estruct
  implicit none

  interface
   subroutine initresp (lBfield, dresp, respBl, qresp)
    use response
    logical :: lBfield 
    class(dp_resp) :: dresp
    type(dp_resp) :: respBl
    type(qp_resp),optional :: qresp
   end subroutine

   subroutine globfac (icubic, algo, mesh, resp, hpresp)
    use params
    use types
    use response
    use estruct, only:vol
    integer :: icubic
    type(algorithm) :: algo
    type(kpointmesh) :: mesh
    class(dp_resp) :: resp
    type(qp_resp), optional ::hpresp
   end subroutine

   subroutine derresp (icubic, algo, resp, hpresp)
    use response
    use types
    integer :: icubic
    type(algorithm) :: algo
    class(dp_resp) :: resp
    type(qp_resp),optional :: hpresp
   end subroutine

   subroutine wrtresp (iT, nalpha, algo, sct, resp, respinter, respBl, hpresp)
    use types
    use response
    integer, intent(in) :: iT, nalpha
    type(algorithm) :: algo
    type(scatrate)  :: sct
    type(dp_resp)   :: resp
    type(dp_respinter) :: respinter
    type(dp_resp)   :: respBl
    type(qp_resp),optional :: hpresp
   end subroutine

  end interface
  
  class(dp_resp), pointer :: pdpresp !local pointer to switch between datatypes
  type(algorithm) :: algo
  type(kpointmesh):: mesh
  type(edisp)     :: ek  
  type(tetramesh) :: thdr
  type(scatrate)  :: sct
  type(dp_resp),target   :: dresp       !intraband response
  type(dp_respinter),target  :: dderesp !derivatives of intraband conductivity 
  type(dp_respinter), target :: dinter  !interband response
  type(dp_resp)   :: respBl !Boltzmann response
  type(qp_resp), optional :: qresp
  real(8) :: mu,fac,facB
  real(8) :: drhodT(sct%nT)
  real(16):: facQ,facBQ
  integer :: iT,ib
  integer :: itet, ik
  integer :: ialpha,ibeta,idiag,icubic,nalpha,ia
  integer :: ktot
    
  complex(8),external  :: wpsipg
  complex(16),external :: wpsipghp

  real(8) :: tmp1, tmp2

  !what happens to the tight binding case? 
  nalpha=3
  if ((mesh%a(1)==mesh%a(2)) .and. (mesh%a(3)==mesh%a(2))) then
     icubic=1 
  else
     icubic=0
  endif
  if (icubic==1) nalpha=1
    
  !initialise the datatype variables
  ! eM: incredibly the qp_resp type seems to be there also when 
  ! it is not passed by the summoning routine. I think this is a 
  ! gfortran compiler glitch in the generation of the 
  ! implicit interface, using ldebug instead of the 
  ! intrinsic fortran present(qresp) solves the issue 
   
  if (algo%ldebug) then
     pdpresp => dresp
     call initresp (algo%lBfield, pdpresp, respBl)
     nullify(pdpresp)
     pdpresp => dderesp
     call initresp (.false., pdpresp, respBl)
     nullify(pdpresp)
     pdpresp => dinter
     call initresp (.false., pdpresp, respBl)
     nullify(pdpresp)
  else
     pdpresp => dresp
     call initresp (algo%lBfield, pdpresp, respBl, qresp)
     nullify(pdpresp)
     pdpresp => dderesp
     call initresp (.false., pdpresp, respBl, qresp)
     nullify(pdpresp)
     pdpresp => dinter
     call initresp (.false., pdpresp, respBl, qresp)
     nullify(pdpresp)
  endif


  ! outer k-loop   
  ! (quantities inside the tetrahedra will be interpolated over)
  ! eM: I decided to treat the tetrahedron and the regular k-mesh cases
  ! independently because later different parallelisation strategies can be 
  ! devised for the two approaches 

  if (algo%ltetra ) then
  !!!!!!!!!! TEST (vltot must add up to 1 -whole BZ-) 
     if (iT == sct%nT) then
        thdr%vltot=0.0d0
        tmp1=0.0d0
        do itet=1,thdr%ntet
           thdr%vltot=thdr%vltot+thdr%vltet(itet)
           do ik =1,4
              if (myid==0) write(33,*) itet, mesh%k_coord(:,thdr%idtet(ik,itet))
           enddo
           tmp1=tmp1+thdr%idtet(0,itet)
           write(34,*) itet, thdr%idtet(0,itet), thdr%vltet(itet)
        enddo
        tmp2=((2*pi)**3)/vol
        if (myid==0) write(*,*) 'tetrahedra volume',thdr%vltot,thdr%vltot*((2*pi)**3)/vol,tmp2
     endif
  !!!!!!!!!! TEST END

     do itet=iqstr,iqend

        !intraband transitions
        call respintet (mu, iT, itet, nalpha, thdr, algo, ek, sct, dresp ) 
        call respintet_Bl (mu, iT, itet, nalpha, thdr, algo, ek, sct, respBl ) 
        if (.not.algo%ldebug) then
           call respintet_qp (mu, iT, itet, nalpha, thdr, algo, ek, sct, qresp)
        endif

        !evaluate the derivatives of the response functions
        !if the chemical potential is fixed use a semplified kernel for the 
        !derivatives (assuming also gamma to be not T dependent)
        !if ((sct%Tstar == 0.0d0) .or. (sct%Tflat == 0.0d0)) then
           if (algo%ltbind .and. (algo%imurestart==2)) then
              !since one has to evaluate the derivative of mu then the first point must be skipped 
              if (iT < sct%nT) call resdertet_symm(mu, iT, itet, thdr, algo, ek, sct, dderesp) 
           else
              !since one has to evaluate the derivative of mu then the first point must be skipped 
              if (iT < sct%nT) call resdertet(iT, itet, nalpha, thdr, algo, ek, sct, dderesp) 
           endif
        !endif 

        !interband transitions
        if (algo%ldebug) then
           !!!!!!!!!!!!! TEST
           !if((iT==sct%nT) .and. (itet==1)) write(*,*)'test for 2-band symmetrical SC, check input parameters!!'
           !call respintert_symm(mu, iT, itet, nalpha, thdr, algo, ek, sct, dinter)
           !!!!!!!!!!!!! TEST END
           call respintert(mu, iT, itet, nalpha, thdr, algo, ek, sct, dinter) 
        else
           call respintert(mu, iT, itet, nalpha, thdr, algo, ek, sct, dinter) 
        endif


        pdpresp => dresp
        !interpolate within tetrahedra intraband response
        call interptra_mu (iT, itet, mu, nalpha, .false., mesh, ek, thdr, sct, pdpresp )   
        call interptra_mu (iT, itet, mu, nalpha, .true. , mesh, ek, thdr, sct, respBl )   
        if (.not.algo%ldebug) then
           call interptra_mu (iT, itet, mu, nalpha, .false., mesh, ek, thdr, sct, pdpresp, qresp )
        endif
        nullify(pdpresp)
        
        pdpresp => dderesp
        !interpolate within tetrahedra intraband response derivatives
        call interptra_mu (iT, itet, mu, nalpha, .false., mesh, ek, thdr, sct, pdpresp )   
        nullify(pdpresp)

        pdpresp => dinter
        !interpolate within tetrahedra interband response
        call interptra_mu (iT, itet, mu, nalpha, .false., mesh, ek, thdr, sct, pdpresp )   
        nullify(pdpresp)

        ! add to tetrahedra-summed response functions (bands have been traced over in INTERPTRA_mu )   
        ! 
        ! functions had to be multiplied by 2 for spin multiplicity
        
        do ialpha=1,nalpha
           do ibeta=ialpha,nalpha
              dresp%s_tot(ialpha,ibeta)=dresp%s_tot(ialpha,ibeta)+dresp%s_tet(ialpha,ibeta)*2.0d0
              dresp%a_tot(ialpha,ibeta)=dresp%a_tot(ialpha,ibeta)+dresp%a_tet(ialpha,ibeta)*2.0d0
              
              dderesp%s_tot(ialpha,ibeta)=dderesp%s_tot(ialpha,ibeta)+dderesp%s_tet(ialpha,ibeta)*2.0d0
              dderesp%a_tot(ialpha,ibeta)=dderesp%a_tot(ialpha,ibeta)+dderesp%a_tet(ialpha,ibeta)*2.0d0
              
              dinter%s_tot(ialpha,ibeta)=dinter%s_tot(ialpha,ibeta)+dinter%s_tet(ialpha,ibeta)*2.0d0
              dinter%a_tot(ialpha,ibeta)=dinter%a_tot(ialpha,ibeta)+dinter%a_tet(ialpha,ibeta)*2.0d0
              
              if (algo%lBfield .and. algo%ltbind ) then
                 dresp%sB_tot(ialpha,ibeta)=dresp%sB_tot(ialpha,ibeta)+dresp%sB_tet(ialpha,ibeta)*2.0d0
                 dresp%aB_tot(ialpha,ibeta)=dresp%aB_tot(ialpha,ibeta)+dresp%aB_tet(ialpha,ibeta)*2.0d0
              endif
              
              respBl%s_tot(ialpha,ibeta)=respBl%s_tot(ialpha,ibeta)+respBl%s_tet(ialpha,ibeta)*2.0d0
              respBl%a_tot(ialpha,ibeta)=respBl%a_tot(ialpha,ibeta)+respBl%a_tet(ialpha,ibeta)*2.0d0
              if (algo%lBfield .and. algo%ltbind ) then
                 respBl%sB_tot(ialpha,ibeta)=respBl%sB_tot(ialpha,ibeta)+respBl%sB_tet(ialpha,ibeta)*2.0d0
                 respBl%aB_tot(ialpha,ibeta)=respBl%aB_tot(ialpha,ibeta)+respBl%aB_tet(ialpha,ibeta)*2.0d0
              endif
              
              if (.not.algo%ldebug) then
                 qresp%s_tot(ialpha,ibeta)=qresp%s_tot(ialpha,ibeta)+qresp%s_tet(ialpha,ibeta)*2.0q0
                 qresp%a_tot(ialpha,ibeta)=qresp%a_tot(ialpha,ibeta)+qresp%a_tet(ialpha,ibeta)*2.0q0
                 if (algo%lBfield .and. algo%ltbind) then
                    qresp%sB_tot(ialpha,ibeta)=qresp%sB_tot(ialpha,ibeta)+qresp%sB_tet(ialpha,ibeta)*2.0q0
                    qresp%aB_tot(ialpha,ibeta)=qresp%aB_tot(ialpha,ibeta)+qresp%aB_tet(ialpha,ibeta)*2.0q0
                 endif
              endif
           enddo !ibeta 
        enddo    !ialpha
     enddo ! loop over tetrahedra 
     
     !Distribute the accumulated variables and sum them up
     if (nproc > 1) then
        call MPI_ALLREDUCE(MPI_IN_PLACE, dresp%s_tot(1,1), 9, &
               MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, mpierr)
        call MPI_ALLREDUCE(MPI_IN_PLACE, dresp%a_tot(1,1), 9, &
               MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, mpierr)
        if (algo%lBfield .and. algo%ltbind ) then
           call MPI_ALLREDUCE(MPI_IN_PLACE, dresp%sB_tot(1,1), 9, &
                  MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, mpierr)
           call MPI_ALLREDUCE(MPI_IN_PLACE, dresp%aB_tot(1,1), 9, &
                  MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, mpierr)
        endif
        !derivative
        call MPI_ALLREDUCE(MPI_IN_PLACE, dderesp%s_tot(1,1), 9, &
               MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, mpierr)
        call MPI_ALLREDUCE(MPI_IN_PLACE, dderesp%a_tot(1,1), 9, &
               MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, mpierr)
        !interband
        call MPI_ALLREDUCE(MPI_IN_PLACE, dinter%s_tot(1,1), 9, &
               MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, mpierr)
        call MPI_ALLREDUCE(MPI_IN_PLACE, dinter%a_tot(1,1), 9, &
               MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, mpierr)
        !Boltzmann
        call MPI_ALLREDUCE(MPI_IN_PLACE, respBl%s_tot(1,1), 9, &
               MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, mpierr)
        call MPI_ALLREDUCE(MPI_IN_PLACE, respBl%a_tot(1,1), 9, &
               MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, mpierr)
        if (algo%lBfield .and. algo%ltbind ) then
           call MPI_ALLREDUCE(MPI_IN_PLACE, respBl%sB_tot(1,1), 9, &
                  MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, mpierr)
           call MPI_ALLREDUCE(MPI_IN_PLACE, respBl%aB_tot(1,1), 9, &
                  MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, mpierr)
        endif
        !intraband QP
        if (.not.algo%ldebug) then
           do ialpha=1,3
              do ibeta=1,3
                 call mpi_reduce_quad(qresp%s_tot(ialpha,ibeta),qresp%s_tot(ialpha,ibeta))
                 call mpi_reduce_quad(qresp%a_tot(ialpha,ibeta),qresp%a_tot(ialpha,ibeta))
              enddo
           enddo
           if (algo%lBfield .and. algo%ltbind ) then
              do ialpha=1,3
                 do ibeta=1,3
                    call mpi_reduce_quad(qresp%sB_tot(ialpha,ibeta),qresp%sB_tot(ialpha,ibeta))
                    call mpi_reduce_quad(qresp%aB_tot(ialpha,ibeta),qresp%aB_tot(ialpha,ibeta))
                 enddo
              enddo
           endif
        endif 
     endif   !nproc>1
        
  else ! no tetrahedron method
        
     if ((nproc > 1) .and.(.not. allocated(dresp%s_local))) then
        allocate (dresp%s_local(ek%nband_max, nalpha, nalpha) )
     endif
     if ( allocated(dresp%s_local)) then
        dresp%s_local(:,:,:) = 0.0d0
     endif
     do ik =iqstr,iqend 
     ! evaluate the trace over bands at each specific k-point of the mesh
        call respinkm (mu, iT, ik, nalpha, algo, ek, sct, dresp)
        call respinkm_Bl (mu, iT, ik, nalpha, algo, ek, sct, respBl)
        if (.not.algo%ldebug) then
           call respinkm_qp (mu, iT, ik, nalpha, algo, ek, sct, qresp)
        endif
        !evaluate the derivatives of the response functions
        !if the chemical potential is fixed use a semplified kernel for the 
        !derivatives (assuming also gamma to be not T dependent)
        !if ((sct%Tstar == 0.0d0) .or. (sct%Tflat == 0.0d0)) then
           if (algo%ltbind .and. (algo%imurestart==2)) then
              if (iT < sct%nT) call resderkm_symm(mu, iT, ik, algo, ek, sct, dderesp) 
           else
              !since one has to evaluate the derivative of mu then the first point must be skipped 
              if (iT < sct%nT) call resderkm(iT, ik, nalpha, algo, ek, sct, dderesp) 
           endif
        !endif

        !interband transitions
        if (algo%ldebug) then
           !!!!!!!!!!!!! TEST
           !if((iT==sct%nT) .and. (itet==1)) write(*,*)'test for 2-band symmetrical SC, check input parameters!!'
           !call respinterkm_symm(mu, iT, ik, nalpha, algo, ek, sct, dinter)
           !!!!!!!!!!!!! TEST END
           call respinterkm(mu, iT, ik, nalpha, algo, ek, sct, dinter) 
        else
           call respinterkm(mu, iT, ik, nalpha, algo, ek, sct, dinter) 
        endif
     enddo

     if (nproc == 1) then
        do ialpha=1,nalpha
           do ibeta=ialpha,nalpha
              do ib=1,ek%nband_max
                 if(ib<ek%nbopt_min) cycle
                 if(ib>ek%nbopt_max) cycle
                 do ik =iqstr,iqend 
                    !multiply by spin multiplicity 
                    dresp%s(ib,ialpha,ibeta)=dresp%s(ib,ialpha,ibeta)+dresp%s_tmp(ik,ib,ialpha,ibeta)*2.0d0
                    dresp%a(ib,ialpha,ibeta)=dresp%a(ib,ialpha,ibeta)+dresp%a_tmp(ik,ib,ialpha,ibeta)*2.0d0
                    dderesp%s(ib,ialpha,ibeta)=dderesp%s(ib,ialpha,ibeta)+dderesp%s_tmp(ik,ib,ialpha,ibeta)*2.0d0
                    dderesp%a(ib,ialpha,ibeta)=dderesp%a(ib,ialpha,ibeta)+dderesp%a_tmp(ik,ib,ialpha,ibeta)*2.0d0
                    dinter%s(ib,ialpha,ibeta)=dinter%s(ib,ialpha,ibeta)+dinter%s_tmp(ik,ib,ialpha,ibeta)*2.0d0
                    dinter%a(ib,ialpha,ibeta)=dinter%a(ib,ialpha,ibeta)+dinter%a_tmp(ik,ib,ialpha,ibeta)*2.0d0
                    if (algo%lBfield .and. algo%ltbind ) then
                       dresp%sB(ib,ialpha,ibeta)=dresp%sB(ib,ialpha,ibeta)+dresp%sB_tmp(ik,ib,ialpha,ibeta)*2.0d0
                       dresp%aB(ib,ialpha,ibeta)=dresp%aB(ib,ialpha,ibeta)+dresp%aB_tmp(ik,ib,ialpha,ibeta)*2.0d0
                    endif
                 
                    respBl%s(ib,ialpha,ibeta)=respBl%s(ib,ialpha,ibeta)+respBl%s_tmp(ik,ib,ialpha,ibeta)*2.0d0
                    respBl%a(ib,ialpha,ibeta)=respBl%a(ib,ialpha,ibeta)+respBl%a_tmp(ik,ib,ialpha,ibeta)*2.0d0
                    if (algo%lBfield .and. algo%ltbind ) then
                       respBl%sB(ib,ialpha,ibeta)=respBl%sB(ib,ialpha,ibeta)+respBl%sB_tmp(ik,ib,ialpha,ibeta)*2.0d0
                       respBl%aB(ib,ialpha,ibeta)=respBl%aB(ib,ialpha,ibeta)+respBl%aB_tmp(ik,ib,ialpha,ibeta)*2.0d0
                    endif
                 
                    if (.not.algo%ldebug) then
                       qresp%s(ib,ialpha,ibeta)=qresp%s(ib,ialpha,ibeta)+qresp%s_tmp(ik,ib,ialpha,ibeta)*2.0q0
                       qresp%a(ib,ialpha,ibeta)=qresp%a(ib,ialpha,ibeta)+qresp%a_tmp(ik,ib,ialpha,ibeta)*2.0q0
                       if (algo%lBfield .and. algo%ltbind ) then
                          qresp%sB(ib,ialpha,ibeta)=qresp%sB(ib,ialpha,ibeta)+qresp%sB_tmp(ik,ib,ialpha,ibeta)*2.0q0
                          qresp%aB(ib,ialpha,ibeta)=qresp%aB(ib,ialpha,ibeta)+qresp%aB_tmp(ik,ib,ialpha,ibeta)*2.0q0
                       endif
                    endif
                 enddo !over kpoints
              enddo !nbands
           enddo !ibeta
        enddo !ialpha

     else !nproc > 1

        do ialpha=1,nalpha
           do ibeta=ialpha,nalpha
              do ib=1,ek%nband_max
                 if(ib<ek%nbopt_min) cycle
                 if(ib>ek%nbopt_max) cycle
                 !accumulate locally over k-points, then reduce over cores
                 do ik=iqstr,iqend 
                    dresp%s_local(ib,ialpha,ibeta)=dresp%s_local(ib,ialpha,ibeta)+dresp%s_tmp(ik,ib,ialpha,ibeta)*2.0d0
                    dresp%a_local(ib,ialpha,ibeta)=dresp%a_local(ib,ialpha,ibeta)+dresp%a_tmp(ik,ib,ialpha,ibeta)*2.0d0
                    dderesp%s_local(ib,ialpha,ibeta)=dderesp%s_local(ib,ialpha,ibeta)+dderesp%s_tmp(ik,ib,ialpha,ibeta)*2.0d0
                    dderesp%a_local(ib,ialpha,ibeta)=dderesp%a_local(ib,ialpha,ibeta)+dderesp%a_tmp(ik,ib,ialpha,ibeta)*2.0d0
                    dinter%s_local(ib,ialpha,ibeta)=dinter%s_local(ib,ialpha,ibeta)+dinter%s_tmp(ik,ib,ialpha,ibeta)*2.0d0
                    dinter%a_local(ib,ialpha,ibeta)=dinter%a_local(ib,ialpha,ibeta)+dinter%a_tmp(ik,ib,ialpha,ibeta)*2.0d0
                    if (algo%lBfield .and. algo%ltbind ) then
                       dresp%sB_local(ib,ialpha,ibeta)=dresp%sB_local(ib,ialpha,ibeta)+dresp%sB_tmp(ik,ib,ialpha,ibeta)*2.0d0
                       dresp%aB_local(ib,ialpha,ibeta)=dresp%aB_local(ib,ialpha,ibeta)+dresp%aB_tmp(ik,ib,ialpha,ibeta)*2.0d0
                    endif
                 
                    respBl%s_local(ib,ialpha,ibeta)=respBl%s_local(ib,ialpha,ibeta)+respBl%s_tmp(ik,ib,ialpha,ibeta)*2.0d0
                    respBl%a_local(ib,ialpha,ibeta)=respBl%a_local(ib,ialpha,ibeta)+respBl%a_tmp(ik,ib,ialpha,ibeta)*2.0d0
                    if (algo%lBfield .and. algo%ltbind ) then
                       respBl%sB_local(ib,ialpha,ibeta)=respBl%sB_local(ib,ialpha,ibeta)+respBl%sB_tmp(ik,ib,ialpha,ibeta)*2.0d0
                       respBl%aB_local(ib,ialpha,ibeta)=respBl%aB_local(ib,ialpha,ibeta)+respBl%aB_tmp(ik,ib,ialpha,ibeta)*2.0d0
                    endif
                 
                    if (.not.algo%ldebug) then
                       qresp%s_local(ib,ialpha,ibeta)=qresp%s_local(ib,ialpha,ibeta)+qresp%s_tmp(ik,ib,ialpha,ibeta)*2.0q0
                       qresp%a_local(ib,ialpha,ibeta)=qresp%a_local(ib,ialpha,ibeta)+qresp%a_tmp(ik,ib,ialpha,ibeta)*2.0q0
                       if (algo%lBfield .and. algo%ltbind ) then
                          qresp%sB_local(ib,ialpha,ibeta)=qresp%sB_local(ib,ialpha,ibeta)+qresp%sB_tmp(ik,ib,ialpha,ibeta)*2.0q0
                          qresp%aB_local(ib,ialpha,ibeta)=qresp%aB_local(ib,ialpha,ibeta)+qresp%aB_tmp(ik,ib,ialpha,ibeta)*2.0q0
                       endif
                    endif
                 enddo
                 !intraband
                 call MPI_REDUCE(dresp%s_local(ib,ialpha,ibeta), dresp%s(ib,ialpha,ibeta), 1, &
                        MPI_REAL8, MPI_SUM, master, MPI_COMM_WORLD, mpierr)
                 call MPI_REDUCE(dresp%a_local(ib,ialpha,ibeta), dresp%a(ib,ialpha,ibeta), 1, &
                        MPI_REAL8, MPI_SUM, master, MPI_COMM_WORLD, mpierr)
                 if (algo%lBfield .and. algo%ltbind ) then
                    call MPI_REDUCE(dresp%sB_local(ib,ialpha,ibeta), dresp%sB(ib,ialpha,ibeta), 1, &
                           MPI_REAL8, MPI_SUM, master, MPI_COMM_WORLD, mpierr)
                    call MPI_REDUCE(dresp%aB_local(ib,ialpha,ibeta), dresp%aB(ib,ialpha,ibeta), 1, &
                           MPI_REAL8, MPI_SUM, master, MPI_COMM_WORLD, mpierr)
                 endif
                 !derivative
                 call MPI_REDUCE(dderesp%s_local(ib,ialpha,ibeta), dderesp%s(ib,ialpha,ibeta), 1, &
                        MPI_REAL8, MPI_SUM, master, MPI_COMM_WORLD, mpierr)
                 call MPI_REDUCE(dderesp%a_local(ib,ialpha,ibeta), dderesp%a(ib,ialpha,ibeta), 1, &
                        MPI_REAL8, MPI_SUM, master, MPI_COMM_WORLD, mpierr)
                 !interband
                 call MPI_REDUCE(dinter%s_local(ib,ialpha,ibeta), dinter%s(ib,ialpha,ibeta), 1, &
                        MPI_REAL8, MPI_SUM, master, MPI_COMM_WORLD, mpierr)
                 call MPI_REDUCE(dinter%a_local(ib,ialpha,ibeta), dinter%a(ib,ialpha,ibeta), 1, &
                        MPI_REAL8, MPI_SUM, master, MPI_COMM_WORLD, mpierr)
                 !Boltzmann
                 call MPI_REDUCE(respBl%s_local(ib,ialpha,ibeta), respBl%s(ib,ialpha,ibeta), 1, &
                        MPI_REAL8, MPI_SUM, master, MPI_COMM_WORLD, mpierr)
                 call MPI_REDUCE(respBl%a_local(ib,ialpha,ibeta), respBl%a(ib,ialpha,ibeta), 1, &
                        MPI_REAL8, MPI_SUM, master, MPI_COMM_WORLD, mpierr)
                 if (algo%lBfield .and. algo%ltbind ) then
                    call MPI_REDUCE(respBl%sB_local(ib,ialpha,ibeta), respBl%sB(ib,ialpha,ibeta), 1, &
                           MPI_REAL8, MPI_SUM, master, MPI_COMM_WORLD, mpierr)
                    call MPI_REDUCE(respBl%aB_local(ib,ialpha,ibeta), respBl%aB(ib,ialpha,ibeta), 1, &
                           MPI_REAL8, MPI_SUM, master, MPI_COMM_WORLD, mpierr)
                 endif
                 !intraband QP
                 if (.not.algo%ldebug) then
                    call mpi_reduce_quad(qresp%s_local(ib,ialpha,ibeta),qresp%s(ib,ialpha,ibeta))  
                    call mpi_reduce_quad(qresp%a_local(ib,ialpha,ibeta),qresp%a(ib,ialpha,ibeta))  
                    if (algo%lBfield .and. algo%ltbind ) then
                       call mpi_reduce_quad(qresp%sB_local(ib,ialpha,ibeta),qresp%sB(ib,ialpha,ibeta))  
                       call mpi_reduce_quad(qresp%aB_local(ib,ialpha,ibeta),qresp%aB(ib,ialpha,ibeta))  
                    endif
                 endif
              enddo !nbands
           enddo !ibeta
        enddo !ialpha
     endif !nproc

     !the loop below has to be protected with myid.eq.master because the band resolved variables have been 
     !accumulated already
     if (myid.eq.master) then  
        ! at this point I need to multiply the factor beta/gamma (or beta/gamma^2 for magnetic responses)
        ! in the tetrahedron method this is done in the INTERPTRA_mu routine
        do ib=1,ek%nband_max
           if (ib < ek%nbopt_min) cycle
           if (ib > ek%nbopt_max) cycle
           dresp%gamma=real(sct%gam(iT,ib)*sct%z,8)
           if (.not.algo%ldebug) qresp%gamma=real(sct%gam(iT,ib)*sct%z,16)
        
           dresp%s(ib,:,:)   = dresp%s(ib,:,:)*beta/dresp%gamma
           dresp%a(ib,:,:)   = dresp%a(ib,:,:)*(beta**2)/dresp%gamma
           dresp%s_tot(:,:)  = dresp%s_tot(:,:) + dresp%s(ib,:,:)
           dresp%a_tot(:,:)  = dresp%a_tot(:,:) + dresp%a(ib,:,:)
           !same treatment for the derivatives...
           dderesp%s_tot(:,:)= dderesp%s_tot(:,:) + dderesp%s(ib,:,:)
           dderesp%a_tot(:,:)= dderesp%a_tot(:,:) + dderesp%a(ib,:,:)
           !and for the interband response
           dinter%s_tot(:,:)= dinter%s_tot(:,:) + dinter%s(ib,:,:)
           dinter%a_tot(:,:)= dinter%a_tot(:,:) + dinter%a(ib,:,:)
           if (algo%lBfield .and. algo%ltbind ) then
              dresp%sB(ib,:,:) = dresp%sB(ib,:,:)*beta/(dresp%gamma**2)
              dresp%aB(ib,:,:) = dresp%aB(ib,:,:)*(beta**2)/(dresp%gamma**2)
              dresp%sB_tot(:,:)= dresp%sB_tot(:,:) + dresp%sB(ib,:,:)
              dresp%aB_tot(:,:)= dresp%aB_tot(:,:) + dresp%aB(ib,:,:)
           endif
        
           respBl%s(ib,:,:) = respBl%s(ib,:,:)/dresp%gamma
           respBl%a(ib,:,:) = respBl%a(ib,:,:)/dresp%gamma
           respBl%s_tot(:,:)= respBl%s_tot(:,:) + respBl%s(ib,:,:)
           respBl%a_tot(:,:)= respBl%a_tot(:,:) + respBl%a(ib,:,:)
           if (algo%lBfield .and. algo%ltbind ) then
              respBl%sB(ib,:,:) = respBl%sB(ib,:,:)/(dresp%gamma**2)
              respBl%aB(ib,:,:) = respBl%aB(ib,:,:)/(dresp%gamma**2)
              respBl%sB_tot(:,:)= respBl%sB_tot(:,:) + dresp%sB(ib,:,:)
              respBl%aB_tot(:,:)= respBl%aB_tot(:,:) + dresp%aB(ib,:,:)
           endif
        
           if (.not.algo%ldebug) then
              qresp%s(ib,:,:) = qresp%s(ib,:,:)*betaQ/qresp%gamma
              qresp%a(ib,:,:) = qresp%a(ib,:,:)*(betaQ**2)/qresp%gamma
              qresp%s_tot(:,:)= qresp%s_tot(:,:) + qresp%s(ib,:,:)
              qresp%a_tot(:,:)= qresp%a_tot(:,:) + qresp%a(ib,:,:)
              if (algo%lBfield .and. algo%ltbind ) then
                 qresp%sB(ib,:,:) = qresp%sB(ib,:,:)*betaQ/(qresp%gamma**2)
                 qresp%aB(ib,:,:) = qresp%aB(ib,:,:)*betaQ/(qresp%gamma**2)
                 qresp%sB_tot(:,:)= qresp%sB_tot(:,:) + qresp%sB(ib,:,:)
                 qresp%aB_tot(:,:)= qresp%aB_tot(:,:) + qresp%aB(ib,:,:)
              endif
           endif
        
        enddo !over bands
     endif

  endif !ltetra

  if (myid.eq.master) then
  ! At this stage solve the equation (d^2 rho)/(d beta^2) =0
  ! since the following multiplications in globfac affect equally the 
  ! conductivity and its derivatives there is  no impact on the 
  ! resulting T* temperature (T at which the resistivity has an inflection)
     !if (sct%Tstar == 0.0d0) then
        if (iT < sct%nT) call findrhoflex (iT, dresp, dderesp, sct)
     !endif
     if (sct%Tflat == 0.0d0) then
        if (iT < sct%nT) call findrhoflat (iT, dresp, dderesp, sct)
     endif
     !if ((iT<sct%nT) .and. (iT>1)) write(*,*)'calling finddrhomax'
     if ((iT<sct%nT) .and. (iT>1)) call finddrhomax (iT, dresp, dderesp, sct, drhodT)
  ! At this point the values in the response datatypes should be consistent (in terms of prefactors/dimensionality)
  ! there are some global factors missing that are taken care of in the following routine 
     pdpresp => dresp
     if (algo%ldebug) then
        call globfac(icubic, algo, mesh, pdpresp)  
     else 
        call globfac(icubic, algo, mesh, pdpresp, qresp)  
     endif
     nullify(pdpresp)
     call globfac(icubic, algo, mesh, respBl) 
     pdpresp => dinter
     call globfac(icubic, algo, mesh, pdpresp)  
     nullify(pdpresp)

  ! The derived variables (Seebeck, Nernst, R_H) are computed in the following routine
     pdpresp => dresp
     if (algo%ldebug) then
        call derresp(icubic, algo, pdpresp) 
     else
        call derresp(icubic, algo, pdpresp, qresp) 
     endif
     nullify(pdpresp)
     call derresp(icubic, algo, respBl) 
     pdpresp => dinter
     call derresp(icubic, algo, pdpresp) 
     nullify(pdpresp)
     
     if (algo%ldebug) then
        call wrtresp(iT, nalpha, algo, sct, dresp, dinter, respBl)
     else
        call wrtresp(iT, nalpha, algo, sct, dresp, dinter, respBl, qresp)
     endif
  endif
end subroutine calc_response

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! INITRESP
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine initialises the datatypes
!
subroutine initresp (lBfield, dresp, respBl, qresp)
  use response
  use types
  implicit none
  logical :: lBfield
  class(dp_resp) :: dresp
  type(dp_resp) :: respBl
  type(qp_resp),optional :: qresp

  dresp%s_tot = 0.0d0
  dresp%s  = 0.d0; dresp%s_tmp  = 0.d0
  dresp%a_tot = 0.0d0
  dresp%a  = 0.d0; dresp%a_tmp  = 0.d0
  if (lBfield) then
     dresp%sB_tot = 0.0d0
     dresp%sB = 0.d0; dresp%sB_tmp = 0.d0
     dresp%aB_tot = 0.0d0
     dresp%aB = 0.d0; dresp%aB_tmp = 0.d0
  endif

  respBl%s_tot = 0.0d0 
  respBl%s  = 0.d0; respBl%s_tmp  = 0.d0
  respBl%a_tot = 0.0d0
  respBl%a  = 0.d0; respBl%a_tmp  = 0.d0
  if (lBfield) then
     respBl%sB_tot = 0.0d0 
     respBl%sB = 0.d0; respBl%sB_tmp = 0.d0
     respBl%aB_tot = 0.0d0
     respBl%aB = 0.d0; respBl%aB_tmp = 0.d0
  endif

  !local (on core) accummulation variables
  !for mpi version
  if (allocated(dresp%s_local)) then
     dresp%s_local  = 0.d0; dresp%a_local  = 0.d0
     respBl%s_local = 0.d0; respBl%a_local = 0.d0
     if (lBfield) then
        dresp%sB_local  = 0.d0; dresp%aB_local  = 0.d0
        respBl%sB_local = 0.d0; respBl%aB_local = 0.d0
     endif
  endif

  if (present(qresp)) then
     qresp%s_tot = 0.0q0 
     qresp%s  = 0.q0; qresp%s_tmp  = 0.q0
     qresp%a_tot = 0.0q0
     qresp%a  = 0.q0; qresp%a_tmp  = 0.q0
     if (lBfield) then
        qresp%sB_tot = 0.0q0 
        qresp%sB = 0.q0; qresp%sB_tmp = 0.q0
        qresp%aB_tot = 0.0q0
        qresp%aB = 0.q0; qresp%aB_tmp = 0.q0
     endif
     !local (on core) accummulation variables
     !for mpi version
     if (allocated(dresp%s_local)) then
        qresp%s_local  = 0.d0; qresp%a_local  = 0.d0
        if (lBfield) then
           qresp%sB_local  = 0.d0; qresp%aB_local  = 0.d0
        endif
     endif
  endif 
end subroutine initresp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! RESPINTET
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This routine (and the following one in QP)
! evaluate the conduction kernel, and the 
! full response functions on the vertices of a
! tetrahedron in the KUBO formalism
!
subroutine respintet(mu, iT, itet, nalpha, thdr, algo, ek, sct, resp)
  use response
  use types
  use params
  implicit none 
  type (dp_resp) :: resp ! dynamical datatype allow only for an inclusion through extension of the parent type,
                         ! I could have declared a dummy class pointer that would be assigned to either dp or qp response type.
                         ! To do so, the varaibles defined in the individual types must have had different names and since I
                         ! REALLY dislike having that extra Q for each varible I decided to stick to static datatypes.
  type (algorithm) :: algo
  type (tetramesh) :: thdr
  type (edisp) :: ek
  type (scatrate) :: sct
  real(8), intent(in) :: mu
  integer, intent(in) :: iT
  integer, intent(in) :: itet
  integer, intent(in) :: nalpha
  integer :: iband, ik, ipg
  integer :: ialpha,ibeta 
  complex(8),external  :: wpsipg
  complex(16),external :: wpsipghp
!local variables
  real(8), allocatable :: s_tmp_tetra(:,:,:,:),  a_tmp_tetra(:,:,:,:)
  real(8), allocatable :: sB_tmp_tetra(:,:,:,:), aB_tmp_tetra(:,:,:,:)

  !allocation  
  if(.not. allocated(s_tmp_tetra)) allocate(s_tmp_tetra(4,ek%nband_max,3,3))
  if(.not. allocated(a_tmp_tetra)) allocate(a_tmp_tetra(4,ek%nband_max,3,3))
  if (algo%lBfield .and. algo%ltbind) then
     if(.not. allocated(sB_tmp_tetra)) allocate(sB_tmp_tetra(4,ek%nband_max,3,3))
     if(.not. allocated(aB_tmp_tetra)) allocate(aB_tmp_tetra(4,ek%nband_max,3,3))
  endif
  !initialisation
  s_tmp_tetra=0.0d0 ; a_tmp_tetra=0.0d0
  if (algo%lBfield .and. algo%ltbind) then
     sB_tmp_tetra=0.0d0 ; aB_tmp_tetra=0.0d0
  endif

   do ik=1,4  !loop over corners of the tetrahedron

      do iband=1,ek%nband_max !loop over bands (these will be traced over)

         ! if the band is not contained in the optical matrices just do nothing
         if (iband < ek%nbopt_min) cycle
         if (iband > ek%nbopt_max) cycle
         resp%z=real(sct%z,8)
         resp%gamma=resp%z*real(sct%gam(iT,iband),8)
         ! pre-compute all needed digamma functions   
         resp%aqp=resp%z*real(ek%band(thdr%idtet(ik,itet),iband)-mu,8)
         resp%zarg=0.5d0+beta2p*(ci*resp%aqp+resp%gamma)
         do ipg=1,3 ! XXX need 0 for alphaxy ????
            resp%ctmp=wpsipg(resp%zarg,ipg)
            resp%RePolyGamma(ipg)=real(resp%ctmp,8)
            resp%ImPolyGamma(ipg)=imag(resp%ctmp)
         enddo
         
         ! compute transport kernels (omega-part)
         ! 
         resp%tmp=resp%z**2 / (4.d0*pi**3) ! missing: beta/gamma (multiplied later to keep numbers reasonable here)        
         resp%s_ker = resp%tmp * (resp%RePolyGamma(1) - resp%gamma*beta2p * resp%RePolyGamma(2) )
         resp%a_ker = resp%tmp * ( resp%aqp * resp%RePolyGamma(1) - resp%aqp*resp%gamma*beta2p*resp%RePolyGamma(2) &
                    - resp%gamma**2.d0 * beta2p * resp%ImPolyGamma(2) )
         
         
         resp%tmp=resp%tmp*3.d0*resp%z/(4.d0*pi) ! additionally missing: 1/gamma (multiplied later) XXX

         if(algo%lBfield .and. algo%ltbind) then        
            resp%sB_ker = resp%tmp*(-resp%RePolyGamma(1)-resp%gamma*beta2p*resp%RePolyGamma(2)-(beta2p*resp%gamma)**2.d0/3.d0 &
                          * resp%RePolyGamma(3))
            
            resp%aB_ker=resp%tmp*(resp%aqp*resp%RePolyGamma(1)-resp%aqp*beta2p*resp%gamma*resp%RePolyGamma(2)+resp%gamma**2/3.d0 &
                 * beta2p * resp%ImPolyGamma(2) &
                 - resp%aqp*resp%gamma**2.d0 / 3.d0 * beta2p**2.d0 * resp%ImPolyGamma(3) &
                 + resp%gamma**3.d0 / 3.d0 * beta2p**2.d0 * resp%RePolyGamma(3) )
          endif 
         
         ! B = 0  
         !tmp=vka(ik,ib,ialpha)*vka(ik,ib,ibeta)
         do ialpha=1,nalpha
         
            if (algo%ltbind) then
               resp%tmp=ek%Mopt(ialpha,thdr%idtet(ik,itet), iband, iband)*ek%Mopt(ialpha,thdr%idtet(ik,itet), iband, iband)
            else
               resp%tmp=ek%Mopt(ialpha,thdr%idtet(ik,itet), iband, iband) !the optical matrix elements given by Wien2k are squared already
            endif
         
            s_tmp_tetra(ik,iband,ialpha,ialpha)=resp%s_ker * resp%tmp  
            a_tmp_tetra(ik,iband,ialpha,ialpha)=resp%a_ker * resp%tmp 
               
            do ibeta=ialpha+1,nalpha
               resp%tmp=ek%Mopt(ialpha+ibeta+1,thdr%idtet(ik,itet), iband, iband) !the optical matrix elements given by Wien2k are squared already
               s_tmp_tetra(ik,iband,ialpha,ibeta)=resp%s_ker * resp%tmp  
               a_tmp_tetra(ik,iband,ialpha,ibeta)=resp%a_ker * resp%tmp 
            enddo !ibeta  
         enddo ! ialpha           
      
         ! B .ne. 0 
         !tmp=vka(ik,ib,ialpha)*( vkab(ik,ib,ibeta,ialpha)*vka(ik,ib,ibeta) - vkab(ik,ib,ibeta,ibeta)*vka(ik,ib,ialpha)    )
         if (algo%lBfield .and. algo%ltbind ) then
            do ialpha=1,nalpha
               do ibeta=ialpha+1,3
                  
                  resp%tmp =ek%Mopt(ialpha,thdr%idtet(ik,itet), iband, iband)* &
                    (ek%M2(ibeta, ialpha, thdr%idtet(ik,itet), iband)*ek%Mopt(ialpha,thdr%idtet(ik,itet), iband, iband) - &
                    ek%M2(ibeta, ibeta, thdr%idtet(ik,itet), iband)* ek%Mopt(ialpha,thdr%idtet(ik,itet), iband, iband) )
                  sB_tmp_tetra(ik,iband,ialpha,ibeta)=resp%sB_ker * resp%tmp 
                  aB_tmp_tetra(ik,iband,ialpha,ibeta)=resp%aB_ker * resp%tmp 
                  
               enddo !ibeta  
            enddo ! ialpha
         endif !lBfield

         ! Now copy the local variable into the datastructure that will be passed to the interptra_mu
         resp%s_tmp(ik,iband,:,:) = s_tmp_tetra(ik,iband,:,:)
         resp%a_tmp(ik,iband,:,:) = a_tmp_tetra(ik,iband,:,:)
         if(algo%lBfield .and. algo%ltbind ) then
            resp%sB_tmp(ik,iband,:,:) = sB_tmp_tetra(ik,iband,:,:)
            resp%aB_tmp(ik,iband,:,:) = aB_tmp_tetra(ik,iband,:,:)
         endif
         
      enddo ! iband

   enddo ! ik   

  !!!!!!!!!!!!! TEST
  ! if ((mod(itet,100) ==0) .and. (sct%TT(iT)==90.00)) then
  !   write(itet+1,*) s_tmp_tetra(1,:,1,1), a_tmp_tetra(1,:,1,1)
  !   write(itet+2,*) s_tmp_tetra(2,:,1,1), a_tmp_tetra(2,:,1,1)
  !   write(itet+3,*) s_tmp_tetra(3,:,1,1), a_tmp_tetra(3,:,1,1)
  !   write(itet+4,*) s_tmp_tetra(4,:,1,1), a_tmp_tetra(4,:,1,1)
  !   write(itet+10,*) resp%s_tmp(1,:,1,1), resp%a_tmp(1,:,1,1)
  !   write(itet+20,*) resp%s_tmp(2,:,1,1), resp%a_tmp(2,:,1,1)
  !   write(itet+30,*) resp%s_tmp(3,:,1,1), resp%a_tmp(3,:,1,1)
  !   write(itet+40,*) resp%s_tmp(4,:,1,1), resp%a_tmp(4,:,1,1)
  ! endif
  !!!!!!!!!!!!! TEST END

  deallocate(s_tmp_tetra)
  deallocate(a_tmp_tetra)
  if (allocated(sB_tmp_tetra)) then
     deallocate(sB_tmp_tetra)
     deallocate(aB_tmp_tetra)
  endif

end subroutine respintet

subroutine respintet_qp(mu, iT, itet, nalpha, thdr, algo, ek, sct, resp)
  use response
  use types
  use params
  implicit none 
  type (qp_resp) :: resp ! dynamical datatype allow only for an inclusion through extension of the parenttype
                         ! I could have declared a dummy class pointer that would be assigned to either dp or qp response type
                         ! to do so the varaibles defined in the individual types must have had different names and since I
                         ! REALLY dislike having that extra Q for each varible I decided to stick to static datatypes
  type (algorithm) :: algo
  type (tetramesh) :: thdr
  type (edisp) :: ek
  type (scatrate) :: sct
  real(8), intent(in) :: mu
  integer, intent(in) :: iT
  integer, intent(in) :: itet
  integer, intent(in) :: nalpha
  integer :: iband, ik, ipg
  integer :: ialpha,ibeta
  complex(8),external  :: wpsipg
  complex(16),external :: wpsipghp
!local variables
  real(16), allocatable :: s_tmp_tetra(:,:,:,:),  a_tmp_tetra(:,:,:,:)
  real(16), allocatable :: sB_tmp_tetra(:,:,:,:), aB_tmp_tetra(:,:,:,:)

  !allocation  
  if(.not. allocated(s_tmp_tetra)) allocate(s_tmp_tetra(4,ek%nband_max,3,3))
  if(.not. allocated(a_tmp_tetra)) allocate(a_tmp_tetra(4,ek%nband_max,3,3))
  if (algo%lBfield .and. algo%ltbind) then
     if(.not. allocated(sB_tmp_tetra)) allocate(sB_tmp_tetra(4,ek%nband_max,3,3))
     if(.not. allocated(aB_tmp_tetra)) allocate(aB_tmp_tetra(4,ek%nband_max,3,3))
  endif
  !initialisation
  s_tmp_tetra=0.0q0 ; a_tmp_tetra=0.0q0
  if (algo%lBfield .and. algo%ltbind) then
     sB_tmp_tetra=0.0q0 ; aB_tmp_tetra=0.0q0
  endif

  betaQ=real(beta,16)
  beta2pQ=betaQ/(2.q0*piQ)


   do ik=1,4  !loop over corners of the tetrahedron

      do iband=1,ek%nband_max !loop over bands (these will be traced over)

         ! if the band is not contained in the optical matrices just do nothing
         if (iband < ek%nbopt_min) cycle
         if (iband > ek%nbopt_max) cycle
         resp%z=real(sct%z,16)
         resp%gamma=resp%z*real(sct%gam(iT,iband),16)
         ! pre-compute all needed digamma functions   
         resp%aqp=resp%z*real(ek%band(thdr%idtet(ik,itet),iband)-mu,16)
         resp%zarg=0.5q0+beta2pQ*(ciQ*resp%aqp+resp%gamma)
         do ipg=1,3 ! XXX need 0 for alphaxy ????
            resp%ctmp=wpsipghp(resp%zarg,ipg)
            resp%RePolyGamma(ipg)=real(resp%ctmp,16)
            resp%ImPolyGamma(ipg)=imag(resp%ctmp)
         enddo
         
         ! compute transport kernels (omega-part)
         ! 
         resp%tmp=resp%z**2 / (4.q0*piQ**3) ! missing: beta/gamma (multiplied later to keep number reasonable here)        
         resp%s_ker = resp%tmp * ( resp%RePolyGamma(1) - resp%gamma*beta2pQ * resp%RePolyGamma(2) )
         resp%a_ker = resp%tmp * ( resp%aqp * resp%RePolyGamma(1) - resp%aqp*resp%gamma*beta2pQ*resp%RePolyGamma(2) &
                    - resp%gamma**2.q0 * beta2pQ * resp%ImPolyGamma(2) )
         
         
         resp%tmp=resp%tmp*3.q0*resp%z/(4.q0*piQ) ! additionally missing: 1/gamma (multiplied later) XXX

         if(algo%lBfield) then        
            resp%sB_ker=resp%tmp*(-resp%RePolyGamma(1)-resp%gamma*beta2pQ*resp%RePolyGamma(2)-(beta2pQ*resp%gamma)**2/3.q0 &
                       * resp%RePolyGamma(3))
            
            
            resp%aB_ker=resp%tmp*(resp%aqp*resp%RePolyGamma(1)-resp%aqp*beta2pQ*resp%gamma*resp%RePolyGamma(2)+resp%gamma**2/3.q0 &
                 * beta2pQ * resp%ImPolyGamma(2) &
                 - resp%aqp*resp%gamma**2.q0 / 3.q0 * beta2pQ**2.q0 * resp%ImPolyGamma(3) &
                 + resp%gamma**3.q0 / 3.q0 * beta2pQ**2.q0 * resp%RePolyGamma(3) )
          endif 
         
         
         ! B = 0  
         !tmp=vka(ik,ib,ialpha)*vka(ik,ib,ibeta)
         do ialpha=1,nalpha
         
            if (algo%ltbind) then
               resp%tmp=real(ek%Mopt(ialpha,thdr%idtet(ik,itet),iband,iband)*ek%Mopt(ialpha,thdr%idtet(ik,itet),iband,iband),16)
            else
               resp%tmp=real(ek%Mopt(ialpha,thdr%idtet(ik,itet), iband, iband),16) !the optical matrix elements given by Wien2k are squared already
            endif
         
            s_tmp_tetra(ik,iband,ialpha,ialpha)=resp%s_ker * resp%tmp  
            a_tmp_tetra(ik,iband,ialpha,ialpha)=resp%a_ker * resp%tmp 
               
            do ibeta=ialpha+1,nalpha
               resp%tmp=real(ek%Mopt(ialpha+ibeta+1,thdr%idtet(ik,itet), iband, iband),16) 
               s_tmp_tetra(ik,iband,ialpha,ibeta)=resp%s_ker * resp%tmp  
               a_tmp_tetra(ik,iband,ialpha,ibeta)=resp%a_ker * resp%tmp 
            enddo !ibeta  
         enddo ! ialpha           

      
         ! B .ne. 0 
         if (algo%lBfield .and. algo%ltbind ) then
            do ialpha=1,nalpha
               do ibeta=ialpha+1,nalpha
                  
                  !tmp=vka(ik,ib,ialpha)*( vkab(ik,ib,ibeta,ialpha)*vka(ik,ib,ibeta) - vkab(ik,ib,ibeta,ibeta)*vka(ik,ib,ialpha)    )
                  resp%tmp = real(ek%Mopt(ialpha,thdr%idtet(ik,itet), iband, iband)* &
                    (ek%M2(ibeta, ialpha, thdr%idtet(ik,itet), iband)*ek%Mopt(ialpha,thdr%idtet(ik,itet), iband, iband) - &
                    ek%M2(ibeta, ibeta, thdr%idtet(ik,itet), iband)* ek%Mopt(ialpha,thdr%idtet(ik,itet), iband, iband) ),16) 
                  sB_tmp_tetra(ik,iband,ialpha,ibeta)=resp%sB_ker * resp%tmp 
                  aB_tmp_tetra(ik,iband,ialpha,ibeta)=resp%aB_ker * resp%tmp 
                  
               enddo !ibeta  
            enddo ! ialpha
         endif !lBfield

         ! Now copy the local variable into the datastructure that will be passed to the interptra_mu
         resp%s_tmp(ik,iband,:,:) = s_tmp_tetra(ik,iband,:,:)
         resp%a_tmp(ik,iband,:,:) = a_tmp_tetra(ik,iband,:,:)
         if(algo%lBfield .and. algo%ltbind ) then
            resp%sB_tmp(ik,iband,:,:) = sB_tmp_tetra(ik,iband,:,:)
            resp%aB_tmp(ik,iband,:,:) = aB_tmp_tetra(ik,iband,:,:)
         endif
         
      enddo ! iband
   enddo ! ik   

end subroutine respintet_qp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! RESPINTET_BL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This routine (and the following one in QP)
! evaluate the conduction kernel, and the 
! full response functions on the vertices of a
! tetrahedron in the BOLTZMANN formalism
!
subroutine respintet_Bl(mu, iT, itet, nalpha, thdr, algo, ek, sct, resp)
  use response
  use types
  use params
  implicit none 
  type (dp_resp) :: resp ! dynamical datatype allow only for an inclusion through extension of the parenttype
                         ! I could have declared a dummy class pointer that would be assigned to either dp or qp response type
                         ! to do so the varaibles defined in the individual types must have had different names and since I
                         ! REALLY dislike having that extra Q for each varible I decided to stick to static datatypes
  type (algorithm) :: algo
  type (tetramesh) :: thdr
  type (edisp) :: ek
  type (scatrate) :: sct
  real(8), intent(in) :: mu
  integer, intent(in) :: iT
  integer, intent(in) :: itet
  integer, intent(in) :: nalpha
  integer :: iband, ik, ipg
  integer :: ialpha,ibeta
!local variables
  real(8), allocatable :: s_tmp_tetra(:,:,:,:),  a_tmp_tetra(:,:,:,:)
  real(8), allocatable :: sB_tmp_tetra(:,:,:,:), aB_tmp_tetra(:,:,:,:)
!external variables
  real(8), external :: dfermi

  !allocation  
  if(.not. allocated(s_tmp_tetra)) allocate(s_tmp_tetra(4,ek%nband_max,3,3))
  if(.not. allocated(a_tmp_tetra)) allocate(a_tmp_tetra(4,ek%nband_max,3,3))
  if (algo%lBfield .and. algo%ltbind) then
     if(.not. allocated(sB_tmp_tetra)) allocate(sB_tmp_tetra(4,ek%nband_max,3,3))
     if(.not. allocated(aB_tmp_tetra)) allocate(aB_tmp_tetra(4,ek%nband_max,3,3))
  endif
  !initialisation
  s_tmp_tetra=0.0d0 ; a_tmp_tetra=0.0d0
  if (algo%lBfield .and. algo%ltbind) then
     sB_tmp_tetra=0.0d0 ; aB_tmp_tetra=0.0d0
  endif


   do ik=1,4  !loop over corners of the tetrahedron

      do iband=1,ek%nband_max !loop over bands (these will be traced over)

         if (iband < ek%nbopt_min) cycle
         if (iband > ek%nbopt_max) cycle
        
         ! if the band is not contained in the optical matrices just do nothing
         resp%z=real(sct%z,8)
         resp%gamma=resp%z*real(sct%gam(iT,iband),8)
         resp%aqp=resp%z*real(ek%band(thdr%idtet(ik,itet),iband)-mu,8)
         
         ! compute transport kernels (omega-part)
         ! 
         resp%tmp=resp%z**2 / (2.d0*pi) ! missing: 1/gamma (multiplied later to keep number reasonable here)        
         resp%s_ker = resp%tmp * dfermi(resp%aqp, beta) 
         resp%a_ker = resp%tmp * resp%aqp * dfermi(resp%aqp, beta)
         
         
         resp%tmp=resp%tmp*3.d0*resp%z/(4.d0*pi) ! additionally missing: 1/gamma (multiplied later) XXX

         if(algo%lBfield .and. algo%ltbind) then        
            resp%sB_ker = resp%tmp * ( -dfermi(resp%aqp, beta)   )
            resp%aB_ker = resp%tmp * ( resp%aqp * dfermi(resp%aqp, beta) )
          endif 
         
         
         ! B = 0  
         !tmp=vka(ik,ib,ialpha)*vka(ik,ib,ibeta)
         do ialpha=1,nalpha

            if (algo%ltbind) then
               resp%tmp=ek%Mopt(ialpha,thdr%idtet(ik,itet), iband, iband)*ek%Mopt(ialpha,thdr%idtet(ik,itet), iband, iband)
            else
               resp%tmp=ek%Mopt(ialpha,thdr%idtet(ik,itet), iband, iband) !the optical matrix elements given by Wien2k are squared already
            endif
            s_tmp_tetra(ik,iband,ialpha,ialpha)=resp%s_ker * resp%tmp  
            a_tmp_tetra(ik,iband,ialpha,ialpha)=resp%a_ker * resp%tmp 
               
            do ibeta=ialpha+1,nalpha
               resp%tmp=ek%Mopt(ialpha+ibeta+1,thdr%idtet(ik,itet), iband, iband) !the optical matrix elements given by Wien2k are squared already
               s_tmp_tetra(ik,iband,ialpha,ibeta)=resp%s_ker * resp%tmp  
               a_tmp_tetra(ik,iband,ialpha,ibeta)=resp%a_ker * resp%tmp 
            enddo !ibeta  
         enddo ! ialpha           
      
         ! B .ne. 0 
         if (algo%lBfield .and. algo%ltbind ) then
            do ialpha=1,nalpha
               do ibeta=ialpha+1,3
                  
                  !tmp=vka(ik,ib,ialpha)*( vkab(ik,ib,ibeta,ialpha)*vka(ik,ib,ibeta) - vkab(ik,ib,ibeta,ibeta)*vka(ik,ib,ialpha)    )
                  resp%tmp =ek%Mopt(ialpha,thdr%idtet(ik,itet), iband, iband)* &
                    (ek%M2(ibeta, ialpha, thdr%idtet(ik,itet), iband)*ek%Mopt(ialpha,thdr%idtet(ik,itet), iband, iband) - &
                    ek%M2(ibeta, ibeta, thdr%idtet(ik,itet), iband)* ek%Mopt(ialpha,thdr%idtet(ik,itet), iband, iband) )

                  sB_tmp_tetra(ik,iband,ialpha,ibeta)=resp%sB_ker * resp%tmp 
                  aB_tmp_tetra(ik,iband,ialpha,ibeta)=resp%aB_ker * resp%tmp 
                  
               enddo !ibeta  
            enddo ! ialpha
         endif !lBfield

         ! Now copy the local variable into the datastructure that will be passed to the interptra_mu
         resp%s_tmp(ik,iband,:,:) = s_tmp_tetra(ik,iband,:,:)
         resp%a_tmp(ik,iband,:,:) = a_tmp_tetra(ik,iband,:,:)
         if(algo%lBfield .and. algo%ltbind ) then
            resp%sB_tmp(ik,iband,:,:) = sB_tmp_tetra(ik,iband,:,:)
            resp%aB_tmp(ik,iband,:,:) = aB_tmp_tetra(ik,iband,:,:)
         endif
         
      enddo ! iband
   enddo ! ik   

end subroutine respintet_Bl

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! RESPINTERT
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This routine (QP counterpart missing)
! evaluates the conduction kernel, and the 
! full response functions on the vertices of a
! tetrahedron in the KUBO formalism
! for interband transitions
! singularities might arise at conical intersections 
! between the two bands
!
subroutine respintert(mu, iT, itet, nalpha, thdr, algo, ek, sct, resp)
  use response
  use types
  use params
  implicit none 
  type (dp_respinter) :: resp 
  type (algorithm) :: algo
  type (tetramesh) :: thdr
  type (edisp) :: ek
  type (scatrate) :: sct
  real(8), intent(in) :: mu
  integer, intent(in) :: iT
  integer, intent(in) :: itet
  integer, intent(in) :: nalpha
  integer :: ib1, ib2, ik, ipg !band1, band2, k-point,  
  integer :: ialpha,ibeta 
  complex(8),external  :: wpsipg
  complex(16),external :: wpsipghp
!local variables
  real(8), allocatable :: s_tmp_tetra(:,:,:,:),  a_tmp_tetra(:,:,:,:)
  real(8) :: Dqp, Dgamma, Ggamma !qp energy difference, scattering rate difference and sum  
  real(8) :: DD1, DD2   !denominators
  real(8) :: ReK, ImK, tmp_s, tmp_a 

  !allocation  
  if(.not. allocated(s_tmp_tetra)) allocate(s_tmp_tetra(4,ek%nband_max,3,3))
  if(.not. allocated(a_tmp_tetra)) allocate(a_tmp_tetra(4,ek%nband_max,3,3))
  !initialisation
  s_tmp_tetra=0.0d0 ; a_tmp_tetra=0.0d0

   do ik=1,4  !loop over corners of the tetrahedron

      do ib1=1,ek%nband_max !loop over bands (these will be traced over)

         ! if the band is not contained in the optical matrices just do nothing
         if (ib1 < ek%nbopt_min) cycle
         if (ib1 > ek%nbopt_max) cycle
        
         resp%z1=real(sct%z,8)
         resp%gamma1=resp%z1*real(sct%gam(iT,ib1),8)
         resp%aqp1=resp%z1*real(ek%band(thdr%idtet(ik,itet),ib1)-mu,8)
         !the first state has to belong to the occupied manifold
         if (resp%aqp1 > 0.0d0) cycle
         resp%zarg=0.5d0+beta2p*(ci*resp%aqp1+resp%gamma1)
         do ipg=1,1 
            resp%ctmp=wpsipg(resp%zarg,ipg)
            resp%RePolyGamma1(ipg)=real(resp%ctmp,8)
            resp%ImPolyGamma1(ipg)=imag(resp%ctmp)
         enddo
         
         ! compute transport kernels (omega-part)
         do ib2=1,ek%nband_max
            if (ib2 < ek%nbopt_min) cycle
            if (ib2 > ek%nbopt_max) cycle
            if (ib2 == ib1 ) cycle
            !singularities might arise if ek%band1 = ek%band2

            !second band variables and derived quantities
            resp%z2=real(sct%z,8)
            resp%gamma2=resp%z2*real(sct%gam(iT,ib2),8)
            resp%aqp2=resp%z2*real(ek%band(thdr%idtet(ik,itet),ib2)-mu,8)
            !the second state has to belong to the unoccupied manifold
            if (resp%aqp2 < 0.0d0) cycle
            resp%zarg=0.5d0+beta2p*(ci*resp%aqp2+resp%gamma2)
            do ipg=1,1 
               resp%ctmp=wpsipg(resp%zarg,ipg)
               resp%RePolyGamma2(ipg)=real(resp%ctmp,8)
               resp%ImPolyGamma2(ipg)=imag(resp%ctmp)
            enddo

            Dqp    = resp%aqp1 - resp%aqp2     !Delta csi in eq
            Dgamma = resp%gamma1 - resp%gamma2 !Delta in eq
            Ggamma = resp%gamma1 + resp%gamma2 !Gamma in eq
            DD1 = 1.0d0/(Dqp**2 + Ggamma**2)
            DD2 = 1.0d0/(Dqp**2 + Dgamma**2)

            ReK = 2.0d0*resp%gamma1*resp%gamma2*DD2*( (resp%gamma2*resp%aqp1) - (resp%gamma1*resp%aqp2) )
            ImK = resp%gamma1*resp%gamma2*DD2*( (Ggamma*Dgamma) + (resp%aqp1+resp%aqp2)*Dqp )

            tmp_s = (resp%z1*resp%z2 * resp%gamma1*resp%gamma2)*DD1*beta/(pi**3)
            tmp_a = 0.5d0*resp%z1*resp%z2*DD1*(beta**2)/(pi**3)
              

            resp%s_ker = tmp_s * ( ((DD2*Dgamma + 0.5d0/resp%gamma2)*resp%RePolyGamma2(1)) &
                       - ((DD2*Dgamma - 0.5d0/resp%gamma1)*resp%RePolyGamma1(1)) &
                       + (DD2*Dqp*(resp%ImPolyGamma2(1) - resp%ImPolyGamma1(1))) )

            resp%a_ker = tmp_a * ( (resp%aqp1*resp%gamma2*resp%RePolyGamma1(1)) + (resp%aqp2*resp%gamma1*resp%RePolyGamma2(1)) &
                       + (ReK*(resp%RePolyGamma1(1) - resp%RePolyGamma2(1))) + (ImK*(resp%ImPolyGamma2(1) - resp%ImPolyGamma1(1))) )
                    
            ! B = 0  
            ! tmp=vka(ik,ib,ialpha)*vka(ik,ib,ibeta)
            do ialpha=1,nalpha
               if (algo%ltbind) then
                  resp%tmp=ek%Mopt(ialpha,thdr%idtet(ik,itet), ib1, ib2)*ek%Mopt(ialpha,thdr%idtet(ik,itet), ib1, ib2)
               else
                  resp%tmp=ek%Mopt(ialpha,thdr%idtet(ik,itet), ib1, ib2) !the optical matrix elements given by Wien2k are squared already
               endif
            
               s_tmp_tetra(ik,ib1,ialpha,ialpha)=s_tmp_tetra(ik,ib1,ialpha,ialpha) + (resp%s_ker * resp%tmp)  
               a_tmp_tetra(ik,ib1,ialpha,ialpha)=a_tmp_tetra(ik,ib1,ialpha,ialpha) + (resp%a_ker * resp%tmp) 
               
               do ibeta=ialpha+1,nalpha
                  resp%tmp=ek%Mopt(ialpha+ibeta+1,thdr%idtet(ik,itet), ib1, ib2) 
                  s_tmp_tetra(ik,ib1,ialpha,ibeta)=s_tmp_tetra(ik,ib1,ialpha,ibeta) + (resp%s_ker * resp%tmp)  
                  a_tmp_tetra(ik,ib1,ialpha,ibeta)=a_tmp_tetra(ik,ib1,ialpha,ibeta) + (resp%a_ker * resp%tmp) 
               enddo !ibeta  
            enddo ! ialpha           
         enddo !ib2

         ! Now copy the local variable into the datastructure that will be passed to the interptra_mu
         resp%s_tmp(ik,ib1,:,:) = s_tmp_tetra(ik,ib1,:,:)
         resp%a_tmp(ik,ib1,:,:) = a_tmp_tetra(ik,ib1,:,:)
         
      enddo ! ib1
   enddo ! ik   

end subroutine respintert

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! RESPINTERT_SYMM
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This routine (QP counterpart missing)
! evaluates the conduction kernel, and the 
! full response functions on the vertices of a
! tetrahedron in the KUBO formalism
! for interband transitions.
! It uses a semplified kernel obtained for a 
! 2 band symmetrical semiconductor with 
! a given band gap !
! Expressions valid only at the Gamma point 
! have been commented out
!
subroutine respintert_symm(mu, iT, itet, nalpha, thdr, algo, ek, sct, resp)
  use response
  use types
  use params
  implicit none 
  type (dp_respinter) :: resp 
  type (algorithm) :: algo
  type (tetramesh) :: thdr
  type (edisp) :: ek
  type (scatrate) :: sct
  real(8), intent(in) :: mu
  !real(8), intent(in) :: gap
  integer, intent(in) :: iT
  integer, intent(in) :: itet
  integer, intent(in) :: nalpha
  integer :: ib1, ib2, ik, ipg !band1, band2, k-point,  
  integer :: ialpha,ibeta 
  complex(8),external  :: wpsipg
  complex(16),external :: wpsipghp
!local variables
  real(8), allocatable :: s_tmp_tetra(:,:,:,:),  a_tmp_tetra(:,:,:,:)
  real(8) :: Dqp, Dgamma, Ggamma !qp energy difference, scattering rate difference and sum  
  real(8) :: DD1, DD2   !denominators
  real(8) :: ReK, ImK, tmp_s, tmp_a 

  !allocation  
  if(.not. allocated(s_tmp_tetra)) allocate(s_tmp_tetra(4,ek%nband_max,3,3))
  if(.not. allocated(a_tmp_tetra)) allocate(a_tmp_tetra(4,ek%nband_max,3,3))
  !initialisation
  s_tmp_tetra=0.0d0 ; a_tmp_tetra=0.0d0

   do ik=1,4  !loop over corners of the tetrahedron

      do ib1=1,ek%nband_max !loop over bands (these will be traced over)

         ! if the band is not contained in the optical matrices just do nothing
         if (ib1 < ek%nbopt_min) cycle
         if (ib1 > ek%nbopt_max) cycle
         resp%z1=real(sct%z,8)
         resp%gamma1=resp%z1*real(sct%gam(iT,ib1),8)
         resp%aqp1=resp%z1*real(ek%band(thdr%idtet(ik,itet),ib1),8) !in a symmetric SC mu=0
         ! if the band is unoccupied cycle
         if(resp%aqp1 > mu) cycle
         resp%zarg=0.5d0+beta2p*((ci*resp%aqp1)+resp%gamma1)
         do ipg=1,1 
            resp%ctmp=wpsipg(resp%zarg,ipg)
            resp%RePolyGamma1(ipg)=real(resp%ctmp,8)
            resp%ImPolyGamma1(ipg)=imag(resp%ctmp)
         enddo
         
         ! compute transport kernels (omega-part)
         ! 
         do ib2=1,ek%nband_max
            if (ib2 < ek%nbopt_min) cycle
            if (ib2 > ek%nbopt_max) cycle
            if (ib2 == ib1 ) cycle

            !second band variables and derived quantities
            resp%z2=real(sct%z,8)
            resp%gamma2=resp%z2*resp%gamma1   !only one gamma required !real(sct%gam(iT,ib2),8)
            resp%aqp2=resp%z2*real(ek%band(thdr%idtet(ik,itet),ib2),8) !in a symmetric SC mu=0
            ! if the second state is occupied cycle (interband contribution)
            if(resp%aqp2 < mu) cycle
            resp%zarg=0.5d0+beta2p*(ci*resp%aqp2+resp%gamma2)
            do ipg=1,1 
               resp%ctmp=wpsipg(resp%zarg,ipg)
               resp%RePolyGamma2(ipg)=real(resp%ctmp,8)
               resp%ImPolyGamma2(ipg)=imag(resp%ctmp)
            enddo

            Dqp    = resp%aqp1 - resp%aqp2     !Delta csi in eq
            !DD1 = 1.0d0/(gap**2 + 4.0d0*(resp%gamma1**2) )
            DD1 = 1.0d0/(Dqp**2 + 4.0d0*(resp%gamma1**2) )

            tmp_s = DD1*((resp%z1 * resp%gamma1)**2)*beta/(pi**3)
            tmp_a = DD1*((resp%z1 * beta)**2)/(2.0d0*(pi**3))
              

            resp%s_ker = tmp_s * ( (resp%RePolyGamma2(1) + resp%RePolyGamma1(1))/(2.0d0*resp%gamma1) &
                       + (resp%ImPolyGamma2(1) - resp%ImPolyGamma1(1))/Dqp )
                       !- (resp%ImPolyGamma2(1) - resp%ImPolyGamma1(1))/gap ) !only at the Gamma point!!
            resp%a_ker = tmp_a * ( resp%gamma1*(resp%aqp1*resp%RePolyGamma1(1) + resp%aqp2*resp%RePolyGamma2(1)) &
                       + (resp%gamma1**2)*(resp%aqp1+resp%aqp2)*(resp%ImPolyGamma2(1)-resp%ImPolyGamma1(1))/Dqp  &
                       + (resp%gamma1**3)*2.0d0*(resp%RePolyGamma1(1) - resp%RePolyGamma2(1))/Dqp )

            !only at the Gamma point!!           
            !resp%a_ker = tmp_a * ( resp%gamma1*abs(resp%aqp1)*(resp%RePolyGamma2(1)-resp%RePolyGamma1(1)) &
            !           + abs(resp%aqp1)*(resp%gamma1**3)*(resp%RePolyGamma2(1)-resp%RePolyGamma1(1))/(resp%aqp1**2) ) 
                       
            ! B = 0  
            !tmp=vka(ik,ib,ialpha)*vka(ik,ib,ibeta)
            do ialpha=1,nalpha
               if (algo%ltbind) then
                  resp%tmp=ek%Mopt(ialpha,thdr%idtet(ik,itet), ib1, ib2)*ek%Mopt(ialpha,thdr%idtet(ik,itet), ib1, ib2)
               else
                  resp%tmp=ek%Mopt(ialpha,thdr%idtet(ik,itet), ib1, ib2) !the optical matrix elements given by Wien2k are squared already
               endif
            
               s_tmp_tetra(ik,ib1,ialpha,ialpha)=s_tmp_tetra(ik,ib1,ialpha,ialpha) + (resp%s_ker * resp%tmp)  
               a_tmp_tetra(ik,ib1,ialpha,ialpha)=a_tmp_tetra(ik,ib1,ialpha,ialpha) + (resp%a_ker * resp%tmp) 
                  
               do ibeta=ialpha+1,nalpha
                  resp%tmp=ek%Mopt(ialpha+ibeta+1,thdr%idtet(ik,itet), ib1, ib2) 
                  s_tmp_tetra(ik,ib1,ialpha,ibeta)=s_tmp_tetra(ik,ib1,ialpha,ibeta) + (resp%s_ker * resp%tmp)  
                  a_tmp_tetra(ik,ib1,ialpha,ibeta)=a_tmp_tetra(ik,ib1,ialpha,ibeta) + (resp%a_ker * resp%tmp) 
               enddo !ibeta  
            enddo ! ialpha           
         enddo !ib2

         ! Now copy the local variable into the datastructure that will be passed to the interptra_mu
         resp%s_tmp(ik,ib1,:,:) = s_tmp_tetra(ik,ib1,:,:)
         resp%a_tmp(ik,ib1,:,:) = a_tmp_tetra(ik,ib1,:,:)
         
      enddo ! ib1
   enddo ! ik   

end subroutine respintert_symm

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! RESPINKM
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This routine evaluates the conduction kernel
! for a given k-point in the KUBO formalism
! for intraband transitions.
!
subroutine respinkm(mu, iT, ik, nalpha, algo, ek, sct, resp)
  use response
  use types
  use params
  implicit none 
  type (dp_resp) :: resp ! dynamical datatype allow only for an inclusion through extension of the parenttype
                         ! I could have declared a dummy class pointer that would be assigned to either dp or qp response type
                         ! to do so the varaibles defined in the individual types must have had different names and since I
                         ! REALLY dislike having that extra Q for each varible I decided to stick to static datatypes
  type (algorithm) :: algo
  type (edisp) :: ek
  type (scatrate) :: sct
  real(8), intent(in) :: mu
  integer, intent(in) :: iT
  integer, intent(in) :: ik
  integer, intent(in) :: nalpha
  integer :: iband, ipg
  integer :: ialpha,ibeta
  complex(8),external  :: wpsipg
  complex(16),external :: wpsipghp

   !loop over k-points is external
   do iband=1,ek%nband_max !loop over bands (these will be traced over)

     ! if the band is not contained in the optical matrices just do nothing
     if (iband < ek%nbopt_min) cycle
     if (iband > ek%nbopt_max) cycle
     resp%z=real(sct%z,8)
     resp%gamma=resp%z*real(sct%gam(iT,iband),8)
     ! pre-compute all needed digamma functions   
     resp%aqp=resp%z*real(ek%band(ik,iband)-mu,8)
     ! TESTED 25.05.2018 tetragonal
     !resp%aqp=real(ek%band(1,iband)-mu,8) !passed loptic=false
     !resp%aqp=real(ek%band(8,iband)-mu,8) !passed loptic=false
     ! TESTED 28.05.2018 tetragonal
     !resp%aqp=real(ek%band(2,iband)-mu,8) !not passed with loptic=true
     !TEST END
     resp%zarg=0.5d0+beta2p*(ci*resp%aqp+resp%gamma)
     do ipg=1,3 ! XXX need 0 for alphaxy ????
        resp%ctmp=wpsipg(resp%zarg,ipg)
        resp%RePolyGamma(ipg)=real(resp%ctmp,8)
        resp%ImPolyGamma(ipg)=imag(resp%ctmp)
     enddo
     
     ! compute transport kernels (omega-part)
     ! 
     resp%tmp=resp%z**2 / (4.d0*pi**3) ! missing: beta/gamma (multiplied later to keep number reasonable here)        
     resp%s_ker = resp%tmp * (resp%RePolyGamma(1) - resp%gamma*beta2p * resp%RePolyGamma(2) )
     resp%a_ker = resp%tmp * (resp%aqp * resp%RePolyGamma(1) - resp%aqp*resp%gamma*beta2p*resp%RePolyGamma(2) &
                - resp%gamma**2.d0 * beta2p * resp%ImPolyGamma(2) )
     
     
     resp%tmp=resp%tmp*3.d0*resp%z/(4.d0*pi) ! additionally missing: 1/gamma (multiplied later) XXX

     if(algo%lBfield .and. algo%ltbind) then        
        resp%sB_ker=resp%tmp*(-resp%RePolyGamma(1)-resp%gamma*beta2p*resp%RePolyGamma(2)-(beta2p*resp%gamma)**2.d0/3.d0 &
                   * resp%RePolyGamma(3)) 
        
        resp%aB_ker=resp%tmp*(resp%aqp*resp%RePolyGamma(1)-resp%aqp*beta2p*resp%gamma*resp%RePolyGamma(2)+resp%gamma**2.d0/3.d0  &
             * beta2p * resp%ImPolyGamma(2) &
             - resp%aqp*resp%gamma**2.d0 / 3.d0 * beta2p**2.d0 * resp%ImPolyGamma(3) &
             + resp%gamma**3.d0 / 3.d0 * beta2p**2.d0 * resp%RePolyGamma(3) )
      endif 
     
     
     ! B = 0  
     !tmp=vka(ik,ib,ialpha)*vka(ik,ib,ibeta)
     do ialpha=1,nalpha
        if (algo%ltbind) then
           resp%tmp=ek%Mopt(ialpha,ik, iband, iband)*ek%Mopt(ialpha,ik, iband, iband)
        else
           !the expression requires only the diagonal of the optical matrix elements because a trace is evaluated 
           resp%tmp=ek%Mopt(ialpha,ik, iband, iband) !the optical matrix elements given by Wien2k are squared already
        endif
     
        resp%s_tmp(ik,iband,ialpha,ialpha)=resp%s_ker * resp%tmp  
        resp%a_tmp(ik,iband,ialpha,ialpha)=resp%a_ker * resp%tmp 
           
     ! TESTED 25.05.2018
        do ibeta=ialpha+1,nalpha
           resp%tmp=ek%Mopt(ialpha+ibeta+1,ik, iband, iband) 
           resp%s_tmp(ik,iband,ialpha,ibeta)=resp%s_ker * resp%tmp  
           resp%a_tmp(ik,iband,ialpha,ibeta)=resp%a_ker * resp%tmp 
        enddo !ibeta  
     enddo ! ialpha           

  
     ! B .ne. 0 
     if (algo%lBfield .and. algo%ltbind ) then
        do ialpha=1,nalpha
           do ibeta=ialpha+1,3
              
              !tmp=vka(ik,ib,ialpha)*( vkab(ik,ib,ibeta,ialpha)*vka(ik,ib,ibeta) - vkab(ik,ib,ibeta,ibeta)*vka(ik,ib,ialpha)    )
              resp%tmp =ek%Mopt(ialpha, ik, iband, iband)*( ek%M2(ibeta, ialpha, ik, iband)*ek%Mopt(ialpha, ik, iband, iband) - &
                 ek%M2(ibeta, ibeta, ik, iband)* ek%Mopt(ialpha, ik, iband, iband) )
                  resp%sB_tmp(ik,iband,ialpha,ibeta)=resp%sB_ker * resp%tmp 
                  resp%aB_tmp(ik,iband,ialpha,ibeta)=resp%aB_ker * resp%tmp 
                 
           enddo !ibeta  
        enddo ! ialpha
     endif !lBfield
     
   enddo ! iband

end subroutine respinkm

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! RESPINKM_QP
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This routine evaluates the conduction kernel
! for a given k-point  in the KUBO formalism
! for intraband transitions. 
! Quadruple precision conterpart of the 
! preceeding routine
!
subroutine respinkm_qp(mu, iT, ik, nalpha, algo, ek, sct, resp)
  use response
  use types
  use params
  implicit none 
  type (qp_resp) :: resp ! dynamical datatype allow only for an inclusion through extension of the parenttype
                         ! I could have declared a dummy class pointer that would be assigned to either dp or qp response type
                         ! to do so the varaibles defined in the individual types must have had different names and since I
                         ! REALLY dislike having that extra Q for each varible I decided to stick to static datatypes
  type (algorithm) :: algo
  type (edisp) :: ek
  type (scatrate) :: sct
  real(8), intent(in) :: mu
  integer, intent(in) :: iT
  integer, intent(in) :: ik
  integer, intent(in) :: nalpha
  integer :: iband, ipg
  integer :: ialpha,ibeta
  complex(8),external  :: wpsipg
  complex(16),external :: wpsipghp

  betaQ=real(beta,16)
  beta2pQ=betaQ/(2.q0*piQ)


   do iband=1,ek%nband_max !loop over bands (these will be traced over)

      if (iband < ek%nbopt_min) cycle
      if (iband > ek%nbopt_max) cycle
      ! if the band is not contained in the optical matrices just do nothing
      resp%z=real(sct%z,16)
      resp%gamma=resp%z*real(sct%gam(iT,iband),16)
      ! pre-compute all needed digamma functions   
      resp%aqp=resp%z*real(ek%band(ik,iband)-mu,16)
      resp%zarg=0.5q0+beta2pQ*(ciQ*resp%aqp+resp%gamma)
      do ipg=1,3 ! XXX need 0 for alphaxy ????
         resp%ctmp=wpsipghp(resp%zarg,ipg)
         resp%RePolyGamma(ipg)=real(resp%ctmp,16)
         resp%ImPolyGamma(ipg)=imag(resp%ctmp)
      enddo
      
      ! compute transport kernels (omega-part)
      ! 
      resp%tmp=resp%z**2 / (4.q0*piQ**3) ! missing: beta/gamma (multiplied later to keep number reasonable here)        
      resp%s_ker=resp%tmp*(resp%RePolyGamma(1)-resp%gamma*beta2pQ*resp%RePolyGamma(2) )
      resp%a_ker=resp%tmp*(resp%aqp*resp%RePolyGamma(1)-resp%aqp*resp%gamma*beta2pQ*resp%RePolyGamma(2) &
                -resp%gamma**2*beta2pQ*resp%ImPolyGamma(2))
      
      
      resp%tmp=resp%tmp*3.q0*resp%z/(4.q0*piQ) ! additionally missing: 1/gamma (multiplied later) XXX

      if(algo%lBfield) then        
         resp%sB_ker=resp%tmp*(-resp%RePolyGamma(1)-resp%gamma*beta2pQ*resp%RePolyGamma(2)-(beta2pQ*resp%gamma)**2/3.q0 &
                    *resp%RePolyGamma(3))
         
         
         resp%aB_ker=resp%tmp*(resp%aqp*resp%RePolyGamma(1)-resp%aqp*beta2pQ*resp%gamma*resp%RePolyGamma(2)+resp%gamma**2.q0/3.q0  &
              * beta2pQ * resp%ImPolyGamma(2) &
              - resp%aqp*resp%gamma**2.q0 / 3.q0 * beta2pQ**2.q0 * resp%ImPolyGamma(3) &
              + resp%gamma**3.q0 / 3.q0 * beta2pQ**2.q0 * resp%RePolyGamma(3) )
       endif 
      
      
      ! B = 0  
      !tmp=vka(ik,ib,ialpha)*vka(ik,ib,ibeta)
      do ialpha=1,nalpha
      
         if (algo%ltbind) then
            resp%tmp=real(ek%Mopt(ialpha,ik, iband, iband)*ek%Mopt(ialpha,ik, iband, iband),16)
         else
            resp%tmp=real(ek%Mopt(ialpha,ik, iband, iband),16) !the optical matrix elements given by Wien2k are squared already
         endif
      
         resp%s_tmp(ik,iband,ialpha,ialpha)=resp%s_ker * resp%tmp  
         resp%a_tmp(ik,iband,ialpha,ialpha)=resp%a_ker * resp%tmp 
            
         do ibeta=ialpha+1,nalpha
            resp%tmp=real(ek%Mopt(ialpha+ibeta+1,ik, iband, iband),16) 
            resp%s_tmp(ik,iband,ialpha,ibeta)=resp%s_ker * resp%tmp  
            resp%a_tmp(ik,iband,ialpha,ibeta)=resp%a_ker * resp%tmp 
         enddo !ibeta  
      enddo ! ialpha           

   
      ! B .ne. 0 
      if (algo%lBfield .and. algo%ltbind ) then
         do ialpha=1,nalpha
            do ibeta=ialpha+1,3
               
               !tmp=vka(ik,ib,ialpha)*( vkab(ik,ib,ibeta,ialpha)*vka(ik,ib,ibeta) - vkab(ik,ib,ibeta,ibeta)*vka(ik,ib,ialpha)    )
               resp%tmp = real(ek%Mopt(ialpha,ik,iband,iband)*(ek%M2(ibeta,ialpha,ik,iband)*ek%Mopt(ialpha,ik,iband,iband) &
                        - ek%M2(ibeta, ibeta, ik, iband)* ek%Mopt(ialpha, ik, iband, iband) ),16) 
               resp%sB_tmp(ik,iband,ialpha,ibeta)=resp%sB_ker * resp%tmp 
               resp%aB_tmp(ik,iband,ialpha,ibeta)=resp%aB_ker * resp%tmp 
               
            enddo !ibeta  
         enddo ! ialpha
      endif !lBfield

   enddo ! iband

end subroutine respinkm_qp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! RESPINKM_BL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This routine evaluates the conduction kernel
! for a given k-point in the BOLTZMANN formalism
! for intraband transitions.
!
subroutine respinkm_Bl(mu, iT, ik, nalpha, algo, ek, sct, resp)
  use response
  use types
  use params
  implicit none 
  type (dp_resp) :: resp ! dynamical datatype allow only for an inclusion through extension of the parenttype
                         ! I could have declared a dummy class pointer that would be assigned to either dp or qp response type
                         ! to do so the varaibles defined in the individual types must have had different names and since I
                         ! REALLY dislike having that extra Q for each varible I decided to stick to static datatypes
  type (algorithm) :: algo
  type (edisp) :: ek
  type (scatrate) :: sct
  real(8), intent(in) :: mu
  integer, intent(in) :: iT
  integer, intent(in) :: ik
  integer, intent(in) :: nalpha
  integer :: iband, ipg
  integer :: ialpha,ibeta
!external variables
  real(8), external :: dfermi

   !loop over k-points is external
   do iband=1,ek%nband_max !loop over bands (these will be traced over)

     ! if the band is not contained in the optical matrices just do nothing
     if (iband < ek%nbopt_min) cycle
     if (iband > ek%nbopt_max) cycle
     resp%z=real(sct%z,8)
     resp%gamma=resp%z*real(sct%gam(iT,iband),8)
     ! pre-compute all needed digamma functions   
     resp%aqp=resp%z*real(ek%band(ik,iband)-mu,8)
     
     ! compute transport kernels (omega-part)
     ! 
     resp%tmp=resp%z**2 / (2.d0*pi) ! missing: 1/gamma (multiplied later to keep number reasonable here)        
     resp%s_ker = resp%tmp * dfermi(resp%aqp, beta) 
     resp%a_ker = resp%tmp * resp%aqp * dfermi(resp%aqp, beta)
     
     
     resp%tmp=resp%tmp*3.d0*resp%z/(4.d0*pi) ! additionally missing: 1/gamma (multiplied later) XXX

     if(algo%lBfield .and. algo%ltbind) then        
        resp%sB_ker = resp%tmp * ( -dfermi(resp%aqp, beta)  )
        resp%aB_ker = resp%tmp * ( resp%aqp * dfermi(resp%aqp, beta) )
      endif 
     
     
     ! B = 0  
     !tmp=vka(ik,ib,ialpha)*vka(ik,ib,ibeta)
     do ialpha=1,nalpha
        if (algo%ltbind) then
           resp%tmp=ek%Mopt(ialpha,ik, iband, iband)*ek%Mopt(ialpha,ik, iband, iband)
        else
           resp%tmp=ek%Mopt(ialpha,ik, iband, iband) !the optical matrix elements given by Wien2k are squared already
        endif
     
        resp%s_tmp(ik,iband,ialpha,ialpha)=resp%s_ker * resp%tmp  
        resp%a_tmp(ik,iband,ialpha,ialpha)=resp%a_ker * resp%tmp 
           
        do ibeta=ialpha+1,nalpha
           resp%tmp=ek%Mopt(ialpha+ibeta+1,ik, iband, iband) !the optical matrix elements given by Wien2k are squared already
           resp%s_tmp(ik,iband,ialpha,ibeta)=resp%s_ker * resp%tmp  
           resp%a_tmp(ik,iband,ialpha,ibeta)=resp%a_ker * resp%tmp 
        enddo !ibeta  
     enddo ! ialpha           

  
     ! B .ne. 0 
     if (algo%lBfield .and. algo%ltbind ) then
        do ialpha=1,nalpha
           do ibeta=ialpha+1,3
              
              !tmp=vka(ik,ib,ialpha)*( vkab(ik,ib,ibeta,ialpha)*vka(ik,ib,ibeta) - vkab(ik,ib,ibeta,ibeta)*vka(ik,ib,ialpha)    )
              resp%tmp =ek%Mopt(ialpha, ik, iband, iband)*( ek%M2(ibeta, ialpha, ik, iband)*ek%Mopt(ialpha, ik, iband, iband) - &
                 ek%M2(ibeta, ibeta, ik, iband)* ek%Mopt(ialpha, ik, iband, iband) )
                  resp%sB_tmp(ik,iband,ialpha,ibeta)=resp%sB_ker * resp%tmp 
                  resp%aB_tmp(ik,iband,ialpha,ibeta)=resp%aB_ker * resp%tmp 
                 
           enddo !ibeta  
        enddo ! ialpha
     endif !lBfield
     
   enddo ! iband

end subroutine respinkm_Bl

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! RESPINTERKM
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This routine (QP counterpart missing)
! evaluates the conduction kernel, and the 
! full response functions for a given k-point
! in the KUBO formalism for interband transitions
! singularities might arise at conical intersections 
! between the two bands
!
subroutine respinterkm(mu, iT, ik, nalpha, algo, ek, sct, resp)
  use response
  use types
  use params
  implicit none 
  type (dp_respinter) :: resp 
  type (algorithm) :: algo
  type (edisp) :: ek
  type (scatrate) :: sct
  real(8), intent(in) :: mu
  integer, intent(in) :: iT
  integer, intent(in) :: ik
  integer, intent(in) :: nalpha
  integer :: ib1, ib2, ipg !band1, band2, degree of Polygamma f'ns  
  integer :: ialpha,ibeta 
  complex(8),external  :: wpsipg
  complex(16),external :: wpsipghp
!local variables
  real(8) :: Dqp, Dgamma, Ggamma !qp energy difference, scattering rate difference and sum  
  real(8) :: DD1, DD2   !denominators
  real(8) :: ReK, ImK, tmp_s, tmp_a 

   do ib1=1,ek%nband_max !loop over bands (these will be traced over)
      ! if the band is not contained in the optical matrices just do nothing
      if (ib1 < ek%nbopt_min) cycle
      if (ib1 > ek%nbopt_max) cycle
      resp%z1=real(sct%z,8)
      resp%gamma1=resp%z1*real(sct%gam(iT,ib1),8)
      resp%aqp1=resp%z1*real(ek%band(ik,ib1)-mu,8)
      !the first state has to belong to the occupied manifold
      if (resp%aqp1 > 0.0d0) cycle
      resp%zarg=0.5d0+beta2p*(ci*resp%aqp1+resp%gamma1)
      do ipg=1,1 
         resp%ctmp=wpsipg(resp%zarg,ipg)
         resp%RePolyGamma1(ipg)=real(resp%ctmp,8)
         resp%ImPolyGamma1(ipg)=imag(resp%ctmp)
      enddo
      
      ! compute transport kernels (omega-part)
      ! 
      !do ib2=ib1+1,ek%nband_max
      do ib2=1,ek%nband_max
         if (ib2 < ek%nbopt_min) cycle
         if (ib2 > ek%nbopt_max) cycle
         if (ib2 == ib1 ) cycle
         !singularities might arise if ek%band1 = ek%band2

         !second band variables and derived quantities
         resp%z2=real(sct%z,8)
         resp%gamma2=resp%z2*real(sct%gam(iT,ib2),8)
         resp%aqp2=resp%z2*real(ek%band(ik,ib2)-mu,8)
         !the second state has to belong to the unoccupied manifold
         if (resp%aqp2 < 0.0d0) cycle
         resp%zarg=0.5d0+beta2p*(ci*resp%aqp2+resp%gamma2)
         do ipg=1,1 
            resp%ctmp=wpsipg(resp%zarg,ipg)
            resp%RePolyGamma2(ipg)=real(resp%ctmp,8)
            resp%ImPolyGamma2(ipg)=imag(resp%ctmp)
         enddo

         Dqp    = resp%aqp1 - resp%aqp2     !Delta csi in eq
         Dgamma = resp%gamma1 - resp%gamma2 !Delta in eq
         Ggamma = resp%gamma1 + resp%gamma2 !Gamma in eq
         DD1 = 1.0d0/(Dqp**2 + Ggamma**2)
         DD2 = 1.0d0/(Dqp**2 + Dgamma**2)

         ReK = 2.0d0*resp%gamma1*resp%gamma2*DD2*( (resp%gamma2*resp%aqp1) - (resp%gamma1*resp%aqp2) )
         ImK = resp%gamma1*resp%gamma2*DD2*( (Ggamma*Dgamma) + (resp%aqp1+resp%aqp2)*Dqp )

         tmp_s = (resp%z1*resp%z2 * resp%gamma1*resp%gamma2)*DD1*beta/(pi**3)
         tmp_a = 0.5d0*resp%z1*resp%z2*DD1*(beta**2)/(pi**3)
           

         resp%s_ker = tmp_s * ( ((DD2*Dgamma + 0.5d0/resp%gamma2)*resp%RePolyGamma2(1)) &
                    - ((DD2*Dgamma - 0.5d0/resp%gamma1)*resp%RePolyGamma1(1)) &
                    + (DD2*Dqp*(resp%ImPolyGamma2(1) - resp%ImPolyGamma1(1))) )

         resp%a_ker = tmp_a * ( (resp%aqp1*resp%gamma2*resp%RePolyGamma1(1)) + (resp%aqp2*resp%gamma1*resp%RePolyGamma2(1)) &
                    + (ReK*(resp%RePolyGamma1(1) - resp%RePolyGamma2(1))) + (ImK*(resp%ImPolyGamma2(1) - resp%ImPolyGamma1(1))) )
                 
         ! B = 0  
         !tmp=vka(ik,ib,ialpha)*vka(ik,ib,ibeta)
         do ialpha=1,nalpha
            if (algo%ltbind) then
               resp%tmp=ek%Mopt(ialpha,ik, ib1, ib2)*ek%Mopt(ialpha,ik, ib1, ib2)
            else
               resp%tmp=ek%Mopt(ialpha,ik, ib1, ib2) !the optical matrix elements given by Wien2k are squared already
            endif
         
            resp%s_tmp(ik,ib1,ialpha,ialpha)=resp%s_tmp(ik,ib1,ialpha,ialpha) + (resp%s_ker * resp%tmp)  
            resp%a_tmp(ik,ib1,ialpha,ialpha)=resp%a_tmp(ik,ib1,ialpha,ialpha) + (resp%a_ker * resp%tmp) 
               
            do ibeta=ialpha+1,nalpha
               resp%tmp=ek%Mopt(ialpha+ibeta+1,ik, ib1, ib2) !the optical matrix elements given by Wien2k are squared already
               resp%s_tmp(ik,ib1,ialpha,ibeta)=resp%s_tmp(ik,ib1,ialpha,ibeta) + (resp%s_ker * resp%tmp)  
               resp%a_tmp(ik,ib1,ialpha,ibeta)=resp%a_tmp(ik,ib1,ialpha,ibeta) + (resp%a_ker * resp%tmp) 
            enddo !ibeta  
         enddo ! ialpha          
 
      enddo !ib2
   enddo ! ib1

end subroutine respinterkm

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! RESPINTERKM_SYMM
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This routine (QP counterpart missing)
! evaluates the conduction kernel, and the 
! full response functions at a given k-point 
! in the KUBO formalism for interband transitions.
! It uses a semplified kernel obtained for a 
! 2 band symmetrical semiconductor with 
! a given band gap !
! Expressions valid only at the Gamma point 
! have been commented out
!
subroutine respinterkm_symm(mu, iT, ik, nalpha, algo, ek, sct, resp)
  use response
  use types
  use params
  implicit none 
  type (dp_respinter) :: resp 
  type (algorithm) :: algo
  type (edisp) :: ek
  type (scatrate) :: sct
  real(8), intent(in) :: mu
  !real(8), intent(in) :: gap
  integer, intent(in) :: iT
  integer, intent(in) :: ik
  integer, intent(in) :: nalpha
  integer :: ib1, ib2, ipg !band1, band2, degree of Polygamma f'ns  
  integer :: ialpha,ibeta 
  complex(8),external  :: wpsipg
  complex(16),external :: wpsipghp
!local variables
  real(8) :: Dqp, Dgamma, Ggamma !qp energy difference, scattering rate difference and sum  
  real(8) :: DD1, DD2   !denominators
  real(8) :: ReK, ImK, tmp_s, tmp_a 

   do ib1=1,ek%nband_max !loop over bands (these will be traced over)
      ! if the band is not contained in the optical matrices just do nothing
      if (ib1 < ek%nbopt_min) cycle
      if (ib1 > ek%nbopt_max) cycle
      resp%z1=real(sct%z,8)
      resp%gamma1=resp%z1*real(sct%gam(iT,ib1),8)
      resp%aqp1=resp%z1*real(ek%band(ik,ib1),8) !in a symmetric SC mu=0
      ! if the band is unoccupied cycle
      if(resp%aqp1 > mu) cycle
      resp%zarg=0.5d0+beta2p*((ci*resp%aqp1)+resp%gamma1)
      do ipg=1,1 
         resp%ctmp=wpsipg(resp%zarg,ipg)
         resp%RePolyGamma1(ipg)=real(resp%ctmp,8)
         resp%ImPolyGamma1(ipg)=imag(resp%ctmp)
      enddo
      
      ! compute transport kernels (omega-part)
      ! 
      do ib2=1,ek%nband_max
         if (ib2 < ek%nbopt_min) cycle
         if (ib2 > ek%nbopt_max) cycle
         if (ib2 == ib1 ) cycle

         !second band variables and derived quantities
         resp%z2=real(sct%z,8)
         resp%gamma2=resp%z2*resp%gamma1   !only one gamma required !real(sct%gam(iT,ib2),8)
         resp%aqp2=resp%z2*real(ek%band(ik,ib2),8) !in a symmetric SC mu=0
         ! if the second state is occupied cycle (interband contribution)
         if(resp%aqp2 < mu) cycle
         resp%zarg=0.5d0+beta2p*(ci*resp%aqp2+resp%gamma2)
         do ipg=1,1 
            resp%ctmp=wpsipg(resp%zarg,ipg)
            resp%RePolyGamma2(ipg)=real(resp%ctmp,8)
            resp%ImPolyGamma2(ipg)=imag(resp%ctmp)
         enddo

         Dqp    = resp%aqp1 - resp%aqp2     !Delta csi in eq
         !DD1 = 1.0d0/(gap**2 + 4.0d0*(resp%gamma1**2) )
         DD1 = 1.0d0/(Dqp**2 + 4.0d0*(resp%gamma1**2) )

         tmp_s = DD1*((resp%z1 * resp%gamma1)**2)*beta/(pi**3)
         tmp_a = DD1*((resp%z1 * beta)**2)/(2.0d0*(pi**3))
           

         resp%s_ker = tmp_s * ( (resp%RePolyGamma2(1) + resp%RePolyGamma1(1))/(2.0d0*resp%gamma1) &
                    + (resp%ImPolyGamma2(1) - resp%ImPolyGamma1(1))/Dqp )
                    !- (resp%ImPolyGamma2(1) - resp%ImPolyGamma1(1))/gap ) !only at the Gamma point!!
         resp%a_ker = tmp_a * ( resp%gamma1*(resp%aqp1*resp%RePolyGamma1(1) + resp%aqp2*resp%RePolyGamma2(1)) &
                    + (resp%gamma1**2)*(resp%aqp1+resp%aqp2)*(resp%ImPolyGamma2(1)-resp%ImPolyGamma1(1))/Dqp  &
                    + (resp%gamma1**3)*2.0d0*(resp%RePolyGamma1(1) - resp%RePolyGamma2(1))/Dqp )

         !only at the Gamma point!!           
         !resp%a_ker = tmp_a * ( resp%gamma1*abs(resp%aqp1)*(resp%RePolyGamma2(1)-resp%RePolyGamma1(1)) &
         !           + abs(resp%aqp1)*(resp%gamma1**3)*(resp%RePolyGamma2(1)-resp%RePolyGamma1(1))/(resp%aqp1**2) ) 
                    
         ! B = 0  
         !tmp=vka(ik,ib,ialpha)*vka(ik,ib,ibeta)
         do ialpha=1,nalpha
            if (algo%ltbind) then
               resp%tmp=ek%Mopt(ialpha,ik, ib1, ib2)*ek%Mopt(ialpha,ik, ib1, ib2)
            else
               resp%tmp=ek%Mopt(ialpha,ik, ib1, ib2) !the optical matrix elements given by Wien2k are squared already
            endif
         
            resp%s_tmp(ik,ib1,ialpha,ialpha)=resp%s_tmp(ik,ib1,ialpha,ialpha) + (resp%s_ker * resp%tmp)  
            resp%a_tmp(ik,ib1,ialpha,ialpha)=resp%a_tmp(ik,ib1,ialpha,ialpha) + (resp%a_ker * resp%tmp) 
               
            do ibeta=ialpha+1,nalpha
               resp%tmp=ek%Mopt(ialpha+ibeta+1,ik, ib1, ib2) 
               resp%s_tmp(ik,ib1,ialpha,ibeta)=resp%s_tmp(ik,ib1,ialpha,ibeta) + (resp%s_ker * resp%tmp)  
               resp%a_tmp(ik,ib1,ialpha,ibeta)=resp%a_tmp(ik,ib1,ialpha,ibeta) + (resp%a_ker * resp%tmp) 
            enddo !ibeta  
         enddo ! ialpha          
 
      enddo !ib2
   enddo ! ib1

end subroutine respinterkm_symm

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! RESDERTET_SYMM
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This routine 
! evaluates the conduction kernel derivatives
! with respect to beta=1/(kB*T) under the assumption
! that both the chemical potential and the 
! scattering rate are temperature independent 
! (hence it can give meaningful results only for
! a symmetric semiconductor).
! the first derivative is saved in resp%s
! the second derivative is saved in resp%a 
!
subroutine resdertet_symm(mu, iT, itet, thdr, algo, ek, sct, resp)
  use response
  use types
  use params
  implicit none 
  type (dp_respinter) :: resp 
  type (algorithm) :: algo
  type (tetramesh) :: thdr
  type (edisp) :: ek
  type (scatrate) :: sct
  real(8), intent(in) :: mu
  integer, intent(in) :: iT
  integer, intent(in) :: itet
  integer :: iband, ik, ipg
  integer :: ialpha,ibeta 
  complex(8),external  :: wpsipg
!local variables
  real(8), allocatable :: s_tmp_tetra(:,:,:,:),  a_tmp_tetra(:,:,:,:)

  !allocation  
  if(.not. allocated(s_tmp_tetra)) allocate(s_tmp_tetra(4,ek%nband_max,3,3))
  if(.not. allocated(a_tmp_tetra)) allocate(a_tmp_tetra(4,ek%nband_max,3,3))
  !initialisation
  s_tmp_tetra=0.0d0 ; a_tmp_tetra=0.0d0

   do ik=1,4  !loop over corners of the tetrahedron
      do iband=1,ek%nband_max !loop over bands (these will be traced over)

         ! if the band is not contained in the optical matrices just do nothing
         if (iband < ek%nbopt_min) cycle
         if (iband > ek%nbopt_max) cycle
         resp%z=real(sct%z,8)
         resp%gamma=resp%z*real(sct%gam(iT,iband),8)
         ! pre-compute all needed digamma functions   
         resp%aqp=resp%z*real(ek%band(thdr%idtet(ik,itet),iband)-mu,8)
         resp%zarg=0.5d0+beta2p*(ci*resp%aqp+resp%gamma)
         do ipg=1,4 
            resp%ctmp=wpsipg(resp%zarg,ipg)
            resp%RePolyGamma(ipg)=real(resp%ctmp,8)
            resp%ImPolyGamma(ipg)=imag(resp%ctmp)
         enddo
         
         ! compute transport kernel derivatives (omega-part)
         ! 
         ! 1st derivative w.r.t. beta
         resp%tmp=resp%z**2 / (4.d0*pi**3) ! for the 2nd derivative there is a factor 1/pi missing   
         resp%s_ker = resp%tmp * ((1.0d0/resp%gamma)*resp%RePolyGamma(1) - beta2p*resp%RePolyGamma(2) &
                    - beta2p*(resp%aqp/resp%gamma)*resp%ImPolyGamma(2) + resp%aqp*(beta2p**2)*resp%ImPolyGamma(3) &
                    - resp%gamma*(beta2p**2)*resp%RePolyGamma(3) ) 

         ! 
         ! 2nd derivative w.r.t. beta
         resp%tmp=resp%z**2 / (4.d0*pi**4)    
         resp%a_ker = resp%tmp * (-(resp%aqp/resp%gamma)*resp%ImPolyGamma(2) &
                    - 0.5d0*resp%gamma*beta2p*(3.0d0+(resp%aqp/resp%gamma)**2 )*resp%RePolyGamma(3) &
                    + beta2p*resp%aqp*resp%ImPolyGamma(3) + (beta2p**2)*resp%aqp*resp%gamma*resp%ImPolyGamma(4) &
                    + 0.5d0*(beta2p**2)*(resp%aqp**2 - resp%gamma**2)*resp%RePolyGamma(4) )
         
         !only the xx component has been evaluated
         do ialpha=1,1
            do ibeta=1,1
         
               !tmp=vka(ik,ib,ialpha)*vka(ik,ib,ibeta)
               if (algo%ltbind) then
                  resp%tmp=ek%Mopt(ialpha,thdr%idtet(ik,itet), iband, iband)*ek%Mopt(ialpha,thdr%idtet(ik,itet), iband, iband)
               else
                  write(*,*) 'resdertet_symm: the expression for the derivatives is only valid for a symmetric SC'
                  STOP
               endif
         
               s_tmp_tetra(ik,iband,ialpha,ibeta)=resp%s_ker * resp%tmp  
               a_tmp_tetra(ik,iband,ialpha,ibeta)=resp%a_ker * resp%tmp 
               
            enddo !ibeta  
         enddo ! ialpha           

         ! Now copy the local variable into the datastructure that will be passed to the interptra_mu
         resp%s_tmp(ik,iband,1,1) = s_tmp_tetra(ik,iband,1,1)
         resp%a_tmp(ik,iband,1,1) = a_tmp_tetra(ik,iband,1,1)
         
      enddo ! iband
   enddo ! ik   

end subroutine resdertet_symm

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! RESDERTET
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This routine 
! evaluates the conduction kernel derivatives
! with respect to beta=1/(kB*T) with a temperature 
! dependent chemical potential and  
! scattering rate. The 2nd derivative of 
! the chemical potential is neglected
! for the scattering rate it is evaluated assuming
! gamma(T) = gc0 + gc2*T^2, so 
! gam2dot = 6*gc2/(kB^2 * beta^4) 
! the first derivative of the conductivity is saved in resp%s
! the second derivative of the conductivity is saved in resp%a 
!
subroutine resdertet(iT, itet, nalpha, thdr, algo, ek, sct, resp)
  use response
  use types
  use params
  implicit none 
  type (dp_respinter) :: resp 
  type (algorithm) :: algo
  type (tetramesh) :: thdr
  type (edisp) :: ek
  type (scatrate) :: sct
  integer, intent(in) :: iT
  integer, intent(in) :: itet
  integer, intent(in) :: nalpha
  complex(8),external  :: wpsipg
!local variables
  real(8), allocatable :: s_tmp_tetra(:,:,:,:),  a_tmp_tetra(:,:,:,:)
  real(8) :: muder, mudot   ! derivative of the chemical potential w.r.t. T, beta 
  real(8) :: gamder, gamdot ! derivative of the scattering rate w.r.t. T, beta 
  real(8) :: dlogg          ! gammadot/gamma
  real(8) :: gam2dot        ! 2nd derivative of gamma w.r.t. beta, assuming gamma(T) = gc0 + gc2*T^2
  real(8) :: csim           ! aqp/beta - mudot
  integer :: iband, ik, ipg
  integer :: ialpha,ibeta 
  
  !allocation  
  if(.not. allocated(s_tmp_tetra)) allocate(s_tmp_tetra(4,ek%nband_max,3,3))
  if(.not. allocated(a_tmp_tetra)) allocate(a_tmp_tetra(4,ek%nband_max,3,3))
  !initialisation
  s_tmp_tetra=0.0d0 ; a_tmp_tetra=0.0d0
  
  !need to call this subroutine for iT<nT
  if(iT == sct%nT) STOP !hopefully this case has been taken care of before calling the routine 
  muder = (sct%mu(iT+1)-sct%mu(iT))/sct%dT 
  mudot = -kB*muder*((sct%TT(iT))**2)

   do ik=1,4  !loop over corners of the tetrahedron
      do iband=1,ek%nband_max !loop over bands (these will be traced over)

         ! if the band is not contained in the optical matrices just do nothing
         if (iband < ek%nbopt_min) cycle
         if (iband > ek%nbopt_max) cycle
         resp%z=real(sct%z,8)
         resp%gamma=resp%z*real(sct%gam(iT,iband),8)
         gamder = (sct%gam(iT+1,iband)-sct%gam(iT,iband))/sct%dT !safeguard condition set above
         gamdot = -kB*gamder*((sct%TT(iT))**2)
         dlogg  = gamdot/resp%gamma
         if (sct%ng ==2) then
            gam2dot = 6.0d0*sct%gc(2)/((kB**2)*(beta**4))
         else
            gam2dot = 0.0d0
         endif

         resp%aqp=resp%z*real(ek%band(thdr%idtet(ik,itet),iband)-sct%mu(iT),8)
         csim = (resp%aqp/beta)-mudot
         ! pre-compute all needed digamma functions   
         resp%zarg=0.5d0+beta2p*(ci*resp%aqp+resp%gamma)
         do ipg=1,4 
            resp%ctmp=wpsipg(resp%zarg,ipg)
            resp%RePolyGamma(ipg)=real(resp%ctmp,8)
            resp%ImPolyGamma(ipg)=imag(resp%ctmp)
         enddo
         
         ! compute transport kernel derivatives 
           
         ! 1st derivative w.r.t. beta (there is a term 2/beta * sigma that needs to be added up)
         resp%tmp=((resp%z*beta)**2) / (8.d0*pi**4)   
         resp%s_ker = resp%tmp*((dlogg+1.0d0/beta)*(resp%RePolyGamma(2)-(1.0d0/(beta2p*resp%gamma))*resp%RePolyGamma(1)) &
                    - (csim/resp%gamma)*resp%ImPolyGamma(2) + beta2p*csim*resp%ImPolyGamma(3) &
                    - (dlogg + 1.0d0/beta)*beta2p*resp%gamma*resp%RePolyGamma(3) )


         ! 2nd derivative w.r.t. beta (there is a term 2/beta^2 * sigma that needs to be added up)
         resp%tmp=((resp%z*beta)**2) / (8.d0*pi**4)   
         resp%a_ker = resp%tmp*( ((2.0d0*dlogg/beta) + (2.0d0/(beta**2)) + (gam2dot/resp%gamma) ) &
                    * (resp%RePolyGamma(2) - (1.0d0/(beta2p*resp%gamma))*resp%RePolyGamma(1)) &
                    + ((csim*(2.0d0*dlogg - 3.0d0/beta)/resp%gamma) + (csim/(beta*resp%gamma)))*resp%ImPolyGamma(2) &
                    - beta2p*(resp%gamma*(3.0d0/(beta**2) + 4.0d0*dlogg/beta + gam2dot/resp%gamma - dlogg**2 ) &
                              + (csim**2)/resp%gamma )*resp%RePolyGamma(3) &
                    - (2.0d0*beta2p*csim*dlogg + (4.0d0*mudot -2.0d0*resp%aqp/beta)/(2.0d0*pi) )*resp%ImPolyGamma(3) &
                    + (beta2p**2)*(csim**2 - (gamdot + resp%gamma/beta)**2)*resp%RePolyGamma(4) &
                    + 2.0d0*(beta2p**2)*(csim*(gamdot + resp%gamma/beta))*resp%ImPolyGamma(4) )

         ! Now add the missing terms:
         resp%tmp=(((resp%z*beta)**2)/(8.d0*pi**4)) * ((1.0d0/(beta2p*resp%gamma))*resp%RePolyGamma(1) - resp%RePolyGamma(2))
         resp%s_ker = resp%s_ker + (2.0d0*resp%tmp/beta)
         resp%a_ker = resp%a_ker + (2.0d0*resp%tmp/(beta**2))
         
         !tmp=vka(ik,ib,ialpha)*vka(ik,ib,ibeta)
         do ialpha=1,nalpha
         
            if (algo%ltbind) then
               resp%tmp=ek%Mopt(ialpha,thdr%idtet(ik,itet), iband, iband)*ek%Mopt(ialpha,thdr%idtet(ik,itet), iband, iband)
            else
               resp%tmp=ek%Mopt(ialpha,thdr%idtet(ik,itet), iband, iband) !the optical matrix elements given by Wien2k are squared already
            endif
         
            s_tmp_tetra(ik,iband,ialpha,ialpha)=resp%s_ker * resp%tmp  
            a_tmp_tetra(ik,iband,ialpha,ialpha)=resp%a_ker * resp%tmp 
               
            do ibeta=ialpha+1,nalpha   
               resp%tmp=ek%Mopt(ialpha+ibeta+1,thdr%idtet(ik,itet), iband, iband) 
               s_tmp_tetra(ik,iband,ialpha,ibeta)=resp%s_ker * resp%tmp  
               a_tmp_tetra(ik,iband,ialpha,ibeta)=resp%a_ker * resp%tmp 
            enddo !ibeta  
         enddo ! ialpha           

         ! Now copy the local variable into the datastructure that will be passed to the interptra_mu
         resp%s_tmp(ik,iband,:,:) = s_tmp_tetra(ik,iband,:,:)
         resp%a_tmp(ik,iband,:,:) = a_tmp_tetra(ik,iband,:,:)
         
      enddo ! iband
   enddo ! ik   

end subroutine resdertet

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! RESDERKM_SYMM
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This routine 
! evaluates the conduction kernel derivatives
! with respect to beta=1/(kB*T) under the assumption
! that both the chemical potential and the 
! scattering rate are temperature independent 
! (hence it can give meaningful results only for
! a symmetric semiconductor).
! the first derivative is saved in resp%s
! the second derivative is saved in resp%a 
!
subroutine resderkm_symm(mu, iT, ik, algo, ek, sct, resp)
  use response
  use types
  use params
  implicit none 
  type (dp_respinter) :: resp 
  type (algorithm) :: algo
  type (edisp) :: ek
  type (scatrate) :: sct
  real(8), intent(in) :: mu
  integer, intent(in) :: iT
  integer, intent(in) :: ik
  integer :: iband, ipg
  integer :: ialpha,ibeta 
  complex(8),external  :: wpsipg

   do iband=1,ek%nband_max !loop over bands (these will be traced over)

     ! if the band is not contained in the optical matrices just do nothing
      if (iband < ek%nbopt_min) cycle
      if (iband > ek%nbopt_max) cycle

      resp%z=real(sct%z,8)
      resp%gamma=resp%z*real(sct%gam(iT,iband),8)
      ! pre-compute all needed digamma functions   
      resp%aqp=resp%z*real(ek%band(ik,iband)-mu,8)
      resp%zarg=0.5d0+beta2p*(ci*resp%aqp+resp%gamma)
      do ipg=1,4 
         resp%ctmp=wpsipg(resp%zarg,ipg)
         resp%RePolyGamma(ipg)=real(resp%ctmp,8)
         resp%ImPolyGamma(ipg)=imag(resp%ctmp)
      enddo
      
      ! compute transport kernel derivatives (omega-part)
        
      ! 1st derivative w.r.t. beta
      resp%tmp=resp%z**2 / (4.d0*pi**3) ! for the 2nd derivative there is a factor 1/pi missing   
      resp%s_ker = resp%tmp * ((1.0d0/resp%gamma)*resp%RePolyGamma(1) - beta2p*resp%RePolyGamma(2) &
                 - beta2p*(resp%aqp/resp%gamma)*resp%ImPolyGamma(2) + resp%aqp*(beta2p**2)*resp%ImPolyGamma(3) &
                 - resp%gamma*(beta2p**2)*resp%RePolyGamma(3) ) 

      ! 2nd derivative w.r.t. beta
      resp%tmp=resp%z**2 / (4.d0*pi**4)    
      resp%a_ker = resp%tmp * (-(resp%aqp/resp%gamma)*resp%ImPolyGamma(2) &
                 - 0.5d0*resp%gamma*beta2p*(3.0d0+(resp%aqp/resp%gamma)**2 )*resp%RePolyGamma(3) &
                 + beta2p*resp%aqp*resp%ImPolyGamma(3) + (beta2p**2)*resp%aqp*resp%gamma*resp%ImPolyGamma(4) &
                 + 0.5d0*(beta2p**2)*(resp%aqp**2 - resp%gamma**2)*resp%RePolyGamma(4) )
      
      !only the xx component has been evaluated
      do ialpha=1,1
         do ibeta=1,1
      
            !tmp=vka(ik,ib,ialpha)*vka(ik,ib,ibeta)
            if (algo%ltbind) then
               resp%tmp=ek%Mopt(ialpha, ik, iband, iband)*ek%Mopt(ialpha, ik, iband, iband)
            else
               write(*,*) 'resderkm_symm: the expression for the derivatives is only valid for a symmetric SC'
               STOP
            endif
      
            resp%s_tmp(ik,iband,ialpha,ibeta)=resp%s_ker * resp%tmp  
            resp%a_tmp(ik,iband,ialpha,ibeta)=resp%a_ker * resp%tmp 
            
         enddo !ibeta  
      enddo ! ialpha           

   enddo ! iband

end subroutine resderkm_symm

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! RESDERKM
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This routine 
! evaluates the conduction kernel derivatives
! with respect to beta=1/(kB*T) with a temperature 
! dependent chemical potential and  
! scattering rate. The 2nd derivatives of both
! the chemical potential and gamma are neglected 
! the first derivative of the conductivity is saved in resp%s
! the second derivative of the conductivity is saved in resp%a 
!
subroutine resderkm(iT, ik, nalpha, algo, ek, sct, resp)
  use response
  use types
  use params
  implicit none 
  type (dp_respinter) :: resp 
  type (algorithm) :: algo
  type (edisp) :: ek
  type (scatrate) :: sct
  integer, intent(in) :: iT
  integer, intent(in) :: ik
  integer, intent(in) :: nalpha
  complex(8),external  :: wpsipg
  !complex(16),external :: wpsipghp
!local variables
  real(8) :: muder, mudot   ! derivative of the chemical potential w.r.t. T, beta 
  real(8) :: gamder, gamdot ! derivative of the scattering rate w.r.t. T, beta 
  real(8) :: dlogg          ! gammadot/gamma
  real(8) :: csim           ! aqp/beta - mudot
  integer :: iband, ipg
  integer :: ialpha,ibeta 
  

  !need to call this subroutine for iT<nT
  if(iT == sct%nT) stop !hopefully this case has been taken care of before calling the routine 
  muder = (sct%mu(iT+1)-sct%mu(iT))/sct%dT 
  mudot = -kB*muder*((sct%TT(iT))**2)

   do iband=1,ek%nband_max !loop over bands (these will be traced over)

      ! if the band is not contained in the optical matrices just do nothing
      if (iband < ek%nbopt_min) cycle
      if (iband > ek%nbopt_max) cycle

      resp%z=real(sct%z,8)
      resp%gamma=resp%z*real(sct%gam(iT,iband),8)
      gamder = (sct%gam(iT+1,iband)-sct%gam(iT,iband))/sct%dT !safeguard condition set above
      gamdot = -kB*gamder*((sct%TT(iT))**2)
      dlogg  = gamdot/resp%gamma

      resp%aqp=resp%z*real(ek%band(ik,iband)-sct%mu(iT),8)
      csim = (resp%aqp/beta)-mudot
      ! pre-compute all needed digamma functions   
      resp%zarg=0.5d0+beta2p*(ci*resp%aqp+resp%gamma)
      do ipg=1,4 
         resp%ctmp=wpsipg(resp%zarg,ipg)
         resp%RePolyGamma(ipg)=real(resp%ctmp,8)
         resp%ImPolyGamma(ipg)=imag(resp%ctmp)
      enddo
      
      ! compute transport kernel derivatives 
        
      ! 1st derivative w.r.t. beta (there is a term 2/beta * sigma that needs to be added up)
      resp%tmp=((resp%z*beta)**2) / (8.d0*pi**4)   
      resp%s_ker = resp%tmp*((dlogg+1.0d0/beta)*(resp%RePolyGamma(2)-(1.0d0/(beta2p*resp%gamma))*resp%RePolyGamma(1)) &
                 - (csim/resp%gamma)*resp%ImPolyGamma(2) + beta2p*csim*resp%ImPolyGamma(3) &
                 - (dlogg + 1.0d0/beta)*beta2p*resp%gamma*resp%RePolyGamma(3) )

      ! 2nd derivative w.r.t. beta (there is a term 2/beta^2 * sigma that needs to be added up)
      resp%tmp=((resp%z*beta)**2) / (8.d0*pi**4)   
      resp%a_ker = resp%tmp*( ((2.0d0*dlogg/beta) + (2.0d0/(beta**2))) &
                 * (resp%RePolyGamma(2) - (1.0d0/(beta2p*resp%gamma))*resp%RePolyGamma(1)) &
                 + ((csim*(2.0d0*dlogg - 3.0d0/beta)/resp%gamma) + (csim/(beta*resp%gamma)))*resp%ImPolyGamma(2) &
                 - beta2p*(resp%gamma*(3.0d0/(beta**2) + 4.0d0*dlogg/beta) + (csim**2)/resp%gamma )*resp%RePolyGamma(3) &
                 - (2.0d0*beta2p*csim*dlogg + (4.0d0*mudot -2.0d0*resp%aqp/beta)/(2.0d0*pi) )*resp%ImPolyGamma(3) &
                 + (beta2p**2)*(csim**2 - (gamdot + resp%gamma/beta)**2)*resp%RePolyGamma(4) &
                 + 2.0d0*(beta2p**2)*(csim*(gamdot + resp%gamma/beta))*resp%ImPolyGamma(4) )


      ! Now add the missing terms:
      resp%tmp=(((resp%z*beta)**2)/(8.d0*pi**4)) * ((1.0d0/(beta2p*resp%gamma))*resp%RePolyGamma(1) - resp%RePolyGamma(2))
      resp%s_ker = resp%s_ker + (2.0d0*resp%tmp/beta)
      resp%a_ker = resp%a_ker + (2.0d0*resp%tmp/(beta**2))
      
      do ialpha=1,nalpha
            if (algo%ltbind) then
               !the expression requires only the diagonal of the optical matrix elements because a trace is evaluated 
               resp%tmp=ek%Mopt(ialpha, ik, iband, iband)*ek%Mopt(ialpha, ik, iband, iband)
            else
               !the expression requires only the diagonal of the optical matrix elements because a trace is evaluated 
               resp%tmp=ek%Mopt(ialpha, ik, iband, iband) !the optical matrix elements given by Wien2k are squared already
            endif
      
            resp%s_tmp(ik,iband,ialpha,ialpha)=resp%s_ker * resp%tmp  
            resp%a_tmp(ik,iband,ialpha,ialpha)=resp%a_ker * resp%tmp 

         do ibeta=ialpha+1,nalpha   
            resp%tmp=ek%Mopt(ialpha+ibeta+1,ik, iband, iband) 
            resp%s_tmp(ik,iband,ialpha,ibeta)=resp%s_ker * resp%tmp  
            resp%a_tmp(ik,iband,ialpha,ibeta)=resp%a_ker * resp%tmp 
         enddo !ibeta  
      enddo ! ialpha           

   enddo ! iband

end subroutine resderkm

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! FINDRHOFLEX
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This routine tries to solve the equation
! 2*(d sigma/d beta)^2 - sigma*(d^2 sigma/d beta^2) = 0
! for a given temperature and stores away the first 
! temperature that fulfills it (in sct%Tstar)
! since the equation is difficult to solve point by point
! (i.e. w/o knowing the global behaviour with T)
! I neglect the term 2*(d sigma/d beta)^2 and find the T
! for which the 2nd derivative changes sign
! the first derivative of the conductivity saved in resp2%s
! the second derivative of the conductivity saved in resp2%a 
!
subroutine findrhoflex(iT, resp1, resp2, sct)
   use types
   use response
   implicit none

   integer, intent(in) :: iT 
   type(dp_resp) :: resp1      !contains the intraband conductivity
   type(dp_respinter) :: resp2 !contains its derivatives
   type(scatrate):: sct
   !local variables
   real(8), parameter :: tol=1.0d-12
   real(8) :: tmp, tmp1, tmp2

   if (.not. allocated(sct%d0)) then 
      allocate(sct%d0(1:sct%nT))
      allocate(sct%d1(1:sct%nT))
      allocate(sct%d2(1:sct%nT))
   endif
 
   tmp1 = 2.0d0*((resp2%s_tot(1,1))**2)
   tmp2 = resp1%s_tot(1,1)*resp2%a_tot(1,1)
   tmp  = tmp1 - tmp2
   sct%d0(iT) = tmp
   sct%d1(iT) = tmp1
   sct%d2(iT) = tmp2
   !if (tmp2<0.0d0 ) then
   !   sct%Tstar = sct%TT(iT)
   !else
   !   sct%Tstar = 0.0d0 
   !endif 

end subroutine findrhoflex

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! FINDRHOFLAT
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This routine finds the value of T for which 
! the derivative  d sigma/d beta changes sign
! (onset of saturation for the resistivity)
! the first derivative of the conductivity saved in resp2%s
! the second derivative of the conductivity saved in resp2%a 
!
subroutine findrhoflat(iT, resp1, resp2, sct)
   use types
   use response
   implicit none

   integer, intent(in) :: iT 
   type(dp_resp) :: resp1      !contains the intraband conductivity
   type(dp_respinter) :: resp2 !contains its derivatives
   type(scatrate):: sct
   !local variables
   real(8), parameter :: tol=1.0d-12
   real(8) :: tmp

   tmp = -resp2%s_tot(1,1)
   if (tmp < tol) then
      sct%Tflat = sct%TT(iT)
   else
      sct%Tflat = 0.0d0 
   endif 

end subroutine findrhoflat

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! FINDDRHOMAX
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This routine evaluates the 1st derivative of the
! resistivity at a given temperature and saves it 
! in the variable sct%drhodT 
! the first derivative of the conductivity saved in resp2%s
! the second derivative of the conductivity saved in resp2%a 
!
subroutine finddrhomax(iT, resp1, resp2, sct, drhodT)
   use params
   use types
   use response
   implicit none

   integer, intent(in) :: iT 
   type(dp_resp) :: resp1      !contains the intraband conductivity
   type(dp_respinter) :: resp2 !contains its derivatives
   type(scatrate):: sct
   real(8), intent(out) :: drhodT(sct%nT)

   !drhodT(iT) = resp2%s_tot(1,1)/(kB*((sct%TT(iT)*resp1%s_tot(1,1))**2))
   drhodT(iT) = -resp2%s_tot(1,1)/(resp1%s_tot(1,1)**2)

end subroutine finddrhomax

subroutine globfac(icubic, algo, mesh, resp, hpresp)
   use params
   use types
   use response
   use estruct, only:vol
   implicit none

   integer :: icubic
   type(algorithm) :: algo
   type(kpointmesh) :: mesh
   class(dp_resp) :: resp
   type(qp_resp), optional ::hpresp
!local variables
   !integer :: ktot
   real(8) :: ktot
   real(8) :: fac,facB
   real(16):: facQ,facBQ

   fac   = 2.d0 * pi * ( echarge / (vol*hbarevs)) * 1.d10 
   facB  = 2.d0 * pi**2 * ( echarge / (vol*hbarevs) ) * (1.d-10 / hbareVs)
   facQ  = 2.q0 * piQ * ( real(echarge,16) / real(vol*hbarevs,16)) * 1.q10
   facBQ = 2.q0 * piQ**2 * ( real(echarge,16) / real(vol*hbarevs,16)) *  (1.q-10 / real(hbareVs,16))
  
   if (algo%ltetra) then
      ktot=1.0d0
   else
      ktot=real(mesh%ktot) 
   endif

   !global symmetries of a cubic system
   if (icubic==1) then
      resp%s(:,2,2)=resp%s(:,1,1)
      resp%s(:,3,3)=resp%s(:,1,1)
      resp%a(:,2,2)=resp%a(:,1,1)
      resp%a(:,3,3)=resp%a(:,1,1)
      resp%s_tot(2,2)=resp%s_tot(1,1)   
      resp%s_tot(3,3)=resp%s_tot(1,1)
      resp%a_tot(2,2)=resp%a_tot(1,1)
      resp%a_tot(3,3)=resp%a_tot(1,1)
      if(algo%lBfield) then
         resp%sB(:,1,3)=resp%sB(:,1,2)
         resp%sB(:,2,3)=resp%sB(:,1,2)
         resp%aB(:,1,3)=resp%aB(:,1,2)
         resp%aB(:,2,3)=resp%aB(:,1,2)
         resp%sB_tot(1,3)=resp%sB_tot(1,2)
         resp%sB_tot(2,3)=resp%sB_tot(1,2)
         resp%aB_tot(1,3)=resp%aB_tot(1,2)
         resp%aB_tot(2,3)=resp%aB_tot(1,2)
      endif
      if (present(hpresp)) then
          hpresp%s(:,2,2)=hpresp%s(:,1,1)
          hpresp%s(:,3,3)=hpresp%s(:,1,1)
          hpresp%a(:,2,2)=hpresp%a(:,1,1)
          hpresp%a(:,3,3)=hpresp%a(:,1,1)
          hpresp%s_tot(2,2)=hpresp%s_tot(1,1)
          hpresp%s_tot(3,3)=hpresp%s_tot(1,1)
          hpresp%a_tot(2,2)=hpresp%a_tot(1,1)
          hpresp%a_tot(3,3)=hpresp%a_tot(1,1)
          if(algo%lBfield) then
             hpresp%sB(:,1,3)=hpresp%sB(:,1,2)
             hpresp%sB(:,2,3)=hpresp%sB(:,1,2)
             hpresp%aB(:,1,3)=hpresp%aB(:,1,2)
             hpresp%aB(:,2,3)=hpresp%aB(:,1,2)
             hpresp%sB_tot(1,3)=hpresp%sB_tot(1,2)
             hpresp%sB_tot(2,3)=hpresp%sB_tot(1,2)
             hpresp%aB_tot(1,3)=hpresp%aB_tot(1,2)
             hpresp%aB_tot(2,3)=hpresp%aB_tot(1,2)
          endif
      endif
   endif

   resp%s  = resp%s/real(ktot,8) * fac ! --> sigma in 1/(Ohm m)     [vk's are in eV*Angstroem]   
   resp%a  = resp%a/real(ktot,8) * fac * ( - beta * kb)  ! --> S=alpha/sigma in units V/K (below conversion to mV/K for output of S) 
   resp%s_tot  = resp%s_tot/real(ktot,8) * fac  
   resp%a_tot  = resp%a_tot/real(ktot,8) * fac * ( - beta * kb)  
   if(algo%lBfield) then
      resp%sB = resp%sB/real(ktot,8) * facB
      resp%aB = resp%aB/real(ktot,8) * facB * ( - beta * kb)
      resp%sB_tot = resp%sB_tot/real(ktot,8) * facB
      resp%aB_tot = resp%aB_tot/real(ktot,8) * facB * ( - beta * kb)
   endif
 
   if (present(hpresp)) then
      hpresp%s  = hpresp%s/real(ktot,16) * facQ ! --> sigma in 1/(Ohm m)     [vk's are in eV*Angstroem]   
      hpresp%a  = hpresp%a/real(ktot,16) * facQ * ( - betaQ * kbQ)  ! --> S=alpha/sigma in units V/K (below conversion to mV/K for output of S) 
      hpresp%s_tot  = hpresp%s_tot/real(ktot,16) * facQ  
      hpresp%a_tot  = hpresp%a_tot/real(ktot,16) * facQ * ( - betaQ * kbQ)  
      if(algo%lBfield) then
         hpresp%sB = hpresp%sB/real(ktot,16) * facBQ
         hpresp%aB = hpresp%aB/real(ktot,16) * facBQ * ( - betaQ * kbQ)
         hpresp%sB_tot = hpresp%sB_tot/real(ktot,16) * facBQ
         hpresp%aB_tot = hpresp%aB_tot/real(ktot,16) * facBQ * ( - betaQ * kbQ)
      endif
   endif

end subroutine globfac

subroutine derresp(icubic, algo, resp, hpresp)
   use response
   use types
   implicit none

   integer :: icubic
   type(algorithm) :: algo
   class(dp_resp)  :: resp
   type(qp_resp),optional :: hpresp
!local variables
   integer :: ix
   
! In Seebeck: *1000 so as to yield [S]=mV/K        
     
     do ix=1,3
        if (icubic==0) then
           resp%Seebeck(ix)=1000.d0*resp%a_tot(ix,ix)/resp%s_tot(ix,ix) 
           if (present(hpresp)) hpresp%Seebeck(ix)=1000.q0*hpresp%a_tot(ix,ix)/hpresp%s_tot(ix,ix)
        else
           resp%Seebeck(ix)=1000.d0*resp%a_tot(1,1)/resp%s_tot(1,1) 
           if (present(hpresp)) hpresp%Seebeck(ix)=1000.q0*hpresp%a_tot(1,1)/hpresp%s_tot(1,1)
        endif
     enddo


!     1 = xy
!     2 = xz
!     3 = yz

     if (algo%lBfield) then
        resp%Nernst(1) = (resp%aB_tot(1,2)*resp%s_tot(1,1)-resp%a_tot(1,1)*resp%sB_tot(1,2))/(resp%s_tot(1,1)**2) 
        resp%Nernst(2) = (resp%aB_tot(1,3)*resp%s_tot(1,1)-resp%a_tot(1,1)*resp%sB_tot(1,3))/(resp%s_tot(1,1)**2) 
        resp%Nernst(3) = (resp%aB_tot(2,3)*resp%s_tot(2,2)-resp%a_tot(2,2)*resp%sB_tot(2,3))/(resp%s_tot(2,2)**2) 
        resp%Nernst = resp%Nernst * 1000.d0 ! V/K --> mV/K
        if (present(hpresp)) then
           hpresp%Nernst(1) = (hpresp%aB_tot(1,2)*hpresp%s_tot(1,1)-hpresp%a_tot(1,1)*hpresp%sB_tot(1,2))/(hpresp%s_tot(1,1)**2) 
           hpresp%Nernst(2) = (hpresp%aB_tot(1,3)*hpresp%s_tot(1,1)-hpresp%a_tot(1,1)*hpresp%sB_tot(1,3))/(hpresp%s_tot(1,1)**2) 
           hpresp%Nernst(3) = (hpresp%aB_tot(2,3)*hpresp%s_tot(2,2)-hpresp%a_tot(2,2)*hpresp%sB_tot(2,3))/(hpresp%s_tot(2,2)**2) 
           hpresp%Nernst = hpresp%Nernst * 1000.q0 ! V/K --> mV/K

           hpresp%Nernstpart(1) = ( hpresp%aB_tot(1,2)*hpresp%s_tot(1,1))/(hpresp%s_tot(1,1)**2)*1000.q0
           hpresp%Nernstpart(2) = (-hpresp%a_tot(1,1)*hpresp%sB_tot(1,2))/(hpresp%s_tot(1,1)**2)*1000.q0
        endif 
        resp%RH(1) = -resp%sB_tot(1,2)/(resp%s_tot(1,1)*resp%s_tot(2,2))
        resp%RH(2) = -resp%sB_tot(1,3)/(resp%s_tot(1,1)*resp%s_tot(3,3))
        resp%RH(3) = -resp%sB_tot(2,3)/(resp%s_tot(3,3)*resp%s_tot(3,3))
        resp%RH = resp%RH * 1.d+7 ! --> 10^-7 m^3/C
        if (present(hpresp)) then
           hpresp%RH(1) = -hpresp%sB_tot(1,2)/(hpresp%s_tot(1,1)*hpresp%s_tot(2,2))
           hpresp%RH(2) = -hpresp%sB_tot(1,3)/(hpresp%s_tot(1,1)*hpresp%s_tot(3,3))
           hpresp%RH(3) = -hpresp%sB_tot(2,3)/(hpresp%s_tot(3,3)*hpresp%s_tot(3,3))
           hpresp%RH = hpresp%RH * 1.q+7 ! --> 10^-7 m^3/C
        endif 
     endif
end subroutine derresp

subroutine wrtresp(iT, nalpha, algo, sct, resp, respinter, respBl, hpresp) 
   use types
   use response
   implicit none

   integer, intent(in) :: iT, nalpha
   type(algorithm) :: algo
   type(scatrate)  :: sct
   type(dp_resp)   :: resp
   type(dp_respinter) :: respinter
   type(dp_resp)   :: respBl
   type(qp_resp),optional :: hpresp
!local variable
   integer :: ia
   real(8) :: T, det, tmp, Mtmp(3,3)
   real(16):: dpdet, dptmp

   T = sct%TT(iT)
   
   write(30,'(100E20.12)') T,(resp%s_tot(ia,ia), ia=1,nalpha )
   write(31,'(100E20.12)') T,(resp%s(:,ia,ia),   ia=1,nalpha )
   write(40,'(100E20.12)') T,(resp%a_tot(ia,ia), ia=1,nalpha )
   write(41,'(100E20.12)') T,(resp%a(:,ia,ia),   ia=1,nalpha )
   if (present(hpresp)) then
      write(130,'(100E20.12)') T,(hpresp%s_tot(ia,ia), ia=1,nalpha )
      write(131,'(100E20.12)') T,(hpresp%s(:,ia,ia),   ia=1,nalpha )
      write(140,'(100E20.12)') T,(hpresp%a_tot(ia,ia), ia=1,nalpha )
      write(141,'(100E20.12)') T,(hpresp%a(:,ia,ia),   ia=1,nalpha )
   endif
   write(230,'(100E20.12)') T,(respBl%s_tot(ia,ia), ia=1,nalpha )
   write(231,'(100E20.12)') T,(respBl%s(:,ia,ia),   ia=1,nalpha )
   write(240,'(100E20.12)') T,(respBl%a_tot(ia,ia), ia=1,nalpha )
   write(241,'(100E20.12)') T,(respBl%a(:,ia,ia),   ia=1,nalpha )

   write(330,'(100E20.12)') T,(respinter%s_tot(ia,ia), ia=1,nalpha )
   write(331,'(100E20.12)') T,(respinter%s(:,ia,ia),   ia=1,nalpha )
   write(340,'(100E20.12)') T,(respinter%a_tot(ia,ia), ia=1,nalpha )
   write(341,'(100E20.12)') T,(respinter%a(:,ia,ia),   ia=1,nalpha )
     
! XXX ACHTUNG not writing out all combinations...
   if (algo%lBfield) then
      write(50,'(100E20.12)') T,resp%sB_tot(1,2)
      write(51,'(100E20.12)') T,resp%sB(:,1,2)
      write(60,'(100E20.12)') T,resp%aB_tot(1,2)
      write(61,'(100E20.12)') T,resp%aB(:,1,2)
      if (present(hpresp)) then
         write(150,'(100E20.12)') T,hpresp%sB_tot(1,2)
         write(151,'(100E20.12)') T,hpresp%sB(:,1,2)
         write(160,'(100E20.12)') T,hpresp%aB_tot(1,2)
         write(161,'(100E20.12)') T,hpresp%aB(:,1,2)
      endif
      write(250,'(100E20.12)') T,respBl%sB_tot(1,2)
      write(251,'(100E20.12)') T,respBl%sB(:,1,2)
      write(260,'(100E20.12)') T,respBl%aB_tot(1,2)
      write(261,'(100E20.12)') T,respBl%aB(:,1,2)
   endif

   write(70,'(100E20.12)') T,resp%Seebeck(:) ! in mV/K
   write(370,'(100E20.12)') T,respinter%Seebeck(:) ! in mV/K
   if (algo%lBfield) then
      write(71,'(100E20.12)') T, resp%Nernst(:)  ! in mV/KT
      write(72,'(100E20.12)') T, resp%RH(:) ! in 10^-7 m^3/C
      write(73,'(100E20.12)') T, resp%sB_tot(1,2) / resp%s_tot(1,1) ! in 1/T
      write(74,'(100E20.12)') T, resp%aB_tot(1,2) / resp%a_tot(1,1) ! in 1/T
   endif
   if (nalpha==1) then
   ! diagonal conductivity tensor
      write(75,'(100E20.12)') T, (1.d0/resp%s_tot(1,1)) ! resistivity in Ohm m
      write(375,'(100E20.12)') T, (1.d0/respinter%s_tot(1,1)) ! resistivity in Ohm m
   else
      if (algo%ldebug) then
         write(75,'(100E20.12)') T, (1.d0/resp%s_tot(ia,ia),ia=1,nalpha) ! resistivity in Ohm m
         write(375,'(100E20.12)') T, (1.d0/respinter%s_tot(ia,ia),ia=1,nalpha) ! resistivity in Ohm m
      else
   ! evaluate the inverse of conductivity tensor
         det = (resp%s_tot(1,1)*resp%s_tot(2,2)*resp%s_tot(3,3)) &
               + (2.0d0*resp%s_tot(1,2)*resp%s_tot(2,3)*resp%s_tot(1,3)) &
               - (resp%s_tot(1,1)*(resp%s_tot(2,3)**2)) &
               - (resp%s_tot(2,2)*(resp%s_tot(1,3)**2)) &
               - (resp%s_tot(3,3)*(resp%s_tot(1,2)**2))
         tmp = ((resp%s_tot(3,3)*resp%s_tot(2,2)) - (resp%s_tot(2,3)**2))/det
         write(75,'(100E20.12)') T, tmp, det, resp%s_tot(1,2),resp%s_tot(2,3),resp%s_tot(1,3) ! resistivity in Ohm m (xx component of the inverted matrix)
         Mtmp(:,:) = resp%s_tot(:,:)+respinter%s_tot(:,:)
         det = (Mtmp(1,1)*Mtmp(2,2)*Mtmp(3,3)) &
               + (2.0d0*Mtmp(1,2)*Mtmp(2,3)*Mtmp(1,3)) &
               - (Mtmp(1,1)*(Mtmp(2,3)**2)) &
               - (Mtmp(2,2)*(Mtmp(1,3)**2)) &
               - (Mtmp(3,3)*(Mtmp(1,2)**2))
         tmp = ((Mtmp(3,3)*Mtmp(2,2)) - (Mtmp(2,3)**2))/det
         write(375,'(100E20.12)') T, tmp ! resistivity in Ohm m (xx component of the inverted matrix)
      endif
   endif

   if (present(hpresp)) then
      write(170,'(100E20.12)') T, real(hpresp%Seebeck(:),8) ! in mV/K
      if (algo%lBfield) then
         write(171,'(100E20.12)') T, real(hpresp%Nernst(:),8)  ! in mV/KT
         write(172,'(100E20.12)') T, real(hpresp%RH(:),8)      ! in 10^-7 m^3/C
         write(173,'(100E20.12)') T, real(hpresp%sB_tot(1,2) / resp%s_tot(1,1) ,8)    ! in 1/T
         write(174,'(100E20.12)') T, real(hpresp%aB_tot(1,2) / resp%a_tot(1,1) ,8)    ! in 1/T

         write(180,'(100E20.12)') T,real(hpresp%Nernstpart(1),8)  ! in mV/K
         write(181,'(100E20.12)') T,real(hpresp%Nernstpart(2),8)  ! in mV/K

      endif
      if (nalpha==1) then
      ! diagonal conductivity tensor
         write(175,'(100E20.12)') T, (real(1.q0/hpresp%s_tot(ia,ia),16),ia=1,nalpha) ! resistivity in Ohm m
      else
         if (algo%ldebug) then
            write(175,'(100E20.12)') T, (real(1.q0/hpresp%s_tot(ia,ia),16),ia=1,nalpha) ! resistivity in Ohm m
         else
      ! evaluate the inverse of conductivity tensor
            dpdet = (hpresp%s_tot(1,1)*hpresp%s_tot(2,2)*hpresp%s_tot(3,3)) &
                    + (2.0q0*hpresp%s_tot(1,2)*hpresp%s_tot(2,3)*hpresp%s_tot(1,3)) &
                    - (hpresp%s_tot(1,1)*(hpresp%s_tot(2,3)**2)) &
                    - (hpresp%s_tot(2,2)*(hpresp%s_tot(1,3)**2)) &
                    - (hpresp%s_tot(3,3)*(hpresp%s_tot(1,2)**2))
            dptmp = ((hpresp%s_tot(3,3)*hpresp%s_tot(2,2)) - (hpresp%s_tot(2,3)**2))/dpdet
            write(175,'(100E20.12)') T, real(dptmp,16),dpdet,hpresp%s_tot(1,2),hpresp%s_tot(2,3),hpresp%s_tot(1,3)! resistivity in Ohm m (xx component of the inverted matrix)
         endif
      endif
   endif 

   write(270,'(100E20.12)') T,respBl%Seebeck(:) ! in mV/K
   if (algo%lBfield) then
      write(271,'(100E20.12)') T, respBl%Nernst(:)  ! in mV/KT
      write(272,'(100E20.12)') T, respBl%RH(:) ! in 10^-7 m^3/C
      write(273,'(100E20.12)') T, respBl%sB_tot(1,2) / resp%s_tot(1,1) ! in 1/T
      write(274,'(100E20.12)') T, respBl%aB_tot(1,2) / resp%a_tot(1,1) ! in 1/T
   endif
   if (nalpha==1) then
   ! diagonal conductivity tensor
      write(275,'(100E20.12)') T, (1.d0/respBl%s_tot(ia,ia),ia=1,nalpha) ! resistivity in Ohm m
   else
      if (algo%ldebug) then
         write(275,'(100E20.12)') T, (1.d0/respBl%s_tot(ia,ia),ia=1,nalpha) ! resistivity in Ohm m
      else
   ! evaluate the inverse of conductivity tensor
         det = (respBl%s_tot(1,1)*respBl%s_tot(2,2)*respBl%s_tot(3,3)) &
               + (2.0d0*respBl%s_tot(1,2)*respBl%s_tot(2,3)*respBl%s_tot(1,3)) &
               - (respBl%s_tot(1,1)*respBl%s_tot(2,3)**2) &
               - (respBl%s_tot(2,2)*respBl%s_tot(1,3)**2) &
               - (respBl%s_tot(3,3)*respBl%s_tot(1,2)**2)
         tmp = ((respBl%s_tot(3,3)*respBl%s_tot(2,2)) - (respBl%s_tot(2,3)**2))/det
         write(275,'(100E20.12)') T, tmp ! resistivity in Ohm m (xx component of the inverted matrix)
      endif
   endif

end subroutine wrtresp

subroutine response_open_files(algo)
   use types
   implicit none
   type(algorithm) :: algo

   open(30,file='sigma_tot.dat',status='unknown')
   open(31,file='sigma_band_xx.dat',status='unknown')
   open(40,file='peltier_tot.dat',status='unknown')
   open(41,file='peltier_band_xx.dat',status='unknown')
! 50 for sxy                                                                                
   if (algo%lBfield) then
   open(50,file='sigmaB_tot.dat',status='unknown')
   open(51,file='sigmaB_band_xy.dat',status='unknown')
! 60 for axy                                                                                
   open(60,file='peltierB_tot.dat',status='unknown')
   open(61,file='peltierB_band_xy.dat',status='unknown')

   open(71,file='Nernst.dat',status='unknown')
   open(72,file='RH.dat',status='unknown')
   open(73,file='muH.dat',status='unknown') ! Hall mobility
   open(74,file='mut.dat',status='unknown') ! thermal counterpart of Hall mobility
   endif

! 70 for aux's: Seebeck, Hall, Nernst...                                                    
   open(70,file='Seebeck.dat',status='unknown')

   open(75,file='resistivity.dat',status='unknown')


! QUAD PRECISION
   open(130,file='sigma_tot.datQ',status='unknown')
   open(131,file='sigma_band_xx.datQ',status='unknown')
   open(140,file='peltier_tot.datQ',status='unknown')
   open(141,file='peltier_band_xx.datQ',status='unknown')
! 50 for sxy                                                                                
   if (algo%lBfield) then
   open(150,file='sigmaB_tot.datQ',status='unknown')
   open(151,file='sigmaB_band_xy.datQ',status='unknown')
! 60 for axy                                                                                
   open(160,file='peltierB_tot.datQ',status='unknown')
   open(161,file='peltierB_band_xy.datQ',status='unknown')
   open(171,file='Nernst.datQ',status='unknown')
   open(172,file='RH.datQ',status='unknown')
   open(173,file='muH.datQ',status='unknown') ! Hall mobility
   open(174,file='mut.datQ',status='unknown') ! thermal counterpart of Hall mobility
   open(180,file='Nernst_part1.datQ',status='unknown')
   open(181,file='Nernst_part2.datQ',status='unknown')
   endif
! 70 for aux's: Seebeck, Hall, Nernst...                                                    
   open(170,file='Seebeck.datQ',status='unknown')

   open(175,file='resistivity.datQ',status='unknown')

! DOUBLE PRECISION, BOLTZMANN MODE
   open(230,file='sigma_tot.datB',status='unknown')
   open(231,file='sigma_band_xx.datB',status='unknown')
   open(240,file='peltier_tot.datB',status='unknown')
   open(241,file='peltier_band_xx.datB',status='unknown')
! 50 for sxy                                                                                
   if (algo%lBfield) then
   open(250,file='sigmaB_tot.datB',status='unknown')
   open(251,file='sigmaB_band_xy.datB',status='unknown')
! 60 for axy                                                                                
   open(260,file='peltierB_tot.datB',status='unknown')
   open(261,file='peltierB_band_xy.datB',status='unknown')
   open(271,file='Nernst.datB',status='unknown')
   open(272,file='RH.datB',status='unknown')
   open(273,file='muH.datB',status='unknown') ! Hall mobility
   open(274,file='mut.datB',status='unknown') ! thermal counterpart of Hall mobility
   !open(280,file='Nernst_part1.datBQ',status='unknown')
   !open(281,file='Nernst_part2.datBQ',status='unknown')
   endif
! 70 for aux's: Seebeck, Hall, Nernst...                                                    
   open(270,file='Seebeck.datB',status='unknown')
   open(275,file='resistivity.datB',status='unknown')

! DOUBLE PRECISION, INTERBAND RESPONSE
   open(330,file='sigma_inter_tot.dat',status='unknown')
   open(331,file='sigma_interband_xx.dat',status='unknown')
   open(340,file='peltier_inter_tot.dat',status='unknown')
   open(341,file='peltier_interband_xx.dat',status='unknown')

! 70 for aux's: Seebeck, Hall, Nernst...                                                    
   open(370,file='Seebeck_inter.dat',status='unknown')

   open(375,file='resistivity_inter.dat',status='unknown')

end subroutine response_open_files

subroutine response_close_files(algo)
   use types
   implicit none
   type(algorithm) :: algo
   
   close(30)
   close(31)
   close(40)
   close(41)
   if (algo%lBfield) then
     close(50)
     close(51)
     close(60)
     close(61)
     close(71)
     close(72)
     close(73)
     close(74)
   endif
   close(70)
   close(75)

   close(130)
   close(131)
   close(140)
   close(141)
   if (algo%lBfield) then
     close(150)
     close(151)
     close(160)
     close(161)
     close(171)
     close(172)
     close(173)
     close(174)
   endif

   close(170)
   close(175)

   close(180)
   close(181)

!Boltzmann
   close(230)
   close(231)
   close(240)
   close(241)
   if (algo%lBfield) then
     close(250)
     close(251)
     close(260)
     close(261)
     close(271)
     close(272)
     close(273)
     close(274)
   endif
   close(270)
   close(275)
   !close(280)
   !close(281)

!Interband
   close(330)
   close(331)
   close(340)
   close(341)
   close(370)
   close(375)

end subroutine response_close_files


subroutine dpresp_alloc(lBfield, dpresp, nk, nband)
use response
implicit none
logical :: lBfield
type(dp_resp)::dpresp 
integer :: nk, nband

! allocate transport variables                                                              
allocate(dpresp%s_tmp(nk,nband,3,3))
allocate(dpresp%a_tmp(nk,nband,3,3))
allocate(dpresp%s(nband,3,3))
allocate(dpresp%a(nband,3,3))
if (nproc > 1) then
   allocate(dpresp%s_local(nband,3,3))
   allocate(dpresp%a_local(nband,3,3))
endif

if (lBfield) then
   allocate(dpresp%sB_tmp(nk,nband,3,3))
   allocate(dpresp%aB_tmp(nk,nband,3,3))
   allocate(dpresp%sB(nband,3,3))
   allocate(dpresp%aB(nband,3,3))
   if (nproc > 1) then
      allocate(dpresp%sB_local(nband,3,3))
      allocate(dpresp%aB_local(nband,3,3))
   endif
endif
end subroutine dpresp_alloc

subroutine qpresp_alloc(lBfield, qpresp, nk, nband)
use response
implicit none
logical :: lBfield
type(qp_resp)::qpresp 
integer :: nk, nband

! allocate transport variables                                                              
allocate(qpresp%s_tmp(nk,nband,3,3))
allocate(qpresp%a_tmp(nk,nband,3,3))
allocate(qpresp%s(nband,3,3))
allocate(qpresp%a(nband,3,3))
if (nproc > 1) then
   allocate(qpresp%s_local(nband,3,3))
   allocate(qpresp%a_local(nband,3,3))
endif

if (lBfield) then
   allocate(qpresp%sB_tmp(nk,nband,3,3))
   allocate(qpresp%aB_tmp(nk,nband,3,3))
   allocate(qpresp%sB(nband,3,3))
   allocate(qpresp%aB(nband,3,3))
   if (nproc > 1) then
      allocate(qpresp%sB_local(nband,3,3))
      allocate(qpresp%aB_local(nband,3,3))
   endif
endif
end subroutine qpresp_alloc


function dfermi(eps,beta)
  implicit none
  real(8) :: dfermi
  real(8) :: eps,beta

!  DFERMIQ=beta / ( QEXP(-beta*eps/2.q0) + QEXP(beta*eps/2.q0) )**2
  dfermi=beta / ( exp(-beta*eps/2.d0) + exp(beta*eps/2.d0) )**2

return
end function dfermi

