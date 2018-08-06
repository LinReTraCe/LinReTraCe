module Mmpi_org

#ifdef MPI
  use mpi
#endif

  integer, parameter   :: master = 0
  integer              :: myid, nproc
  integer, allocatable :: displs(:),rcounts(:)
  integer              :: iqstr, iqend

  integer              :: mpierr, nkthis
  character(len=5)     :: chmyid

  contains

  ! initialize mpi environment
  subroutine mpi_initialize()
    implicit none
#ifdef MPI
    call MPI_INIT(mpierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD,myid,mpierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,mpierr)
#else
    myid = 0
    nproc = 1
#endif
  end subroutine mpi_initialize

  ! close mpi environment
  subroutine mpi_close()
    implicit none
#ifdef MPI
    if (allocated(rcounts)) deallocate(rcounts)
    if (allocated(displs))  deallocate(displs)
    call MPI_FINALIZE( mpierr )
#endif
  end subroutine mpi_close

  subroutine prepare_chmyid(myid,chmyid)
    implicit none
    integer, intent(in) :: myid
    character(len=*)    :: chmyid
    write(chmyid,'(I5.5)') myid
    return
  end subroutine prepare_chmyid

  subroutine mpi_genkstep(nk)
    implicit none
    integer, intent(in) :: nk
    integer             :: i,nktmp
#ifdef MPI
    allocate(displs(nproc),rcounts(nproc))
    displs = 0
    rcounts = 0
    do i=1,nproc-1
       rcounts(i) = (nk - displs(i)) / (nproc+1-i) ! integer division
       displs(i+1) = displs(i) + rcounts(i)
    enddo
    rcounts(nproc) = (nk - displs(nproc))

    iqstr = displs(myid+1) + 1
    iqend = displs(myid+1) + rcounts(myid+1)
#else
    iqstr = 1
    iqend = nk
#endif
  end subroutine mpi_genkstep


!!! COMMUNICATION ROUTINES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! SUBROUTINE MPI_DVECSCAT
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! scatters a double precision vector
! across processors

#ifdef MPI
  subroutine mpi_dvecscat(dvec, n, sub_dvec, nk)
   implicit none
   integer :: n, nk
   real(8) :: dvec(n)
   real(8) :: sub_dvec(nk)
   integer :: ierror
   integer :: sendcount, recvcount

   sendcount=int(n/nproc)
   recvcount=sendcount
   call MPI_SCATTER(dvec(1),sendcount,MPI_double_precision, &
        sub_dvec,recvcount,MPI_double_precision,master,MPI_COMM_WORLD, ierror)

   if (ierror /= 0) then
      write (*,*) 'mpi_dvecscat returns error'
      STOP
   endif

  end subroutine mpi_dvecscat

!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! SUBROUTINE MPI_QVECSCAT
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! scatters a quad precision vector
! across processors
  subroutine mpi_Qvecscat(qvec, n, sub_qvec, nk)
   implicit none
   integer :: n, nk
   real(16) :: qvec(n)
   real(16) :: sub_qvec(nk)
   integer :: ierror
   integer :: sendcount, recvcount

   sendcount=int(n/nproc)
   recvcount=sendcount
   call MPI_SCATTER(qvec(1),sendcount,MPI_REAL16, &
        sub_qvec,recvcount,MPI_REAL16,master,MPI_COMM_WORLD, ierror)
   if (ierror /= 0) then
      write (*,*) 'mpi_Qvecscat returns error'
      STOP
   endif

  end subroutine mpi_Qvecscat

!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! SUBROUTINE MPI_GENKSTEP
!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  SUBROUTINE MPI_ALGOBCAST
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine mpi_algobcast(algo)
    use Mtypes
    implicit none
    type(algorithm) :: algo
    integer :: ierror

    ierror=0
    !write(100+myid,*)'algo%ldebug',algo%ldebug
    call MPI_BCAST(algo%ldebug ,1,MPI_LOGICAL,master,MPI_COMM_WORLD,ierror)
    if (ierror /= 0) then
       write (*,*) 'mpi_algobcast returns error/1'
       STOP
    endif
    !write(100+myid,*)'algo%ltbind',algo%ltbind
    call MPI_BCAST(algo%ltbind ,1,MPI_LOGICAL,master,MPI_COMM_WORLD,ierror)
    if (ierror /= 0) then
       write (*,*) 'mpi_algobcast returns error/2'
       STOP
    endif
    !write(100+myid,*)'algo%ltetra',algo%ltetra
    call MPI_BCAST(algo%ltetra ,1,MPI_LOGICAL,master,MPI_COMM_WORLD,ierror)
    if (ierror /= 0) then
       write (*,*) 'mpi_algobcast returns error/3'
       STOP
    endif
    !write(100+myid,*)'algo%lw2k',algo%lw2k
    call MPI_BCAST(algo%lw2k   ,1,MPI_LOGICAL,master,MPI_COMM_WORLD,ierror)
    if (ierror /= 0) then
       write (*,*) 'mpi_algobcast returns error/4'
       STOP
    endif
    !write(100+myid,*)'algo%loptic',algo%loptic
    call MPI_BCAST(algo%loptic ,1,MPI_LOGICAL,master,MPI_COMM_WORLD,ierror)
    if (ierror /= 0) then
       write (*,*) 'mpi_algobcast returns error/5'
       STOP
    endif
    !write(100+myid,*)'algo%lBfield',algo%lBfield
    call MPI_BCAST(algo%lBfield,1,MPI_LOGICAL,master,MPI_COMM_WORLD,ierror)
    if (ierror /= 0) then
       write (*,*) 'mpi_algobcast returns error/6'
       STOP
    endif
    !write(100+myid,*)'algo%lsymm',algo%lsymm
    call MPI_BCAST(algo%lsymm  ,1,MPI_LOGICAL,master,MPI_COMM_WORLD,ierror)
    if (ierror /= 0) then
       write (*,*) 'mpi_algobcast returns error/7'
       STOP
    endif
    !write(100+myid,*)'algo%imurestart',algo%imurestart
    call MPI_BCAST(algo%imurestart ,1,MPI_INTEGER,master,MPI_COMM_WORLD,ierror)
    if (ierror /= 0) then
       write (*,*) 'mpi_algobcast returns error/8'
       STOP
    endif

  end subroutine mpi_algobcast

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  SUBROUTINE MPI_KMESHBCAST
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine mpi_kmeshbcast(kmesh)
    use Mtypes
    implicit none
    type (kpointmesh) :: kmesh
    integer :: ierror

    !write(200+myid,*)'kmesh%k_id(3,3,3)',kmesh%k_id(3,3,3)
    call MPI_BCAST(kmesh%k_id(1,1,1),size(kmesh%k_id),MPI_INTEGER,master,MPI_COMM_WORLD,ierror)
    if (ierror /= 0) then
       write (*,*) 'mpi_kmeshbcast failed 1'
       STOP
    endif
    !write(200+myid,*)'kmesh%ktot',kmesh%ktot,kmesh%kx,kmesh%ky,kmesh%kz
    call MPI_BCAST(kmesh%ktot,1,MPI_INTEGER,master,MPI_COMM_WORLD,ierror)
    if (ierror /= 0) then
       write (*,*) 'mpi_kmeshbcast failed 2'
       STOP
    endif
    call MPI_BCAST(kmesh%kx,1,MPI_INTEGER,master,MPI_COMM_WORLD,ierror)
    if (ierror /= 0) then
       write (*,*) 'mpi_kmeshbcast failed 3'
       STOP
    endif
    call MPI_BCAST(kmesh%ky,1,MPI_INTEGER,master,MPI_COMM_WORLD,ierror)
    if (ierror /= 0) then
       write (*,*) 'mpi_kmeshbcast failed 4'
       STOP
    endif
    call MPI_BCAST(kmesh%kz,1,MPI_INTEGER,master,MPI_COMM_WORLD,ierror)
    if (ierror /= 0) then
       write (*,*) 'mpi_kmeshbcast failed 5'
       STOP
    endif
    !write(200+myid,*)'kmesh%alat',kmesh%alat,kmesh%a(:)
    call MPI_BCAST(kmesh%alat,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierror)
    if (ierror /= 0) then
       write (*,*) 'mpi_kmeshbcast failed 6'
       STOP
    endif
    call MPI_BCAST(kmesh%a(1),3,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierror)
    if (ierror /= 0) then
       write (*,*) 'mpi_kmeshbcast failed 7'
       STOP
    endif

  end subroutine mpi_kmeshbcast

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  SUBROUTINE MPI_EDISPBCAST
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine mpi_edispbcast(ek)
    use Mtypes
    implicit none
    type (edisp) :: ek
    integer :: ierror

    !write(300+myid,*)'ek%nband_max',ek%nband_max
    call MPI_BCAST(ek%nband_max,1,MPI_INTEGER,master,MPI_COMM_WORLD,ierror)
    if (ierror /= 0) then
       write (*,*) 'mpi_edispbcast failed 1'
       STOP
    endif
    !write(300+myid,*)'ek%nbopt_min',ek%nbopt_min,'ek%nbopt_max',ek%nbopt_max
    call MPI_BCAST(ek%nbopt_min,1,MPI_INTEGER,master,MPI_COMM_WORLD,ierror)
    if (ierror /= 0) then
       write (*,*) 'mpi_edispbcast failed 2'
       STOP
    endif
    call MPI_BCAST(ek%nbopt_max,1,MPI_INTEGER,master,MPI_COMM_WORLD,ierror)
    if (ierror /= 0) then
       write (*,*) 'mpi_edispbcast failed 3'
       STOP
    endif
    !write(300+myid,*)'ek%efer',ek%efer,'ek%nelect',ek%nelect
    call MPI_BCAST(ek%efer,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierror)
    if (ierror /= 0) then
       write (*,*) 'mpi_edispbcast failed 4'
       STOP
    endif
    !call MPI_BCAST(ek%nelect,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierror)
    if (ierror /= 0) then
       write (*,*) 'mpi_edispbcast failed 5'
       STOP
    endif
    !write(300+myid,*)'ek%band(1,1)',ek%band(1,1),'ek%band(1,10)',ek%band(1,10)
    call MPI_BCAST(ek%band(1,1),size(ek%band),MPI_REAL8,master,MPI_COMM_WORLD,ierror)
    if (ierror /= 0) then
       write (*,*) 'mpi_edispbcast failed 6'
       STOP
    endif
    !write(300+myid,*)'ek%Mopt(1,1,1,1)',ek%Mopt(1,1,1,1),'ek%Mopt(1,2,1,10)',ek%Mopt(1,2,1,10)
    call MPI_BCAST(ek%Mopt(1,1,1,1),size(ek%Mopt),MPI_REAL8,master,MPI_COMM_WORLD,ierror)
    if (ierror /= 0) then
       write (*,*) 'mpi_edispbcast failed 7'
       STOP
    endif

  end subroutine mpi_edispbcast

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  SUBROUTINE MPI_THDRBCAST
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine mpi_thdrbcast(thdr)
    use Mtypes
    implicit none
    type (tetramesh) :: thdr
    integer :: ierror

    call MPI_BCAST(thdr%ntet,1,MPI_INTEGER,master,MPI_COMM_WORLD,ierror)
    if (ierror /= 0) then
       write (*,*) 'mpi_thdrbcast failed 1'
       STOP
    endif
    call MPI_BCAST(thdr%idtet(1,1),size(thdr%idtet),MPI_INTEGER,master,MPI_COMM_WORLD,ierror)
    if (ierror /= 0) then
       write (*,*) 'mpi_thdrbcast failed 2'
       STOP
    endif
    call MPI_BCAST(thdr%vltot,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierror)
    if (ierror /= 0) then
       write (*,*) 'mpi_thdrbcast failed 3'
       STOP
    endif
    call MPI_BCAST(thdr%vltet,size(thdr%vltet),MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierror)
    if (ierror /= 0) then
       write (*,*) 'mpi_thdrbcast failed 4'
       STOP
    endif

  end subroutine mpi_thdrbcast

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  SUBROUTINE MPI_SCTRBCAST
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine mpi_sctrbcast(sct)
    use Mtypes
    implicit none
    type (scatrate) :: sct
    integer :: ierror, n, m

    ierror=0
    n=size(sct%TT)
    m=size(sct%gam)
    !write(400+myid,'(F8.5/)') sct%TT(:)
    !write(410+myid,'(F8.5/)') sct%gam(:,1)
    call MPI_BCAST(sct%nT,1,MPI_INTEGER,master,MPI_COMM_WORLD,ierror)
    call MPI_BCAST(sct%TT(1),n,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierror)
    if (ierror /= 0) then
       write (*,*) 'mpi_sctrbcast failed 1'
       STOP
    endif
    !call MPI_BCAST(sct%gam(1,1),m,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierror)
    if (ierror /= 0) then
       write (*,*) 'mpi_sctrbcast failed 2'
       STOP
    endif

  end subroutine mpi_sctrbcast

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  SUBROUTINE MPI_REDUCE_QUAD
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine mpi_reduce_quad(qp0,qp1)

    implicit none

    !input/output
    real(16) :: qp0,qp1
    !parameters
    real(16) :: QCUT
    parameter(QCUT=1.Q14)  ! relevant digits  14?
    !tmps
    real(16) :: qtmp,qtmp2,qtmp3
    !Transferrables
    real(8) :: DP0
    integer(8) :: INT0,IEXP
    !Transferred into
    real(8) :: DP0_T(nproc)
    integer(8) :: INT0_T(nproc),IEXP_T(nproc)
    !mpi
    integer :: iproc

    ! TEST NUMBER
    !qp0=1.23456789123456789123456q-3
!    write(*,*)'I ', myid,qp0

    !Get integer exponent  IEXP of qp0
    if (qp0.ne.0.q0) then
       IEXP=int(log10(abs(qp0)),8)
    else
       IEXP=0
    endif

!    write(*,*)'IEXP ',IEXP

    !Slide qp0 by IEXP digits --> leading digit
    !Then slide by number of DP-relevant digits via QCUT (usually QCUT=1.Q12)
    qtmp=(qp0/(10.q0**iEXP))*QCUT
!    write(*,*)'SLIDED  ',qtmp

    !Define integer INT0 that contains all DP-relevant digits, and the <0 difference, cast into DP0
    INT0=int(qtmp,8)
    qtmp2=qtmp-real(INT0,16)
    DP0=real(qtmp2,8)

!    write(*,*)'INT  ',INT0
!    write(*,*)'DIFF ',qtmp2
!    write(*,*)'DIFF-DP ',DP0

    ! transfer IEXP, INT0 and DP0
    ! then lrecombine

!    write(*,*)' Transfer : ', iexp,int0,dp0


!    write(*,*)'============================='

!    call MPI_BARRIER( MPI_COMM_WORLD, mpierr )
!    write(*,*)'B ', myid,INT0
!    write(*,*)'B ', myid,IEXP
!    write(*,*)'B ', myid,DP0
!    write(*,*)

!    call MPI_BARRIER( MPI_COMM_WORLD, mpierr ) ! do I need this here...?
!    if (myid.eq.master) then
!       write(*,*) 'I     ',qp0
!       write(*,*) 'INT0  ',INT0
!       write(*,*) 'IEXP0 ',IEXP
!       write(*,*) 'DP0   ',DP0
!    endif

! GATHER INT0 ... could do gather instead of allgather... but dont wanna recommunicate the result...
    call MPI_BARRIER( MPI_COMM_WORLD, mpierr )
    call MPI_ALLGATHER(INT0,1,MPI_INTEGER8,INT0_T,1,MPI_INTEGER8,MPI_COMM_WORLD,mpierr)

! GATHER IEXP
    call MPI_BARRIER( MPI_COMM_WORLD, mpierr )
    call MPI_ALLGATHER(IEXP,1,MPI_INTEGER8,IEXP_T,1,MPI_INTEGER8,MPI_COMM_WORLD,mpierr)

! GATHER DP0
    call MPI_BARRIER( MPI_COMM_WORLD, mpierr )
    call MPI_ALLGATHER(DP0,1,MPI_DOUBLE_PRECISION,DP0_T,1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,mpierr)

!    call MPI_BARRIER( MPI_COMM_WORLD, mpierr )
 !   call MPI_BARRIER( MPI_COMM_WORLD, mpierr ) ! do I need this here...?
 !   if (myid.eq.master) then
 !      write(*,*) 'INT0_T  ',INT0_T
 !      write(*,*) 'IEXP0_T ',IEXP_T
 !      write(*,*) 'DP0_T   ',DP0_T
 !      write(*,*)
 !   endif


!    call MPI_REDUCE(dp1,dp1_tot,1,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_WORLD,mpierr
!    call MPI_BARRIER( MPI_COMM_WORLD, mpierr ) ! do I need this here...?
!    call MPI_REDUCE(dp2,dp2_tot,1,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_WORLD,mpierr)
!    call MPI_BARRIER( MPI_COMM_WORLD, mpierr ) ! do I need this here...?

!    write(*,*)'============================='


    qp1=0.q0
    do iproc=1,nproc

       qtmp3=real(INT0_T(iproc),16)
       qtmp3=qtmp3*(10.q0**iEXP_T(iproc) / QCUT )
       qtmp3=qtmp3+real(DP0_T(iproc),16)*(10.q0**iEXP_T(iproc) / QCUT )

       qp1=qp1+qtmp3

    enddo

!    call MPI_BARRIER( MPI_COMM_WORLD, mpierr ) ! do I need this here...?
!    if (myid.eq.master) then
!    write(*,*)'F ',qp1
!    write(*,*)'DIFF ', qp0-qp1

!    open(99,file='debug_'//trim(chmyid),status='unknown',position='append')
!    write(99,'(100E20.12)')qp0-qtrans(myid+1)
!    close(99)

!       do iproc=1,nproc
!          write(*,*) 'F ',qtrans(iproc)
!       enddo

!    write(*,*)
!    write(*,*)
!    write(*,*)
!    endif

!    pause

    return
  end subroutine mpi_reduce_quad


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  SUBROUTINE MPI_BCAST_QUAD_OLD
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine mpi_bcast_quad_old(qp0)!,qp1)
    implicit none
    real(16) :: qptmp,qp0!,qp1
    real(8) :: dp1,dp2,cutQ
    parameter(cutQ=1Q12)

    qptmp=real(int(qp0*cutQ,8),16)/cutQ
    dp1=real(qptmp,8)
!    dp2=real(qp0-qptmp,8)
    dp2=real(qp0-real(dp1,16),8)

    call MPI_BCAST(dp1,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpierr)
    call MPI_BARRIER( MPI_COMM_WORLD, mpierr ) ! do I need this here...?
    call MPI_BCAST(dp2,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpierr)
    call MPI_BARRIER( MPI_COMM_WORLD, mpierr ) ! do I need this here...?

    qp0=real(dp1,16)+real(dp2,16)

    return
  end subroutine mpi_bcast_quad_old


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  SUBROUTINE MPI_BCAST_QUAD
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine mpi_bcast_quad(qp0)
    implicit none

    !input/output
    real(16) :: qp0
    !parameters
    real(16) :: QCUT
    parameter(QCUT=1.Q14)  ! relevant digits
    !tmps
    real(16) :: qtmp,qtmp2,qp1
    !Transferrables
    real(8) :: DP0
    integer(8) :: INT0,IEXP

    if (myid.eq.master) then
!       write(*,*)'I ',qp0

       !Get integer exponent  IEXP of qp0
       !careful here... if qp0=0 this will overflow...
       if (qp0.ne.0.q0) then
          IEXP=int(log10(abs(qp0)),8)
       else
          IEXP=0
       endif
!       write(*,*)'IEXP ',IEXP

       !Slide qp0 by IEXP digits --> leading digit
       !Then slide by number of DP-relevant digits via QCUT (usually QCUT=1.Q12)
       qtmp=(qp0/(10.q0**iEXP))*QCUT
!       write(*,*)'SLIDED  ',qtmp


       !Define integer INT0 that contains all DP-relevant digits, and the <0 difference, cast into DP0
       INT0=int(qtmp,8)
       qtmp2=qtmp-real(INT0,16)
       DP0=real(qtmp2,8)

!       write(*,*)'INT  ',INT0
!       write(*,*)'DIFF ',qtmp2
!       write(*,*)'DIFF-DP ',DP0

    endif !master

! transfer IEXP, INT0 and DP0
! then recombine

    call MPI_BARRIER( MPI_COMM_WORLD, mpierr )
    call MPI_BCAST(dp0,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpierr)
    call MPI_BARRIER( MPI_COMM_WORLD, mpierr )
    call MPI_BCAST(INT0,1,MPI_INTEGER8,master,MPI_COMM_WORLD,mpierr)
    call MPI_BARRIER( MPI_COMM_WORLD, mpierr )
    call MPI_BCAST(IEXP,1,MPI_INTEGER8,master,MPI_COMM_WORLD,mpierr)

!    if (myid.eq.master) then
!       write(*,*)' Transfer : ', iexp,int0,dp0
!       write(*,*)'============================='
!    endif

    qp1=real(INT0,16)
!    write(*,*)qp0
    qp1=qp1*(10.q0**iEXP / QCUT )
!    write(*,*)qp0

    qp1=qp1+real(DP0,16)*(10.q0**iEXP / QCUT )

!    write(*,*)'E ',myid, qp1
!    if (myid.eq.master) write(*,*)'DIFF ', qp0-qp1
!    write(*,*)

!pause

!    call MPI_BCAST(dp1,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpierr)
!    call MPI_BARRIER( MPI_COMM_WORLD, mpierr ) ! do I need this here...?
!    call MPI_BCAST(dp2,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpierr)
!    call MPI_BARRIER( MPI_COMM_WORLD, mpierr ) ! do I need this here...?


    qp0=qp1

    return
  end subroutine mpi_bcast_quad

! mpi communication tools
#endif

end module Mmpi_org
