module Mresponse
  use Mmpi_org
  use Mtypes
  use Mparams
  use Mfermi
  implicit none

  interface allocate_response
    module procedure dpresp_alloc, qpresp_alloc
  end interface

  interface calc_polygamma
    module procedure calc_polygamma_D, calc_polygamma_Q
  end interface calc_polygamma

  contains

subroutine calc_response(PolyGamma, mu, edisp, sct, kmesh, algo, info, dresp, respBl, dinter, qresp)
  implicit none
  real(8), intent(in) :: mu

  type(energydisp)    :: edisp
  type(scattering)    :: sct
  type(kpointmesh)    :: kmesh
  type(algorithm)     :: algo
  type(runinfo)       :: info

  type(response_dp) :: dresp   !intraband response
  type(response_dp) :: dinter  !interband response
  type(response_dp) :: respBl  !intrabond response in the Boltzman regime
  type(response_qp) :: qresp   !intraband response in quadruple precision

  complex(8), intent(in) :: PolyGamma(3,edisp%nbopt_min:edisp%nbopt_max, ikstr:ikend, edisp%ispin)

  algo%lInterbandQuantities = .true.

  call response_intra_km(dresp,  PolyGamma, mu, edisp, sct, kmesh, algo, info)
  if (algo%lInterbandquantities) then
    call response_inter_km(dinter, PolyGamma, mu, edisp, sct, kmesh, algo, info)
  endif
  call response_intra_Boltzmann_km(respBl, mu, edisp, sct, kmesh, algo, info)

end subroutine calc_response

subroutine initresp(algo, dresp)
  implicit none
  type(algorithm)   :: algo
  type(response_dp) :: dresp

  dresp%s_full = 0.d0
  dresp%a_full = 0.d0
  dresp%s_sum = 0.d0
  dresp%a_sum = 0.d0

  if (algo%lBfield) then
     dresp%sB_full = 0.d0
     dresp%aB_full = 0.d0
     dresp%sB_sum = 0.d0
     dresp%aB_sum = 0.d0
  endif
end subroutine initresp

subroutine initresp_qp (algo, qresp)
  implicit none
  type(algorithm)   :: algo
  type(response_qp) :: qresp

  qresp%s_full = 0.q0
  qresp%a_full = 0.q0
  qresp%s_sum = 0.q0
  qresp%a_sum = 0.q0

  if (algo%lBfield) then
     qresp%sB_full = 0.q0
     qresp%aB_full = 0.q0
     qresp%sB_sum = 0.q0
     qresp%aB_sum = 0.q0
  endif
end subroutine initresp_qp

subroutine response_intra_km(resp, PolyGamma, mu, edisp, sct, kmesh, algo, info)
  implicit none
  real(8), intent(in) :: mu
  type (response_dp)      :: resp

  type(energydisp)    :: edisp
  type(scattering)    :: sct
  type(kpointmesh)    :: kmesh
  type(algorithm)     :: algo
  type(runinfo)       :: info

  complex(8)          :: PolyGamma(3,edisp%nbopt_min:edisp%nbopt_max, ikstr:ikend, edisp%ispin)

  real(8) :: zqp
  real(8) :: gam
  real(8), allocatable :: enrgy(:,:)

  integer :: i,j
  integer :: iband


  allocate(enrgy(edisp%nbopt_min:edisp%nbopt_max,edisp%ispin))
  ! first we write the kernel into the 1 1 component
  if (algo%lScatteringFile) then
    enrgy = sct%zqp(:,info%ik,:) * (edisp%band(edisp%nbopt_min:edisp%nbopt_max,info%ik,:) - mu)

    resp%s_full(1,1,:,:,info%ik) = real(PolyGamma(1,:,info%ik,:)) &
                                    - info%beta2p*sct%gam(:,info%ik,:)*real(PolyGamma(2,:,info%ik,:))
    resp%s_full(1,1,:,:,info%ik) = resp%s_full(1,1,:,:,info%ik) &
                                    * sct%zqp(:,info%ik,:)**2 * info%beta &
                                    / (4.d0 * pi**3 * sct%gam(:,info%ik,:))

    resp%a_full(1,1,:,:,info%ik) = enrgy * real(PolyGamma(1,:,info%ik,:)) &
                                      - enrgy * sct%gam(:,info%ik,:) * info%beta2p * real(PolyGamma(2,:,info%ik,:)) &
                                      - sct%gam(:,info%ik,:)**2 * info%beta2p * aimag(PolyGamma(2,:,info%ik,:))

    resp%a_full(1,1,:,:,info%ik) = resp%a_full(1,1,:,:,info%ik) &
                                    * sct%zqp(:,info%ik,:)**2 * info%beta**2 &
                                    / (4.d0 * pi**3 * sct%gam(:,info%ik,:))

    if (algo%lBfield) then

      resp%sB_full(1,1,:,:,info%ik) = -sct%zqp(:,info%ik,:)**2 * info%beta**2 / (4.d0 * pi**2) &
                                      * real(PolyGamma(3,:,info%ik,:)) &
                                      + 3.d0 * sct%gam(:,info%ik,:) * info%beta2p &
                                      * real(PolyGamma(2,:,info%ik,:)) &
                                      - 3.d0 * real(PolyGamma(1,:,info%ik,:))
      resp%sB_full(1,1,:,:,info%ik) = resp%sB_full(1,1,:,:,info%ik) * &
                                      (-sct%zqp(:,info%ik,:)**3 * info%beta / (16.d0 * pi**4 * sct%gam(:,info%ik,:)))

      resp%aB_full(1,1,:,:,info%ik) = -3.d0 * enrgy * sct%gam(:,info%ik,:) * info%beta2p &
                                      * real(PolyGamma(2,:,info%ik,:)) &
                                      - sct%gam(:,info%ik,:) * info%beta2p &
                                      * aimag(PolyGamma(2,:,info%ik,:)) &
                                      + 3.d0 * enrgy * real(PolyGamma(1,:,info%ik,:))

    endif
  else
    enrgy = (sct%zqpscalar * (edisp%band(edisp%nbopt_min:edisp%nbopt_max,info%ik,:) - mu))

    resp%s_full(1,1,:,:,info%ik) = real(PolyGamma(1,:,info%ik,:)) &
                                    - info%beta2p*sct%gamscalar*real(PolyGamma(2,:,info%ik,:))
    resp%s_full(1,1,:,:,info%ik) = resp%s_full(1,1,:,:,info%ik) &
                                    * sct%zqpscalar**2 * info%beta &
                                    / (4.d0 * pi**3 * sct%gamscalar)

    resp%a_full(1,1,:,:,info%ik) = enrgy * real(PolyGamma(1,:,info%ik,:)) &
                                      - enrgy * sct%gamscalar * info%beta2p * real(PolyGamma(2,:,info%ik,:)) &
                                      - sct%gamscalar**2 * info%beta2p * aimag(PolyGamma(2,:,info%ik,:))

    resp%a_full(1,1,:,:,info%ik) = resp%a_full(1,1,:,:,info%ik) &
                                    * sct%zqpscalar**2 * info%beta**2 &
                                    / (4.d0 * pi**3 * sct%gamscalar)

    if (algo%lBfield) then

      resp%sB_full(1,1,:,:,info%ik) = -sct%zqpscalar**2 * info%beta**2 / (4.d0 * pi**2) &
                                      * real(PolyGamma(3,:,info%ik,:)) &
                                      + 3.d0 * sct%gamscalar * info%beta2p &
                                      * real(PolyGamma(2,:,info%ik,:)) &
                                      - 3.d0 * real(PolyGamma(1,:,info%ik,:))
      resp%sB_full(1,1,:,:,info%ik) = resp%sB_full(1,1,:,:,info%ik) * &
                                      (-sct%zqpscalar**3 * info%beta / (16.d0 * pi**4 * sct%gamscalar))

      resp%aB_full(1,1,:,:,info%ik) = -3.d0 * enrgy * sct%gamscalar * info%beta2p &
                                      * real(PolyGamma(2,:,info%ik,:)) &
                                      - sct%gamscalar * info%beta2p &
                                      * aimag(PolyGamma(2,:,info%ik,:)) &
                                      + 3.d0 * enrgy * real(PolyGamma(1,:,info%ik,:))

    endif
  endif
  deallocate(enrgy)

  call response_optical_weights(resp, edisp, info)
  if (algo%lBfield) then
    call response_peierls_weights(resp, edisp, info)
  endif


end subroutine response_intra_km

subroutine response_inter_km(resp, PolyGamma, mu, edisp, sct, kmesh, algo, info)
  implicit none
  real(8), intent(in) :: mu
  type (response_dp)  :: resp

  type(energydisp)    :: edisp
  type(scattering)    :: sct
  type(kpointmesh)    :: kmesh
  type(algorithm)     :: algo
  type(runinfo)       :: info

  complex(8)          :: PolyGamma(3,edisp%nbopt_min:edisp%nbopt_max, ikstr:ikend, edisp%ispin)

  real(8) :: zqp
  real(8) :: gam
  real(8), allocatable :: enrgy(:,:)
  real(8), allocatable :: enrgydiff(:)

  integer :: i,j
  integer :: iband1, iband2, iband



  allocate(enrgy(edisp%nbopt_min:edisp%nbopt_max,edisp%ispin))
  allocate(enrgydiff(edisp%ispin))

  ! first we write the kernel into the 1 1 component
  if (algo%lScatteringFile) then
    enrgy = sct%zqp(:,info%ik,:) * (edisp%band(edisp%nbopt_min:edisp%nbopt_max,info%ik,:) - mu)

    do iband1 = edisp%nbopt_min, edisp%nbopt_max
      do iband2 = edisp%nbopt_min, edisp%nbopt_max
        if (iband1 == iband2) cycle
        enrgydiff = enrgy(iband1,:) - enrgy(iband2,:)

        resp%s_full(1,1,iband1,:,info%ik) = resp%s_full(1,1,iband1,:,info%ik) &
            + real(PolyGamma(1,iband1,info%ik,:)) &
              * ( enrgydiff**2 + sct%gam(iband2,info%ik,:)**2 - sct%gam(iband1,info%ik,:)**2) &
              / sct%gam(iband1,info%ik,:)

        resp%s_full(1,1,iband1,:,info%ik) = resp%s_full(1,1,iband1,:,info%ik) &
            + aimag(PolyGamma(1,iband1,info%ik,:)) &
              * ( 2.d0* (-enrgydiff) )

        resp%s_full(1,1,iband1,:,info%ik) = resp%s_full(1,1,iband1,:,info%ik) &
            + real(PolyGamma(1,iband2,info%ik,:)) &
              * ( enrgydiff**2 + sct%gam(iband1,info%ik,:)**2 - sct%gam(iband2,info%ik,:)**2) &
              / sct%gam(iband2,info%ik,:)

        resp%s_full(1,1,iband1,:,info%ik) = resp%s_full(1,1,iband1,:,info%ik) &
            + aimag(PolyGamma(1,iband1,info%ik,:)) &
              * ( 2.d0* enrgydiff )

        resp%s_full(1,1,iband1,:,info%ik) = resp%s_full(1,1,iband1,:,info%ik) &
            * sct%gam(iband1,info%ik,:) * sct%gam(iband2,info%ik,:) &
            * sct%zqp(iband1,info%ik,:) * sct%zqp(iband2,info%ik,:) &
            * info%beta

        resp%s_full(1,1,iband1,:,info%ik) = resp%s_full(1,1,iband1,:,info%ik) &
            / (2.d0 * pi**2 * ( enrgydiff**2 + (sct%gam(iband1,info%ik,:) - sct%gam(iband2,info%ik,:))**2)) &
            / ( enrgydiff**2 * (sct%gam(iband1,info%ik,:) + sct%gam(iband2,info%ik,:))**2)


        resp%a_full(1,1,iband1,:,info%ik) = resp%a_full(1,1,iband1,:,info%ik) &
            + real(PolyGamma(1,iband1,info%ik,:)) &
              * (enrgy(iband1,:) * (enrgydiff**2 + sct%gam(iband2,info%ik,:)**2 - sct%gam(iband1,info%ik,:)**2) &
              / sct%gam(iband1,info%ik,:) + 2.d0*sct%gam(iband1,info%ik,:)*enrgydiff)

        resp%a_full(1,1,iband1,:,info%ik) = resp%a_full(1,1,iband1,:,info%ik) &
            + aimag(PolyGamma(1,iband1,info%ik,:)) &
              * (enrgydiff**2 + sct%gam(iband2,info%ik,:)**2 - sct%gam(iband1,info%ik,:)**2 &
              - 2.d0*enrgy(iband1,:)*enrgydiff)

        resp%a_full(1,1,iband1,:,info%ik) = resp%a_full(1,1,iband1,:,info%ik) &
            + real(PolyGamma(1,iband2,info%ik,:)) &
              * (enrgy(iband2,:) * (enrgydiff**2 + sct%gam(iband1,info%ik,:)**2 - sct%gam(iband2,info%ik,:)**2) &
              / sct%gam(iband2,info%ik,:) + 2.d0*sct%gam(iband2,info%ik,:)*enrgydiff*(-1.d0))

        resp%a_full(1,1,iband1,:,info%ik) = resp%a_full(1,1,iband1,:,info%ik) &
            + aimag(PolyGamma(1,iband2,info%ik,:)) &
              * (enrgydiff**2 + sct%gam(iband1,info%ik,:)**2 - sct%gam(iband2,info%ik,:)**2 &
              - 2.d0*enrgy(iband2,:)*enrgydiff*(-1.d0))

        resp%a_full(1,1,iband1,:,info%ik) = resp%a_full(1,1,iband1,:,info%ik) &
            * sct%gam(iband1,info%ik,:) * sct%gam(iband2,info%ik,:) &
            * sct%zqp(iband1,info%ik,:) * sct%zqp(iband2,info%ik,:) &
            * info%beta

        resp%a_full(1,1,iband1,:,info%ik) = resp%a_full(1,1,iband1,:,info%ik) &
            / (2.d0 * pi**3 * ( enrgydiff**2 + (sct%gam(iband1,info%ik,:) - sct%gam(iband2,info%ik,:))**2)) &
            / ( enrgydiff**2 * (sct%gam(iband1,info%ik,:) + sct%gam(iband2,info%ik,:))**2)
      enddo
    enddo
  else
    enrgy = sct%zqpscalar * (edisp%band(edisp%nbopt_min:edisp%nbopt_max,info%ik,:) - mu)

    do iband1 = edisp%nbopt_min, edisp%nbopt_max
      do iband2 = edisp%nbopt_min, edisp%nbopt_max
        if (iband1 == iband2) cycle
        enrgydiff = enrgy(iband1,:) - enrgy(iband2,:)

        resp%s_full(1,1,iband1,:,info%ik) = resp%s_full(1,1,iband1,:,info%ik) &
            + real(PolyGamma(1,iband1,info%ik,:)) &
              * enrgydiff**2  / sct%gamscalar

        resp%s_full(1,1,iband1,:,info%ik) = resp%s_full(1,1,iband1,:,info%ik) &
            + aimag(PolyGamma(1,iband1,info%ik,:)) &
              * ( 2.d0* (-enrgydiff) )

        resp%s_full(1,1,iband1,:,info%ik) = resp%s_full(1,1,iband1,:,info%ik) &
            + real(PolyGamma(1,iband2,info%ik,:)) &
              * enrgydiff**2 / sct%gamscalar

        resp%s_full(1,1,iband1,:,info%ik) = resp%s_full(1,1,iband1,:,info%ik) &
            + aimag(PolyGamma(1,iband1,info%ik,:)) &
              * ( 2.d0* enrgydiff )

        resp%s_full(1,1,iband1,:,info%ik) = resp%s_full(1,1,iband1,:,info%ik) &
            * sct%gamscalar**2 &
            * sct%zqpscalar**2 &
            * info%beta

        resp%s_full(1,1,iband1,:,info%ik) = resp%s_full(1,1,iband1,:,info%ik) &
            / (2.d0 * pi**2 * enrgydiff**2) &
            / ( enrgydiff**2 * (2.d0 * sct%gamscalar)**2)


        resp%a_full(1,1,iband1,:,info%ik) = resp%a_full(1,1,iband1,:,info%ik) &
            + real(PolyGamma(1,iband1,info%ik,:)) &
              * (enrgy(iband1,:) * enrgydiff**2 &
              / sct%gamscalar + 2.d0*sct%gamscalar*enrgydiff)

        resp%a_full(1,1,iband1,:,info%ik) = resp%a_full(1,1,iband1,:,info%ik) &
            + aimag(PolyGamma(1,iband1,info%ik,:)) &
              * (enrgydiff**2 - 2.d0*enrgy(iband1,:)*enrgydiff)

        resp%a_full(1,1,iband1,:,info%ik) = resp%a_full(1,1,iband1,:,info%ik) &
            + real(PolyGamma(1,iband2,info%ik,:)) &
              * (enrgy(iband2,:) * enrgydiff**2 &
              / sct%gamscalar - 2.d0*sct%gamscalar*enrgydiff)

        resp%a_full(1,1,iband1,:,info%ik) = resp%a_full(1,1,iband1,:,info%ik) &
            + aimag(PolyGamma(1,iband2,info%ik,:)) &
              * (enrgydiff**2 + 2.d0*enrgy(iband2,:)*enrgydiff)

        resp%a_full(1,1,iband1,:,info%ik) = resp%a_full(1,1,iband1,:,info%ik) &
            * sct%gamscalar**2 &
            * sct%zqpscalar**2 * sct%zqp(iband2,info%ik,:) &
            * info%beta

        resp%a_full(1,1,iband1,:,info%ik) = resp%a_full(1,1,iband1,:,info%ik) &
            / (2.d0 * pi**3 * enrgydiff**2) &
            / ( enrgydiff**2 * (2.d0 * sct%gamscalar)**2)
      enddo
    enddo


  endif
  deallocate(enrgy)
  deallocate(enrgydiff)

  call response_optical_weights(resp, edisp, info)

end subroutine response_inter_km

subroutine response_intra_Boltzmann_km(resp, mu, edisp, sct, kmesh, algo, info)
  implicit none
  real(8), intent(in) :: mu
  type (response_dp)      :: resp

  type(energydisp)    :: edisp
  type(scattering)    :: sct
  type(kpointmesh)    :: kmesh
  type(algorithm)     :: algo
  type(runinfo)       :: info

  real(8) :: zqp
  real(8) :: gam
  real(8), allocatable :: enrgy(:,:)

  integer :: i,j
  integer :: iband


  allocate(enrgy(edisp%nbopt_min:edisp%nbopt_max,edisp%ispin))
  ! first we write the kernel into the 1 1 component
  if (algo%lScatteringFile) then
    enrgy = sct%zqp(:,info%ik,:) * (edisp%band(edisp%nbopt_min:edisp%nbopt_max,info%ik,:) - mu)

    ! asymptotic term in the Gamma-> 0 limit
    resp%s_full(1,1,:,:,info%ik) = dfermi(enrgy,info%beta) &
                                    * sct%zqp(:,info%ik,:)**2 * info%beta &
                                    / (4.d0 * pi**3 * sct%gam(:,info%ik,:))

    resp%a_full(1,1,:,:,info%ik) = resp%s_full(1,1,:,:,info%ik) &
                                    * info%beta * enrgy

    if (algo%lBfield) then

      resp%sB_full(1,1,:,:,info%ik) = dfermi(enrgy,info%beta) &
                                      * 3.d0 * sct%zqp(:,info%ik,:)**3 * info%beta &
                                      / (16.d0 * pi**4 * sct%gam(:,info%ik,:)**2)

      resp%a_full(1,1,:,:,info%ik) = resp%sB_full(1,1,:,:,info%ik) &
                                      * info%beta * enrgy

    endif
  else
    enrgy = (sct%zqpscalar * (edisp%band(edisp%nbopt_min:edisp%nbopt_max,info%ik,:) - mu))

    ! asymptotic term in the Gamma-> 0 limit
    resp%s_full(1,1,:,:,info%ik) = dfermi(enrgy,info%beta) &
                                    * sct%zqpscalar**2 * info%beta &
                                    / (4.d0 * pi**3 * sct%gamscalar)

    resp%a_full(1,1,:,:,info%ik) = resp%s_full(1,1,:,:,info%ik) &
                                    * info%beta * enrgy

    if (algo%lBfield) then

      resp%sB_full(1,1,:,:,info%ik) = dfermi(enrgy,info%beta) &
                                      * 3.d0 * sct%zqpscalar**3 * info%beta &
                                      / (16.d0 * pi**4 * sct%gamscalar**2)

      resp%a_full(1,1,:,:,info%ik) = resp%sB_full(1,1,:,:,info%ik) &
                                      * info%beta * enrgy

    endif
  endif

  deallocate(enrgy)

  call response_optical_weights(resp, edisp, info)
  if (algo%lBfield) then
    call response_peierls_weights(resp, edisp, info)
  endif

end subroutine response_intra_Boltzmann_km


! multiply optical elements onto quantities without B-Field
subroutine response_optical_weights(resp, edisp, info)
  implicit none
  type (response_dp) :: resp

  type(energydisp)   :: edisp
  type(runinfo)      :: info

  integer :: index1(9), index2(9)
  integer :: iband, idir

  index1 = (/1,2,3,1,1,2,1,1,2/)
  index2 = (/1,2,3,2,3,3,2,3,3/)

  !( 1 4+i7 5+i8 )
  !( - 2    6+i9 )
  !( - -    3    )

  do iband = edisp%nbopt_min, edisp%nbopt_max
    ! the kernels are saved in the 1 1 directions
    ! so we calculate the 1 1 component at the end
    do idir = 2,edisp%iOptical
      if (idir <= 6) then
        resp%s_full(index1(idir),index2(idir),iband,:,info%ik) = resp%s_full(1,1,iband,:,info%ik) &
                                                               * edisp%Mopt(idir,iband,iband,:)
        resp%a_full(index1(idir),index2(idir),iband,:,info%ik) = resp%a_full(1,1,iband,:,info%ik) &
                                                               * edisp%Mopt(idir,iband,iband,:)
      else
        ! here we ADD the complex part to the response
        resp%s_full(index1(idir),index2(idir),iband,:,info%ik) = resp%s_full(index1(idir),index2(idir),iband,:,info%ik) &
                                                               * edisp%Mopt(idir,iband,iband,:) * ci
        resp%a_full(index1(idir),index2(idir),iband,:,info%ik) = resp%a_full(index1(idir),index2(idir),iband,:,info%ik) &
                                                               * edisp%Mopt(idir,iband,iband,:) * ci
      endif
    enddo

    resp%s_full(1,1,iband,:,info%ik) = resp%s_full(1,1,iband,:,info%ik) &
                                       * edisp%Mopt(1,iband,iband,:)
    resp%a_full(1,1,iband,:,info%ik) = resp%a_full(1,1,iband,:,info%ik) &
                                       * edisp%Mopt(1,iband,iband,:)

    ! symmetrize
    do idir=6,9
      resp%s_full(index2(idir),index1(idir),iband,:,info%ik) = conjg(resp%s_full(index1(idir),index2(idir),iband,:,info%ik))
    enddo
  enddo
end subroutine response_optical_weights

subroutine response_peierls_weights(resp, edisp, info)
  implicit none
  type (response_dp) :: resp

  type(energydisp)   :: edisp
  type(runinfo)      :: info

  integer :: index1(9), index2(9)
  integer :: iband, idir

  index1 = (/1,2,3,1,1,2,1,1,2/)
  index2 = (/1,2,3,2,3,3,2,3,3/)

  ! for the time being
  ! this is only with Bfield in z-direction
  ! TODO: complete ... arbitrary direction

  do iband = edisp%nbopt_min, edisp%nbopt_max
    resp%sB_full(1,2,iband,:,info%ik) = resp%sB_full(1,1,iband,:,info%ik) &
              * (edisp%band_dk(1,iband,info%ik,:)*edisp%band_dk(2,iband,info%ik,:) &
                 *edisp%band_d2k(4,iband,info%ik,:) &
                 -edisp%band_dk(1,iband,info%ik,:)*edisp%band_dk(1,iband,info%ik,:) &
                 *edisp%band_d2k(2,iband,info%ik,:))

    resp%sB_full(2,1,iband,:,info%ik) = resp%sB_full(1,1,iband,:,info%ik) &
              * (edisp%band_dk(2,iband,info%ik,:)*edisp%band_dk(1,iband,info%ik,:) &
                 *edisp%band_d2k(4,iband,info%ik,:) &
                 -edisp%band_dk(2,iband,info%ik,:)*edisp%band_dk(2,iband,info%ik,:) &
                 *edisp%band_d2k(1,iband,info%ik,:))

    resp%sB_full(1,1,iband,:,info%ik) = resp%sB_full(1,1,iband,:,info%ik) &
              * (edisp%band_dk(1,iband,info%ik,:)*edisp%band_dk(1,iband,info%ik,:) &
                 *edisp%band_d2k(4,iband,info%ik,:) &
                 -edisp%band_dk(1,iband,info%ik,:)*edisp%band_dk(2,iband,info%ik,:) &
                 *edisp%band_d2k(1,iband,info%ik,:))


    resp%aB_full(1,2,iband,:,info%ik) = resp%aB_full(1,1,iband,:,info%ik) &
              * (edisp%band_dk(1,iband,info%ik,:)*edisp%band_dk(2,iband,info%ik,:) &
                 *edisp%band_d2k(4,iband,info%ik,:) &
                 -edisp%band_dk(1,iband,info%ik,:)*edisp%band_dk(1,iband,info%ik,:) &
                 *edisp%band_d2k(2,iband,info%ik,:))

    resp%aB_full(2,1,iband,:,info%ik) = resp%aB_full(1,1,iband,:,info%ik) &
              * (edisp%band_dk(2,iband,info%ik,:)*edisp%band_dk(1,iband,info%ik,:) &
                 *edisp%band_d2k(4,iband,info%ik,:) &
                 -edisp%band_dk(2,iband,info%ik,:)*edisp%band_dk(2,iband,info%ik,:) &
                 *edisp%band_d2k(1,iband,info%ik,:))

    resp%aB_full(1,1,iband,:,info%ik) = resp%aB_full(1,1,iband,:,info%ik) &
              * (edisp%band_dk(1,iband,info%ik,:)*edisp%band_dk(1,iband,info%ik,:) &
                 *edisp%band_d2k(4,iband,info%ik,:) &
                 -edisp%band_dk(1,iband,info%ik,:)*edisp%band_dk(2,iband,info%ik,:) &
                 *edisp%band_d2k(1,iband,info%ik,:))
  enddo
end subroutine

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
!   type(response_qp), optional ::hpresp
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
!   type(response_qp),optional :: hpresp
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
!   type(response_qp),optional :: hpresp
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
                  (sct%gam(:,ikstr:ikend,:) + ci*sct%zqp(:,ikstr:ikend,:)*(edisp%band(:,ikstr:ikend,:) - mu))
  else
    to_evaluate = 0.5d0 + info%beta2p * &
                  (sct%gamscalar + ci*sct%zqpscalar*(edisp%band(:,ikstr:ikend,:) - mu))
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

subroutine calc_total_energy(mu, energy_tot, edisp, sct, kmesh, algo, info)
  real(8), intent(in)  :: mu
  real(8), intent(out) :: energy_tot

  type(energydisp) :: edisp
  type(scattering) :: sct
  type(kpointmesh) :: kmesh
  type(algorithm)  :: algo
  type(runinfo)    :: info
  !local variables

  real(8) :: energy_loc
  integer :: is, ik, iband
  complex(8), allocatable :: to_evaluate(:,:,:)
  real(8), allocatable    :: energy_post_factor(:,:,:)
  real(8), allocatable    :: energy(:,:,:)
  !external variables
  complex(8), external :: wpsipg

  allocate(to_evaluate(edisp%nband_max, ikstr:ikend, edisp%ispin))
  allocate(energy_post_factor(edisp%nband_max, ikstr:ikend, edisp%ispin))
  allocate(energy(edisp%nband_max, ikstr:ikend, edisp%ispin))

  if (algo%lScatteringFile) then
    to_evaluate = 0.5d0 + info%beta2p * &
                  (sct%gam(:,ikstr:ikend,:) - ci*(sct%zqp(:,ikstr:ikend,:)*edisp%band(:,ikstr:ikend,:) - mu))
    energy_post_factor = sct%zqp(:,ikstr:ikend,:) * edisp%band(:,ikstr:ikend,:) - mu
  else
    to_evaluate = 0.5d0 + info%beta2p * &
                  (sct%gamscalar - ci*(sct%zqpscalar*edisp%band(:,ikstr:ikend,:) - mu))
    energy_post_factor = sct%zqpscalar * edisp%band(:,ikstr:ikend,:) - mu
  endif

  ! evaluate the function
  do is = 1,edisp%ispin
    do ik = ikstr, ikend
      do iband=1,edisp%nband_max
        energy(iband,ik,is) = 0.5d0 + aimag(wpsipg(to_evaluate(iband,ik,is),0))/pi ! this is the occupation
        energy(iband,ik,is) = energy(iband,ik,is) * kmesh%weight(ik) * energy_post_factor(iband,ik,is) ! multiplied with weight and energy gives the energy
      enddo
    enddo
  enddo

  deallocate(to_evaluate)
  deallocate(energy_post_factor)
  energy_loc = sum(energy)
  deallocate(energy)

#ifdef MPI
  call MPI_ALLREDUCE(energy_loc, energy_tot, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, mpierr)
#else
  energy_tot = energy_loc
#endif

end subroutine

subroutine intldos(mu, dos, edisp, sct, kmesh, algo, info)
  !passed variables
  real(8), intent(in) :: mu
  type(dosgrid)    :: dos
  type(algorithm)  :: algo
  type(energydisp) :: edisp
  type(kpointmesh) :: kmesh
  type(scattering) :: sct
  type(runinfo)    :: info

  !local variables
  integer :: ee, ik, iband
  real(8) :: eps, eta
  real(8) :: s, dee
  real(8), allocatable :: AA(:)
  complex(8) :: iu, G0

  real(8) :: dw = 1.d-3
  integer :: nw

  nw = int(10.d0/dw)

  allocate(AA(-nw:nw))



  AA = 0.0d0
  do ee=-nw,nw
    do ik=1,kmesh%nkp
       do iband= 1,edisp%nband_max
           eps =  mu - edisp%band(iband,ik,1)
           G0 = sct%zqp(iband,ik,1) / (dw*ee + sct%zqp(iband,ik,1) * (eps + ci*sct%gam(iband,ik,1)))
           ! A = -1/pi  Im G0
           AA(ee) = AA(ee) - aimag(G0)*kmesh%weight(ik)
        enddo
     enddo
  enddo !ee
  AA = AA/pi

  !!use trapezoidal rule to evaluate the number of electrons
  !s=0.5d0*(AA(1)+AA(dos%nnrg))
  !do ee=2,dos%nnrg-1
  !   s=s+AA(ee)
  !enddo

  do ee=-nw,nw
     write(120,'(E12.6,3x,E12.6)') dw*ee,AA(ee)
  enddo

  deallocate(AA)

end subroutine intldos


end module Mresponse
