module Mresponse
  use Mmpi_org
  use Mtypes
  use Mparams
  use Mfermi
  use hdf5
  use hdf5_wrapper

  implicit none

  interface allocate_response
    module procedure dpresp_alloc, qpresp_alloc
  end interface

  interface calc_polygamma
    module procedure calc_polygamma_D, calc_polygamma_Q
  end interface calc_polygamma

  contains

subroutine calc_response(PolyGamma, mu, edisp, sct, kmesh, algo, info, &
           resp_intra, resp_intra_Boltzmann, &
           resp_inter, resp_inter_Boltzmann)
  implicit none
  real(8), intent(in) :: mu

  type(energydisp)    :: edisp
  type(scattering)    :: sct
  type(kpointmesh)    :: kmesh
  type(algorithm)     :: algo
  type(runinfo)       :: info

  type(response_dp) :: resp_intra
  type(response_dp) :: resp_intra_Boltzmann
  type(response_dp) :: resp_inter
  type(response_dp) :: resp_inter_Boltzmann

  complex(8), intent(in) :: PolyGamma(3,edisp%nbopt_min:edisp%nbopt_max, ikstr:ikend, edisp%ispin)

  call response_intra_km(resp_intra,  PolyGamma, mu, edisp, sct, kmesh, algo, info)
  if (algo%lBoltzmann) then
    call response_intra_Boltzmann_km(resp_intra_Boltzmann, mu, edisp, sct, kmesh, algo, info)
  endif

  if (algo%lInterbandquantities) then
    call response_inter_km(resp_inter, PolyGamma, mu, edisp, sct, kmesh, algo, info)
    if (algo%lBoltzmann) then
      call response_inter_Boltzmann_km(resp_inter_Boltzmann, mu, edisp, sct, kmesh, algo, info)
    endif
  endif

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

  call response_intra_optical_weights(resp, edisp, info)
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

  complex(8), allocatable :: calc_cond(:)
  complex(8), allocatable :: calc_seeb(:)

  integer :: index1(9), index2(9)
  integer :: i,j,idir
  integer :: iband1, iband2, iband

  index1 = (/1,2,3,1,1,2,1,1,2/)
  index2 = (/1,2,3,2,3,3,2,3,3/)


  allocate(enrgy(edisp%nbopt_min:edisp%nbopt_max,edisp%ispin))
  allocate(enrgydiff(edisp%ispin))
  allocate(calc_cond(edisp%ispin))
  allocate(calc_seeb(edisp%ispin))

  ! first we write the kernel into the 1 1 component
  if (algo%lScatteringFile) then
    enrgy = sct%zqp(:,info%ik,:) * (edisp%band(edisp%nbopt_min:edisp%nbopt_max,info%ik,:) - mu)

    do iband1 = edisp%nbopt_min, edisp%nbopt_max
      do iband2 = edisp%nbopt_min, edisp%nbopt_max
        if (iband1 == iband2) cycle
        enrgydiff = enrgy(iband1,:) - enrgy(iband2,:)

        calc_cond = 0.d0

        calc_cond = calc_cond &
            + real(PolyGamma(1,iband1,info%ik,:)) &
              * ( enrgydiff**2 + sct%gam(iband2,info%ik,:)**2 - sct%gam(iband1,info%ik,:)**2) &
              / sct%gam(iband1,info%ik,:)

        calc_cond = calc_cond &
            + aimag(PolyGamma(1,iband1,info%ik,:)) &
              * ( 2.d0* (-enrgydiff) )

        calc_cond = calc_cond &
            + real(PolyGamma(1,iband2,info%ik,:)) &
              * ( enrgydiff**2 + sct%gam(iband1,info%ik,:)**2 - sct%gam(iband2,info%ik,:)**2) &
              / sct%gam(iband2,info%ik,:)

        calc_cond = calc_cond &
            + aimag(PolyGamma(1,iband2,info%ik,:)) &
              * ( 2.d0* enrgydiff )

        calc_cond = calc_cond &
            * sct%gam(iband1,info%ik,:) * sct%gam(iband2,info%ik,:) &
            * sct%zqp(iband1,info%ik,:) * sct%zqp(iband2,info%ik,:) &
            * info%beta

        calc_cond = calc_cond &
            / (2.d0 * pi**3 * ( enrgydiff**2 + (sct%gam(iband1,info%ik,:) - sct%gam(iband2,info%ik,:))**2)) &
            / ( enrgydiff**2 + (sct%gam(iband1,info%ik,:) + sct%gam(iband2,info%ik,:))**2)



        calc_seeb = 0.d0

        calc_seeb = calc_seeb &
            + real(PolyGamma(1,iband1,info%ik,:)) &
              * (enrgy(iband1,:) * (enrgydiff**2 + sct%gam(iband2,info%ik,:)**2 - sct%gam(iband1,info%ik,:)**2) &
              / sct%gam(iband1,info%ik,:) + 2.d0*sct%gam(iband1,info%ik,:)*enrgydiff)

        calc_seeb = calc_seeb &
            + aimag(PolyGamma(1,iband1,info%ik,:)) &
              * (enrgydiff**2 + sct%gam(iband2,info%ik,:)**2 - sct%gam(iband1,info%ik,:)**2 &
              - 2.d0*enrgy(iband1,:)*enrgydiff)

        calc_seeb = calc_seeb &
            + real(PolyGamma(1,iband2,info%ik,:)) &
              * (enrgy(iband2,:) * (enrgydiff**2 + sct%gam(iband1,info%ik,:)**2 - sct%gam(iband2,info%ik,:)**2) &
              / sct%gam(iband2,info%ik,:) + 2.d0*sct%gam(iband2,info%ik,:)*enrgydiff*(-1.d0))

        calc_seeb = calc_seeb &
            + aimag(PolyGamma(1,iband2,info%ik,:)) &
              * (enrgydiff**2 + sct%gam(iband1,info%ik,:)**2 - sct%gam(iband2,info%ik,:)**2 &
              - 2.d0*enrgy(iband2,:)*enrgydiff*(-1.d0))

        calc_seeb = calc_seeb &
            * sct%gam(iband1,info%ik,:) * sct%gam(iband2,info%ik,:) &
            * sct%zqp(iband1,info%ik,:) * sct%zqp(iband2,info%ik,:) &
            * info%beta

        calc_seeb = calc_seeb &
            / (2.d0 * pi**3 * ( enrgydiff**2 + (sct%gam(iband1,info%ik,:) - sct%gam(iband2,info%ik,:))**2)) &
            / ( enrgydiff**2 + (sct%gam(iband1,info%ik,:) + sct%gam(iband2,info%ik,:))**2)


        ! multiply optical elements
        do idir = 1,edisp%iOptical
          if (idir <= 6) then
            ! ATTENTION
            ! we read the wien2k files via 1 - 2
            ! we save those in hdf5
            ! and read them via Fortran -> implicit transposition
            ! we have to use 2 - 1 here
            resp%s_full(index1(idir),index2(idir),iband1,:,info%ik) = &
            resp%s_full(index1(idir),index2(idir),iband1,:,info%ik) + calc_cond * edisp%Mopt(idir,iband2,iband1,:)

            resp%a_full(index1(idir),index2(idir),iband1,:,info%ik) = &
            resp%a_full(index1(idir),index2(idir),iband1,:,info%ik) + calc_seeb * edisp%Mopt(idir,iband2,iband1,:)
          else
            ! here we ADD the complex part to the response
            resp%s_full(index1(idir),index2(idir),iband1,:,info%ik) = &
            resp%s_full(index1(idir),index2(idir),iband1,:,info%ik) + calc_cond * edisp%Mopt(idir,iband2,iband1,:) * ci

            resp%a_full(index1(idir),index2(idir),iband1,:,info%ik) = &
            resp%a_full(index1(idir),index2(idir),iband1,:,info%ik) + calc_seeb * edisp%Mopt(idir,iband2,iband1,:) * ci
          endif
        enddo

      enddo
    enddo
  else

    enrgy = sct%zqpscalar * (edisp%band(edisp%nbopt_min:edisp%nbopt_max,info%ik,:) - mu)

    do iband1 = edisp%nbopt_min, edisp%nbopt_max
      do iband2 = edisp%nbopt_min, edisp%nbopt_max
        if (iband1 == iband2) cycle
        enrgydiff = enrgy(iband1,:) - enrgy(iband2,:)

        calc_cond = 0.d0

        calc_cond = calc_cond &
            + real(PolyGamma(1,iband1,info%ik,:)) &
              * enrgydiff**2 &
              / sct%gamscalar

        calc_cond = calc_cond &
            + aimag(PolyGamma(1,iband1,info%ik,:)) &
              * ( 2.d0* (-enrgydiff) )

        calc_cond = calc_cond &
            + real(PolyGamma(1,iband2,info%ik,:)) &
              * enrgydiff**2 &
              / sct%gamscalar

        calc_cond = calc_cond &
            + aimag(PolyGamma(1,iband2,info%ik,:)) &
              * ( 2.d0* enrgydiff )

        calc_cond = calc_cond &
            * sct%gamscalar**2 &
            * sct%zqpscalar**2 &
            * info%beta

        calc_cond = calc_cond &
            / (2.d0 * pi**2 * enrgydiff**2 ) &
            / ( enrgydiff**2 + (2.d0 * sct%gamscalar)**2)



        calc_seeb = 0.d0

        calc_seeb = calc_seeb &
            + real(PolyGamma(1,iband1,info%ik,:)) &
              * (enrgy(iband1,:) * enrgydiff**2 &
              / sct%gamscalar + 2.d0*sct%gamscalar*enrgydiff)

        calc_seeb = calc_seeb &
            + aimag(PolyGamma(1,iband1,info%ik,:)) &
              * (enrgydiff**2 - 2.d0*enrgy(iband1,:)*enrgydiff)

        calc_seeb = calc_seeb &
            + real(PolyGamma(1,iband2,info%ik,:)) &
              * (enrgy(iband2,:) * enrgydiff**2 &
              / sct%gamscalar + 2.d0*sct%gamscalar*enrgydiff*(-1.d0))

        calc_seeb = calc_seeb &
            + aimag(PolyGamma(1,iband2,info%ik,:)) &
              * (enrgydiff**2 - 2.d0*enrgy(iband2,:)*enrgydiff*(-1.d0))

        calc_seeb = calc_seeb &
            * sct%gamscalar**2 &
            * sct%zqpscalar**2 &
            * info%beta

        calc_seeb = calc_seeb &
            / (2.d0 * pi**3 * enrgydiff**2 ) &
            / ( enrgydiff**2 + (2.d0 * sct%gamscalar)**2)


        ! multiply optical elements
        do idir = 1,edisp%iOptical
          if (idir <= 6) then
            ! ATTENTION
            ! we read the wien2k files via 1 - 2
            ! we save those in hdf5
            ! and read them via Fortran -> implicit transposition
            ! we have to use 2 - 1 here
            resp%s_full(index1(idir),index2(idir),iband1,:,info%ik) = &
            resp%s_full(index1(idir),index2(idir),iband1,:,info%ik) + calc_cond * edisp%Mopt(idir,iband2,iband1,:)

            resp%a_full(index1(idir),index2(idir),iband1,:,info%ik) = &
            resp%a_full(index1(idir),index2(idir),iband1,:,info%ik) + calc_seeb * edisp%Mopt(idir,iband2,iband1,:)
          else
            ! here we ADD the complex part to the response
            resp%s_full(index1(idir),index2(idir),iband1,:,info%ik) = &
            resp%s_full(index1(idir),index2(idir),iband1,:,info%ik) + calc_cond * edisp%Mopt(idir,iband2,iband1,:) * ci

            resp%a_full(index1(idir),index2(idir),iband1,:,info%ik) = &
            resp%a_full(index1(idir),index2(idir),iband1,:,info%ik) + calc_seeb * edisp%Mopt(idir,iband2,iband1,:) * ci
          endif
        enddo

      enddo
    enddo

  endif


  deallocate(enrgy)
  deallocate(enrgydiff)
  deallocate(calc_cond, calc_seeb)

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
    resp%s_full(1,1,:,:,info%ik) = polygamma2fermi(enrgy,info%beta) &
                                    * sct%zqp(:,info%ik,:)**2 * info%beta &
                                    / (4.d0 * pi**3 * sct%gam(:,info%ik,:))

    resp%a_full(1,1,:,:,info%ik) = resp%s_full(1,1,:,:,info%ik) &
                                    * info%beta * enrgy

    if (algo%lBfield) then

      resp%sB_full(1,1,:,:,info%ik) = polygamma2fermi(enrgy,info%beta) &
                                      * 3.d0 * sct%zqp(:,info%ik,:)**3 * info%beta &
                                      / (16.d0 * pi**4 * sct%gam(:,info%ik,:)**2)

      resp%a_full(1,1,:,:,info%ik) = resp%sB_full(1,1,:,:,info%ik) &
                                      * info%beta * enrgy

    endif
  else
    enrgy = (sct%zqpscalar * (edisp%band(edisp%nbopt_min:edisp%nbopt_max,info%ik,:) - mu))

    ! asymptotic term in the Gamma-> 0 limit
    resp%s_full(1,1,:,:,info%ik) = polygamma2fermi(enrgy,info%beta) &
                                    * sct%zqpscalar**2 * info%beta &
                                    / (4.d0 * pi**3 * sct%gamscalar)

    resp%a_full(1,1,:,:,info%ik) = resp%s_full(1,1,:,:,info%ik) &
                                    * info%beta * enrgy

    if (algo%lBfield) then

      resp%sB_full(1,1,:,:,info%ik) = polygamma2fermi(enrgy,info%beta) &
                                      * 3.d0 * sct%zqpscalar**3 * info%beta &
                                      / (16.d0 * pi**4 * sct%gamscalar**2)

      resp%a_full(1,1,:,:,info%ik) = resp%sB_full(1,1,:,:,info%ik) &
                                      * info%beta * enrgy

    endif
  endif

  deallocate(enrgy)

  call response_intra_optical_weights(resp, edisp, info)
  if (algo%lBfield) then
    call response_peierls_weights(resp, edisp, info)
  endif

end subroutine response_intra_Boltzmann_km

subroutine response_inter_Boltzmann_km(resp, mu, edisp, sct, kmesh, algo, info)
  implicit none
  real(8), intent(in) :: mu
  type (response_dp)  :: resp

  type(energydisp)    :: edisp
  type(scattering)    :: sct
  type(kpointmesh)    :: kmesh
  type(algorithm)     :: algo
  type(runinfo)       :: info

  real(8) :: zqp
  real(8) :: gam
  real(8), allocatable :: enrgy(:,:)
  real(8), allocatable :: enrgydiff(:)

  complex(8), allocatable :: calc_cond(:)
  complex(8), allocatable :: calc_seeb(:)

  integer :: index1(9), index2(9)
  integer :: i,j,idir
  integer :: iband1, iband2, iband

  index1 = (/1,2,3,1,1,2,1,1,2/)
  index2 = (/1,2,3,2,3,3,2,3,3/)


  allocate(enrgy(edisp%nbopt_min:edisp%nbopt_max,edisp%ispin))
  allocate(enrgydiff(edisp%ispin))
  allocate(calc_cond(edisp%ispin))
  allocate(calc_seeb(edisp%ispin))

  ! first we write the kernel into the 1 1 component
  if (algo%lScatteringFile) then
    enrgy = sct%zqp(:,info%ik,:) * (edisp%band(edisp%nbopt_min:edisp%nbopt_max,info%ik,:) - mu)

    do iband1 = edisp%nbopt_min, edisp%nbopt_max
      do iband2 = edisp%nbopt_min, edisp%nbopt_max
        if (iband1 == iband2) cycle
        enrgydiff = enrgy(iband1,:) - enrgy(iband2,:)

        calc_cond = 0.d0

        calc_cond = calc_cond &
            + polygamma2fermi(enrgy(iband1,:), info%beta) &
              *  enrgydiff**2 / sct%gam(iband1,info%ik,:)

        calc_cond = calc_cond &
            + polygamma2fermi(enrgy(iband2,:), info%beta) &
              * enrgydiff**2 / sct%gam(iband2,info%ik,:)

        calc_cond = calc_cond &
            * sct%gam(iband1,info%ik,:) * sct%gam(iband2,info%ik,:) &
            * sct%zqp(iband1,info%ik,:) * sct%zqp(iband2,info%ik,:) &
            * info%beta

        calc_cond = calc_cond &
            / (2.d0 * pi**3 * ( enrgydiff**2 + (sct%gam(iband1,info%ik,:) - sct%gam(iband2,info%ik,:))**2)) &
            / ( enrgydiff**2 + (sct%gam(iband1,info%ik,:) + sct%gam(iband2,info%ik,:))**2)


        calc_seeb = 0.d0

        calc_seeb = calc_seeb &
            + polygamma2fermi(enrgy(iband1,:), info%beta) &
              * enrgy(iband1,:) * enrgydiff**2 / sct%gam(iband1,info%ik,:)

        calc_seeb = calc_seeb &
            + polygamma2fermi(enrgy(iband2,:), info%beta) &
              * enrgy(iband2,:) * enrgydiff**2 / sct%gam(iband2,info%ik,:)

        calc_seeb = calc_seeb &
            * sct%gam(iband1,info%ik,:) * sct%gam(iband2,info%ik,:) &
            * sct%zqp(iband1,info%ik,:) * sct%zqp(iband2,info%ik,:) &
            * info%beta

        calc_seeb = calc_seeb &
            / (2.d0 * pi**3 * ( enrgydiff**2 + (sct%gam(iband1,info%ik,:) - sct%gam(iband2,info%ik,:))**2)) &
            / ( enrgydiff**2 + (sct%gam(iband1,info%ik,:) + sct%gam(iband2,info%ik,:))**2)


        ! multiply optical elements
        do idir = 1,edisp%iOptical
          if (idir <= 6) then
            ! ATTENTION
            ! we read the wien2k files via 1 - 2
            ! we save those in hdf5
            ! and read them via Fortran -> implicit transposition
            ! we have to use 2 - 1 here
            resp%s_full(index1(idir),index2(idir),iband1,:,info%ik) = &
            resp%s_full(index1(idir),index2(idir),iband1,:,info%ik) + calc_cond * edisp%Mopt(idir,iband2,iband1,:)

            resp%a_full(index1(idir),index2(idir),iband1,:,info%ik) = &
            resp%a_full(index1(idir),index2(idir),iband1,:,info%ik) + calc_seeb * edisp%Mopt(idir,iband2,iband1,:)
          else
            ! here we ADD the complex part to the response
            resp%s_full(index1(idir),index2(idir),iband1,:,info%ik) = &
            resp%s_full(index1(idir),index2(idir),iband1,:,info%ik) + calc_cond * edisp%Mopt(idir,iband2,iband1,:) * ci

            resp%a_full(index1(idir),index2(idir),iband1,:,info%ik) = &
            resp%a_full(index1(idir),index2(idir),iband1,:,info%ik) + calc_seeb * edisp%Mopt(idir,iband2,iband1,:) * ci
          endif
        enddo

      enddo
    enddo
  else

    enrgy = sct%zqpscalar * (edisp%band(edisp%nbopt_min:edisp%nbopt_max,info%ik,:) - mu)

    do iband1 = edisp%nbopt_min, edisp%nbopt_max
      do iband2 = edisp%nbopt_min, edisp%nbopt_max
        if (iband1 == iband2) cycle
        enrgydiff = enrgy(iband1,:) - enrgy(iband2,:)

        calc_cond = 0.d0

        calc_cond = calc_cond &
            + polygamma2fermi(enrgy(iband1,:), info%beta) &
              *  enrgydiff**2 / sct%gamscalar

        calc_cond = calc_cond &
            + polygamma2fermi(enrgy(iband2,:), info%beta) &
              * enrgydiff**2 / sct%gamscalar

        calc_cond = calc_cond &
            * sct%gamscalar**2 &
            * sct%zqpscalar**2 &
            * info%beta

        calc_cond = calc_cond &
            / (2.d0 * pi**2 * enrgydiff**2 ) &
            / ( enrgydiff**2 + (2.d0 * sct%gamscalar)**2)


        calc_seeb = 0.d0

        calc_seeb = calc_seeb &
            + polygamma2fermi(enrgy(iband1,:), info%beta) &
              * enrgy(iband1,:) * enrgydiff**2 / sct%gamscalar

        calc_seeb = calc_seeb &
            + polygamma2fermi(enrgy(iband2,:), info%beta) &
              * enrgy(iband2,:) * enrgydiff**2 / sct%gamscalar

        calc_seeb = calc_seeb &
            * sct%gamscalar**2 &
            * sct%zqpscalar**2 &
            * info%beta

        calc_seeb = calc_seeb &
            / (2.d0 * pi**3 * enrgydiff**2 ) &
            / ( enrgydiff**2 + (2.d0 * sct%gamscalar)**2)


        ! multiply optical elements
        do idir = 1,edisp%iOptical
          if (idir <= 6) then
            ! ATTENTION
            ! we read the wien2k files via 1 - 2
            ! we save those in hdf5
            ! and read them via Fortran -> implicit transposition
            ! we have to use 2 - 1 here
            resp%s_full(index1(idir),index2(idir),iband1,:,info%ik) = &
            resp%s_full(index1(idir),index2(idir),iband1,:,info%ik) + calc_cond * edisp%Mopt(idir,iband2,iband1,:)

            resp%a_full(index1(idir),index2(idir),iband1,:,info%ik) = &
            resp%a_full(index1(idir),index2(idir),iband1,:,info%ik) + calc_seeb * edisp%Mopt(idir,iband2,iband1,:)
          else
            ! here we ADD the complex part to the response
            resp%s_full(index1(idir),index2(idir),iband1,:,info%ik) = &
            resp%s_full(index1(idir),index2(idir),iband1,:,info%ik) + calc_cond * edisp%Mopt(idir,iband2,iband1,:) * ci

            resp%a_full(index1(idir),index2(idir),iband1,:,info%ik) = &
            resp%a_full(index1(idir),index2(idir),iband1,:,info%ik) + calc_seeb * edisp%Mopt(idir,iband2,iband1,:) * ci
          endif
        enddo

      enddo
    enddo

  endif


  deallocate(enrgy)
  deallocate(enrgydiff)
  deallocate(calc_cond, calc_seeb)

end subroutine response_inter_Boltzmann_km


! multiply optical elements onto quantities without B-Field
subroutine response_intra_optical_weights(resp, edisp, info)
  implicit none
  type (response_dp) :: resp

  type(energydisp)   :: edisp
  type(runinfo)      :: info

  integer :: index1(9), index2(9)
  integer :: iband, idir

  !( 1 4+i7 5+i8 )
  !( - 2    6+i9 )
  !( - -    3    )

  ! we use these two index lists to move along the described order in the 3x3 matrix
  index1 = (/1,2,3,1,1,2,1,1,2/)
  index2 = (/1,2,3,2,3,3,2,3,3/)

  do iband = edisp%nbopt_min, edisp%nbopt_max
    ! the kernels are saved in the 1 1 directions
    ! so we calculate the 1 1 component at the end
    do idir = 2,edisp%iOptical
      if (idir <= 6) then
        resp%s_full(index1(idir),index2(idir),iband,:,info%ik) = &
        resp%s_full(1,1,iband,:,info%ik) * edisp%Mopt(idir,iband,iband,:)

        resp%a_full(index1(idir),index2(idir),iband,:,info%ik) = &
        resp%a_full(1,1,iband,:,info%ik) * edisp%Mopt(idir,iband,iband,:)
      else
        ! here we ADD the complex part to the response
        resp%s_full(index1(idir),index2(idir),iband,:,info%ik) = &
        resp%s_full(index1(idir),index2(idir),iband,:,info%ik) * edisp%Mopt(idir,iband,iband,:) * ci

        resp%a_full(index1(idir),index2(idir),iband,:,info%ik) = &
        resp%a_full(index1(idir),index2(idir),iband,:,info%ik) * edisp%Mopt(idir,iband,iband,:) * ci
      endif
    enddo

    resp%s_full(1,1,iband,:,info%ik) = &
    resp%s_full(1,1,iband,:,info%ik) * edisp%Mopt(1,iband,iband,:)

    resp%a_full(1,1,iband,:,info%ik) = &
    resp%a_full(1,1,iband,:,info%ik) * edisp%Mopt(1,iband,iband,:)

    ! symmetrize
    ! do idir=6,9
    !   resp%s_full(index2(idir),index1(idir),iband,:,info%ik) = conjg(resp%s_full(index1(idir),index2(idir),iband,:,info%ik))
    ! enddo
  enddo
end subroutine response_intra_optical_weights

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

subroutine response_h5_output(resp, gname, edisp, algo, info, temp, kmesh, lBfield)
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
  integer(hid_t)     :: ifile

  integer :: iband, ik


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
  if (algo%lFullOutput) then
    ! we gather all the data at the master node and write it to hdf5
    if (myid .eq. master) then
      allocate(resp%s_gather(3,3,edisp%nbopt_min:edisp%nbopt_max,edisp%ispin,kmesh%nkp))
    else
      allocate(resp%s_gather(1,1,1,1,1))
    endif
#ifdef MPI
    call MPI_gatherv(resp%s_full,(ikend-ikstr+1)*9*(edisp%nbopt_max-edisp%nbopt_min+1)*edisp%ispin, &
                     MPI_DOUBLE_COMPLEX, resp%s_gather, rcounts*9*(edisp%nbopt_max-edisp%nbopt_min+1)*edisp%ispin, &
                     displs*9*(edisp%nbopt_max-edisp%nbopt_min+1)*edisp%ispin, MPI_DOUBLE_COMPLEX, master, &
                     MPI_COMM_WORLD, mpierr)
#else
    resp%s_gather = resp%s_full
#endif

    if (myid .eq. master) then
      allocate(resp%a_gather(3,3,edisp%nbopt_min:edisp%nbopt_max,edisp%ispin,kmesh%nkp))
    else
      allocate(resp%a_gather(1,1,1,1,1))
    endif
#ifdef MPI
    call MPI_gatherv(resp%a_full,(ikend-ikstr+1)*9*(edisp%nbopt_max-edisp%nbopt_min+1)*edisp%ispin, &
                     MPI_DOUBLE_COMPLEX, resp%a_gather, rcounts*9*(edisp%nbopt_max-edisp%nbopt_min+1)*edisp%ispin, &
                     displs*9*(edisp%nbopt_max-edisp%nbopt_min+1)*edisp%ispin, MPI_DOUBLE_COMPLEX, master, &
                     MPI_COMM_WORLD, mpierr)
#else
    resp%a_gather = resp%a_full
#endif

    resp%s_gather = resp%s_gather * 2.d0 * pi * ( echarge / (kmesh%vol*hbarevs)) * 1.d10
    resp%a_gather = resp%a_gather * 2.d0 * pi * ( echarge / (kmesh%vol*hbarevs)) * 1.d10 * kB

    if (myid .eq. master) then
      write(string,'(I6.6)') info%iT
      string = trim(string) // "/conductivity/" // trim(adjustl(gname)) // "/full"
      call hdf5_write_data(ifile, string, resp%s_gather)

      write(string,'(I6.6)') info%iT
      string = trim(string) // "/peltier/" // trim(adjustl(gname)) // "/full"
      call hdf5_write_data(ifile, string, resp%a_gather)

    endif

    deallocate(resp%s_gather)
    deallocate(resp%a_gather)
  endif ! full output

  ! perform a local summation
  do ik = ikstr,ikend
    do iband = edisp%nbopt_min,edisp%nbopt_max
      resp%s_sum(:,:,:) = resp%s_sum(:,:,:) + resp%s_full(:,:,iband,:,ik) * kmesh%weight(ik)
      resp%a_sum(:,:,:) = resp%a_sum(:,:,:) + resp%a_full(:,:,iband,:,ik) * kmesh%weight(ik)
    enddo
  enddo

  ! perform MPI summation
#ifdef MPI
  if (myid.eq.master) then
    call MPI_REDUCE(MPI_IN_PLACE, resp%s_sum, 9*edisp%ispin, MPI_DOUBLE_COMPLEX, MPI_SUM, master, MPI_COMM_WORLD, mpierr)
  else
    call MPI_REDUCE(resp%s_sum, resp%s_sum, 9*edisp%ispin, MPI_DOUBLE_COMPLEX, MPI_SUM, master, MPI_COMM_WORLD, mpierr)
  endif

  if (myid.eq.master) then
    call MPI_REDUCE(MPI_IN_PLACE, resp%a_sum, 9*edisp%ispin, MPI_DOUBLE_COMPLEX, MPI_SUM, master, MPI_COMM_WORLD, mpierr)
  else
    call MPI_REDUCE(resp%a_sum, resp%a_sum, 9*edisp%ispin, MPI_DOUBLE_COMPLEX, MPI_SUM, master, MPI_COMM_WORLD, mpierr)
  endif
#endif

  resp%s_sum = resp%s_sum * 2.d0 * pi * ( echarge / (kmesh%vol*hbarevs)) * 1.d10
  resp%a_sum = resp%a_sum * 2.d0 * pi * ( echarge / (kmesh%vol*hbarevs)) * 1.d10 * kB

  if (myid .eq. master) then
    write(string,'(I6.6)') info%iT
    string = trim(string) // "/conductivity/" // trim(adjustl(gname)) // "/sum"
    call hdf5_write_data(ifile, string, resp%s_sum)

    write(string,'(I6.6)') info%iT
    string = trim(string) // "/peltier/" // trim(adjustl(gname)) // "/sum"
    call hdf5_write_data(ifile, string, resp%a_sum)

  endif


  if (lBoutput) then
    ! conductivity and seebeck coefficient without B-field
    if (algo%lFullOutput) then
      if (myid .eq. master) then
        allocate(resp%sB_gather(3,3,edisp%nbopt_min:edisp%nbopt_max,edisp%ispin,kmesh%nkp))
      else
        allocate(resp%sB_gather(1,1,1,1,1))
      endif
#ifdef MPI
      call MPI_gatherv(resp%sB_full,(ikend-ikstr+1)*9*(edisp%nbopt_max-edisp%nbopt_min+1)*edisp%ispin, &
                       MPI_DOUBLE_COMPLEX, resp%sB_gather, rcounts*9*(edisp%nbopt_max-edisp%nbopt_min+1)*edisp%ispin, &
                       displs*9*(edisp%nbopt_max-edisp%nbopt_min+1)*edisp%ispin, MPI_DOUBLE_COMPLEX, master, &
                       MPI_COMM_WORLD, mpierr)
#else
      resp%sB_gather = resp%sB_full
#endif

      if (myid .eq. master) then
        allocate(resp%aB_gather(3,3,edisp%nbopt_min:edisp%nbopt_max,edisp%ispin,kmesh%nkp))
      else
        allocate(resp%aB_gather(1,1,1,1,1))
      endif
#ifdef MPI
      call MPI_gatherv(resp%aB_full,(ikend-ikstr+1)*9*(edisp%nbopt_max-edisp%nbopt_min+1)*edisp%ispin, &
                       MPI_DOUBLE_COMPLEX, resp%aB_gather, rcounts*9*(edisp%nbopt_max-edisp%nbopt_min+1)*edisp%ispin, &
                       displs*9*(edisp%nbopt_max-edisp%nbopt_min+1)*edisp%ispin, MPI_DOUBLE_COMPLEX, master, &
                       MPI_COMM_WORLD, mpierr)
#else
      resp%aB_gather = resp%aB_full
#endif

      resp%sB_gather = resp%sB_gather * 2.d0 * pi**2 * ( echarge / (kmesh%vol*hbarevs)) * (1.d-10/hbarevs)
      resp%aB_gather = resp%aB_gather * 2.d0 * pi**2 * ( echarge / (kmesh%vol*hbarevs)) * (1.d-10/hbarevs) * kB

      if (myid .eq. master) then
        write(string,'(I6.6)') info%iT
        string = trim(string) // "/conductivity/" // trim(adjustl(gname)) // "/fullB"
        call hdf5_write_data(ifile, string, resp%sB_gather)

        write(string,'(I6.6)') info%iT
        string = trim(string) // "/peltier/" // trim(adjustl(gname)) // "/fullB"
        call hdf5_write_data(ifile, string, resp%aB_gather)

      endif

      deallocate(resp%sB_gather)
      deallocate(resp%aB_gather)
    endif ! full output

    ! perform a local summation
    do ik = ikstr,ikend
      do iband = edisp%nbopt_min,edisp%nbopt_max
        resp%sB_sum(:,:,:) = resp%sB_sum(:,:,:) + resp%sB_full(:,:,iband,:,ik) * kmesh%weight(ik)
        resp%aB_sum(:,:,:) = resp%aB_sum(:,:,:) + resp%aB_full(:,:,iband,:,ik) * kmesh%weight(ik)
      enddo
    enddo

  ! perform MPI summation
#ifdef MPI
  if (myid.eq.master) then
    call MPI_REDUCE(MPI_IN_PLACE, resp%sB_sum, 9*edisp%ispin, MPI_DOUBLE_COMPLEX, MPI_SUM, master, MPI_COMM_WORLD, mpierr)
  else
    call MPI_REDUCE(resp%sB_sum, resp%sB_sum, 9*edisp%ispin, MPI_DOUBLE_COMPLEX, MPI_SUM, master, MPI_COMM_WORLD, mpierr)
  endif

  if (myid.eq.master) then
    call MPI_REDUCE(MPI_IN_PLACE, resp%aB_sum, 9*edisp%ispin, MPI_DOUBLE_COMPLEX, MPI_SUM, master, MPI_COMM_WORLD, mpierr)
  else
    call MPI_REDUCE(resp%aB_sum, resp%aB_sum, 9*edisp%ispin, MPI_DOUBLE_COMPLEX, MPI_SUM, master, MPI_COMM_WORLD, mpierr)
  endif
#endif

    resp%sB_sum = resp%sB_sum * 2.d0 * pi**2 * ( echarge / (kmesh%vol*hbarevs)) * (1.d-10/hbarevs)
    resp%aB_sum = resp%aB_sum * 2.d0 * pi**2 * ( echarge / (kmesh%vol*hbarevs)) * (1.d-10/hbarevs) * kB

    if (myid .eq. master) then
      write(string,'(I6.6)') info%iT
      string = trim(string) // "/conductivity/" // trim(adjustl(gname)) // "/sumB"
      call hdf5_write_data(ifile, string, resp%sB_sum)

      write(string,'(I6.6)') info%iT
      string = trim(string) // "/peltier/" // trim(adjustl(gname)) // "/sumB"
      call hdf5_write_data(ifile, string, resp%aB_sum)

    endif
  endif ! Boutput

  if (myid.eq.master) then
    call hdf5_close_file(ifile)
  endif

end subroutine


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
