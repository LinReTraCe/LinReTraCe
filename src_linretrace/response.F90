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
  dresp%x_full = 0.d0
  dresp%s_sum = 0.d0
  dresp%a_sum = 0.d0
  dresp%x_sum = 0.d0

  if (algo%lBfield) then
     dresp%sB_full = 0.d0
     dresp%aB_full = 0.d0
     dresp%xB_full = 0.d0
     dresp%sB_sum = 0.d0
     dresp%aB_sum = 0.d0
     dresp%xB_sum = 0.d0
  endif
end subroutine initresp

subroutine initresp_qp (algo, qresp)
  implicit none
  type(algorithm)   :: algo
  type(response_qp) :: qresp

  qresp%s_full = 0.q0
  qresp%a_full = 0.q0
  qresp%x_full = 0.q0
  qresp%s_sum = 0.q0
  qresp%a_sum = 0.q0
  qresp%x_sum = 0.q0

  if (algo%lBfield) then
     qresp%sB_full = 0.q0
     qresp%aB_full = 0.q0
     qresp%xB_full = 0.q0
     qresp%sB_sum = 0.q0
     qresp%aB_sum = 0.q0
     qresp%xB_sum = 0.q0
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

  real(8), allocatable :: enrgy(:,:)

  allocate(enrgy(edisp%nbopt_min:edisp%nbopt_max,edisp%ispin))
  ! first we write the kernel into the 1 1 component
  enrgy = sct%zqp(edisp%nbopt_min:edisp%nbopt_max,info%ik,:) * (edisp%band(edisp%nbopt_min:edisp%nbopt_max,info%ik,:) - mu)

  resp%s_full(1,1,:,:,info%ik) = real(PolyGamma(1,:,info%ik,:)) &
                               - real(PolyGamma(2,:,info%ik,:)) * info%beta2p*sct%gam(edisp%nbopt_min:edisp%nbopt_max,info%ik,:)

  resp%s_full(1,1,:,:,info%ik) = resp%s_full(1,1,:,:,info%ik) &
                               * sct%zqp(edisp%nbopt_min:edisp%nbopt_max,info%ik,:)**2 * info%beta &
                               / (4.d0 * pi**3 * sct%gam(edisp%nbopt_min:edisp%nbopt_max,info%ik,:))


  resp%a_full(1,1,:,:,info%ik) =  real(PolyGamma(1,:,info%ik,:)) * enrgy &
                               -  real(PolyGamma(2,:,info%ik,:)) * enrgy &
                                  * sct%gam(edisp%nbopt_min:edisp%nbopt_max,info%ik,:) * info%beta2p &
                               - aimag(PolyGamma(2,:,info%ik,:)) * sct%gam(edisp%nbopt_min:edisp%nbopt_max,info%ik,:)**2 &
                                 * info%beta2p

  resp%a_full(1,1,:,:,info%ik) = resp%a_full(1,1,:,:,info%ik) &
                               * sct%zqp(edisp%nbopt_min:edisp%nbopt_max,info%ik,:)**2 * info%beta &
                               / (4.d0 * pi**3 * sct%gam(edisp%nbopt_min:edisp%nbopt_max,info%ik,:))


  resp%x_full(1,1,:,:,info%ik) =  real(PolyGamma(1,:,info%ik,:)) &
                                  * (enrgy**2 + sct%gam(edisp%nbopt_min:edisp%nbopt_max,info%ik,:)**2) &
                               +  real(PolyGamma(2,:,info%ik,:)) &
                                  * info%beta2p * sct%gam(edisp%nbopt_min:edisp%nbopt_max,info%ik,:) &
                                  * (sct%gam(edisp%nbopt_min:edisp%nbopt_max,info%ik,:)**2 - enrgy**2) &
                               - aimag(PolyGamma(2,:,info%ik,:)) &
                                  * info%beta / pi * enrgy * sct%gam(edisp%nbopt_min:edisp%nbopt_max,info%ik,:)**2

  resp%x_full(1,1,:,:,info%ik) = resp%x_full(1,1,:,:,info%ik) &
                               * sct%zqp(edisp%nbopt_min:edisp%nbopt_max,info%ik,:)**2 * info%beta &
                               / (4.d0 * pi**3 * sct%gam(edisp%nbopt_min:edisp%nbopt_max,info%ik,:))

  if (algo%lBfield) then

    resp%sB_full(1,1,:,:,info%ik) = real(PolyGamma(3,:,info%ik,:)) &
                                      * sct%gam(edisp%nbopt_min:edisp%nbopt_max,info%ik,:)**2 * info%beta**2 / (4.d0 * pi**2) &
                                  - real(PolyGamma(2,:,info%ik,:)) &
                                    * 3.d0 * sct%gam(edisp%nbopt_min:edisp%nbopt_max,info%ik,:) * info%beta2p &
                                  - real(PolyGamma(1,:,info%ik,:)) * 3.d0

    resp%sB_full(1,1,:,:,info%ik) = resp%sB_full(1,1,:,:,info%ik) &
                                  * sct%zqp(edisp%nbopt_min:edisp%nbopt_max,info%ik,:)**3 * info%beta &
                                  / (16.d0 * pi**4 * sct%gam(edisp%nbopt_min:edisp%nbopt_max,info%ik,:)**2)


    resp%aB_full(1,1,:,:,info%ik) =  real(PolyGamma(3,:,info%ik,:)) &
                                       * enrgy * sct%gam(edisp%nbopt_min:edisp%nbopt_max,info%ik,:)**2 &
                                       * info%beta**2 / (4.d0 * pi**2) &
                                  + aimag(PolyGamma(3,:,info%ik,:)) &
                                       * sct%gam(edisp%nbopt_min:edisp%nbopt_max,info%ik,:)**3 * info%beta**3 / (4.d0 * pi**2) &
                                  -  real(PolyGamma(2,:,info%ik,:)) &
                                       * 3.d0 * enrgy * sct%gam(edisp%nbopt_min:edisp%nbopt_max,info%ik,:) * info%beta2p &
                                  - aimag(PolyGamma(2,:,info%ik,:)) &
                                       * sct%gam(edisp%nbopt_min:edisp%nbopt_max,info%ik,:)**2 * info%beta2p &
                                  +  real(PolyGamma(1,:,info%ik,:)) * 3.d0 * enrgy

    resp%aB_full(1,1,:,:,info%ik) = resp%aB_full(1,1,:,:,info%ik) &
                                  * sct%zqp(edisp%nbopt_min:edisp%nbopt_max,info%ik,:)**3 * info%beta / (16.d0 * pi**4 &
                                  * sct%gam(edisp%nbopt_min:edisp%nbopt_max,info%ik,:)**2)


    resp%xB_full(1,1,:,:,info%ik) =  real(PolyGamma(3,:,info%ik,:)) &
                                       * info%beta**2 * sct%gam(edisp%nbopt_min:edisp%nbopt_max,info%ik,:)**2 &
                                       * (sct%gam(edisp%nbopt_min:edisp%nbopt_max,info%ik,:)**2 - enrgy**2) &
                                       / (4.d0 * pi**2) &
                                  - aimag(PolyGamma(3,:,info%ik,:)) &
                                       * info%beta**2 * enrgy * sct%gam(edisp%nbopt_min:edisp%nbopt_max,info%ik,:)**3 &
                                       / (2.d0 * pi) &
                                  -  real(PolyGamma(2,:,info%ik,:)) &
                                       * info%beta * sct%gam(edisp%nbopt_min:edisp%nbopt_max,info%ik,:) &
                                       * (sct%gam(edisp%nbopt_min:edisp%nbopt_max,info%ik,:)**2 + 3.d0*enrgy**2) &
                                       / (2.d0 * pi) &
                                  + aimag(PolyGamma(2,:,info%ik,:)) &
                                       * info%beta * enrgy * sct%gam(edisp%nbopt_min:edisp%nbopt_max,info%ik,:)**2 / pi &
                                  -  real(PolyGamma(1,:,info%ik,:)) &
                                       * (sct%gam(edisp%nbopt_min:edisp%nbopt_max,info%ik,:)**2 + 3.d0*enrgy**2)

    resp%xB_full(1,1,:,:,info%ik) = resp%xB_full(1,1,:,:,info%ik) &
                                  * (-1.d0) * sct%zqp(edisp%nbopt_min:edisp%nbopt_max,info%ik,:)**3 * info%beta &
                                  / (16.d0 * pi**4 * sct%gam(edisp%nbopt_min:edisp%nbopt_max,info%ik,:)**2)

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
  real(8), allocatable :: gamdiff(:)

  complex(8) :: calc_sigma
  complex(8) :: calc_alpha
  complex(8) :: calc_xi

  integer :: index1(9), index2(9)
  integer :: i,j,idir
  integer :: iband1, iband2, iband, is

  ! index1 = (/1,2,3,1,1,2,1,1,2/)
  ! index2 = (/1,2,3,2,3,3,2,3,3/)

  ! NOTE: here we transpose it internally in Fortran
  ! so the output (hdf5 is in the correct order)
  index1 = (/1,2,3,2,3,3,2,3,3/)
  index2 = (/1,2,3,1,1,2,1,1,2/)


  allocate(enrgy(edisp%nbopt_min:edisp%nbopt_max,edisp%ispin))
  allocate(enrgydiff(edisp%ispin))

  ! first we write the kernel into the 1 1 component
  enrgy = sct%zqp(edisp%nbopt_min:edisp%nbopt_max,info%ik,:) &
          * (edisp%band(edisp%nbopt_min:edisp%nbopt_max,info%ik,:) - mu)

  do iband1 = edisp%nbopt_min, edisp%nbopt_max
    do iband2 = edisp%nbopt_min, edisp%nbopt_max
      if (iband1 == iband2) cycle
      enrgydiff = enrgy(iband1,:) - enrgy(iband2,:)
      gamdiff   = sct%gam(iband1,info%ik,:) - sct%gam(iband2,info%ik,:)

      do is = 1,edisp%ispin

        if ((abs(enrgydiff(is)) .lt. 1d-6) .and. (abs(gamdiff(is)) .lt. 1d-6)) then

          ! use the intra-band limit .....
          calc_sigma  = real(PolyGamma(1,iband1,info%ik,is)) &
                      - real(PolyGamma(2,iband1,info%ik,is)) * info%beta2p*sct%gam(iband1,info%ik,is)

          calc_sigma = calc_sigma &
                     * sct%zqp(iband1,info%ik,is)**2 * info%beta &
                     / (4.d0 * pi**3 * sct%gam(iband1,info%ik,is))

          calc_alpha =  real(PolyGamma(1,iband1,info%ik,is)) * enrgy(iband1,is) &
                     -  real(PolyGamma(2,iband1,info%ik,is)) * enrgy(iband1,is) * sct%gam(iband1,info%ik,is) * info%beta2p &
                     - aimag(PolyGamma(2,iband1,info%ik,is)) * sct%gam(iband1,info%ik,is)**2 * info%beta2p

          calc_alpha = calc_alpha &
                     * sct%zqp(iband1,info%ik,is)**2 * info%beta &
                     / (4.d0 * pi**3 * sct%gam(iband1,info%ik,is))

          calc_xi    =  real(PolyGamma(1,iband1,info%ik,is)) &
                        * (enrgy(iband1,is)**2 + sct%gam(iband1,info%ik,is)**2) &
                     +  real(PolyGamma(2,iband1,info%ik,is)) &
                        * info%beta2p * sct%gam(iband1,info%ik,is) * (sct%gam(iband1,info%ik,is)**2 - enrgy(iband1,is)**2) &
                     - aimag(PolyGamma(2,iband1,info%ik,is)) &
                        * info%beta / pi * enrgy(iband1,is) * sct%gam(iband1,info%ik,is)**2

          calc_xi    = calc_xi &
                       * sct%zqp(iband1,info%ik,is)**2 * info%beta &
                       / (4.d0 * pi**3 * sct%gam(iband1,info%ik,is))

        else
          calc_sigma = real(PolyGamma(1,iband1,info%ik,is)) &
                     * (enrgydiff(is)**2 + sct%gam(iband2,info%ik,is)**2 - sct%gam(iband1,info%ik,is)**2) &
                     / sct%gam(iband1,info%ik,is)

          calc_sigma = calc_sigma &
                     + aimag(PolyGamma(1,iband1,info%ik,is)) &
                     * (2.d0* (-enrgydiff(is)))

          calc_sigma = calc_sigma &
                     + real(PolyGamma(1,iband2,info%ik,is)) &
                     * (enrgydiff(is)**2 + sct%gam(iband1,info%ik,is)**2 - sct%gam(iband2,info%ik,is)**2) &
                     / sct%gam(iband2,info%ik,is)

          calc_sigma = calc_sigma &
                     + aimag(PolyGamma(1,iband2,info%ik,is)) &
                     * (2.d0* enrgydiff(is))

          calc_sigma = calc_sigma &
                     * sct%gam(iband1,info%ik,is) * sct%gam(iband2,info%ik,is) &
                     * sct%zqp(iband1,info%ik,is) * sct%zqp(iband2,info%ik,is) &
                     * info%beta

          calc_sigma = calc_sigma &
                     / (2.d0 * pi**3 * ( enrgydiff(is)**2 + (sct%gam(iband1,info%ik,is) - sct%gam(iband2,info%ik,is))**2)) &
                     / ( enrgydiff(is)**2 + (sct%gam(iband1,info%ik,is) + sct%gam(iband2,info%ik,is))**2)



          calc_alpha = real(PolyGamma(1,iband1,info%ik,is)) &
                     * (enrgy(iband1,is) * (enrgydiff(is)**2 + sct%gam(iband2,info%ik,is)**2 - sct%gam(iband1,info%ik,is)**2) &
                     / sct%gam(iband1,info%ik,is) + 2.d0*sct%gam(iband1,info%ik,is)*enrgydiff(is))

          calc_alpha = calc_alpha &
                     + aimag(PolyGamma(1,iband1,info%ik,is)) &
                     * (enrgydiff(is)**2 + sct%gam(iband2,info%ik,is)**2 - sct%gam(iband1,info%ik,is)**2 &
                     - 2.d0*enrgy(iband1,is)*enrgydiff(is))

          calc_alpha = calc_alpha &
                     + real(PolyGamma(1,iband2,info%ik,is)) &
                     * (enrgy(iband2,is) * (enrgydiff(is)**2 + sct%gam(iband1,info%ik,is)**2 - sct%gam(iband2,info%ik,is)**2) &
                     / sct%gam(iband2,info%ik,is) + 2.d0*sct%gam(iband2,info%ik,is)*enrgydiff(is)*(-1.d0))

          calc_alpha = calc_alpha &
                     + aimag(PolyGamma(1,iband2,info%ik,is)) &
                     * (enrgydiff(is)**2 + sct%gam(iband1,info%ik,is)**2 - sct%gam(iband2,info%ik,is)**2 &
                     - 2.d0*enrgy(iband2,is)*enrgydiff(is)*(-1.d0))

          calc_alpha = calc_alpha &
                     * sct%gam(iband1,info%ik,is) * sct%gam(iband2,info%ik,is) &
                     * sct%zqp(iband1,info%ik,is) * sct%zqp(iband2,info%ik,is) &
                     * info%beta

          calc_alpha = calc_alpha &
                     / (2.d0 * pi**3 * ( enrgydiff(is)**2 + (sct%gam(iband1,info%ik,is) - sct%gam(iband2,info%ik,is))**2)) &
                     / ( enrgydiff(is)**2 + (sct%gam(iband1,info%ik,is) + sct%gam(iband2,info%ik,is))**2)


          calc_xi    = real(PolyGamma(1,iband1,info%ik,is) &
                     * (enrgy(iband1,is) - ci*sct%gam(iband1,info%ik,is))**2 &
                     * (enrgydiff(is) + ci*gamdiff(is)) &
                     * (enrgydiff(is) + ci*(sct%gam(iband1,info%ik,is) + sct%gam(iband2,info%ik,is)))) &
                     * sct%zqp(iband1,info%ik,is) * sct%zqp(iband2,info%ik,is) * sct%gam(iband2,info%ik,is) &
                     * info%beta / (2.d0 * pi**3)

          calc_xi    = calc_xi &
                     + real(PolyGamma(1,iband2,info%ik,is) &
                     * (enrgy(iband2,is) - ci*sct%gam(iband2,info%ik,is))**2 &
                     * (-enrgydiff(is) - ci*gamdiff(is)) &
                     * (-enrgydiff(is) + ci*(sct%gam(iband1,info%ik,is) + sct%gam(iband2,info%ik,is)))) &
                     * sct%zqp(iband1,info%ik,is) * sct%zqp(iband2,info%ik,is) * sct%gam(iband1,info%ik,is) &
                     * info%beta / (2.d0 * pi**3)

          calc_xi    = calc_xi &
                     / ((enrgydiff(is)**2 + gamdiff(is)**2) &
                      * (enrgydiff(is)**2 + (sct%gam(iband1,info%ik,is) + sct%gam(iband2,info%ik,is))**2))
        endif


        ! multiply optical elements
        do idir = 1,edisp%iOptical
          if (idir <= 6) then
            ! ATTENTION
            ! we read the wien2k files via 1 - 2
            ! we save those in hdf5
            ! and read them via Fortran -> implicit transposition
            ! we have to use 2 - 1 here
            resp%s_full(index1(idir),index2(idir),iband1,is,info%ik) = &
            resp%s_full(index1(idir),index2(idir),iband1,is,info%ik) + calc_sigma * edisp%Mopt(idir,iband2,iband1,is,info%ik)

            resp%a_full(index1(idir),index2(idir),iband1,is,info%ik) = &
            resp%a_full(index1(idir),index2(idir),iband1,is,info%ik) + calc_alpha * edisp%Mopt(idir,iband2,iband1,is,info%ik)

            resp%x_full(index1(idir),index2(idir),iband1,is,info%ik) = &
            resp%x_full(index1(idir),index2(idir),iband1,is,info%ik) + calc_xi    * edisp%Mopt(idir,iband2,iband1,is,info%ik)
          else
            ! here we ADD the complex part to the response
            resp%s_full(index1(idir),index2(idir),iband1,is,info%ik) = &
            resp%s_full(index1(idir),index2(idir),iband1,is,info%ik) + calc_sigma * edisp%Mopt(idir,iband2,iband1,is,info%ik) * ci

            resp%a_full(index1(idir),index2(idir),iband1,is,info%ik) = &
            resp%a_full(index1(idir),index2(idir),iband1,is,info%ik) + calc_alpha * edisp%Mopt(idir,iband2,iband1,is,info%ik) * ci

            resp%x_full(index1(idir),index2(idir),iband1,is,info%ik) = &
            resp%x_full(index1(idir),index2(idir),iband1,is,info%ik) + calc_xi     * edisp%Mopt(idir,iband2,iband1,is,info%ik) * ci
          endif
        enddo
      enddo
    enddo
  enddo

  ! TODO: FIX THIS
  ! i.e. we save all this data to python
  ! that means an implicit transposition takes place
  ! i.e. we either hav eto treat all of this the other way around
  if (edisp%iOptical > 3) then
    do iband = edisp%nbopt_min, edisp%nbopt_max
      resp%s_full(2,1,iband,:,info%ik) = conjg(resp%s_full(1,2,iband,:,info%ik))
      resp%s_full(3,1,iband,:,info%ik) = conjg(resp%s_full(1,3,iband,:,info%ik))
      resp%s_full(3,2,iband,:,info%ik) = conjg(resp%s_full(2,3,iband,:,info%ik))

      resp%a_full(2,1,iband,:,info%ik) = conjg(resp%a_full(1,2,iband,:,info%ik))
      resp%a_full(3,1,iband,:,info%ik) = conjg(resp%a_full(1,3,iband,:,info%ik))
      resp%a_full(3,2,iband,:,info%ik) = conjg(resp%a_full(2,3,iband,:,info%ik))

      resp%x_full(2,1,iband,:,info%ik) = conjg(resp%x_full(1,2,iband,:,info%ik))
      resp%x_full(3,1,iband,:,info%ik) = conjg(resp%x_full(1,3,iband,:,info%ik))
      resp%x_full(3,2,iband,:,info%ik) = conjg(resp%x_full(2,3,iband,:,info%ik))
    enddo
  endif

  deallocate(enrgy)
  deallocate(enrgydiff)
  deallocate(gamdiff)

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
  enrgy = sct%zqp(edisp%nbopt_min:edisp%nbopt_max,info%ik,:) &
          * (edisp%band(edisp%nbopt_min:edisp%nbopt_max,info%ik,:) - mu)

  ! asymptotic term in the Gamma-> 0 limit
  ! the polygamma2fermi already contains the pi^2 / 2
  resp%s_full(1,1,:,:,info%ik) = polygamma2fermi(enrgy,info%beta) &
                                  * sct%zqp(edisp%nbopt_min:edisp%nbopt_max,info%ik,:)**2 * info%beta &
                                  / (4.d0 * pi**3 * sct%gam(edisp%nbopt_min:edisp%nbopt_max,info%ik,:))

  resp%a_full(1,1,:,:,info%ik) = resp%s_full(1,1,:,:,info%ik) * enrgy

  resp%x_full(1,1,:,:,info%ik) = resp%s_full(1,1,:,:,info%ik) * enrgy**2

  if (algo%lBfield) then

    resp%sB_full(1,1,:,:,info%ik) = polygamma2fermi(enrgy,info%beta) &
                                    * 3.d0 * sct%zqp(edisp%nbopt_min:edisp%nbopt_max,info%ik,:)**3 * info%beta &
                                    / (16.d0 * pi**4 * sct%gam(edisp%nbopt_min:edisp%nbopt_max,info%ik,:)**2)

    resp%aB_full(1,1,:,:,info%ik) = resp%sB_full(1,1,:,:,info%ik) * enrgy

    resp%xB_full(1,1,:,:,info%ik) = resp%sB_full(1,1,:,:,info%ik) * enrgy**2

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
  real(8), allocatable :: gamdiff(:)

  complex(8) :: calc_sigma
  complex(8) :: calc_alpha
  complex(8) :: calc_xi

  integer :: index1(9), index2(9)
  integer :: i,j,idir
  integer :: iband1, iband2, iband, is

  ! NOTE: here we transpose it internally in Fortran
  ! so the output (hdf5 is in the correct order)
  index1 = (/1,2,3,2,3,3,2,3,3/)
  index2 = (/1,2,3,1,1,2,1,1,2/)


  allocate(enrgy(edisp%nbopt_min:edisp%nbopt_max,edisp%ispin))
  allocate(enrgydiff(edisp%ispin))
  allocate(gamdiff(edisp%ispin))

  ! first we write the kernel into the 1 1 component
  enrgy = sct%zqp(edisp%nbopt_min:edisp%nbopt_max,info%ik,:) &
          * (edisp%band(edisp%nbopt_min:edisp%nbopt_max,info%ik,:) - mu)

  do iband1 = edisp%nbopt_min, edisp%nbopt_max
    do iband2 = edisp%nbopt_min, edisp%nbopt_max
      if (iband1 == iband2) cycle
      enrgydiff = enrgy(iband1,:) - enrgy(iband2,:)
      gamdiff   = sct%gam(iband1,info%ik,:) - sct%gam(iband2,info%ik,:)

      do is = 1,edisp%ispin
        if ((abs(enrgydiff(is)) .lt. 1d-6) .and. (abs(gamdiff(is)) .lt. 1d-6)) then
        ! use the intra-band limit .....
          calc_sigma = polygamma2fermi(enrgy(iband1,is),info%beta) &
                     * sct%zqp(iband1,info%ik,is)**2 * info%beta &
                     / (4.d0 * pi**3 * sct%gam(iband1,info%ik,is))

          calc_alpha = calc_sigma * enrgy(iband1,is)

          calc_xi    = calc_sigma * enrgy(iband1,is)**2

        else

          calc_sigma = polygamma2fermi(enrgy(iband1,is), info%beta) &
                     * enrgydiff(is)**2 / sct%gam(iband1,info%ik,is)

          calc_sigma = calc_sigma &
                     + polygamma2fermi(enrgy(iband2,is), info%beta) &
                     * enrgydiff(is)**2 / sct%gam(iband2,info%ik,is)

          calc_sigma = calc_sigma &
                     * sct%gam(iband1,info%ik,is) * sct%gam(iband2,info%ik,is) &
                     * sct%zqp(iband1,info%ik,is) * sct%zqp(iband2,info%ik,is) &
                     * info%beta

          calc_sigma = calc_sigma &
                     / (2.d0 * pi**3 * ( enrgydiff(is)**2 + (sct%gam(iband1,info%ik,is) - sct%gam(iband2,info%ik,is))**2)) &
                     / ( enrgydiff(is)**2 + (sct%gam(iband1,info%ik,is) + sct%gam(iband2,info%ik,is))**2)


          calc_alpha = polygamma2fermi(enrgy(iband1,is), info%beta) &
                     * enrgy(iband1,is) * enrgydiff(is)**2 / sct%gam(iband1,info%ik,is)

          calc_alpha = calc_alpha &
                     + polygamma2fermi(enrgy(iband2,is), info%beta) &
                     * enrgy(iband2,is) * enrgydiff(is)**2 / sct%gam(iband2,info%ik,is)

          calc_alpha = calc_alpha &
                     * sct%gam(iband1,info%ik,is) * sct%gam(iband2,info%ik,is) &
                     * sct%zqp(iband1,info%ik,is) * sct%zqp(iband2,info%ik,is) &
                     * info%beta

          calc_alpha = calc_alpha &
                     / (2.d0 * pi**3 * ( enrgydiff(is)**2 + (sct%gam(iband1,info%ik,is) - sct%gam(iband2,info%ik,is))**2)) &
                     / ( enrgydiff(is)**2 + (sct%gam(iband1,info%ik,is) + sct%gam(iband2,info%ik,is))**2)

          calc_xi    = polygamma2fermi(enrgy(iband1,is), info%beta) &
                     * enrgy(iband1,is)**2 * enrgydiff(is)**2 / sct%gam(iband1,info%ik,is)

          calc_xi    = calc_xi &
                     + polygamma2fermi(enrgy(iband2,is), info%beta) &
                     * enrgy(iband2,is)**2 * enrgydiff(is)**2 / sct%gam(iband2,info%ik,is)

          calc_xi    = calc_xi &
                     * sct%gam(iband1,info%ik,is) * sct%gam(iband2,info%ik,is) &
                     * sct%zqp(iband1,info%ik,is) * sct%zqp(iband2,info%ik,is) &
                     * info%beta

          calc_xi    = calc_xi &
                     / (2.d0 * pi**3 * ( enrgydiff(is)**2 + (sct%gam(iband1,info%ik,is) - sct%gam(iband2,info%ik,is))**2)) &
                     / ( enrgydiff(is)**2 + (sct%gam(iband1,info%ik,is) + sct%gam(iband2,info%ik,is))**2)
        endif


        ! multiply optical elements
        do idir = 1,edisp%iOptical
          if (idir <= 6) then
            ! ATTENTION
            ! we read the wien2k files via 1 - 2
            ! we save those in hdf5
            ! and read them via Fortran -> implicit transposition
            ! we have to use 2 - 1 here
            resp%s_full(index1(idir),index2(idir),iband1,:,info%ik) = &
            resp%s_full(index1(idir),index2(idir),iband1,:,info%ik) + calc_sigma * edisp%Mopt(idir,iband2,iband1,:,info%ik)

            resp%a_full(index1(idir),index2(idir),iband1,:,info%ik) = &
            resp%a_full(index1(idir),index2(idir),iband1,:,info%ik) + calc_alpha * edisp%Mopt(idir,iband2,iband1,:,info%ik)

            resp%x_full(index1(idir),index2(idir),iband1,:,info%ik) = &
            resp%x_full(index1(idir),index2(idir),iband1,:,info%ik) + calc_xi    * edisp%Mopt(idir,iband2,iband1,:,info%ik)
          else
            ! here we ADD the complex part to the response
            resp%s_full(index1(idir),index2(idir),iband1,:,info%ik) = &
            resp%s_full(index1(idir),index2(idir),iband1,:,info%ik) + calc_sigma * edisp%Mopt(idir,iband2,iband1,:,info%ik) * ci

            resp%a_full(index1(idir),index2(idir),iband1,:,info%ik) = &
            resp%a_full(index1(idir),index2(idir),iband1,:,info%ik) + calc_alpha * edisp%Mopt(idir,iband2,iband1,:,info%ik) * ci

            resp%x_full(index1(idir),index2(idir),iband1,:,info%ik) = &
            resp%x_full(index1(idir),index2(idir),iband1,:,info%ik) + calc_xi    * edisp%Mopt(idir,iband2,iband1,:,info%ik) * ci
          endif
        enddo !idir

      enddo !is
    enddo !iband2
  enddo ! iband1

  if (edisp%iOptical > 3) then
    do iband = edisp%nbopt_min, edisp%nbopt_max
      resp%s_full(2,1,iband,:,info%ik) = conjg(resp%s_full(1,2,iband,:,info%ik))
      resp%s_full(3,1,iband,:,info%ik) = conjg(resp%s_full(1,3,iband,:,info%ik))
      resp%s_full(3,2,iband,:,info%ik) = conjg(resp%s_full(2,3,iband,:,info%ik))

      resp%a_full(2,1,iband,:,info%ik) = conjg(resp%a_full(1,2,iband,:,info%ik))
      resp%a_full(3,1,iband,:,info%ik) = conjg(resp%a_full(1,3,iband,:,info%ik))
      resp%a_full(3,2,iband,:,info%ik) = conjg(resp%a_full(2,3,iband,:,info%ik))

      resp%x_full(2,1,iband,:,info%ik) = conjg(resp%x_full(1,2,iband,:,info%ik))
      resp%x_full(3,1,iband,:,info%ik) = conjg(resp%x_full(1,3,iband,:,info%ik))
      resp%x_full(3,2,iband,:,info%ik) = conjg(resp%x_full(2,3,iband,:,info%ik))
    enddo
  endif

  deallocate(enrgy)
  deallocate(enrgydiff)
  deallocate(gamdiff)

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

  ! NOTE: here we transpose it internally in Fortran
  ! so the output (hdf5 is in the correct order)
  index1 = (/1,2,3,2,3,3,2,3,3/)
  index2 = (/1,2,3,1,1,2,1,1,2/)

  do iband = edisp%nbopt_min, edisp%nbopt_max
    ! the kernels are saved in the 1 1 directions
    ! so we calculate the 1 1 component at the end
    do idir = 2,edisp%iOptical
      if (idir <= 6) then
        resp%s_full(index1(idir),index2(idir),iband,:,info%ik) = &
        resp%s_full(1,1,iband,:,info%ik) * edisp%MoptDiag(idir,iband,:,info%ik)

        resp%a_full(index1(idir),index2(idir),iband,:,info%ik) = &
        resp%a_full(1,1,iband,:,info%ik) * edisp%MoptDiag(idir,iband,:,info%ik)

        resp%x_full(index1(idir),index2(idir),iband,:,info%ik) = &
        resp%x_full(1,1,iband,:,info%ik) * edisp%MoptDiag(idir,iband,:,info%ik)
      else
        ! here we ADD the complex part to the response
        resp%s_full(index1(idir),index2(idir),iband,:,info%ik) = &
        resp%s_full(index1(idir),index2(idir),iband,:,info%ik) + &
        resp%s_full(1,1,iband,:,info%ik) * edisp%MoptDiag(idir,iband,:,info%ik) * ci

        resp%a_full(index1(idir),index2(idir),iband,:,info%ik) = &
        resp%a_full(index1(idir),index2(idir),iband,:,info%ik) + &
        resp%a_full(1,1,iband,:,info%ik) * edisp%MoptDiag(idir,iband,:,info%ik) * ci

        resp%x_full(index1(idir),index2(idir),iband,:,info%ik) = &
        resp%x_full(index1(idir),index2(idir),iband,:,info%ik) + &
        resp%x_full(1,1,iband,:,info%ik) * edisp%MoptDiag(idir,iband,:,info%ik) * ci
      endif
    enddo

    if (edisp%iOptical > 3) then ! 'symmetrize' the whole thing
      resp%s_full(1,2,iband,:,info%ik) = conjg(resp%s_full(2,1,iband,:,info%ik))
      resp%s_full(1,3,iband,:,info%ik) = conjg(resp%s_full(3,1,iband,:,info%ik))
      resp%s_full(2,3,iband,:,info%ik) = conjg(resp%s_full(3,2,iband,:,info%ik))

      resp%a_full(1,2,iband,:,info%ik) = conjg(resp%a_full(2,1,iband,:,info%ik))
      resp%a_full(1,3,iband,:,info%ik) = conjg(resp%a_full(3,1,iband,:,info%ik))
      resp%a_full(2,3,iband,:,info%ik) = conjg(resp%a_full(3,2,iband,:,info%ik))

      resp%x_full(1,2,iband,:,info%ik) = conjg(resp%x_full(2,1,iband,:,info%ik))
      resp%x_full(1,3,iband,:,info%ik) = conjg(resp%x_full(3,1,iband,:,info%ik))
      resp%x_full(2,3,iband,:,info%ik) = conjg(resp%x_full(3,2,iband,:,info%ik))
    endif


    resp%s_full(1,1,iband,:,info%ik) = &
    resp%s_full(1,1,iband,:,info%ik) * edisp%MoptDiag(1,iband,:,info%ik)

    resp%a_full(1,1,iband,:,info%ik) = &
    resp%a_full(1,1,iband,:,info%ik) * edisp%MoptDiag(1,iband,:,info%ik)

    resp%x_full(1,1,iband,:,info%ik) = &
    resp%x_full(1,1,iband,:,info%ik) * edisp%MoptDiag(1,iband,:,info%ik)
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
  ! TODO: FIX THIS

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

    if (myid .eq. master) then
      allocate(resp%x_gather(3,3,edisp%nbopt_min:edisp%nbopt_max,edisp%ispin,kmesh%nkp))
    else
      allocate(resp%x_gather(1,1,1,1,1))
    endif
#ifdef MPI
    call MPI_gatherv(resp%x_full,(ikend-ikstr+1)*9*(edisp%nbopt_max-edisp%nbopt_min+1)*edisp%ispin, &
                     MPI_DOUBLE_COMPLEX, resp%x_gather, rcounts*9*(edisp%nbopt_max-edisp%nbopt_min+1)*edisp%ispin, &
                     displs*9*(edisp%nbopt_max-edisp%nbopt_min+1)*edisp%ispin, MPI_DOUBLE_COMPLEX, master, &
                     MPI_COMM_WORLD, mpierr)
#else
    resp%x_gather = resp%x_full
#endif

    ! in order to keep the numerical values in a good range
    ! we do not apply the final multiplication with e
    ! this way we essentially get L0, L1/q, L2/q**2
    ! with which the derived quantities are properly defined
    ! e.g. peltier = L1/(q*L0)
    resp%s_gather = resp%s_gather * pi * ( echarge / (kmesh%vol*hbarevs)) * 1.d10 ! -> 1/(Ohm*m) = A / (V * m)
    resp%a_gather = resp%a_gather * pi * ( echarge / (kmesh%vol*hbarevs)) * 1.d10 ! -> A**2 * s / m    (we do not multiply with e here)
    resp%x_gather = resp%x_gather * pi * ( echarge / (kmesh%vol*hbarevs)) * 1.d10 ! -> V * A**3 * s**2 / m (we do not multiply with e**2 here)

    if (myid .eq. master) then
      write(string,'(I6.6)') info%iT
      string = trim(string) // "/L0/" // trim(adjustl(gname)) // "/full"
      call hdf5_write_data(ifile, string, resp%s_gather)

      write(string,'(I6.6)') info%iT
      string = trim(string) // "/L1/" // trim(adjustl(gname)) // "/full"
      call hdf5_write_data(ifile, string, resp%a_gather)

      write(string,'(I6.6)') info%iT
      string = trim(string) // "/L2/" // trim(adjustl(gname)) // "/full"
      call hdf5_write_data(ifile, string, resp%x_gather)
    endif

    deallocate(resp%s_gather)
    deallocate(resp%a_gather)
    deallocate(resp%x_gather)

  ! full output end ... this is always for each T-point
  ! arrays would be too large for this
  endif

  ! perform a local summation
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
    call MPI_REDUCE(MPI_IN_PLACE, resp%s_sum, 9*edisp%ispin, MPI_DOUBLE_COMPLEX, MPI_SUM, master, MPI_COMM_WORLD, mpierr)
  else
    call MPI_REDUCE(resp%s_sum, resp%s_sum, 9*edisp%ispin, MPI_DOUBLE_COMPLEX, MPI_SUM, master, MPI_COMM_WORLD, mpierr)
  endif

  if (myid.eq.master) then
    call MPI_REDUCE(MPI_IN_PLACE, resp%a_sum, 9*edisp%ispin, MPI_DOUBLE_COMPLEX, MPI_SUM, master, MPI_COMM_WORLD, mpierr)
  else
    call MPI_REDUCE(resp%a_sum, resp%a_sum, 9*edisp%ispin, MPI_DOUBLE_COMPLEX, MPI_SUM, master, MPI_COMM_WORLD, mpierr)
  endif
  if (myid.eq.master) then
    call MPI_REDUCE(MPI_IN_PLACE, resp%x_sum, 9*edisp%ispin, MPI_DOUBLE_COMPLEX, MPI_SUM, master, MPI_COMM_WORLD, mpierr)
  else
    call MPI_REDUCE(resp%x_sum, resp%x_sum, 9*edisp%ispin, MPI_DOUBLE_COMPLEX, MPI_SUM, master, MPI_COMM_WORLD, mpierr)
  endif
#endif

  resp%s_sum = resp%s_sum * pi * ( echarge / (kmesh%vol*hbarevs)) * 1.d10
  resp%a_sum = resp%a_sum * pi * ( echarge / (kmesh%vol*hbarevs)) * 1.d10
  resp%x_sum = resp%x_sum * pi * ( echarge / (kmesh%vol*hbarevs)) * 1.d10

  if (myid .eq. master) then
    if (algo%lDebug .and. (index(algo%dbgstr,"ReduceIO") .ne. 0)) then
      ! gather the data in the arrays
      resp%s_sum_temp(:,:,:,info%iT) = resp%s_sum
      resp%a_sum_temp(:,:,:,info%iT) = resp%a_sum
      resp%x_sum_temp(:,:,:,info%iT) = resp%x_sum

      ! output at the last temperature step
      if ((temp%Tstep==1 .and. info%iT==temp%nT) .or. (temp%Tstep==-1 .and. info%iT==1)) then
        string = "/L0/" // trim(adjustl(gname)) // "/sum"
        call hdf5_write_data(ifile, string, resp%s_sum_temp)
        string = "/L1/" // trim(adjustl(gname)) // "/sum"
        call hdf5_write_data(ifile, string, resp%a_sum_temp)
        string = "/L2/" // trim(adjustl(gname)) // "/sum"
        call hdf5_write_data(ifile, string, resp%x_sum_temp)
      endif
    else
        ! output it for each temperature point
      write(string,'(I6.6)') info%iT
      string = trim(string) // "/L0/" // trim(adjustl(gname)) // "/sum"
      call hdf5_write_data(ifile, string, resp%s_sum)

      write(string,'(I6.6)') info%iT
      string = trim(string) // "/L1/" // trim(adjustl(gname)) // "/sum"
      call hdf5_write_data(ifile, string, resp%a_sum)

      write(string,'(I6.6)') info%iT
      string = trim(string) // "/L2/" // trim(adjustl(gname)) // "/sum"
      call hdf5_write_data(ifile, string, resp%x_sum)
    endif
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

      if (myid .eq. master) then
        allocate(resp%xB_gather(3,3,edisp%nbopt_min:edisp%nbopt_max,edisp%ispin,kmesh%nkp))
      else
        allocate(resp%xB_gather(1,1,1,1,1))
      endif
#ifdef MPI
      call MPI_gatherv(resp%xB_full,(ikend-ikstr+1)*9*(edisp%nbopt_max-edisp%nbopt_min+1)*edisp%ispin, &
                       MPI_DOUBLE_COMPLEX, resp%xB_gather, rcounts*9*(edisp%nbopt_max-edisp%nbopt_min+1)*edisp%ispin, &
                       displs*9*(edisp%nbopt_max-edisp%nbopt_min+1)*edisp%ispin, MPI_DOUBLE_COMPLEX, master, &
                       MPI_COMM_WORLD, mpierr)
#else
      resp%xB_gather = resp%xB_full
#endif

      ! same story here ... do not multiply with e
      resp%sB_gather = resp%sB_gather * pi**2 * ( echarge / (kmesh%vol*hbarevs)) * (1.d-10 / hbarevs) ! -> A * m / (V**2 * s)
      resp%aB_gather = resp%aB_gather * pi**2 * ( echarge / (kmesh%vol*hbarevs)) * (1.d-10 / hbarevs) ! -> A**2 * m / V
      resp%xB_gather = resp%xB_gather * pi**2 * ( echarge / (kmesh%vol*hbarevs)) * (1.d-10 / hbarevs) ! -> A**3 * m * s

      if (myid .eq. master) then
        write(string,'(I6.6)') info%iT
        string = trim(string) // "/L0/" // trim(adjustl(gname)) // "/fullM"
        call hdf5_write_data(ifile, string, resp%sB_gather)

        write(string,'(I6.6)') info%iT
        string = trim(string) // "/L1/" // trim(adjustl(gname)) // "/fullM"
        call hdf5_write_data(ifile, string, resp%aB_gather)

        write(string,'(I6.6)') info%iT
        string = trim(string) // "/L2/" // trim(adjustl(gname)) // "/fullM"
        call hdf5_write_data(ifile, string, resp%xB_gather)

      endif

      deallocate(resp%sB_gather)
      deallocate(resp%aB_gather)
      deallocate(resp%xB_gather)
    endif ! full output

    ! perform a local summation
    ! these are already initialized to 0
    do ik = ikstr,ikend
      do iband = edisp%nbopt_min,edisp%nbopt_max
        resp%sB_sum(:,:,:) = resp%sB_sum(:,:,:) + resp%sB_full(:,:,iband,:,ik) * kmesh%weight(ik)
        resp%aB_sum(:,:,:) = resp%aB_sum(:,:,:) + resp%aB_full(:,:,iband,:,ik) * kmesh%weight(ik)
        resp%xB_sum(:,:,:) = resp%xB_sum(:,:,:) + resp%xB_full(:,:,iband,:,ik) * kmesh%weight(ik)
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

  if (myid.eq.master) then
    call MPI_REDUCE(MPI_IN_PLACE, resp%xB_sum, 9*edisp%ispin, MPI_DOUBLE_COMPLEX, MPI_SUM, master, MPI_COMM_WORLD, mpierr)
  else
    call MPI_REDUCE(resp%xB_sum, resp%xB_sum, 9*edisp%ispin, MPI_DOUBLE_COMPLEX, MPI_SUM, master, MPI_COMM_WORLD, mpierr)
  endif
#endif

    resp%sB_sum = resp%sB_sum * pi**2 * ( echarge / (kmesh%vol*hbarevs)) * (1.d-10 / hbarevs)
    resp%aB_sum = resp%aB_sum * pi**2 * ( echarge / (kmesh%vol*hbarevs)) * (1.d-10 * echarge / hbarevs)
    resp%xB_sum = resp%xB_sum * pi**2 * ( echarge / (kmesh%vol*hbarevs)) * (1.d-10 * echarge**2 / hbarevs)

    if (myid .eq. master) then

      if (algo%lDebug .and. (index(algo%dbgstr,"ReduceIO") .ne. 0)) then
        ! gather the data in the arrays
        resp%sB_sum_temp(:,:,:,info%iT) = resp%sB_sum
        resp%aB_sum_temp(:,:,:,info%iT) = resp%aB_sum
        resp%xB_sum_temp(:,:,:,info%iT) = resp%xB_sum

        ! output at the last temperature step
        if ((temp%Tstep==1 .and. info%iT==temp%nT) .or. (temp%Tstep==-1 .and. info%iT==1)) then
          string = "/L0/" // trim(adjustl(gname)) // "/sumM"
          call hdf5_write_data(ifile, string, resp%sB_sum_temp)
          string = "/L1/" // trim(adjustl(gname)) // "/sumM"
          call hdf5_write_data(ifile, string, resp%aB_sum_temp)
          string = "/L2/" // trim(adjustl(gname)) // "/sumM"
          call hdf5_write_data(ifile, string, resp%xB_sum_temp)
        endif
      else
        write(string,'(I6.6)') info%iT
        string = trim(string) // "/L0/" // trim(adjustl(gname)) // "/sumM"
        call hdf5_write_data(ifile, string, resp%sB_sum)

        write(string,'(I6.6)') info%iT
        string = trim(string) // "/L1/" // trim(adjustl(gname)) // "/sumM"
        call hdf5_write_data(ifile, string, resp%aB_sum)

        write(string,'(I6.6)') info%iT
        string = trim(string) // "/L2/" // trim(adjustl(gname)) // "/sumM"
        call hdf5_write_data(ifile, string, resp%xB_sum)
      endif
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

subroutine dpresp_alloc(algo, edisp, temp, dpresp)
  implicit none
  type(algorithm)   :: algo
  type(energydisp)  :: edisp
  type(temperature) :: temp
  type(response_dp) :: dpresp

  ! allocate transport variables
  allocate(dpresp%s_full(3,3,edisp%nbopt_min:edisp%nbopt_max,edisp%ispin,ikstr:ikend))
  allocate(dpresp%a_full(3,3,edisp%nbopt_min:edisp%nbopt_max,edisp%ispin,ikstr:ikend))
  allocate(dpresp%x_full(3,3,edisp%nbopt_min:edisp%nbopt_max,edisp%ispin,ikstr:ikend))
  allocate(dpresp%s_sum(3,3,edisp%iSpin))
  allocate(dpresp%a_sum(3,3,edisp%iSpin))
  allocate(dpresp%x_sum(3,3,edisp%iSpin))

  if (myid.eq.master .and. algo%lDebug .and. (index(algo%dbgstr,"ReduceIO") .ne. 0)) then
    allocate(dpresp%s_sum_temp(3,3,edisp%iSpin,temp%nT))
    allocate(dpresp%a_sum_temp(3,3,edisp%iSpin,temp%nT))
    allocate(dpresp%x_sum_temp(3,3,edisp%iSpin,temp%nT))
  endif


  if (algo%lBfield) then
    allocate(dpresp%sB_full(3,3,edisp%nbopt_min:edisp%nbopt_max,edisp%ispin,ikstr:ikend))
    allocate(dpresp%aB_full(3,3,edisp%nbopt_min:edisp%nbopt_max,edisp%ispin,ikstr:ikend))
    allocate(dpresp%xB_full(3,3,edisp%nbopt_min:edisp%nbopt_max,edisp%ispin,ikstr:ikend))
    allocate(dpresp%sB_sum(3,3,edisp%iSpin))
    allocate(dpresp%aB_sum(3,3,edisp%iSpin))
    allocate(dpresp%xB_sum(3,3,edisp%iSpin))

    if (myid.eq.master .and. algo%lDebug .and. (index(algo%dbgstr,"ReduceIO") .ne. 0)) then
      allocate(dpresp%sB_sum_temp(3,3,edisp%iSpin,temp%nT))
      allocate(dpresp%aB_sum_temp(3,3,edisp%iSpin,temp%nT))
      allocate(dpresp%xB_sum_temp(3,3,edisp%iSpin,temp%nT))
    endif
  endif

end subroutine dpresp_alloc

subroutine qpresp_alloc(algo, edisp, temp, qpresp)
  implicit none
  type(algorithm)   :: algo
  type(energydisp)  :: edisp
  type(temperature) :: temp
  type(response_qp) :: qpresp

  ! allocate transport variables
  allocate(qpresp%s_full(3,3,edisp%nbopt_min:edisp%nbopt_max,edisp%ispin,ikstr:ikend))
  allocate(qpresp%a_full(3,3,edisp%nbopt_min:edisp%nbopt_max,edisp%ispin,ikstr:ikend))
  allocate(qpresp%x_full(3,3,edisp%nbopt_min:edisp%nbopt_max,edisp%ispin,ikstr:ikend))
  allocate(qpresp%s_sum(3,3,edisp%iSpin))
  allocate(qpresp%a_sum(3,3,edisp%iSpin))
  allocate(qpresp%x_sum(3,3,edisp%iSpin))

  if (myid .eq. master .and. algo%lDebug .and. (index(algo%dbgstr,"ReduceIO") .ne. 0)) then
    allocate(qpresp%s_sum_temp(3,3,edisp%iSpin,temp%nT))
    allocate(qpresp%a_sum_temp(3,3,edisp%iSpin,temp%nT))
    allocate(qpresp%x_sum_temp(3,3,edisp%iSpin,temp%nT))
  endif

  if (algo%lBfield) then
    allocate(qpresp%sB_full(3,3,edisp%nbopt_min:edisp%nbopt_max,edisp%ispin,ikstr:ikend))
    allocate(qpresp%aB_full(3,3,edisp%nbopt_min:edisp%nbopt_max,edisp%ispin,ikstr:ikend))
    allocate(qpresp%xB_full(3,3,edisp%nbopt_min:edisp%nbopt_max,edisp%ispin,ikstr:ikend))
    allocate(qpresp%sB_sum(3,3,edisp%iSpin))
    allocate(qpresp%aB_sum(3,3,edisp%iSpin))
    allocate(qpresp%xB_sum(3,3,edisp%iSpin))

    if (myid .eq. master .and. algo%lDebug .and. (index(algo%dbgstr,"ReduceIO") .ne. 0)) then
      allocate(qpresp%sB_sum_temp(3,3,edisp%iSpin,temp%nT))
      allocate(qpresp%aB_sum_temp(3,3,edisp%iSpin,temp%nT))
      allocate(qpresp%xB_sum_temp(3,3,edisp%iSpin,temp%nT))
    endif
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

  to_evaluate = 0.5d0 + info%beta2p * &
                (sct%gam(edisp%nbopt_min:edisp%nbopt_max,ikstr:ikend,:) &
                 + ci*sct%zqp(edisp%nbopt_min:edisp%nbopt_max,ikstr:ikend,:) &
                 * (edisp%band(edisp%nbopt_min:edisp%nbopt_max,ikstr:ikend,:) - mu))

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

  to_evaluate = 0.5q0 + info%beta2pQ * &
                (sct%gam(edisp%nbopt_min:edisp%nbopt_max,ikstr:ikend,:) &
                 + ciQ*sct%zqp(edisp%nbopt_min:edisp%nbopt_max,ikstr:ikend,:) &
                 * (edisp%band(edisp%nbopt_min:edisp%nbopt_max,ikstr:ikend,:) - mu))

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

subroutine calc_total_energy_digamma(mu, energy_tot, edisp, sct, kmesh, imp, algo, info)
  real(8), intent(in)  :: mu
  real(8), intent(out) :: energy_tot

  type(energydisp) :: edisp
  type(scattering) :: sct
  type(kpointmesh) :: kmesh
  type(impurity)   :: imp
  type(algorithm)  :: algo
  type(runinfo)    :: info
  !local variables

  real(8) :: energy_loc
  integer :: is, ik, iband, iimp
  complex(8), allocatable :: to_evaluate(:,:,:)
  real(8), allocatable    :: energy_post_factor(:,:,:)
  real(8), allocatable    :: energy(:,:,:)
  !external variables
  complex(8), external :: wpsipg

  allocate(to_evaluate(edisp%nband_max, ikstr:ikend, edisp%ispin))
  allocate(energy_post_factor(edisp%nband_max, ikstr:ikend, edisp%ispin))
  allocate(energy(edisp%nband_max, ikstr:ikend, edisp%ispin))

  ! for the occupation
  to_evaluate = 0.5d0 + info%beta2p * &
                (sct%gam(:,ikstr:ikend,:) - ci*sct%zqp(:,ikstr:ikend,:)*(edisp%band(:,ikstr:ikend,:) - mu))
  ! energy contribution : occupation * energy of this occupation
  energy_post_factor = sct%zqp(:,ikstr:ikend,:) * (edisp%band(:,ikstr:ikend,:) - mu)

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

  ! if we have impurity levels add their energy contribution here
  if (algo%lImpurities) then
    do iimp = 1,imp%nimp
      energy_tot = energy_tot - imp%Dopant(iimp)*imp%Density(iimp) &
        / (1.d0 + imp%Degeneracy(iimp) * exp(info%beta*imp%Dopant(iimp)*(mu-imp%Energy(iimp)))) &
        * imp%Energy(iimp)
    enddo
  endif

end subroutine

subroutine calc_total_energy_fermi(mu, energy_tot, edisp, sct, kmesh, imp, algo, info)
  real(8), intent(in)  :: mu
  real(8), intent(out) :: energy_tot

  type(energydisp) :: edisp
  type(scattering) :: sct
  type(kpointmesh) :: kmesh
  type(impurity)   :: imp
  type(algorithm)  :: algo
  type(runinfo)    :: info
  !local variables

  real(8) :: energy_loc
  integer :: is, ik, iband, iimp
  real(8), allocatable    :: energy_post_factor(:,:,:)
  real(8), allocatable    :: energy(:,:,:)
  !external variables
  complex(8), external :: wpsipg

  allocate(energy_post_factor(edisp%nband_max, ikstr:ikend, edisp%ispin))
  allocate(energy(edisp%nband_max, ikstr:ikend, edisp%ispin))

  ! energy contribution : occupation * energy of this occupation
  energy_post_factor = sct%zqp(:,ikstr:ikend,:) * (edisp%band(:,ikstr:ikend,:) - mu)

  ! evaluate the function
  do is = 1,edisp%ispin
    do ik = ikstr, ikend
      do iband=1,edisp%nband_max
        energy(iband,ik,is) = fermi_dp(sct%zqp(iband,ik,is)*(edisp%band(iband,ik,is)-mu), info%beta)   ! occupation
        energy(iband,ik,is) = energy(iband,ik,is) * kmesh%weight(ik) * energy_post_factor(iband,ik,is) ! multiplied with weight and energy gives the energy
      enddo
    enddo
  enddo

  deallocate(energy_post_factor)
  energy_loc = sum(energy)
  deallocate(energy)

#ifdef MPI
  call MPI_ALLREDUCE(energy_loc, energy_tot, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, mpierr)
#else
  energy_tot = energy_loc
#endif

  ! if we have impurity levels add their energy contribution here
  if (algo%lImpurities) then
    do iimp = 1,imp%nimp
      energy_tot = energy_tot - imp%Dopant(iimp)*imp%Density(iimp) &
        / (1.d0 + imp%Degeneracy(iimp) * exp(info%beta*imp%Dopant(iimp)*(mu-imp%Energy(iimp)))) &
        * imp%Energy(iimp)
    enddo
  endif

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

subroutine response_intra_km_Q(resp, PolyGamma, mu, edisp, sct, kmesh, algo, info)
  implicit none
  real(8), intent(in) :: mu
  type (response_qp)  :: resp

  type(energydisp)    :: edisp
  type(scattering)    :: sct
  type(kpointmesh)    :: kmesh
  type(algorithm)     :: algo
  type(runinfo)       :: info

  complex(16)          :: PolyGamma(3,edisp%nbopt_min:edisp%nbopt_max, ikstr:ikend, edisp%ispin)

  real(8), allocatable :: enrgy(:,:)

  allocate(enrgy(edisp%nbopt_min:edisp%nbopt_max,edisp%ispin))
  ! first we write the kernel into the 1 1 component
  enrgy = sct%zqp(edisp%nbopt_min:edisp%nbopt_max,info%ik,:) * (edisp%band(edisp%nbopt_min:edisp%nbopt_max,info%ik,:) - mu)

  resp%s_full(1,1,:,:,info%ik) = real(PolyGamma(1,:,info%ik,:)) &
                               - real(PolyGamma(2,:,info%ik,:)) * info%beta2pQ*sct%gam(edisp%nbopt_min:edisp%nbopt_max,info%ik,:)

  resp%s_full(1,1,:,:,info%ik) = resp%s_full(1,1,:,:,info%ik) &
                               * sct%zqp(edisp%nbopt_min:edisp%nbopt_max,info%ik,:)**2 * info%betaQ &
                               / (4.d0 * piQ**3 * sct%gam(edisp%nbopt_min:edisp%nbopt_max,info%ik,:))


  resp%a_full(1,1,:,:,info%ik) =  real(PolyGamma(1,:,info%ik,:)) * enrgy &
                               -  real(PolyGamma(2,:,info%ik,:)) * enrgy &
                                  * sct%gam(edisp%nbopt_min:edisp%nbopt_max,info%ik,:) * info%beta2pQ &
                               - aimag(PolyGamma(2,:,info%ik,:)) * sct%gam(edisp%nbopt_min:edisp%nbopt_max,info%ik,:)**2 &
                                 * info%beta2pQ

  resp%a_full(1,1,:,:,info%ik) = resp%a_full(1,1,:,:,info%ik) &
                               * sct%zqp(edisp%nbopt_min:edisp%nbopt_max,info%ik,:)**2 * info%betaQ &
                               / (4.d0 * piQ**3 * sct%gam(edisp%nbopt_min:edisp%nbopt_max,info%ik,:))


  resp%x_full(1,1,:,:,info%ik) =  real(PolyGamma(1,:,info%ik,:)) &
                                  * (enrgy**2 + sct%gam(edisp%nbopt_min:edisp%nbopt_max,info%ik,:)**2) &
                               +  real(PolyGamma(2,:,info%ik,:)) &
                                  * info%beta2pQ * sct%gam(edisp%nbopt_min:edisp%nbopt_max,info%ik,:) &
                                  * (sct%gam(edisp%nbopt_min:edisp%nbopt_max,info%ik,:)**2 - enrgy**2) &
                               - aimag(PolyGamma(2,:,info%ik,:)) &
                                  * info%betaQ / piQ * enrgy * sct%gam(edisp%nbopt_min:edisp%nbopt_max,info%ik,:)**2

  resp%x_full(1,1,:,:,info%ik) = resp%x_full(1,1,:,:,info%ik) &
                               * sct%zqp(edisp%nbopt_min:edisp%nbopt_max,info%ik,:)**2 * info%betaQ &
                               / (4.d0 * piQ**3 * sct%gam(edisp%nbopt_min:edisp%nbopt_max,info%ik,:))

  if (algo%lBfield) then

    resp%sB_full(1,1,:,:,info%ik) = real(PolyGamma(3,:,info%ik,:)) &
                                      * sct%gam(edisp%nbopt_min:edisp%nbopt_max,info%ik,:)**2 * info%betaQ**2 / (4.d0 * piQ**2) &
                                  - real(PolyGamma(2,:,info%ik,:)) &
                                    * 3.d0 * sct%gam(edisp%nbopt_min:edisp%nbopt_max,info%ik,:) * info%beta2pQ &
                                  - real(PolyGamma(1,:,info%ik,:)) * 3.d0

    resp%sB_full(1,1,:,:,info%ik) = resp%sB_full(1,1,:,:,info%ik) &
                                  * sct%zqp(edisp%nbopt_min:edisp%nbopt_max,info%ik,:)**3 * info%betaQ &
                                  / (16.d0 * piQ**4 * sct%gam(edisp%nbopt_min:edisp%nbopt_max,info%ik,:)**2)


    resp%aB_full(1,1,:,:,info%ik) =  real(PolyGamma(3,:,info%ik,:)) &
                                       * enrgy * sct%gam(edisp%nbopt_min:edisp%nbopt_max,info%ik,:)**2 &
                                       * info%betaQ**2 / (4.d0 * piQ**2) &
                                  + aimag(PolyGamma(3,:,info%ik,:)) &
                                       * sct%gam(edisp%nbopt_min:edisp%nbopt_max,info%ik,:)**3 * info%betaQ**3 / (4.d0 * piQ**2) &
                                  -  real(PolyGamma(2,:,info%ik,:)) &
                                       * 3.d0 * enrgy * sct%gam(edisp%nbopt_min:edisp%nbopt_max,info%ik,:) * info%beta2pQ &
                                  - aimag(PolyGamma(2,:,info%ik,:)) &
                                       * sct%gam(edisp%nbopt_min:edisp%nbopt_max,info%ik,:)**2 * info%beta2pQ &
                                  +  real(PolyGamma(1,:,info%ik,:)) * 3.d0 * enrgy

    resp%aB_full(1,1,:,:,info%ik) = resp%aB_full(1,1,:,:,info%ik) &
                                  * sct%zqp(edisp%nbopt_min:edisp%nbopt_max,info%ik,:)**3 * info%betaQ / (16.d0 * piQ**4 &
                                  * sct%gam(edisp%nbopt_min:edisp%nbopt_max,info%ik,:)**2)


    resp%xB_full(1,1,:,:,info%ik) =  real(PolyGamma(3,:,info%ik,:)) &
                                       * info%betaQ**2 * sct%gam(edisp%nbopt_min:edisp%nbopt_max,info%ik,:)**2 &
                                       * (sct%gam(edisp%nbopt_min:edisp%nbopt_max,info%ik,:)**2 - enrgy**2) &
                                       / (4.d0 * piQ**2) &
                                  - aimag(PolyGamma(3,:,info%ik,:)) &
                                       * info%betaQ**2 * enrgy * sct%gam(edisp%nbopt_min:edisp%nbopt_max,info%ik,:)**3 &
                                       / (2.d0 * piQ) &
                                  -  real(PolyGamma(2,:,info%ik,:)) &
                                       * info%betaQ * sct%gam(edisp%nbopt_min:edisp%nbopt_max,info%ik,:) &
                                       * (sct%gam(edisp%nbopt_min:edisp%nbopt_max,info%ik,:)**2 + 3.d0*enrgy**2) &
                                       / (2.d0 * piQ) &
                                  + aimag(PolyGamma(2,:,info%ik,:)) &
                                       * info%betaQ * enrgy * sct%gam(edisp%nbopt_min:edisp%nbopt_max,info%ik,:)**2 / piQ &
                                  -  real(PolyGamma(1,:,info%ik,:)) &
                                       * (sct%gam(edisp%nbopt_min:edisp%nbopt_max,info%ik,:)**2 + 3.d0*enrgy**2)

    resp%xB_full(1,1,:,:,info%ik) = resp%xB_full(1,1,:,:,info%ik) &
                                  * (-1.d0) * sct%zqp(edisp%nbopt_min:edisp%nbopt_max,info%ik,:)**3 * info%betaQ &
                                  / (16.d0 * piQ**4 * sct%gam(edisp%nbopt_min:edisp%nbopt_max,info%ik,:)**2)

  endif

  deallocate(enrgy)

  call response_intra_optical_weights_Q(resp, edisp, info)
  ! if (algo%lBfield) then
  !   call response_peierls_weights(resp, edisp, info)
  ! endif


end subroutine response_intra_km_Q

subroutine response_inter_km_Q(resp, PolyGamma, mu, edisp, sct, kmesh, algo, info)
  implicit none
  real(8), intent(in) :: mu
  type (response_qp)  :: resp

  type(energydisp)    :: edisp
  type(scattering)    :: sct
  type(kpointmesh)    :: kmesh
  type(algorithm)     :: algo
  type(runinfo)       :: info

  complex(16)         :: PolyGamma(3,edisp%nbopt_min:edisp%nbopt_max, ikstr:ikend, edisp%ispin)

  ! these stay double precision .. because we read them from files
  real(8), allocatable :: enrgy(:,:)
  real(8), allocatable :: enrgydiff(:)
  real(8), allocatable :: gamdiff(:)

  complex(16) :: calc_sigma
  complex(16) :: calc_alpha
  complex(16) :: calc_xi

  integer :: index1(9), index2(9)
  integer :: i,j,idir
  integer :: iband1, iband2, iband, is

  ! index1 = (/1,2,3,1,1,2,1,1,2/)
  ! index2 = (/1,2,3,2,3,3,2,3,3/)

  ! NOTE: here we transpose it internally in Fortran
  ! so the output (hdf5 is in the correct order)
  index1 = (/1,2,3,2,3,3,2,3,3/)
  index2 = (/1,2,3,1,1,2,1,1,2/)


  allocate(enrgy(edisp%nbopt_min:edisp%nbopt_max,edisp%ispin))
  allocate(enrgydiff(edisp%ispin))

  ! first we write the kernel into the 1 1 component
  enrgy = sct%zqp(edisp%nbopt_min:edisp%nbopt_max,info%ik,:) &
          * (edisp%band(edisp%nbopt_min:edisp%nbopt_max,info%ik,:) - mu)

  do iband1 = edisp%nbopt_min, edisp%nbopt_max
    do iband2 = edisp%nbopt_min, edisp%nbopt_max
      if (iband1 == iband2) cycle
      enrgydiff = enrgy(iband1,:) - enrgy(iband2,:)
      gamdiff   = sct%gam(iband1,info%ik,:) - sct%gam(iband2,info%ik,:)

      do is = 1,edisp%ispin

        if ((abs(enrgydiff(is)) .lt. 1d-6) .and. (abs(gamdiff(is)) .lt. 1d-6)) then

          ! use the intra-band limit .....
          calc_sigma  = real(PolyGamma(1,iband1,info%ik,is)) &
                      - real(PolyGamma(2,iband1,info%ik,is)) * info%beta2pQ*sct%gam(iband1,info%ik,is)

          calc_sigma = calc_sigma &
                     * sct%zqp(iband1,info%ik,is)**2 * info%betaQ &
                     / (4.d0 * piQ**3 * sct%gam(iband1,info%ik,is))

          calc_alpha =  real(PolyGamma(1,iband1,info%ik,is)) * enrgy(iband1,is) &
                     -  real(PolyGamma(2,iband1,info%ik,is)) * enrgy(iband1,is) * sct%gam(iband1,info%ik,is) * info%beta2pQ &
                     - aimag(PolyGamma(2,iband1,info%ik,is)) * sct%gam(iband1,info%ik,is)**2 * info%beta2pQ

          calc_alpha = calc_alpha &
                     * sct%zqp(iband1,info%ik,is)**2 * info%betaQ &
                     / (4.d0 * piQ**3 * sct%gam(iband1,info%ik,is))

          calc_xi    =  real(PolyGamma(1,iband1,info%ik,is)) &
                        * (enrgy(iband1,is)**2 + sct%gam(iband1,info%ik,is)**2) &
                     +  real(PolyGamma(2,iband1,info%ik,is)) &
                        * info%beta2pQ * sct%gam(iband1,info%ik,is) * (sct%gam(iband1,info%ik,is)**2 - enrgy(iband1,is)**2) &
                     - aimag(PolyGamma(2,iband1,info%ik,is)) &
                        * info%betaQ / piQ * enrgy(iband1,is) * sct%gam(iband1,info%ik,is)**2

          calc_xi    = calc_xi &
                       * sct%zqp(iband1,info%ik,is)**2 * info%betaQ &
                       / (4.d0 * piQ**3 * sct%gam(iband1,info%ik,is))

        else
          calc_sigma = real(PolyGamma(1,iband1,info%ik,is)) &
                     * (enrgydiff(is)**2 + sct%gam(iband2,info%ik,is)**2 - sct%gam(iband1,info%ik,is)**2) &
                     / sct%gam(iband1,info%ik,is)

          calc_sigma = calc_sigma &
                     + aimag(PolyGamma(1,iband1,info%ik,is)) &
                     * (2.d0* (-enrgydiff(is)))

          calc_sigma = calc_sigma &
                     + real(PolyGamma(1,iband2,info%ik,is)) &
                     * (enrgydiff(is)**2 + sct%gam(iband1,info%ik,is)**2 - sct%gam(iband2,info%ik,is)**2) &
                     / sct%gam(iband2,info%ik,is)

          calc_sigma = calc_sigma &
                     + aimag(PolyGamma(1,iband2,info%ik,is)) &
                     * (2.d0* enrgydiff(is))

          calc_sigma = calc_sigma &
                     * sct%gam(iband1,info%ik,is) * sct%gam(iband2,info%ik,is) &
                     * sct%zqp(iband1,info%ik,is) * sct%zqp(iband2,info%ik,is) &
                     * info%betaQ

          calc_sigma = calc_sigma &
                     / (2.d0 * piQ**3 * ( enrgydiff(is)**2 + (sct%gam(iband1,info%ik,is) - sct%gam(iband2,info%ik,is))**2)) &
                     / ( enrgydiff(is)**2 + (sct%gam(iband1,info%ik,is) + sct%gam(iband2,info%ik,is))**2)



          calc_alpha = real(PolyGamma(1,iband1,info%ik,is)) &
                     * (enrgy(iband1,is) * (enrgydiff(is)**2 + sct%gam(iband2,info%ik,is)**2 - sct%gam(iband1,info%ik,is)**2) &
                     / sct%gam(iband1,info%ik,is) + 2.d0*sct%gam(iband1,info%ik,is)*enrgydiff(is))

          calc_alpha = calc_alpha &
                     + aimag(PolyGamma(1,iband1,info%ik,is)) &
                     * (enrgydiff(is)**2 + sct%gam(iband2,info%ik,is)**2 - sct%gam(iband1,info%ik,is)**2 &
                     - 2.d0*enrgy(iband1,is)*enrgydiff(is))

          calc_alpha = calc_alpha &
                     + real(PolyGamma(1,iband2,info%ik,is)) &
                     * (enrgy(iband2,is) * (enrgydiff(is)**2 + sct%gam(iband1,info%ik,is)**2 - sct%gam(iband2,info%ik,is)**2) &
                     / sct%gam(iband2,info%ik,is) + 2.d0*sct%gam(iband2,info%ik,is)*enrgydiff(is)*(-1.d0))

          calc_alpha = calc_alpha &
                     + aimag(PolyGamma(1,iband2,info%ik,is)) &
                     * (enrgydiff(is)**2 + sct%gam(iband1,info%ik,is)**2 - sct%gam(iband2,info%ik,is)**2 &
                     - 2.d0*enrgy(iband2,is)*enrgydiff(is)*(-1.d0))

          calc_alpha = calc_alpha &
                     * sct%gam(iband1,info%ik,is) * sct%gam(iband2,info%ik,is) &
                     * sct%zqp(iband1,info%ik,is) * sct%zqp(iband2,info%ik,is) &
                     * info%betaQ

          calc_alpha = calc_alpha &
                     / (2.d0 * piQ**3 * ( enrgydiff(is)**2 + (sct%gam(iband1,info%ik,is) - sct%gam(iband2,info%ik,is))**2)) &
                     / ( enrgydiff(is)**2 + (sct%gam(iband1,info%ik,is) + sct%gam(iband2,info%ik,is))**2)


          calc_xi    = real(PolyGamma(1,iband1,info%ik,is) &
                     * (enrgy(iband1,is) - ciQ*sct%gam(iband1,info%ik,is))**2 &
                     * (enrgydiff(is) + ciQ*gamdiff(is)) &
                     * (enrgydiff(is) + ciQ*(sct%gam(iband1,info%ik,is) + sct%gam(iband2,info%ik,is)))) &
                     * sct%zqp(iband1,info%ik,is) * sct%zqp(iband2,info%ik,is) * sct%gam(iband2,info%ik,is) &
                     * info%betaQ / (2.d0 * piQ**3)

          calc_xi    = calc_xi &
                     + real(PolyGamma(1,iband2,info%ik,is) &
                     * (enrgy(iband2,is) - ciQ*sct%gam(iband2,info%ik,is))**2 &
                     * (-enrgydiff(is) - ciQ*gamdiff(is)) &
                     * (-enrgydiff(is) + ciQ*(sct%gam(iband1,info%ik,is) + sct%gam(iband2,info%ik,is)))) &
                     * sct%zqp(iband1,info%ik,is) * sct%zqp(iband2,info%ik,is) * sct%gam(iband1,info%ik,is) &
                     * info%betaQ / (2.d0 * piQ**3)

          calc_xi    = calc_xi &
                     / ((enrgydiff(is)**2 + gamdiff(is)**2) &
                      * (enrgydiff(is)**2 + (sct%gam(iband1,info%ik,is) + sct%gam(iband2,info%ik,is))**2))
        endif


        ! multiply optical elements
        do idir = 1,edisp%iOptical
          if (idir <= 6) then
            ! ATTENTION
            ! we read the wien2k files via 1 - 2
            ! we save those in hdf5
            ! and read them via Fortran -> implicit transposition
            ! we have to use 2 - 1 here
            resp%s_full(index1(idir),index2(idir),iband1,is,info%ik) = &
            resp%s_full(index1(idir),index2(idir),iband1,is,info%ik) + calc_sigma * edisp%Mopt(idir,iband2,iband1,is,info%ik)

            resp%a_full(index1(idir),index2(idir),iband1,is,info%ik) = &
            resp%a_full(index1(idir),index2(idir),iband1,is,info%ik) + calc_alpha * edisp%Mopt(idir,iband2,iband1,is,info%ik)

            resp%x_full(index1(idir),index2(idir),iband1,is,info%ik) = &
            resp%x_full(index1(idir),index2(idir),iband1,is,info%ik) + calc_xi    * edisp%Mopt(idir,iband2,iband1,is,info%ik)
          else
            ! here we ADD the complex part to the response
            resp%s_full(index1(idir),index2(idir),iband1,is,info%ik) = &
            resp%s_full(index1(idir),index2(idir),iband1,is,info%ik) + calc_sigma * edisp%Mopt(idir,iband2,iband1,is,info%ik) * ciQ

            resp%a_full(index1(idir),index2(idir),iband1,is,info%ik) = &
            resp%a_full(index1(idir),index2(idir),iband1,is,info%ik) + calc_alpha * edisp%Mopt(idir,iband2,iband1,is,info%ik) * ciQ

            resp%x_full(index1(idir),index2(idir),iband1,is,info%ik) = &
            resp%x_full(index1(idir),index2(idir),iband1,is,info%ik) + calc_xi     * edisp%Mopt(idir,iband2,iband1,is,info%ik) * ciQ
          endif
        enddo
      enddo
    enddo
  enddo

  ! TODO: FIX THIS
  ! i.e. we save all this data to python
  ! that means an implicit transposition takes place
  ! i.e. we either hav eto treat all of this the other way around
  if (edisp%iOptical > 3) then
    do iband = edisp%nbopt_min, edisp%nbopt_max
      resp%s_full(2,1,iband,:,info%ik) = conjg(resp%s_full(1,2,iband,:,info%ik))
      resp%s_full(3,1,iband,:,info%ik) = conjg(resp%s_full(1,3,iband,:,info%ik))
      resp%s_full(3,2,iband,:,info%ik) = conjg(resp%s_full(2,3,iband,:,info%ik))

      resp%a_full(2,1,iband,:,info%ik) = conjg(resp%a_full(1,2,iband,:,info%ik))
      resp%a_full(3,1,iband,:,info%ik) = conjg(resp%a_full(1,3,iband,:,info%ik))
      resp%a_full(3,2,iband,:,info%ik) = conjg(resp%a_full(2,3,iband,:,info%ik))

      resp%x_full(2,1,iband,:,info%ik) = conjg(resp%x_full(1,2,iband,:,info%ik))
      resp%x_full(3,1,iband,:,info%ik) = conjg(resp%x_full(1,3,iband,:,info%ik))
      resp%x_full(3,2,iband,:,info%ik) = conjg(resp%x_full(2,3,iband,:,info%ik))
    enddo
  endif

  deallocate(enrgy)
  deallocate(enrgydiff)
  deallocate(gamdiff)

end subroutine response_inter_km_Q

subroutine response_intra_optical_weights_Q(resp, edisp, info)
  implicit none
  type (response_qp) :: resp

  type(energydisp)   :: edisp
  type(runinfo)      :: info

  integer :: index1(9), index2(9)
  integer :: iband, idir

  !( 1 4+i7 5+i8 )
  !( - 2    6+i9 )
  !( - -    3    )

  ! we use these two index lists to move along the described order in the 3x3 matrix

  ! NOTE: here we transpose it internally in Fortran
  ! so the output (hdf5 is in the correct order)
  index1 = (/1,2,3,2,3,3,2,3,3/)
  index2 = (/1,2,3,1,1,2,1,1,2/)

  do iband = edisp%nbopt_min, edisp%nbopt_max
    ! the kernels are saved in the 1 1 directions
    ! so we calculate the 1 1 component at the end
    do idir = 2,edisp%iOptical
      if (idir <= 6) then
        resp%s_full(index1(idir),index2(idir),iband,:,info%ik) = &
        resp%s_full(1,1,iband,:,info%ik) * edisp%MoptDiag(idir,iband,:,info%ik)

        resp%a_full(index1(idir),index2(idir),iband,:,info%ik) = &
        resp%a_full(1,1,iband,:,info%ik) * edisp%MoptDiag(idir,iband,:,info%ik)

        resp%x_full(index1(idir),index2(idir),iband,:,info%ik) = &
        resp%x_full(1,1,iband,:,info%ik) * edisp%MoptDiag(idir,iband,:,info%ik)
      else
        ! here we ADD the complex part to the response
        resp%s_full(index1(idir),index2(idir),iband,:,info%ik) = &
        resp%s_full(index1(idir),index2(idir),iband,:,info%ik) + &
        resp%s_full(1,1,iband,:,info%ik) * edisp%MoptDiag(idir,iband,:,info%ik) * ciQ

        resp%a_full(index1(idir),index2(idir),iband,:,info%ik) = &
        resp%a_full(index1(idir),index2(idir),iband,:,info%ik) + &
        resp%a_full(1,1,iband,:,info%ik) * edisp%MoptDiag(idir,iband,:,info%ik) * ciQ

        resp%x_full(index1(idir),index2(idir),iband,:,info%ik) = &
        resp%x_full(index1(idir),index2(idir),iband,:,info%ik) + &
        resp%x_full(1,1,iband,:,info%ik) * edisp%MoptDiag(idir,iband,:,info%ik) * ciQ
      endif
    enddo

    if (edisp%iOptical > 3) then ! 'symmetrize' the whole thing
      resp%s_full(1,2,iband,:,info%ik) = conjg(resp%s_full(2,1,iband,:,info%ik))
      resp%s_full(1,3,iband,:,info%ik) = conjg(resp%s_full(3,1,iband,:,info%ik))
      resp%s_full(2,3,iband,:,info%ik) = conjg(resp%s_full(3,2,iband,:,info%ik))

      resp%a_full(1,2,iband,:,info%ik) = conjg(resp%a_full(2,1,iband,:,info%ik))
      resp%a_full(1,3,iband,:,info%ik) = conjg(resp%a_full(3,1,iband,:,info%ik))
      resp%a_full(2,3,iband,:,info%ik) = conjg(resp%a_full(3,2,iband,:,info%ik))

      resp%x_full(1,2,iband,:,info%ik) = conjg(resp%x_full(2,1,iband,:,info%ik))
      resp%x_full(1,3,iband,:,info%ik) = conjg(resp%x_full(3,1,iband,:,info%ik))
      resp%x_full(2,3,iband,:,info%ik) = conjg(resp%x_full(3,2,iband,:,info%ik))
    endif


    resp%s_full(1,1,iband,:,info%ik) = &
    resp%s_full(1,1,iband,:,info%ik) * edisp%MoptDiag(1,iband,:,info%ik)

    resp%a_full(1,1,iband,:,info%ik) = &
    resp%a_full(1,1,iband,:,info%ik) * edisp%MoptDiag(1,iband,:,info%ik)

    resp%x_full(1,1,iband,:,info%ik) = &
    resp%x_full(1,1,iband,:,info%ik) * edisp%MoptDiag(1,iband,:,info%ik)
  enddo

end subroutine response_intra_optical_weights_Q

subroutine response_h5_output_Q(resp, gname, edisp, algo, info, temp, kmesh, lBfield)
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
  integer(hid_t)     :: ifile

  ! this sucks
  ! quadruple response s/a/x array
  real(16), allocatable :: qrsarr(:,:,:) ! to collect the data
  real(16), allocatable :: qraarr(:,:,:)
  real(16), allocatable :: qrxarr(:,:,:)
  real(16), allocatable :: qisarr(:,:,:)
  real(16), allocatable :: qiaarr(:,:,:)
  real(16), allocatable :: qixarr(:,:,:)

  ! double complex (z) array
  complex(8),  allocatable :: zdarr(:,:,:) ! for output

  integer :: iband, ik, is
  integer :: ii, ij

  allocate(qrsarr(3,3,edisp%ispin))
  allocate(qraarr(3,3,edisp%ispin))
  allocate(qrxarr(3,3,edisp%ispin))
  allocate(qisarr(3,3,edisp%ispin))
  allocate(qiaarr(3,3,edisp%ispin))
  allocate(qixarr(3,3,edisp%ispin))

  allocate(zdarr(3,3,edisp%ispin))

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
  zdarr = 0.d0
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

  qrsarr = qrsarr * pi * ( echarge / (kmesh%vol*hbarevs)) * 1.q10
  qisarr = qisarr * pi * ( echarge / (kmesh%vol*hbarevs)) * 1.q10
  qraarr = qraarr * pi * ( echarge / (kmesh%vol*hbarevs)) * 1.q10
  qiaarr = qiaarr * pi * ( echarge / (kmesh%vol*hbarevs)) * 1.q10
  qrxarr = qrxarr * pi * ( echarge / (kmesh%vol*hbarevs)) * 1.q10
  qixarr = qixarr * pi * ( echarge / (kmesh%vol*hbarevs)) * 1.q10

  ! should work=?
  if (myid .eq. master) then
    if (algo%lDebug .and. (index(algo%dbgstr,"ReduceIO") .ne. 0)) then
      ! gather the data in the arrays
      zdarr = cmplx(real(qrsarr,8),real(qisarr,8))
      resp%s_sum_temp(:,:,:,info%iT) = zdarr
      zdarr = cmplx(real(qrsarr,8),real(qisarr,8))
      resp%a_sum_temp(:,:,:,info%iT) = resp%a_sum
      zdarr = cmplx(real(qrsarr,8),real(qisarr,8))
      resp%x_sum_temp(:,:,:,info%iT) = zdarr

      ! output at the last temperature step
      if ((temp%Tstep==1 .and. info%iT==temp%nT) .or. (temp%Tstep==-1 .and. info%iT==1)) then
        string = "/L0/" // trim(adjustl(gname)) // "/sum"
        call hdf5_write_data(ifile, string, resp%s_sum_temp)
        string = "/L1/" // trim(adjustl(gname)) // "/sum"
        call hdf5_write_data(ifile, string, resp%a_sum_temp)
        string = "/L2/" // trim(adjustl(gname)) // "/sum"
        call hdf5_write_data(ifile, string, resp%x_sum_temp)
      endif

    else
      write(string,'(I6.6)') info%iT
      string = trim(string) // "/L0/" // trim(adjustl(gname)) // "/sum"
      zdarr = cmplx(real(qrsarr,8),real(qisarr,8))
      call hdf5_write_data(ifile, string, zdarr)

      write(string,'(I6.6)') info%iT
      string = trim(string) // "/L1/" // trim(adjustl(gname)) // "/sum"
      zdarr = cmplx(real(qraarr,8),real(qiaarr,8))
      call hdf5_write_data(ifile, string, zdarr)

      write(string,'(I6.6)') info%iT
      string = trim(string) // "/L2/" // trim(adjustl(gname)) // "/sum"
      zdarr = cmplx(real(qrxarr,8),real(qixarr,8))
      call hdf5_write_data(ifile, string, zdarr)
    endif
  endif

  ! if (lBoutput) then
  !   ! perform a local summation
  !   ! these are already initialized to 0
  !   do ik = ikstr,ikend
  !     do iband = edisp%nbopt_min,edisp%nbopt_max
  !       resp%sB_sum(:,:,:) = resp%sB_sum(:,:,:) + resp%sB_full(:,:,iband,:,ik) * kmesh%weightQ(ik)
  !       resp%aB_sum(:,:,:) = resp%aB_sum(:,:,:) + resp%aB_full(:,:,iband,:,ik) * kmesh%weightQ(ik)
  !       resp%xB_sum(:,:,:) = resp%xB_sum(:,:,:) + resp%xB_full(:,:,iband,:,ik) * kmesh%weightQ(ik)
  !     enddo
  !   enddo

  ! ! perform MPI summation
! #ifdef MPI
  ! if (myid.eq.master) then
  !   call MPI_REDUCE(MPI_IN_PLACE, resp%sB_sum, 9*edisp%ispin, MPI_DOUBLE_COMPLEX, MPI_SUM, master, MPI_COMM_WORLD, mpierr)
  ! else
  !   call MPI_REDUCE(resp%sB_sum, resp%sB_sum, 9*edisp%ispin, MPI_DOUBLE_COMPLEX, MPI_SUM, master, MPI_COMM_WORLD, mpierr)
  ! endif

  ! if (myid.eq.master) then
  !   call MPI_REDUCE(MPI_IN_PLACE, resp%aB_sum, 9*edisp%ispin, MPI_DOUBLE_COMPLEX, MPI_SUM, master, MPI_COMM_WORLD, mpierr)
  ! else
  !   call MPI_REDUCE(resp%aB_sum, resp%aB_sum, 9*edisp%ispin, MPI_DOUBLE_COMPLEX, MPI_SUM, master, MPI_COMM_WORLD, mpierr)
  ! endif

  ! if (myid.eq.master) then
  !   call MPI_REDUCE(MPI_IN_PLACE, resp%xB_sum, 9*edisp%ispin, MPI_DOUBLE_COMPLEX, MPI_SUM, master, MPI_COMM_WORLD, mpierr)
  ! else
  !   call MPI_REDUCE(resp%xB_sum, resp%xB_sum, 9*edisp%ispin, MPI_DOUBLE_COMPLEX, MPI_SUM, master, MPI_COMM_WORLD, mpierr)
  ! endif
! #endif

  !   resp%sB_sum = resp%sB_sum * pi**2 * ( echarge / (kmesh%vol*hbarevs)) * (1.d-10 / hbarevs)
  !   resp%aB_sum = resp%aB_sum * pi**2 * ( echarge / (kmesh%vol*hbarevs)) * (1.d-10 * echarge / hbarevs)
  !   resp%xB_sum = resp%xB_sum * pi**2 * ( echarge / (kmesh%vol*hbarevs)) * (1.d-10 * echarge**2 / hbarevs)

  !   if (myid .eq. master) then
  !     write(string,'(I6.6)') info%iT
  !     string = trim(string) // "/L0/" // trim(adjustl(gname)) // "/sumM"
  !     call hdf5_write_data(ifile, string, resp%sB_sum)

  !     write(string,'(I6.6)') info%iT
  !     string = trim(string) // "/L1/" // trim(adjustl(gname)) // "/sumM"
  !     call hdf5_write_data(ifile, string, resp%aB_sum)

  !     write(string,'(I6.6)') info%iT
  !     string = trim(string) // "/L2/" // trim(adjustl(gname)) // "/sumM"
  !     call hdf5_write_data(ifile, string, resp%xB_sum)

  !   endif
  ! endif ! Boutput

  if (myid.eq.master) then
    call hdf5_close_file(ifile)
  endif

  deallocate(qrsarr)
  deallocate(qraarr)
  deallocate(qrxarr)
  deallocate(qisarr)
  deallocate(qiaarr)
  deallocate(qixarr)

  deallocate(zdarr)

end subroutine


end module Mresponse
