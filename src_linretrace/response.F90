module Mresponse
  use Mmpi_org
  use Mtypes
  use Mparams
  use Mfermi

  implicit none

  interface initialize_response
    module procedure initresp, initresp_qp
  end interface

  interface allocate_response
    module procedure dpresp_alloc, qpresp_alloc
  end interface

  interface calc_polygamma
    module procedure calc_polygamma_D, calc_polygamma_Q
  end interface calc_polygamma

  contains

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

subroutine response_intra_km(resp, PolyGamma, edisp, sct, kmesh, algo, info)
  implicit none
  type (response_dp)   :: resp

  type(energydisp)     :: edisp
  type(scattering)     :: sct
  type(kpointmesh)     :: kmesh
  type(algorithm)      :: algo
  type(runinfo)        :: info

  complex(8)           :: PolyGamma(3,edisp%nbopt_min:edisp%nbopt_max, ikstr:ikend, edisp%ispin)
  real(8), allocatable :: enrgy(:,:)

  allocate(enrgy(edisp%nbopt_min:edisp%nbopt_max,edisp%ispin))
  ! first we write the kernel into the 1 1 component
  enrgy = sct%zqp(edisp%nbopt_min:edisp%nbopt_max,info%ik,:) * (edisp%band(edisp%nbopt_min:edisp%nbopt_max,info%ik,:) - info%mu)

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

    resp%sB_full(1,1,1,:,:,info%ik) = real(PolyGamma(3,:,info%ik,:)) &
                                      * sct%gam(edisp%nbopt_min:edisp%nbopt_max,info%ik,:)**2 * info%beta**2 / (4.d0 * pi**2) &
                                  - real(PolyGamma(2,:,info%ik,:)) &
                                    * 3.d0 * sct%gam(edisp%nbopt_min:edisp%nbopt_max,info%ik,:) * info%beta2p &
                                  + real(PolyGamma(1,:,info%ik,:)) * 3.d0

    resp%sB_full(1,1,1,:,:,info%ik) = resp%sB_full(1,1,1,:,:,info%ik) &
                                  * sct%zqp(edisp%nbopt_min:edisp%nbopt_max,info%ik,:)**3 * info%beta &
                                  / (16.d0 * pi**4 * sct%gam(edisp%nbopt_min:edisp%nbopt_max,info%ik,:)**2)


    resp%aB_full(1,1,1,:,:,info%ik) =  real(PolyGamma(3,:,info%ik,:)) &
                                       * enrgy * sct%gam(edisp%nbopt_min:edisp%nbopt_max,info%ik,:)**2 &
                                       * info%beta**2 / (4.d0 * pi**2) &
                                  + aimag(PolyGamma(3,:,info%ik,:)) &
                                       * sct%gam(edisp%nbopt_min:edisp%nbopt_max,info%ik,:)**3 * info%beta**2 / (4.d0 * pi**2) &
                                  -  real(PolyGamma(2,:,info%ik,:)) &
                                       * 3.d0 * enrgy * sct%gam(edisp%nbopt_min:edisp%nbopt_max,info%ik,:) * info%beta2p &
                                  - aimag(PolyGamma(2,:,info%ik,:)) &
                                       * sct%gam(edisp%nbopt_min:edisp%nbopt_max,info%ik,:)**2 * info%beta2p &
                                  +  real(PolyGamma(1,:,info%ik,:)) * 3.d0 * enrgy

    resp%aB_full(1,1,1,:,:,info%ik) = resp%aB_full(1,1,1,:,:,info%ik) &
                                  * sct%zqp(edisp%nbopt_min:edisp%nbopt_max,info%ik,:)**3 * info%beta / (16.d0 * pi**4 &
                                  * sct%gam(edisp%nbopt_min:edisp%nbopt_max,info%ik,:)**2)


    resp%xB_full(1,1,1,:,:,info%ik) =  real(PolyGamma(3,:,info%ik,:)) &
                                       * info%beta**2 * sct%gam(edisp%nbopt_min:edisp%nbopt_max,info%ik,:)**2 &
                                       * (sct%gam(edisp%nbopt_min:edisp%nbopt_max,info%ik,:)**2 - enrgy**2) &
                                       / (4.d0 * pi**2) &
                                  - aimag(PolyGamma(3,:,info%ik,:)) &
                                       * info%beta**2 * enrgy * sct%gam(edisp%nbopt_min:edisp%nbopt_max,info%ik,:)**3 &
                                       / (2.d0 * pi**2) &
                                  +  real(PolyGamma(2,:,info%ik,:)) &
                                       * info%beta * sct%gam(edisp%nbopt_min:edisp%nbopt_max,info%ik,:) &
                                       * (sct%gam(edisp%nbopt_min:edisp%nbopt_max,info%ik,:)**2 + 3.d0*enrgy**2) &
                                       / (2.d0 * pi) &
                                  + aimag(PolyGamma(2,:,info%ik,:)) &
                                       * info%beta * enrgy * sct%gam(edisp%nbopt_min:edisp%nbopt_max,info%ik,:)**2 / pi &
                                  -  real(PolyGamma(1,:,info%ik,:)) &
                                       * (sct%gam(edisp%nbopt_min:edisp%nbopt_max,info%ik,:)**2 + 3.d0*enrgy**2)

    resp%xB_full(1,1,1,:,:,info%ik) = resp%xB_full(1,1,1,:,:,info%ik) &
                                  * (-1.d0) * sct%zqp(edisp%nbopt_min:edisp%nbopt_max,info%ik,:)**3 * info%beta &
                                  / (16.d0 * pi**4 * sct%gam(edisp%nbopt_min:edisp%nbopt_max,info%ik,:)**2)

  endif

  deallocate(enrgy)

  call response_intra_optical_weights(algo, resp, edisp, info)
  if (algo%lBfield) then
    call response_peierls_weights(algo, resp, edisp, info)
  endif


end subroutine response_intra_km

subroutine response_inter_km(resp, PolyGamma, edisp, sct, kmesh, algo, info)
  implicit none
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
          * (edisp%band(edisp%nbopt_min:edisp%nbopt_max,info%ik,:) - info%mu)

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

          calc_sigma = real((enrgydiff(is)**2 + sct%gam(iband2,info%ik,is)**2 - sct%gam(iband1,info%ik,is)**2 &
                             - 2*ci*sct%gam(iband1,info%ik,is)*(-enrgydiff(is))*sct%gam(iband1,info%ik,is)) &
                           * sct%gam(iband2,info%ik,is) * PolyGamma(1,iband1,info%ik,is)) &
                     + real((enrgydiff(is)**2 + sct%gam(iband1,info%ik,is)**2 - sct%gam(iband2,info%ik,is)**2 &
                             - 2*ci*sct%gam(iband2,info%ik,is)*enrgydiff(is)*sct%gam(iband2,info%ik,is)) &
                           * sct%gam(iband1,info%ik,is) * PolyGamma(1,iband2,info%ik,is))

          calc_sigma = calc_sigma * sct%zqp(iband1,info%ik,is) * sct%zqp(iband2,info%ik,is) * info%beta &
                     / (2.d0 * pi**3 * ( enrgydiff(is)**2 + (sct%gam(iband1,info%ik,is) - sct%gam(iband2,info%ik,is))**2)) &
                     / ( enrgydiff(is)**2 + (sct%gam(iband1,info%ik,is) + sct%gam(iband2,info%ik,is))**2)


          calc_alpha = real((enrgydiff(is)**2 + sct%gam(iband2,info%ik,is)**2 - sct%gam(iband1,info%ik,is)**2 &
                             - 2*ci*sct%gam(iband1,info%ik,is)*(-enrgydiff(is))*sct%gam(iband1,info%ik,is)) &
                           * (enrgy(iband1,is)-ci*sct%gam(iband1,info%ik,is)) &
                           * sct%gam(iband2,info%ik,is) * PolyGamma(1,iband1,info%ik,is)) &
                     + real((enrgydiff(is)**2 + sct%gam(iband1,info%ik,is)**2 - sct%gam(iband2,info%ik,is)**2 &
                             - 2*ci*sct%gam(iband2,info%ik,is)*enrgydiff(is)*sct%gam(iband2,info%ik,is)) &
                           * (enrgy(iband2,is)-ci*sct%gam(iband2,info%ik,is)) &
                           * sct%gam(iband1,info%ik,is) * PolyGamma(1,iband2,info%ik,is))

          calc_alpha = calc_alpha * sct%zqp(iband1,info%ik,is) * sct%zqp(iband2,info%ik,is) * info%beta &
                     / (2.d0 * pi**3 * ( enrgydiff(is)**2 + (sct%gam(iband1,info%ik,is) - sct%gam(iband2,info%ik,is))**2)) &
                     / ( enrgydiff(is)**2 + (sct%gam(iband1,info%ik,is) + sct%gam(iband2,info%ik,is))**2)

          calc_xi    = real((enrgydiff(is)**2 + sct%gam(iband2,info%ik,is)**2 - sct%gam(iband1,info%ik,is)**2 &
                             - 2*ci*sct%gam(iband1,info%ik,is)*(-enrgydiff(is))*sct%gam(iband1,info%ik,is)) &
                           * (enrgy(iband1,is)-ci*sct%gam(iband1,info%ik,is))**2 &
                           * sct%gam(iband2,info%ik,is) * PolyGamma(1,iband1,info%ik,is)) &
                     + real((enrgydiff(is)**2 + sct%gam(iband1,info%ik,is)**2 - sct%gam(iband2,info%ik,is)**2 &
                             - 2*ci*sct%gam(iband2,info%ik,is)*enrgydiff(is)*sct%gam(iband2,info%ik,is)) &
                           * (enrgy(iband2,is)-ci*sct%gam(iband2,info%ik,is))**2 &
                           * sct%gam(iband1,info%ik,is) * PolyGamma(1,iband2,info%ik,is))

          calc_xi    = calc_xi * sct%zqp(iband1,info%ik,is) * sct%zqp(iband2,info%ik,is) * info%beta &
                     / (2.d0 * pi**3 * ( enrgydiff(is)**2 + (sct%gam(iband1,info%ik,is) - sct%gam(iband2,info%ik,is))**2)) &
                     / ( enrgydiff(is)**2 + (sct%gam(iband1,info%ik,is) + sct%gam(iband2,info%ik,is))**2)

        endif

        if (algo%ldebug .and. (index(algo%dbgstr,"KernelsOnly") .ne. 0)) then
          cycle
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

subroutine response_intra_Boltzmann_km(resp, edisp, sct, kmesh, algo, info)
  implicit none
  type (response_dp)      :: resp

  type(energydisp)    :: edisp
  type(scattering)    :: sct
  type(kpointmesh)    :: kmesh
  type(algorithm)     :: algo
  type(runinfo)       :: info

  real(8) :: zqp
  real(8) :: gam
  real(8), allocatable :: enrgy(:,:)

  integer :: i,j,is
  integer :: iband


  allocate(enrgy(edisp%nbopt_min:edisp%nbopt_max,edisp%ispin))
  ! first we write the kernel into the 1 1 component
  enrgy = sct%zqp(edisp%nbopt_min:edisp%nbopt_max,info%ik,:) &
          * (edisp%band(edisp%nbopt_min:edisp%nbopt_max,info%ik,:) - info%mu)

  ! asymptotic term in the Gamma-> 0 limit
  ! the polygamma2fermi already contains the pi^2 / 2
  if (algo%lBoltzmannFermi) then
    resp%s_full(1,1,:,:,info%ik) = polygamma2fermi(enrgy,info%beta) &
                                    * sct%zqp(edisp%nbopt_min:edisp%nbopt_max,info%ik,:)**2 * info%beta &
                                    / (4.d0 * pi**3 * sct%gam(edisp%nbopt_min:edisp%nbopt_max,info%ik,:))
  else
    do is=1,edisp%ispin
      do iband=edisp%nbopt_min,edisp%nbopt_max
        resp%s_full(1,1,iband,is,info%ik) = polygamma2psi1(sct%gam(iband,info%ik,is),enrgy(iband,is),info%beta) &
                                        * sct%zqp(iband,info%ik,is)**2 * info%beta &
                                        / (4.d0 * pi**3 * sct%gam(iband,info%ik,is))
      enddo
    enddo
  endif

  resp%a_full(1,1,:,:,info%ik) = resp%s_full(1,1,:,:,info%ik) * enrgy

  resp%x_full(1,1,:,:,info%ik) = resp%s_full(1,1,:,:,info%ik) * enrgy**2

  if (algo%lBfield) then

    if (algo%lBoltzmannFermi) then
      resp%sB_full(1,1,1,:,:,info%ik) = polygamma2fermi(enrgy,info%beta) &
                                      * 3.d0 * sct%zqp(edisp%nbopt_min:edisp%nbopt_max,info%ik,:)**3 * info%beta &
                                      / (16.d0 * pi**4 * sct%gam(edisp%nbopt_min:edisp%nbopt_max,info%ik,:)**2)
    else
      do is=1,edisp%ispin
        do iband=edisp%nbopt_min,edisp%nbopt_max
          resp%sB_full(1,1,1,iband,is,info%ik) = polygamma2psi1(sct%gam(iband,info%ik,is),enrgy(iband,is),info%beta) &
                                          * 3.d0 * sct%zqp(iband,info%ik,is)**3 * info%beta &
                                          / (16.d0 * pi**4 * sct%gam(iband,info%ik,is)**2)
        enddo
      enddo
    endif

    resp%aB_full(1,1,1,:,:,info%ik) = resp%sB_full(1,1,1,:,:,info%ik) * enrgy

    resp%xB_full(1,1,1,:,:,info%ik) = resp%sB_full(1,1,1,:,:,info%ik) * enrgy**2

  endif

  deallocate(enrgy)

  call response_intra_optical_weights(algo, resp, edisp, info)
  if (algo%lBfield) then
    call response_peierls_weights(algo, resp, edisp, info)
  endif

end subroutine response_intra_Boltzmann_km

subroutine response_intra_Boltzmann_km_Q(resp, edisp, sct, kmesh, algo, info)
  implicit none
  type (response_qp)    :: resp

  type(energydisp)      :: edisp
  type(scattering)      :: sct
  type(kpointmesh)      :: kmesh
  type(algorithm)       :: algo
  type(runinfo)         :: info

  real(16), allocatable :: enrgy(:,:)

  integer :: i,j,is
  integer :: iband


  allocate(enrgy(edisp%nbopt_min:edisp%nbopt_max,edisp%ispin))
  ! first we write the kernel into the 1 1 component
  enrgy = sct%zqp(edisp%nbopt_min:edisp%nbopt_max,info%ik,:) &
          * (edisp%band(edisp%nbopt_min:edisp%nbopt_max,info%ik,:) - info%muQ)

  ! asymptotic term in the Gamma-> 0 limit
  ! the polygamma2fermi already contains the pi^2 / 2
  if (algo%lBoltzmannFermi) then
    resp%s_full(1,1,:,:,info%ik) = polygamma2fermi(enrgy,info%betaQ) &
                                    * sct%zqp(edisp%nbopt_min:edisp%nbopt_max,info%ik,:)**2 * info%betaQ &
                                    / (4.q0 * piQ**3 * sct%gam(edisp%nbopt_min:edisp%nbopt_max,info%ik,:))
  else
    do is=1,edisp%ispin
      do iband=edisp%nbopt_min,edisp%nbopt_max
        resp%s_full(1,1,iband,is,info%ik) = polygamma2psi1(sct%gam(iband,info%ik,is),enrgy(iband,is),info%betaQ) &
                                        * sct%zqp(iband,info%ik,is)**2 * info%beta &
                                        / (4.q0 * piQ**3 * sct%gam(iband,info%ik,is))
      enddo
    enddo
  endif

  resp%a_full(1,1,:,:,info%ik) = resp%s_full(1,1,:,:,info%ik) * enrgy

  resp%x_full(1,1,:,:,info%ik) = resp%s_full(1,1,:,:,info%ik) * enrgy**2

  if (algo%lBfield) then

    if (algo%lBoltzmannFermi) then
      resp%sB_full(1,1,1,:,:,info%ik) = polygamma2fermi(enrgy,info%betaQ) &
                                      * 3.q0 * sct%zqp(edisp%nbopt_min:edisp%nbopt_max,info%ik,:)**3 * info%betaQ &
                                      / (16.d0 * piQ**4 * sct%gam(edisp%nbopt_min:edisp%nbopt_max,info%ik,:)**2)
    else
      do is=1,edisp%ispin
        do iband=edisp%nbopt_min,edisp%nbopt_max
          resp%sB_full(1,1,1,iband,is,info%ik) = polygamma2psi1(sct%gam(iband,info%ik,is),enrgy(iband,is),info%betaQ) &
                                          * 3.q0 * sct%zqp(iband,info%ik,is)**3 * info%betaQ &
                                          / (16.q0 * piQ**4 * sct%gam(iband,info%ik,is)**2)
        enddo
      enddo
    endif

    resp%aB_full(1,1,1,:,:,info%ik) = resp%sB_full(1,1,1,:,:,info%ik) * enrgy

    resp%xB_full(1,1,1,:,:,info%ik) = resp%sB_full(1,1,1,:,:,info%ik) * enrgy**2

  endif

  deallocate(enrgy)

  call response_intra_optical_weights_Q(algo, resp, edisp, info)
  if (algo%lBfield) then
    call response_peierls_weights_Q(algo, resp, edisp, info)
  endif

end subroutine response_intra_Boltzmann_km_Q

subroutine response_inter_Boltzmann_km(resp, edisp, sct, kmesh, algo, info)
  implicit none
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
          * (edisp%band(edisp%nbopt_min:edisp%nbopt_max,info%ik,:) - info%mu)

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

        if (algo%ldebug .and. (index(algo%dbgstr,"KernelsOnly") .ne. 0)) then
          cycle
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

subroutine response_inter_Boltzmann_km_Q(resp, edisp, sct, kmesh, algo, info)
  implicit none
  type (response_qp)   :: resp

  type(energydisp)     :: edisp
  type(scattering)     :: sct
  type(kpointmesh)     :: kmesh
  type(algorithm)      :: algo
  type(runinfo)        :: info

  real(16), allocatable :: enrgy(:,:)
  real(16), allocatable :: enrgydiff(:)
  real(16), allocatable :: gamdiff(:)

  complex(16) :: calc_sigma
  complex(16) :: calc_alpha
  complex(16) :: calc_xi

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
          * (edisp%band(edisp%nbopt_min:edisp%nbopt_max,info%ik,:) - info%muQ)

  do iband1 = edisp%nbopt_min, edisp%nbopt_max
    do iband2 = edisp%nbopt_min, edisp%nbopt_max
      if (iband1 == iband2) cycle
      enrgydiff = enrgy(iband1,:) - enrgy(iband2,:)
      gamdiff   = sct%gam(iband1,info%ik,:) - sct%gam(iband2,info%ik,:)

      do is = 1,edisp%ispin
        if ((abs(enrgydiff(is)) .lt. 1d-6) .and. (abs(gamdiff(is)) .lt. 1d-6)) then
        ! use the intra-band limit .....
          calc_sigma = polygamma2fermi(enrgy(iband1,is),info%betaQ) &
                     * sct%zqp(iband1,info%ik,is)**2 * info%betaQ &
                     / (4.d0 * piQ**3 * sct%gam(iband1,info%ik,is))

          calc_alpha = calc_sigma * enrgy(iband1,is)

          calc_xi    = calc_sigma * enrgy(iband1,is)**2

        else

          calc_sigma = polygamma2fermi(enrgy(iband1,is), info%betaQ) &
                     * enrgydiff(is)**2 / sct%gam(iband1,info%ik,is)

          calc_sigma = calc_sigma &
                     + polygamma2fermi(enrgy(iband2,is), info%betaQ) &
                     * enrgydiff(is)**2 / sct%gam(iband2,info%ik,is)

          calc_sigma = calc_sigma &
                     * sct%gam(iband1,info%ik,is) * sct%gam(iband2,info%ik,is) &
                     * sct%zqp(iband1,info%ik,is) * sct%zqp(iband2,info%ik,is) &
                     * info%betaQ

          calc_sigma = calc_sigma &
                     / (2.d0 * piQ**3 * ( enrgydiff(is)**2 + (sct%gam(iband1,info%ik,is) - sct%gam(iband2,info%ik,is))**2)) &
                     / ( enrgydiff(is)**2 + (sct%gam(iband1,info%ik,is) + sct%gam(iband2,info%ik,is))**2)


          calc_alpha = polygamma2fermi(enrgy(iband1,is), info%betaQ) &
                     * enrgy(iband1,is) * enrgydiff(is)**2 / sct%gam(iband1,info%ik,is)

          calc_alpha = calc_alpha &
                     + polygamma2fermi(enrgy(iband2,is), info%betaQ) &
                     * enrgy(iband2,is) * enrgydiff(is)**2 / sct%gam(iband2,info%ik,is)

          calc_alpha = calc_alpha &
                     * sct%gam(iband1,info%ik,is) * sct%gam(iband2,info%ik,is) &
                     * sct%zqp(iband1,info%ik,is) * sct%zqp(iband2,info%ik,is) &
                     * info%betaQ

          calc_alpha = calc_alpha &
                     / (2.d0 * piQ**3 * ( enrgydiff(is)**2 + (sct%gam(iband1,info%ik,is) - sct%gam(iband2,info%ik,is))**2)) &
                     / ( enrgydiff(is)**2 + (sct%gam(iband1,info%ik,is) + sct%gam(iband2,info%ik,is))**2)

          calc_xi    = polygamma2fermi(enrgy(iband1,is), info%betaQ) &
                     * enrgy(iband1,is)**2 * enrgydiff(is)**2 / sct%gam(iband1,info%ik,is)

          calc_xi    = calc_xi &
                     + polygamma2fermi(enrgy(iband2,is), info%betaQ) &
                     * enrgy(iband2,is)**2 * enrgydiff(is)**2 / sct%gam(iband2,info%ik,is)

          calc_xi    = calc_xi &
                     * sct%gam(iband1,info%ik,is) * sct%gam(iband2,info%ik,is) &
                     * sct%zqp(iband1,info%ik,is) * sct%zqp(iband2,info%ik,is) &
                     * info%betaQ

          calc_xi    = calc_xi &
                     / (2.d0 * piQ**3 * ( enrgydiff(is)**2 + (sct%gam(iband1,info%ik,is) - sct%gam(iband2,info%ik,is))**2)) &
                     / ( enrgydiff(is)**2 + (sct%gam(iband1,info%ik,is) + sct%gam(iband2,info%ik,is))**2)
        endif

        if (algo%ldebug .and. (index(algo%dbgstr,"KernelsOnly") .ne. 0)) then
          cycle
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

end subroutine response_inter_Boltzmann_km_Q


! multiply optical elements onto quantities without B-Field
subroutine response_intra_optical_weights(algo, resp, edisp, info)
  implicit none
  type(algorithm)    :: algo
  type(response_dp)  :: resp
  type(energydisp)   :: edisp
  type(runinfo)      :: info

  integer :: index1(9), index2(9)
  integer :: iband, idir

  if (algo%ldebug .and. (index(algo%dbgstr,"KernelsOnly") .ne. 0)) then
    return
  endif

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

subroutine response_peierls_weights(algo, resp, edisp, info)
  implicit none
  type(algorithm)    :: algo
  type(response_dp)  :: resp
  type(energydisp)   :: edisp
  type(runinfo)      :: info

  integer :: iband
  integer :: i,j,k

  if (algo%ldebug .and. (index(algo%dbgstr,"KernelsOnly") .ne. 0)) then
    return
  endif

  do iband = edisp%nbopt_min, edisp%nbopt_max
    do i=1,3
      do j=1,3
        do k=1,3
          if (i==1 .and. j==1 .and. k==1) cycle ! we need to keep the kernel saved

          resp%sB_full(k,j,i,iband,:,info%ik) = resp%sB_full(1,1,1,iband,:,info%ik) &
          * edisp%MBoptDiag(k,j,i,iband,:,info%ik)

          resp%aB_full(k,j,i,iband,:,info%ik) = resp%aB_full(1,1,1,iband,:,info%ik) &
          * edisp%MBoptDiag(k,j,i,iband,:,info%ik)

          resp%xB_full(k,j,i,iband,:,info%ik) = resp%xB_full(1,1,1,iband,:,info%ik) &
          * edisp%MBoptDiag(k,j,i,iband,:,info%ik)

        enddo
      enddo
    enddo

    resp%sB_full(1,1,1,iband,:,info%ik) = resp%sB_full(1,1,1,iband,:,info%ik) &
    * edisp%MBoptDiag(1,1,1,iband,:,info%ik)

    resp%aB_full(1,1,1,iband,:,info%ik) = resp%aB_full(1,1,1,iband,:,info%ik) &
    * edisp%MBoptDiag(1,1,1,iband,:,info%ik)

    resp%xB_full(1,1,1,iband,:,info%ik) = resp%xB_full(1,1,1,iband,:,info%ik) &
    * edisp%MBoptDiag(1,1,1,iband,:,info%ik)

  enddo

end subroutine

subroutine response_peierls_weights_Q(algo, resp, edisp, info)
  implicit none
  type(algorithm)    :: algo
  type(response_qp)  :: resp
  type(energydisp)   :: edisp
  type(runinfo)      :: info

  integer :: iband
  integer :: i,j,k

  if (algo%ldebug .and. (index(algo%dbgstr,"KernelsOnly") .ne. 0)) then
    return
  endif

  do iband = edisp%nbopt_min, edisp%nbopt_max
    do i=1,3
      do j=1,3
        do k=1,3
          if (i==1 .and. j==1 .and. k==1) cycle ! we need to keep the kernel saved

          resp%sB_full(k,j,i,iband,:,info%ik) = resp%sB_full(1,1,1,iband,:,info%ik) &
          * edisp%MBoptDiag(k,j,i,iband,:,info%ik)

          resp%aB_full(k,j,i,iband,:,info%ik) = resp%aB_full(1,1,1,iband,:,info%ik) &
          * edisp%MBoptDiag(k,j,i,iband,:,info%ik)

          resp%xB_full(k,j,i,iband,:,info%ik) = resp%xB_full(1,1,1,iband,:,info%ik) &
          * edisp%MBoptDiag(k,j,i,iband,:,info%ik)

        enddo
      enddo
    enddo

    resp%sB_full(1,1,1,iband,:,info%ik) = resp%sB_full(1,1,1,iband,:,info%ik) &
    * edisp%MBoptDiag(1,1,1,iband,:,info%ik)

    resp%aB_full(1,1,1,iband,:,info%ik) = resp%aB_full(1,1,1,iband,:,info%ik) &
    * edisp%MBoptDiag(1,1,1,iband,:,info%ik)

    resp%xB_full(1,1,1,iband,:,info%ik) = resp%xB_full(1,1,1,iband,:,info%ik) &
    * edisp%MBoptDiag(1,1,1,iband,:,info%ik)

  enddo

end subroutine

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

  if (myid.eq.master) then
    allocate(dpresp%s_sum_range(3,3,edisp%iSpin,algo%steps))
    allocate(dpresp%a_sum_range(3,3,edisp%iSpin,algo%steps))
    allocate(dpresp%x_sum_range(3,3,edisp%iSpin,algo%steps))
  endif


  if (algo%lBfield) then
    allocate(dpresp%sB_full(3,3,3,edisp%nbopt_min:edisp%nbopt_max,edisp%ispin,ikstr:ikend))
    allocate(dpresp%aB_full(3,3,3,edisp%nbopt_min:edisp%nbopt_max,edisp%ispin,ikstr:ikend))
    allocate(dpresp%xB_full(3,3,3,edisp%nbopt_min:edisp%nbopt_max,edisp%ispin,ikstr:ikend))
    allocate(dpresp%sB_sum(3,3,3,edisp%iSpin))
    allocate(dpresp%aB_sum(3,3,3,edisp%iSpin))
    allocate(dpresp%xB_sum(3,3,3,edisp%iSpin))

    if (myid.eq.master) then
      allocate(dpresp%sB_sum_range(3,3,3,edisp%iSpin,algo%steps))
      allocate(dpresp%aB_sum_range(3,3,3,edisp%iSpin,algo%steps))
      allocate(dpresp%xB_sum_range(3,3,3,edisp%iSpin,algo%steps))
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

  if (myid .eq. master) then
    allocate(qpresp%s_sum_range(3,3,edisp%iSpin,algo%steps))
    allocate(qpresp%a_sum_range(3,3,edisp%iSpin,algo%steps))
    allocate(qpresp%x_sum_range(3,3,edisp%iSpin,algo%steps))
  endif

  if (algo%lBfield) then
    allocate(qpresp%sB_full(3,3,3,edisp%nbopt_min:edisp%nbopt_max,edisp%ispin,ikstr:ikend))
    allocate(qpresp%aB_full(3,3,3,edisp%nbopt_min:edisp%nbopt_max,edisp%ispin,ikstr:ikend))
    allocate(qpresp%xB_full(3,3,3,edisp%nbopt_min:edisp%nbopt_max,edisp%ispin,ikstr:ikend))
    allocate(qpresp%sB_sum(3,3,3,edisp%iSpin))
    allocate(qpresp%aB_sum(3,3,3,edisp%iSpin))
    allocate(qpresp%xB_sum(3,3,3,edisp%iSpin))

    if (myid .eq. master) then
      allocate(qpresp%sB_sum_range(3,3,3,edisp%iSpin,algo%steps))
      allocate(qpresp%aB_sum_range(3,3,3,edisp%iSpin,algo%steps))
      allocate(qpresp%xB_sum_range(3,3,3,edisp%iSpin,algo%steps))
    endif
  endif
end subroutine qpresp_alloc

subroutine calc_polygamma_D(PolyGamma, edisp, sct, kmesh, algo, info)
  implicit none
  type(algorithm)  :: algo
  type(energydisp) :: edisp
  type(kpointmesh) :: kmesh
  type(scattering) :: sct
  type(runinfo)    :: info

  complex(8) :: PolyGamma(3,edisp%nbopt_min:edisp%nbopt_max, ikstr:ikend, edisp%ispin)

  complex(8), external :: wpsipg
  complex(8), allocatable :: to_evaluate(:,:,:)
  integer :: ipg, iband, ik, is

  allocate(to_evaluate(edisp%nbopt_min:edisp%nbopt_max, ikstr:ikend, edisp%ispin))

  to_evaluate = 0.5d0 + info%beta2p * &
                (sct%gam(edisp%nbopt_min:edisp%nbopt_max,ikstr:ikend,:) &
                 + ci*sct%zqp(edisp%nbopt_min:edisp%nbopt_max,ikstr:ikend,:) &
                 * (edisp%band(edisp%nbopt_min:edisp%nbopt_max,ikstr:ikend,:) - info%mu))

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

subroutine calc_polygamma_Q(PolyGamma, edisp, sct, kmesh, algo, info)
  implicit none
  type(algorithm)  :: algo
  type(energydisp) :: edisp
  type(kpointmesh) :: kmesh
  type(scattering) :: sct
  type(runinfo)    :: info

  complex(16) :: PolyGamma(3,edisp%nbopt_min:edisp%nbopt_max, ikstr:ikend, edisp%ispin)

  complex(16), external :: wpsipghp
  complex(16), allocatable :: to_evaluate(:,:,:)
  integer :: ipg, iband, ik, is

  allocate(to_evaluate(edisp%nbopt_min:edisp%nbopt_max, ikstr:ikend, edisp%ispin))

  to_evaluate = 0.5q0 + info%beta2pQ * &
                (sct%gam(edisp%nbopt_min:edisp%nbopt_max,ikstr:ikend,:) &
                 + ciQ*sct%zqp(edisp%nbopt_min:edisp%nbopt_max,ikstr:ikend,:) &
                 * (edisp%band(edisp%nbopt_min:edisp%nbopt_max,ikstr:ikend,:) - info%muQ))

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

subroutine calc_total_energy_digamma(energy_tot, edisp, sct, kmesh, imp, algo, info)
  real(8), intent(out) :: energy_tot

  type(energydisp) :: edisp
  type(scattering) :: sct
  type(kpointmesh) :: kmesh
  type(impurity)   :: imp
  type(algorithm)  :: algo
  type(runinfo)    :: info
  !local variables

  real(16) :: energy_loc, energy_sum
  integer :: is, ik, iband, iimp
  complex(16), allocatable :: to_evaluate(:,:,:)
  real(16), allocatable    :: energy_post_factor(:,:,:)
  real(16), allocatable    :: energy(:,:,:)
  !external variables
  complex(16), external :: wpsipghp

  allocate(to_evaluate(edisp%nband_max, ikstr:ikend, edisp%ispin))
  allocate(energy_post_factor(edisp%nband_max, ikstr:ikend, edisp%ispin))
  allocate(energy(edisp%nband_max, ikstr:ikend, edisp%ispin))

  ! for the occupation
  to_evaluate = 0.5q0 + info%beta2pQ * &
                (sct%gam(:,ikstr:ikend,:) - ciQ*sct%zqp(:,ikstr:ikend,:)*(edisp%band(:,ikstr:ikend,:) - info%muQ))
  ! energy contribution : occupation * energy of this occupation
  energy_post_factor = sct%zqp(:,ikstr:ikend,:) * (edisp%band(:,ikstr:ikend,:) - info%muQ)

  ! evaluate the function
  do is = 1,edisp%ispin
    do ik = ikstr, ikend
      do iband=1,edisp%nband_max
        energy(iband,ik,is) = 0.5q0 + aimag(wpsipghp(to_evaluate(iband,ik,is),0))/piQ ! this is the occupation
        energy(iband,ik,is) = energy(iband,ik,is) * kmesh%weightQ(ik) * energy_post_factor(iband,ik,is) ! multiplied with weight and energy gives the energy
      enddo
    enddo
  enddo

  deallocate(to_evaluate)
  deallocate(energy_post_factor)
  energy_loc = sum(energy)
  deallocate(energy)
  energy_sum = 0.q0

#ifdef MPI
  call mpi_reduce_quad(energy_loc, energy_sum)
#else
  energy_sum = energy_loc
#endif

  energy_tot = real(energy_sum,8)
  return

end subroutine

subroutine calc_total_energy_fermi(energy_tot, edisp, sct, kmesh, imp, algo, info)
  real(8), intent(out) :: energy_tot

  type(energydisp) :: edisp
  type(scattering) :: sct
  type(kpointmesh) :: kmesh
  type(impurity)   :: imp
  type(algorithm)  :: algo
  type(runinfo)    :: info
  !local variables

  real(16) :: energy_loc, energy_sum
  integer :: is, ik, iband, iimp
  real(16), allocatable    :: energy_post_factor(:,:,:)
  real(16), allocatable    :: energy(:,:,:)

  allocate(energy_post_factor(edisp%nband_max, ikstr:ikend, edisp%ispin))
  allocate(energy(edisp%nband_max, ikstr:ikend, edisp%ispin))

  ! energy contribution : occupation * energy of this occupation
  energy_post_factor = sct%zqp(:,ikstr:ikend,:) * (edisp%band(:,ikstr:ikend,:) - info%muQ)

  ! evaluate the function
  do is = 1,edisp%ispin
    do ik = ikstr, ikend
      do iband=1,edisp%nband_max
        energy(iband,ik,is) = fermi_qp(sct%zqp(iband,ik,is)*(edisp%band(iband,ik,is)-info%muQ), info%betaQ)   ! occupation
        energy(iband,ik,is) = energy(iband,ik,is) * kmesh%weightQ(ik) * energy_post_factor(iband,ik,is) ! multiplied with weight and energy gives the energy
      enddo
    enddo
  enddo

  deallocate(energy_post_factor)
  energy_loc = sum(energy)
  deallocate(energy)

#ifdef MPI
  call mpi_reduce_quad(energy_loc, energy_sum)
#else
  energy_sum = energy_loc
#endif

  energy_tot = real(energy_sum,8)
  return

end subroutine

subroutine calc_elecholes_digamma(electrons_total, holes_total, edisp, sct, kmesh, imp, algo, info)
  real(8), intent(out) :: electrons_total
  real(8), intent(out) :: holes_total

  type(energydisp) :: edisp
  type(scattering) :: sct
  type(kpointmesh) :: kmesh
  type(impurity)   :: imp
  type(algorithm)  :: algo
  type(runinfo)    :: info
  !local variables

  real(8) :: energy_loc
  integer :: is, ik, iband, iimp


  real(16) :: elecs
  real(16) :: holes
  real(16) :: elecssum
  real(16) :: holessum
  real(16) :: elecsmpi
  real(16) :: holesmpi

  !external variables
  complex(16), external :: wpsipghp

  ! evaluate the function
  elecssum = 0.q0
  holessum = 0.q0

  do is = 1,edisp%ispin
    do ik = ikstr, ikend
      do iband=1,edisp%nband_max
        elecs = 0.5q0 - aimag(wpsipghp(0.5q0 + info%beta2pQ &
          * (sct%gam(iband,ik,is) + ciQ*sct%zqp(iband,ik,is)*(edisp%band(iband,ik,is) - info%muQ)),0))/piQ ! this is the occupation
        holes = 1.q0 - elecs ! should be enough accuracy

        if (algo%lTMODE .and. edisp%gapped(is)) then
          if (iband <= edisp%valenceBand(is)) then
            holessum = holessum + holes * kmesh%weightQ(ik)
          else if (iband >= edisp%conductionBand(is)) then
            elecssum = elecssum + elecs * kmesh%weightQ(ik)
          endif
        else
          if (elecs < holes) then
            elecssum = elecssum + elecs * kmesh%weightQ(ik)
          else
            holessum = holessum + holes * kmesh%weightQ(ik)
          endif
        endif

      enddo
    enddo
  enddo

#ifdef MPI
  call mpi_reduce_quad(elecssum, elecsmpi) ! custom quad reduction
#else
  elecsmpi = elecssum
#endif

#ifdef MPI
  call mpi_reduce_quad(holessum, holesmpi) ! custom quad reduction
#else
  holesmpi = holessum
#endif

  ! backtransform to double precision
  electrons_total = real(elecsmpi,8)
  holes_total = real(holesmpi,8)

end subroutine

subroutine calc_elecholes_fermi(electrons_total, holes_total, edisp, sct, kmesh, imp, algo, info)
  real(8), intent(out) :: electrons_total
  real(8), intent(out) :: holes_total

  type(energydisp) :: edisp
  type(scattering) :: sct
  type(kpointmesh) :: kmesh
  type(impurity)   :: imp
  type(algorithm)  :: algo
  type(runinfo)    :: info
  !local variables

  real(8) :: energy_loc
  integer :: is, ik, iband, iimp


  real(16) :: elecs
  real(16) :: holes
  real(16) :: elecssum
  real(16) :: holessum
  real(16) :: elecsmpi
  real(16) :: holesmpi

  !external variables
  complex(8), external :: wpsipg

  ! evaluate the function
  elecssum = 0.q0
  holessum = 0.q0

  do is = 1,edisp%ispin
    do ik = ikstr, ikend
      do iband=1,edisp%nband_max
        elecs = fermi_qp(sct%zqp(iband,ik,is)*(edisp%band(iband,ik,is)-info%muQ), info%betaQ)
        holes = omfermi_qp(sct%zqp(iband,ik,is)*(edisp%band(iband,ik,is)-info%muQ), info%betaQ)

        if (algo%lTMODE .and. edisp%gapped(is)) then
          if (iband <= edisp%valenceBand(is)) then
            holessum = holessum + holes * kmesh%weightQ(ik)
          else if (iband >= edisp%conductionBand(is)) then
            elecssum = elecssum + elecs * kmesh%weightQ(ik)
          endif
        else
          if (elecs < holes) then
            elecssum = elecssum + elecs * kmesh%weightQ(ik)
          else
            holessum = holessum + holes * kmesh%weightQ(ik)
          endif
        endif
      enddo
    enddo
  enddo

#ifdef MPI
  call mpi_reduce_quad(elecssum, elecsmpi) ! custom quad reduction
#else
  elecsmpi = elecssum
#endif

#ifdef MPI
  call mpi_reduce_quad(holessum, holesmpi) ! custom quad reduction
#else
  holesmpi = holessum
#endif

  ! backtransform to double precision
  electrons_total = real(elecsmpi,8)
  holes_total = real(holesmpi,8)

end subroutine

subroutine response_intra_km_Q(resp, PolyGamma, edisp, sct, kmesh, algo, info)
  implicit none
  type (response_qp)    :: resp

  type(energydisp)      :: edisp
  type(scattering)      :: sct
  type(kpointmesh)      :: kmesh
  type(algorithm)       :: algo
  type(runinfo)         :: info

  complex(16)           :: PolyGamma(3,edisp%nbopt_min:edisp%nbopt_max, ikstr:ikend, edisp%ispin)
  real(16), allocatable :: enrgy(:,:)

  allocate(enrgy(edisp%nbopt_min:edisp%nbopt_max,edisp%ispin))
  ! first we write the kernel into the 1 1 component
  enrgy = sct%zqp(edisp%nbopt_min:edisp%nbopt_max,info%ik,:) * (edisp%band(edisp%nbopt_min:edisp%nbopt_max,info%ik,:) - info%muQ)

  resp%s_full(1,1,:,:,info%ik) = real(PolyGamma(1,:,info%ik,:)) &
                               - real(PolyGamma(2,:,info%ik,:)) * info%beta2pQ*sct%gam(edisp%nbopt_min:edisp%nbopt_max,info%ik,:)

  resp%s_full(1,1,:,:,info%ik) = resp%s_full(1,1,:,:,info%ik) &
                               * sct%zqp(edisp%nbopt_min:edisp%nbopt_max,info%ik,:)**2 * info%betaQ &
                               / (4.q0 * piQ**3 * sct%gam(edisp%nbopt_min:edisp%nbopt_max,info%ik,:))


  resp%a_full(1,1,:,:,info%ik) =  real(PolyGamma(1,:,info%ik,:)) * enrgy &
                               -  real(PolyGamma(2,:,info%ik,:)) * enrgy &
                                  * sct%gam(edisp%nbopt_min:edisp%nbopt_max,info%ik,:) * info%beta2pQ &
                               - aimag(PolyGamma(2,:,info%ik,:)) * sct%gam(edisp%nbopt_min:edisp%nbopt_max,info%ik,:)**2 &
                                 * info%beta2pQ

  resp%a_full(1,1,:,:,info%ik) = resp%a_full(1,1,:,:,info%ik) &
                               * sct%zqp(edisp%nbopt_min:edisp%nbopt_max,info%ik,:)**2 * info%betaQ &
                               / (4.q0 * piQ**3 * sct%gam(edisp%nbopt_min:edisp%nbopt_max,info%ik,:))


  resp%x_full(1,1,:,:,info%ik) =  real(PolyGamma(1,:,info%ik,:)) &
                                  * (enrgy**2 + sct%gam(edisp%nbopt_min:edisp%nbopt_max,info%ik,:)**2) &
                               +  real(PolyGamma(2,:,info%ik,:)) &
                                  * info%beta2pQ * sct%gam(edisp%nbopt_min:edisp%nbopt_max,info%ik,:) &
                                  * (sct%gam(edisp%nbopt_min:edisp%nbopt_max,info%ik,:)**2 - enrgy**2) &
                               - aimag(PolyGamma(2,:,info%ik,:)) &
                                  * info%betaQ / piQ * enrgy * sct%gam(edisp%nbopt_min:edisp%nbopt_max,info%ik,:)**2

  resp%x_full(1,1,:,:,info%ik) = resp%x_full(1,1,:,:,info%ik) &
                               * sct%zqp(edisp%nbopt_min:edisp%nbopt_max,info%ik,:)**2 * info%betaQ &
                               / (4.q0 * piQ**3 * sct%gam(edisp%nbopt_min:edisp%nbopt_max,info%ik,:))

  if (algo%lBfield) then

    resp%sB_full(1,1,1,:,:,info%ik) = real(PolyGamma(3,:,info%ik,:)) &
                                      * sct%gam(edisp%nbopt_min:edisp%nbopt_max,info%ik,:)**2 * info%betaQ**2 / (4.q0 * piQ**2) &
                                  - real(PolyGamma(2,:,info%ik,:)) &
                                    * 3.q0 * sct%gam(edisp%nbopt_min:edisp%nbopt_max,info%ik,:) * info%beta2pQ &
                                  + real(PolyGamma(1,:,info%ik,:)) * 3.q0

    resp%sB_full(1,1,1,:,:,info%ik) = resp%sB_full(1,1,1,:,:,info%ik) &
                                  * sct%zqp(edisp%nbopt_min:edisp%nbopt_max,info%ik,:)**3 * info%betaQ &
                                  / (16.q0 * piQ**4 * sct%gam(edisp%nbopt_min:edisp%nbopt_max,info%ik,:)**2)


    resp%aB_full(1,1,1,:,:,info%ik) =  real(PolyGamma(3,:,info%ik,:)) &
                                       * enrgy * sct%gam(edisp%nbopt_min:edisp%nbopt_max,info%ik,:)**2 &
                                       * info%betaQ**2 / (4.q0 * piQ**2) &
                                  + aimag(PolyGamma(3,:,info%ik,:)) &
                                       * sct%gam(edisp%nbopt_min:edisp%nbopt_max,info%ik,:)**3 * info%betaQ**2 / (4.q0 * piQ**2) &
                                  -  real(PolyGamma(2,:,info%ik,:)) &
                                       * 3.q0 * enrgy * sct%gam(edisp%nbopt_min:edisp%nbopt_max,info%ik,:) * info%beta2pQ &
                                  - aimag(PolyGamma(2,:,info%ik,:)) &
                                       * sct%gam(edisp%nbopt_min:edisp%nbopt_max,info%ik,:)**2 * info%beta2pQ &
                                  +  real(PolyGamma(1,:,info%ik,:)) * 3.q0 * enrgy

    resp%aB_full(1,1,1,:,:,info%ik) = resp%aB_full(1,1,1,:,:,info%ik) &
                                  * sct%zqp(edisp%nbopt_min:edisp%nbopt_max,info%ik,:)**3 * info%betaQ / (16.q0 * piQ**4 &
                                  * sct%gam(edisp%nbopt_min:edisp%nbopt_max,info%ik,:)**2)


    resp%xB_full(1,1,1,:,:,info%ik) =  real(PolyGamma(3,:,info%ik,:)) &
                                       * info%betaQ**2 * sct%gam(edisp%nbopt_min:edisp%nbopt_max,info%ik,:)**2 &
                                       * (sct%gam(edisp%nbopt_min:edisp%nbopt_max,info%ik,:)**2 - enrgy**2) &
                                       / (4.q0 * piQ**2) &
                                  - aimag(PolyGamma(3,:,info%ik,:)) &
                                       * info%betaQ**2 * enrgy * sct%gam(edisp%nbopt_min:edisp%nbopt_max,info%ik,:)**3 &
                                       / (2.q0 * piQ**2) &
                                  +  real(PolyGamma(2,:,info%ik,:)) &
                                       * info%betaQ * sct%gam(edisp%nbopt_min:edisp%nbopt_max,info%ik,:) &
                                       * (sct%gam(edisp%nbopt_min:edisp%nbopt_max,info%ik,:)**2 + 3.q0*enrgy**2) &
                                       / (2.q0 * piQ) &
                                  + aimag(PolyGamma(2,:,info%ik,:)) &
                                       * info%betaQ * enrgy * sct%gam(edisp%nbopt_min:edisp%nbopt_max,info%ik,:)**2 / piQ &
                                  -  real(PolyGamma(1,:,info%ik,:)) &
                                       * (sct%gam(edisp%nbopt_min:edisp%nbopt_max,info%ik,:)**2 + 3.q0*enrgy**2)

    resp%xB_full(1,1,1,:,:,info%ik) = resp%xB_full(1,1,1,:,:,info%ik) &
                                  * (-1.q0) * sct%zqp(edisp%nbopt_min:edisp%nbopt_max,info%ik,:)**3 * info%betaQ &
                                  / (16.q0 * piQ**4 * sct%gam(edisp%nbopt_min:edisp%nbopt_max,info%ik,:)**2)

  endif

  deallocate(enrgy)

  call response_intra_optical_weights_Q(algo, resp, edisp, info)
  if (algo%lBfield) then
    call response_peierls_weights_Q(algo, resp, edisp, info)
  endif

end subroutine response_intra_km_Q

subroutine response_inter_km_Q(resp, PolyGamma, edisp, sct, kmesh, algo, info)
  implicit none
  type (response_qp)   :: resp

  type(energydisp)     :: edisp
  type(scattering)     :: sct
  type(kpointmesh)     :: kmesh
  type(algorithm)      :: algo
  type(runinfo)        :: info

  complex(16)          :: PolyGamma(3,edisp%nbopt_min:edisp%nbopt_max, ikstr:ikend, edisp%ispin)

  ! these stay double precision .. because we read them from files
  real(16), allocatable :: enrgy(:,:)
  real(16), allocatable :: enrgydiff(:)
  real(16), allocatable :: gamdiff(:)

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
          * (edisp%band(edisp%nbopt_min:edisp%nbopt_max,info%ik,:) - info%muQ)

  do iband1 = edisp%nbopt_min, edisp%nbopt_max
    do iband2 = edisp%nbopt_min, edisp%nbopt_max
      if (iband1 == iband2) cycle
      enrgydiff = enrgy(iband1,:) - enrgy(iband2,:)
      gamdiff   = sct%gam(iband1,info%ik,:) - sct%gam(iband2,info%ik,:)

      do is = 1,edisp%ispin

        if ((abs(enrgydiff(is)) .lt. 1q-6) .and. (abs(gamdiff(is)) .lt. 1q-6)) then

          ! use the intra-band limit .....
          calc_sigma  = real(PolyGamma(1,iband1,info%ik,is)) &
                      - real(PolyGamma(2,iband1,info%ik,is)) * info%beta2pQ*sct%gam(iband1,info%ik,is)

          calc_sigma = calc_sigma &
                     * sct%zqp(iband1,info%ik,is)**2 * info%betaQ &
                     / (4.q0 * piQ**3 * sct%gam(iband1,info%ik,is))

          calc_alpha =  real(PolyGamma(1,iband1,info%ik,is)) * enrgy(iband1,is) &
                     -  real(PolyGamma(2,iband1,info%ik,is)) * enrgy(iband1,is) * sct%gam(iband1,info%ik,is) * info%beta2pQ &
                     - aimag(PolyGamma(2,iband1,info%ik,is)) * sct%gam(iband1,info%ik,is)**2 * info%beta2pQ

          calc_alpha = calc_alpha &
                     * sct%zqp(iband1,info%ik,is)**2 * info%betaQ &
                     / (4.q0 * piQ**3 * sct%gam(iband1,info%ik,is))

          calc_xi    =  real(PolyGamma(1,iband1,info%ik,is)) &
                        * (enrgy(iband1,is)**2 + sct%gam(iband1,info%ik,is)**2) &
                     +  real(PolyGamma(2,iband1,info%ik,is)) &
                        * info%beta2pQ * sct%gam(iband1,info%ik,is) * (sct%gam(iband1,info%ik,is)**2 - enrgy(iband1,is)**2) &
                     - aimag(PolyGamma(2,iband1,info%ik,is)) &
                        * info%betaQ / piQ * enrgy(iband1,is) * sct%gam(iband1,info%ik,is)**2

          calc_xi    = calc_xi &
                       * sct%zqp(iband1,info%ik,is)**2 * info%betaQ &
                       / (4.q0 * piQ**3 * sct%gam(iband1,info%ik,is))

        else

          calc_sigma = real((enrgydiff(is)**2 + sct%gam(iband2,info%ik,is)**2 - sct%gam(iband1,info%ik,is)**2 &
                             - 2*ciQ*sct%gam(iband1,info%ik,is)*(-enrgydiff(is))*sct%gam(iband1,info%ik,is)) &
                           * sct%gam(iband2,info%ik,is) * PolyGamma(1,iband1,info%ik,is)) &
                     + real((enrgydiff(is)**2 + sct%gam(iband1,info%ik,is)**2 - sct%gam(iband2,info%ik,is)**2 &
                             - 2*ciQ*sct%gam(iband2,info%ik,is)*enrgydiff(is)*sct%gam(iband2,info%ik,is)) &
                           * sct%gam(iband1,info%ik,is) * PolyGamma(1,iband2,info%ik,is))

          calc_sigma = calc_sigma * sct%zqp(iband1,info%ik,is) * sct%zqp(iband2,info%ik,is) * info%betaQ &
                     / (2.d0 * piQ**3 * ( enrgydiff(is)**2 + (sct%gam(iband1,info%ik,is) - sct%gam(iband2,info%ik,is))**2)) &
                     / ( enrgydiff(is)**2 + (sct%gam(iband1,info%ik,is) + sct%gam(iband2,info%ik,is))**2)


          calc_alpha = real((enrgydiff(is)**2 + sct%gam(iband2,info%ik,is)**2 - sct%gam(iband1,info%ik,is)**2 &
                             - 2*ciQ*sct%gam(iband1,info%ik,is)*(-enrgydiff(is))*sct%gam(iband1,info%ik,is)) &
                           * (enrgy(iband1,is)-ciQ*sct%gam(iband1,info%ik,is)) &
                           * sct%gam(iband2,info%ik,is) * PolyGamma(1,iband1,info%ik,is)) &
                     + real((enrgydiff(is)**2 + sct%gam(iband1,info%ik,is)**2 - sct%gam(iband2,info%ik,is)**2 &
                             - 2*ciQ*sct%gam(iband2,info%ik,is)*enrgydiff(is)*sct%gam(iband2,info%ik,is)) &
                           * (enrgy(iband2,is)-ciQ*sct%gam(iband2,info%ik,is)) &
                           * sct%gam(iband1,info%ik,is) * PolyGamma(1,iband2,info%ik,is))

          calc_alpha = calc_alpha * sct%zqp(iband1,info%ik,is) * sct%zqp(iband2,info%ik,is) * info%betaQ &
                     / (2.d0 * piQ**3 * ( enrgydiff(is)**2 + (sct%gam(iband1,info%ik,is) - sct%gam(iband2,info%ik,is))**2)) &
                     / ( enrgydiff(is)**2 + (sct%gam(iband1,info%ik,is) + sct%gam(iband2,info%ik,is))**2)

          calc_xi    = real((enrgydiff(is)**2 + sct%gam(iband2,info%ik,is)**2 - sct%gam(iband1,info%ik,is)**2 &
                             - 2*ciQ*sct%gam(iband1,info%ik,is)*(-enrgydiff(is))*sct%gam(iband1,info%ik,is)) &
                           * (enrgy(iband1,is)-ciQ*sct%gam(iband1,info%ik,is))**2 &
                           * sct%gam(iband2,info%ik,is) * PolyGamma(1,iband1,info%ik,is)) &
                     + real((enrgydiff(is)**2 + sct%gam(iband1,info%ik,is)**2 - sct%gam(iband2,info%ik,is)**2 &
                             - 2*ciQ*sct%gam(iband2,info%ik,is)*enrgydiff(is)*sct%gam(iband2,info%ik,is)) &
                           * (enrgy(iband2,is)-ciQ*sct%gam(iband2,info%ik,is))**2 &
                           * sct%gam(iband1,info%ik,is) * PolyGamma(1,iband2,info%ik,is))

          calc_xi    = calc_xi * sct%zqp(iband1,info%ik,is) * sct%zqp(iband2,info%ik,is) * info%beta &
                     / (2.d0 * piQ**3 * ( enrgydiff(is)**2 + (sct%gam(iband1,info%ik,is) - sct%gam(iband2,info%ik,is))**2)) &
                     / ( enrgydiff(is)**2 + (sct%gam(iband1,info%ik,is) + sct%gam(iband2,info%ik,is))**2)
        endif

        if (algo%ldebug .and. (index(algo%dbgstr,"KernelsOnly") .ne. 0)) then
          cycle
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

subroutine response_intra_optical_weights_Q(algo, resp, edisp, info)
  implicit none
  type(algorithm)    :: algo
  type(response_qp)  :: resp
  type(energydisp)   :: edisp
  type(runinfo)      :: info

  integer :: index1(9), index2(9)
  integer :: iband, idir

  if (algo%ldebug .and. (index(algo%dbgstr,"KernelsOnly") .ne. 0)) then
    return
  endif

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

subroutine levicivita_peierls(dir1,dir2,dir3, sign1,vdir1,mdir1,sign2,vdir2,mdir2)
  ! for a given B-field quantities identified by
  ! dir1 dir2 dir
  ! return the required directions for the band derivatives / curvatures

  implicit none
  integer, intent(in)  :: dir1,dir2,dir3
  integer, intent(out) :: sign1,vdir1,mdir1, sign2,vdir2,mdir2

  integer :: i,j,eps

  sign1 = 0
  sign2 = 0

  ! sigma _abc = epsilon cij * va * vj * M^-1 bi

  do i=1,3
    do j=1,3
      eps = levicivita(dir3,i,j)
      if (eps == 0) then
        continue
      else
        if (sign1 == 0) then
          sign1 = eps
          vdir1 = index2compound(dir1,j)
          mdir1 = index2compound(dir2,i)
        else
          sign2 = eps
          vdir2 = index2compound(dir1,j)
          mdir2 = index2compound(dir2,i)

          goto 100 ! completely break out of all loops
        endif
      endif
    enddo
  enddo

100 return

end subroutine levicivita_peierls

integer function levicivita(dir1,dir2,dir3)
  implicit none
  integer, intent(in) :: dir1, dir2, dir3
  ! dir1 dir2 dir3 in the range [0,2]
  if ((dir1==dir2) .or. (dir2==dir3) .or. (dir1==dir3)) then
    levicivita =  0
    return
  else
    if (dir1 == 1) then
      if (dir2 == 2) then
        levicivita =  1  ! 123
        return
      else
        levicivita =  -1 ! 132
        return
      endif
    else if (dir1 == 2) then
      if (dir2 == 3) then
        levicivita =  1  ! 231
        return
      else
        levicivita =  -1 ! 213
        return
      endif
    else if (dir1 == 3) then
      if (dir2 == 1) then
        levicivita =  1  ! 312
        return
      else
        levicivita =  -1 ! 321
        return
      endif
    endif
  endif
end function levicivita

integer function index2compound(dir1,dir2)
  implicit none
  integer, intent(in) :: dir1,dir2

  ! maps out
  ! | 1 4 5 |
  ! | 4 2 6 |
  ! | 5 6 3 |

  if (dir1==dir2) then
    index2compound = dir1
  else if ((dir1 == 1 .and. dir2 == 2) .or. (dir1 == 2 .and. dir2 == 1)) then
    index2compound = 4
  else if ((dir1 == 1 .and. dir2 == 3) .or. (dir1 == 3 .and. dir2 == 1)) then
    index2compound = 5
  else if ((dir1 == 2 .and. dir2 == 3) .or. (dir1 == 3 .and. dir2 == 2)) then
    index2compound = 6
  endif

end function index2compound

end module Mresponse
