module Mantiresponse
  use Mmpi_org
  use Mtypes
  use Mparams
  use Mfermi


implicit none

  interface calc_digamma
    module procedure calc_digamma_D, calc_digamma_Q
  end interface calc_digamma

  interface calc_digamma_lim
    module procedure calc_digamma_lim_D, calc_digamma_Q
  end interface calc_digamma_lim

  interface calc_polygamma_lim
    module procedure calc_polygamma_lim_D, calc_polygamma_lim_Q
  end interface calc_polygamma_lim

  contains


subroutine calc_digamma_D(DiGamma, edisp, sct, kmesh, algo, info) !digamma (psi_0) function double precision 
  implicit none
  type(algorithm)  :: algo
  type(energydisp) :: edisp
  type(kpointmesh) :: kmesh
  type(scattering) :: sct
  type(runinfo)    :: info

  complex(8) :: DiGamma(1,edisp%nbopt_min:edisp%nbopt_max, ikstr:ikend, edisp%ispin)

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
        DiGamma(1,iband,ik,is) = wpsipg(to_evaluate(iband,ik,is),0)
      enddo
    enddo
  enddo

  deallocate(to_evaluate)

end subroutine calc_digamma_D

subroutine calc_digamma_Q(DiGammaQ, edisp, sct, kmesh, algo, info) !digamma (psi_0) function quadruple precision 
  implicit none
  type(algorithm)  :: algo
  type(energydisp) :: edisp
  type(kpointmesh) :: kmesh
  type(scattering) :: sct
  type(runinfo)    :: info

  complex(16) :: DiGammaQ(1,edisp%nbopt_min:edisp%nbopt_max, ikstr:ikend, edisp%ispin)

  complex(16), external :: wpsipghp
  complex(16), allocatable :: to_evaluate(:,:,:)
  integer :: iband, ik, is

  allocate(to_evaluate(edisp%nbopt_min:edisp%nbopt_max, ikstr:ikend, edisp%ispin))

  to_evaluate = 0.5q0 + info%beta2pQ * &
                (sct%gam(edisp%nbopt_min:edisp%nbopt_max,ikstr:ikend,:) &
                 + ciQ*sct%zqp(edisp%nbopt_min:edisp%nbopt_max,ikstr:ikend,:) &
                 * (edisp%band(edisp%nbopt_min:edisp%nbopt_max,ikstr:ikend,:) - info%muQ))

  do is = 1,edisp%ispin
    do ik = ikstr, ikend
      do iband = edisp%nbopt_min,edisp%nbopt_max
        DiGammaQ(1,iband,ik,is) = wpsipghp(to_evaluate(iband,ik,is),0)
      enddo
    enddo
  enddo

  deallocate(to_evaluate)

end subroutine calc_digamma_Q

subroutine calc_digamma_lim_D(DiGamLim, edisp, sct, info) !T=0 lim of digamma function double precision
  implicit none
  type(energydisp) :: edisp
  type(scattering) :: sct
  type(runinfo)    :: info

  complex(8) :: DiGamLim(1,edisp%nbopt_min:edisp%nbopt_max, ikstr:ikend, edisp%ispin)

  complex(8), allocatable :: limenergy(:,:,:)
  complex(8), allocatable :: limgam(:,:,:)
  integer :: iband, ik, is

  allocate(limenergy(edisp%nbopt_min:edisp%nbopt_max, ikstr:ikend, edisp%ispin))
  allocate(limgam(edisp%nbopt_min:edisp%nbopt_max, ikstr:ikend, edisp%ispin))

  limenergy = sct%zqp(edisp%nbopt_min:edisp%nbopt_max,ikstr:ikend,:) &
                 * (edisp%band(edisp%nbopt_min:edisp%nbopt_max,ikstr:ikend,:) - info%mu)
  limgam=sct%gam(edisp%nbopt_min:edisp%nbopt_max,ikstr:ikend,:)

  do is = 1,edisp%ispin
    do ik = ikstr, ikend
      do iband = edisp%nbopt_min,edisp%nbopt_max
        DiGamLim(1,iband,ik,is) = log(sqrt(real(limgam(iband,ik,is)**2+limenergy(iband,ik,is)**2)))&
                                  +ci*atan(real(limenergy(iband,ik,is)/limgam(iband,ik,is)))
      enddo
    enddo
  enddo
  deallocate(limenergy)
  deallocate(limgam)

end subroutine calc_digamma_lim_D

subroutine calc_digamma_lim_Q(DiGamLim, edisp, sct, info) !T=0 lim of digamma function quadruple precision
  implicit none
  type(energydisp) :: edisp
  type(scattering) :: sct
  type(runinfo)    :: info

  complex(16) :: DiGamLim(1,edisp%nbopt_min:edisp%nbopt_max, ikstr:ikend, edisp%ispin)

  complex(16), allocatable :: limenergy(:,:,:)
  complex(16), allocatable :: limgam(:,:,:)
  integer :: iband, ik, is

  allocate(limenergy(edisp%nbopt_min:edisp%nbopt_max, ikstr:ikend, edisp%ispin))
  allocate(limgam(edisp%nbopt_min:edisp%nbopt_max, ikstr:ikend, edisp%ispin))

  limenergy = sct%zqp(edisp%nbopt_min:edisp%nbopt_max,ikstr:ikend,:) &
                 * (edisp%band(edisp%nbopt_min:edisp%nbopt_max,ikstr:ikend,:) - info%muQ)
  limgam=sct%gam(edisp%nbopt_min:edisp%nbopt_max,ikstr:ikend,:)

  do is = 1,edisp%ispin
    do ik = ikstr, ikend
      do iband = edisp%nbopt_min,edisp%nbopt_max
        DiGamLim(1,iband,ik,is) = log(sqrt(real(limgam(iband,ik,is)**2+limenergy(iband,ik,is)**2)))&
                                  +ci*atan(real(limenergy(iband,ik,is)/limgam(iband,ik,is)))
      enddo
    enddo
  enddo
  deallocate(limenergy)
  deallocate(limgam)

end subroutine calc_digamma_lim_Q

subroutine calc_polygamma_lim_D(PolyGamLim, edisp, sct, info) !T=0 lim of Psi_1 function double precision
  implicit none
  type(energydisp) :: edisp
  type(scattering) :: sct
  type(runinfo)    :: info

  complex(8) :: PolyGamLim(1,edisp%nbopt_min:edisp%nbopt_max, ikstr:ikend, edisp%ispin)

  complex(8), allocatable :: limenergy(:,:,:)
  complex(8), allocatable :: limgam(:,:,:)
  integer :: iband, ik, is

  allocate(limenergy(edisp%nbopt_min:edisp%nbopt_max, ikstr:ikend, edisp%ispin))
  allocate(limgam(edisp%nbopt_min:edisp%nbopt_max, ikstr:ikend, edisp%ispin))

  limenergy = sct%zqp(edisp%nbopt_min:edisp%nbopt_max,ikstr:ikend,:) &
                 * (edisp%band(edisp%nbopt_min:edisp%nbopt_max,ikstr:ikend,:) - info%mu)
  limgam=sct%gam(edisp%nbopt_min:edisp%nbopt_max,ikstr:ikend,:)
  
  do is = 1,edisp%ispin
    do ik = ikstr, ikend
      do iband = edisp%nbopt_min,edisp%nbopt_max
        PolyGamLim(1,iband,ik,is) = 1/(limgam(iband,ik,is)**2+limenergy(iband,ik,is)**2)&
                                    *(limgam(iband,ik,is)-ci*limenergy(iband,ik,is))
      enddo
    enddo
  enddo
  deallocate(limenergy)
  deallocate(limgam)

end subroutine calc_polygamma_lim_D

subroutine calc_polygamma_lim_Q(PolyGamLim, edisp, sct, info) !T=0 lim of Psi_1 function quadruple precision
  implicit none
  type(energydisp) :: edisp
  type(scattering) :: sct
  type(runinfo)    :: info

  complex(16) :: PolyGamLim(1,edisp%nbopt_min:edisp%nbopt_max, ikstr:ikend, edisp%ispin)

  complex(16), allocatable :: limenergy(:,:,:)
  complex(16), allocatable :: limgam(:,:,:)
  integer :: iband, ik, is

  allocate(limenergy(edisp%nbopt_min:edisp%nbopt_max, ikstr:ikend, edisp%ispin))
  allocate(limgam(edisp%nbopt_min:edisp%nbopt_max, ikstr:ikend, edisp%ispin))

  limenergy = sct%zqp(edisp%nbopt_min:edisp%nbopt_max,ikstr:ikend,:) &
                 * (edisp%band(edisp%nbopt_min:edisp%nbopt_max,ikstr:ikend,:) - info%muQ)
  limgam=sct%gam(edisp%nbopt_min:edisp%nbopt_max,ikstr:ikend,:)
  
  do is = 1,edisp%ispin
    do ik = ikstr, ikend
      do iband = edisp%nbopt_min,edisp%nbopt_max
        PolyGamLim(1,iband,ik,is) = 1/(limgam(iband,ik,is)**2+limenergy(iband,ik,is)**2)&
                                    *(limgam(iband,ik,is)-ci*limenergy(iband,ik,is))
      enddo
    enddo
  enddo
  deallocate(limenergy)
  deallocate(limgam)

end subroutine calc_polygamma_lim_Q


subroutine response_inter_anti_km(resp, DiGamma, PolyGamma, DiGamLim, PolyGamLim, edisp, sct, kmesh, algo, info)
  implicit none
  type (response_dp)  :: resp

  type(energydisp)    :: edisp
  type(scattering)    :: sct
  type(kpointmesh)    :: kmesh
  type(algorithm)     :: algo
  type(runinfo)       :: info

  complex(8)          :: PolyGamma(3,edisp%nbopt_min:edisp%nbopt_max, ikstr:ikend, edisp%ispin)
  complex(8)          :: DiGamma(1,edisp%nbopt_min:edisp%nbopt_max, ikstr:ikend, edisp%ispin)
  complex(8)          :: DiGamLim(1,edisp%nbopt_min:edisp%nbopt_max, ikstr:ikend, edisp%ispin)
  complex(8)          :: PolyGamLim(1,edisp%nbopt_min:edisp%nbopt_max, ikstr:ikend, edisp%ispin)

  real(8) :: zqp
  real(8) :: gam
  real(8), allocatable :: enrgy(:,:)
  real(8), allocatable :: enrgydiff(:)
  real(8), allocatable :: gamdiff(:)
  real(8), allocatable :: gamplus(:)

  real(8), parameter :: tol = 1e-16

  complex(8) :: M0_11
  complex(8) :: M1_11
  complex(8) :: M2_11
  complex(8) :: M3_11

  complex(8) :: M0_12
  complex(8) :: M1_12
  complex(8) :: M2_12
  complex(8) :: M3_12

  complex(8) :: M0_22
  complex(8) :: M1_22
  complex(8) :: M2_22
  complex(8) :: M3_22

  complex(8) :: calc_sigma
  complex(8) :: calc_alpha
  complex(8) :: calc_xi

  integer :: index1(3), index2(3)
  integer :: i,j,idir,idir1,idir2,idir3
  integer :: iband1, iband2, iband, is


  allocate(enrgy(edisp%nbopt_min:edisp%nbopt_max,edisp%ispin))
  allocate(enrgydiff(edisp%ispin))
  allocate(gamdiff(edisp%ispin))
  allocate(gamplus(edisp%ispin))

  enrgy = sct%zqp(edisp%nbopt_min:edisp%nbopt_max,info%ik,:) &
          * (edisp%band(edisp%nbopt_min:edisp%nbopt_max,info%ik,:) - info%mu)

  do iband1 = edisp%nbopt_min, edisp%nbopt_max
    do iband2 = edisp%nbopt_min, edisp%nbopt_max
      if (iband1 == iband2) cycle
      enrgydiff = enrgy(iband1,:) - enrgy(iband2,:)
      gamdiff   = sct%gam(iband1,info%ik,:) - sct%gam(iband2,info%ik,:)
      gamplus   = sct%gam(iband1,info%ik,:) + sct%gam(iband2,info%ik,:)

      do is = 1,edisp%ispin

        !L11-anti kernel calculation 

          M0_11 = ci/sct%gam(iband1,info%ik,is)-1.d0/(enrgydiff(is)+ci*gamdiff(is)) &
                        -1.d0/(enrgydiff(is)+ci*gamplus(is))

          M0_11=M0_11*1.d0/(2.d0*ci*sct%gam(iband1,info%ik,is)**2)*1.d0/(enrgydiff(is) &
                        +ci*gamplus(is))

          M0_11=M0_11+1.d0/(sct%gam(iband2,info%ik,is)*(enrgydiff(is)+ci*gamdiff(is)) &
                        *(enrgydiff(is)-ci*gamplus(is))**2)

          M0_11=M0_11*(-pi/(4.d0*(enrgydiff(is)+ci*gamdiff(is))))

          M1_11=1.d0/((enrgydiff(is)-ci*gamdiff(is))**2*(enrgydiff(is)+ci*gamplus(is))**2 &
                        *2.d0*sct%gam(iband2,info%ik,is)*ci)

          M2_11=1.d0/(ci*sct%gam(iband1,info%ik,is))-1.d0/(enrgydiff(is)-ci*gamplus(is)) &
                        -1.d0/(enrgydiff(is)-ci*gamdiff(is))

          M2_11=M2_11*1.d0/(4.d0*sct%gam(iband1,info%ik,is)**2*(enrgydiff(is)-ci*gamdiff(is)) &
                         * (enrgydiff(is)-ci*gamplus(is)))

          M3_11=1.d0/(4.d0*sct%gam(iband1,info%ik,is)**2*(enrgydiff(is)-ci*gamplus(is)) &
                          *(enrgydiff(is)-ci*gamdiff(is)))

          
          if (abs(info%Temp) <= abs(tol)) then !zero temepature lowest order expansion of kernel function
            calc_sigma = real(M1_11*DiGamLim(1,iband2,info%ik,is))+real(M2_11*DiGamLim(1,iband1,info%ik,is)) &
                              -aimag(M3_11*PolyGamLim(1,iband1,info%ik,is))
          else
            calc_sigma = real(M1_11*DiGamma(1,iband2,info%ik,is))+real(M2_11*DiGamma(1,iband1,info%ik,is)) &
                            -info%beta2p*aimag(M3_11*PolyGamma(1,iband1,info%ik,is))
          endif

          calc_sigma=2.d0*pi*2.d0*sct%zqp(iband1,info%ik,is)**2 * sct%zqp(iband2,info%ik,is) * &
                            sct%gam(iband1,info%ik,is)**2*sct%gam(iband2,info%ik,is)/pi**3* &
                            (M0_11+calc_sigma)
                            


        !L12-anti kernel calculation

          M0_12 = ci/sct%gam(iband1,info%ik,is)-1.d0/(enrgydiff(is)+ci*gamdiff(is)) &
                        -1.d0/(enrgydiff(is)+ci*gamplus(is))+1.d0/(enrgy(iband1,is) &
                        +ci*sct%gam(iband1,info%ik,is))

          M0_12=M0_12*(enrgy(iband1,is)+ci*sct%gam(iband1,info%ik,is))/(2.d0*ci*sct%gam(iband1,info%ik,is)**2)&
                        *1.d0/(enrgydiff(is) +ci*gamplus(is))

          M0_12=M0_12+(enrgy(iband2,is)+ci*sct%gam(iband2,info%ik,is))/(sct%gam(iband2,info%ik,is)* &
                        (enrgydiff(is)+ci*gamdiff(is))*(enrgydiff(is)-ci*gamplus(is))**2)

          M0_12=M0_12*(-pi/(4.d0*(enrgydiff(is)+ci*gamdiff(is))))

          M1_12=(enrgy(iband2,is)-ci*sct%gam(iband2,info%ik,is))/((enrgydiff(is)-ci*gamdiff(is))**2 &
                       *(enrgydiff(is)+ci*gamplus(is))**2*2.d0*sct%gam(iband2,info%ik,is)*ci)

          M2_12=1.d0/(ci*sct%gam(iband1,info%ik,is))-1.d0/(enrgydiff(is)-ci*gamdiff(is)) &
                        -1.d0/(enrgydiff(is)-ci*gamplus(is))+1.d0/(enrgy(iband1,is)-ci*sct%gam(iband1,info%ik,is))

          M2_12=M2_12*(enrgy(iband1,is)-ci*sct%gam(iband1,info%ik,is))/(4.d0*sct%gam(iband1,info%ik,is)**2 &
                         *(enrgydiff(is)-ci*gamdiff(is))*(enrgydiff(is)-ci*gamplus(is)))

          M3_12=(enrgy(iband1,is)-ci*sct%gam(iband1,info%ik,is))/(4.d0*sct%gam(iband1,info%ik,is)**2 &
                          *(enrgydiff(is)-ci*gamplus(is))*(enrgydiff(is)-ci*gamdiff(is)))

          
          if (abs(info%Temp) <= abs(tol)) then !zero temepature lowest order expansion of kernel function
            calc_alpha = real(M1_12*DiGamLim(1,iband2,info%ik,is))+real(M2_12*DiGamLim(1,iband1,info%ik,is)) &
                            -aimag(M3_12*PolyGamLim(1,iband1,info%ik,is))
          else  
            calc_alpha = real(M1_12*DiGamma(1,iband2,info%ik,is))+real(M2_12*DiGamma(1,iband1,info%ik,is)) &
                            -info%beta2p*aimag(M3_12*PolyGamma(1,iband1,info%ik,is))
          endif
          
          calc_alpha=2.d0*pi*2.d0*sct%zqp(iband1,info%ik,is)**2 * sct%zqp(iband2,info%ik,is) * &
                            sct%gam(iband1,info%ik,is)**2*sct%gam(iband2,info%ik,is)/pi**3* &
                            (M0_12+calc_alpha)

        !L22-anti kernel calculation

          M0_22 = ci/sct%gam(iband1,info%ik,is)-1.d0/(enrgydiff(is)+ci*gamdiff(is)) &
                        -1.d0/(enrgydiff(is)+ci*gamplus(is))+2.d0/(enrgy(iband1,is)&
                        +ci*sct%gam(iband1,info%ik,is))

          M0_22=M0_22*(enrgy(iband1,is)+ci*sct%gam(iband1,info%ik,is))**2/ &
                        (2.d0*ci*sct%gam(iband1,info%ik,is)**2) &
                        *1.d0/(enrgydiff(is)+ci*gamplus(is))

          M0_22=M0_22+(enrgy(iband2,is)+ci*sct%gam(iband2,info%ik,is))**2/(sct%gam(iband2,info%ik,is) &
                        *(enrgydiff(is)+ci*gamdiff(is))*(enrgydiff(is)-ci*gamplus(is))**2)

          M0_22=M0_22*(-pi/(4.d0*(enrgydiff(is)+ci*gamdiff(is))))

          M1_22=(enrgy(iband2,is)-ci*sct%gam(iband2,info%ik,is))**2/((enrgydiff(is)-ci*gamdiff(is))**2 &
                        *(enrgydiff(is)+ci*gamplus(is))**2*2.d0*sct%gam(iband2,info%ik,is)*ci)

          M2_22=1.d0/(ci*sct%gam(iband1,info%ik,is))-1.d0/(enrgydiff(is)-ci*gamdiff(is)) &
                        -1.d0/(enrgydiff(is)-ci*gamplus(is))+2.d0/(enrgy(iband1,is)-ci*sct%gam(iband1,info%ik,is))

          M2_22=M2_22*(enrgy(iband1,is)-ci*sct%gam(iband1,info%ik,is))**2/(4.d0*sct%gam(iband1,info%ik,is)**2 &
                        *(enrgydiff(is)-ci*gamdiff(is))*(enrgydiff(is)-ci*gamplus(is)))

          M3_22=(enrgy(iband1,is)-ci*sct%gam(iband1,info%ik,is))**2/(4.d0*sct%gam(iband1,info%ik,is)**2 &
                          *(enrgydiff(is)-ci*gamplus(is))*(enrgydiff(is)-ci*gamdiff(is)))

          
          if (abs(info%Temp) <= abs(tol)) then !zero temepature lowest order expansion of kernel function
            calc_xi = real(M1_22*DiGamLim(1,iband2,info%ik,is))+real(M2_22*DiGamLim(1,iband1,info%ik,is)) &
                            -aimag(M3_22*PolyGamLim(1,iband1,info%ik,is))
          else
            calc_xi = real(M1_22*DiGamma(1,iband2,info%ik,is))+real(M2_22*DiGamma(1,iband1,info%ik,is)) &
                            -info%beta2p*aimag(M3_22*PolyGamma(1,iband1,info%ik,is))
          endif

          calc_xi=2.d0*pi*2.d0*sct%zqp(iband1,info%ik,is)**2 * sct%zqp(iband2,info%ik,is) * &
                            sct%gam(iband1,info%ik,is)**2*sct%gam(iband2,info%ik,is)/pi**3* &
                            (M0_22+calc_xi)


        if (algo%ldebug .and. (index(algo%dbgstr,"KernelsOnly") .ne. 0)) then
          cycle
        endif

        !only have three non-zero Berry curvature components to read in in the order

        ! |- 1 2 |
        ! |- - 3 |
        ! |- - - | 

        index1 = (/1,1,2/)
        index2 = (/2,3,3/)
        
        !(xy,xz,yz)

        ! multiply by Berry curvature element

        do idir = 1,3

          resp%s_full(index1(idir),index2(idir),iband1,is,info%ik) = &
          resp%s_full(index1(idir),index2(idir),iband1,is,info%ik) + calc_sigma &
              * edisp%BerryCurv(idir,iband2,iband1,is,info%ik)

          resp%a_full(index1(idir),index2(idir),iband1,is,info%ik) = &
          resp%a_full(index1(idir),index2(idir),iband1,is,info%ik) + calc_alpha &
              * edisp%BerryCurv(idir,iband2,iband1,is,info%ik)

          resp%x_full(index1(idir),index2(idir),iband1,is,info%ik) = &
          resp%x_full(index1(idir),index2(idir),iband1,is,info%ik) + calc_xi &
              * edisp%BerryCurv(idir,iband2,iband1,is,info%ik)
        enddo
      enddo
    enddo
  enddo


  do iband = edisp%nbopt_min, edisp%nbopt_max ! this is anti-symmetric transport so we need to anti-symmetrize this 

    resp%s_full(2,1,iband,:,info%ik) = -1.d0*resp%s_full(1,2,iband,:,info%ik)
    resp%s_full(3,1,iband,:,info%ik) = -1.d0*resp%s_full(1,3,iband,:,info%ik)
    resp%s_full(3,2,iband,:,info%ik) = -1.d0*resp%s_full(2,3,iband,:,info%ik)

    resp%a_full(2,1,iband,:,info%ik) = -1.d0*resp%a_full(1,2,iband,:,info%ik)
    resp%a_full(3,1,iband,:,info%ik) = -1.d0*resp%a_full(1,3,iband,:,info%ik)
    resp%a_full(3,2,iband,:,info%ik) = -1.d0*resp%a_full(2,3,iband,:,info%ik)

    resp%x_full(2,1,iband,:,info%ik) = -1.d0*resp%x_full(1,2,iband,:,info%ik)
    resp%x_full(3,1,iband,:,info%ik) = -1.d0*resp%x_full(1,3,iband,:,info%ik)
    resp%x_full(3,2,iband,:,info%ik) = -1.d0*resp%x_full(2,3,iband,:,info%ik)
  enddo
  
  deallocate(enrgy)
  deallocate(enrgydiff)
  deallocate(gamdiff)
  deallocate(gamplus)

end subroutine response_inter_anti_km

subroutine response_inter_anti_km_Q(resp, DiGamma, PolyGamma,DiGamLim, PolyGamLim, edisp, sct, kmesh, algo, info)
  implicit none
  type (response_qp)  :: resp

  type(energydisp)    :: edisp
  type(scattering)    :: sct
  type(kpointmesh)    :: kmesh
  type(algorithm)     :: algo
  type(runinfo)       :: info

  complex(16)          :: PolyGamma(3,edisp%nbopt_min:edisp%nbopt_max, ikstr:ikend, edisp%ispin)
  complex(16)          :: DiGamma(1,edisp%nbopt_min:edisp%nbopt_max, ikstr:ikend, edisp%ispin)

  complex(16)          :: DiGamLim(1,edisp%nbopt_min:edisp%nbopt_max, ikstr:ikend, edisp%ispin)
  complex(16)          :: PolyGamLim(1,edisp%nbopt_min:edisp%nbopt_max, ikstr:ikend, edisp%ispin)

  real(16) :: zqp
  real(16) :: gam
  real(16), allocatable :: enrgy(:,:)
  real(16), allocatable :: enrgydiff(:)
  real(16), allocatable :: gamdiff(:)
  real(16), allocatable :: gamplus(:)

  complex(16) :: M0_11
  complex(16) :: M1_11
  complex(16) :: M2_11
  complex(16) :: M3_11

  complex(16) :: M0_12
  complex(16) :: M1_12
  complex(16) :: M2_12
  complex(16) :: M3_12

  complex(16) :: M0_22
  complex(16) :: M1_22
  complex(16) :: M2_22
  complex(16) :: M3_22

  complex(16) :: calc_sigma
  complex(16) :: calc_alpha
  complex(16) :: calc_xi

  integer :: index1(3), index2(3)
  integer :: i,j,idir,idir1,idir2,idir3
  integer :: iband1, iband2, iband, is

  real(16), parameter :: tol = 1e-16


  allocate(enrgy(edisp%nbopt_min:edisp%nbopt_max,edisp%ispin))
  allocate(enrgydiff(edisp%ispin))
  allocate(gamdiff(edisp%ispin))
  allocate(gamplus(edisp%ispin))

  enrgy = sct%zqp(edisp%nbopt_min:edisp%nbopt_max,info%ik,:) &
          * (edisp%band(edisp%nbopt_min:edisp%nbopt_max,info%ik,:) - info%muQ)

  do iband1 = edisp%nbopt_min, edisp%nbopt_max
    do iband2 = edisp%nbopt_min, edisp%nbopt_max
      if (iband1 == iband2) cycle
      enrgydiff = enrgy(iband1,:) - enrgy(iband2,:)
      gamdiff   = sct%gam(iband1,info%ik,:) - sct%gam(iband2,info%ik,:)
      gamplus   = sct%gam(iband1,info%ik,:) + sct%gam(iband2,info%ik,:)

      do is = 1,edisp%ispin

        !L11-anti kernel calculation 

          M0_11 = ci/sct%gam(iband1,info%ik,is)-1.q0/(enrgydiff(is)+ci*gamdiff(is)) &
                        -1.q0/(enrgydiff(is)+ci*gamplus(is))

          M0_11=M0_11*1.q0/(2.q0*ci*sct%gam(iband1,info%ik,is)**2)*1.q0/(enrgydiff(is) &
                        +ci*gamplus(is))

          M0_11=M0_11+1.q0/(sct%gam(iband2,info%ik,is)*(enrgydiff(is)+ci*gamdiff(is)) &
                        *(enrgydiff(is)-ci*gamplus(is))**2)

          M0_11=M0_11*(-piQ/(4.q0*(enrgydiff(is)+ci*gamdiff(is))))

          M1_11=1.q0/((enrgydiff(is)-ci*gamdiff(is))**2*(enrgydiff(is)+ci*gamplus(is))**2 &
                        *2.q0*sct%gam(iband2,info%ik,is)*ci)

          M2_11=1.q0/(ci*sct%gam(iband1,info%ik,is))-1.q0/(enrgydiff(is)-ci*gamdiff(is)) &
                        -1.q0/(enrgydiff(is)-ci*gamdiff(is))

          M2_11=M2_11*1.q0/(4.q0*sct%gam(iband1,info%ik,is)**2*(enrgydiff(is)-ci*gamdiff(is)) &
                         * (enrgydiff(is)-ci*gamplus(is)))

          M3_11=1.q0/(4.q0*sct%gam(iband1,info%ik,is)**2*(enrgydiff(is)-ci*gamplus(is)) &
                          *(enrgydiff(is)-ci*gamdiff(is)))

          
          if (abs(info%Temp) <= abs(tol)) then !zero temepature lowest order expansion of kernel function
            calc_sigma = real(M1_11*DiGamLim(1,iband2,info%ik,is))+real(M2_11*DiGamLim(1,iband1,info%ik,is)) &
                              -aimag(M3_11*PolyGamLim(1,iband1,info%ik,is))
          else
            calc_sigma = real(M1_11*DiGamma(1,iband2,info%ik,is))+real(M2_11*DiGamma(1,iband1,info%ik,is)) &
                            -info%beta2pQ*aimag(M3_11*PolyGamma(1,iband1,info%ik,is))
          endif
          
          calc_sigma=2.q0*piQ*2.q0*sct%zqp(iband1,info%ik,is)**2 * sct%zqp(iband2,info%ik,is) * &
                            sct%gam(iband1,info%ik,is)**2*sct%gam(iband2,info%ik,is)/piQ**3* &
                            (M0_11+calc_sigma)
                            


        !L12-anti kernel calculation

          M0_12 = ci/sct%gam(iband1,info%ik,is)-1.q0/(enrgydiff(is)+ci*gamdiff(is)) &
                        -1.q0/(enrgydiff(is)+ci*gamplus(is))+1.q0/(enrgy(iband1,is) &
                        +ci*sct%gam(iband1,info%ik,is))

          M0_12=M0_12*(enrgy(iband1,is)+ci*sct%gam(iband1,info%ik,is))/ &
                        (2.q0*ci*sct%gam(iband1,info%ik,is)**2) &
                        *1.q0/(enrgydiff(is)+ci*gamplus(is))

          M0_12=M0_12+(enrgy(iband2,is)+ci*sct%gam(iband2,info%ik,is))/(sct%gam(iband2,info%ik,is) &
                       *(enrgydiff(is)+ci*gamdiff(is))*(enrgydiff(is)-ci*gamplus(is))**2)

          M0_12=M0_12*(-piQ/(4.q0*(enrgydiff(is)+ci*gamdiff(is))))

          M1_12=(enrgy(iband2,is)-ci*sct%gam(iband2,info%ik,is))/((enrgydiff(is)-ci*gamdiff(is))**2 &
                        *(enrgydiff(is)+ci*gamplus(is))**2*2.q0*sct%gam(iband2,info%ik,is)*ci)

          M2_12=1.q0/(ci*sct%gam(iband1,info%ik,is))-1.q0/(enrgydiff(is)-ci*gamdiff(is)) &
                        -1.q0/(enrgydiff(is)-ci*gamdiff(is))+1.q0/(enrgy(iband1,is) &
                        -ci*sct%gam(iband1,info%ik,is))

          M2_12=M2_12*(enrgy(iband1,is)-ci*sct%gam(iband1,info%ik,is))/(4.q0*sct%gam(iband1,info%ik,is)**2 &
                         *(enrgydiff(is)-ci*gamdiff(is))*(enrgydiff(is)-ci*gamplus(is)))

          M3_12=(enrgy(iband1,is)-ci*sct%gam(iband1,info%ik,is))/(4.q0*sct%gam(iband1,info%ik,is)**2 &
                          *(enrgydiff(is)-ci*gamplus(is))*(enrgydiff(is)-ci*gamdiff(is)))

          
          if (abs(info%Temp) <= abs(tol)) then !zero temepature lowest order expansion of kernel function
            calc_alpha = real(M1_12*DiGamLim(1,iband2,info%ik,is))+real(M2_12*DiGamLim(1,iband1,info%ik,is)) &
                            -aimag(M3_12*PolyGamLim(1,iband1,info%ik,is))
          else  
            calc_alpha = real(M1_12*DiGamma(1,iband2,info%ik,is))+real(M2_12*DiGamma(1,iband1,info%ik,is)) &
                            -info%beta2pQ*aimag(M3_12*PolyGamma(1,iband1,info%ik,is))
          endif
          
          calc_alpha=2.q0*piQ*2.q0*sct%zqp(iband1,info%ik,is)**2 * sct%zqp(iband2,info%ik,is) * &
                            sct%gam(iband1,info%ik,is)**2*sct%gam(iband2,info%ik,is)/piQ**3* &
                            (M0_11+calc_alpha)

        !L22-anti kernel calculation

          M0_22 = ci/sct%gam(iband1,info%ik,is)-1.q0/(enrgydiff(is)+ci*gamdiff(is)) &
                        -1.q0/(enrgydiff(is)+ci*gamplus(is))+2.q0/(enrgy(iband1,is) &
                        +ci*sct%gam(iband1,info%ik,is))

          M0_22=M0_22*(enrgy(iband1,is)+ci*sct%gam(iband1,info%ik,is))**2/(2.q0*ci*sct%gam(iband1,info%ik,is)**2) &
                        *1.q0/(enrgydiff(is)+ci*gamplus(is))

          M0_22=M0_22+(enrgy(iband2,is)+ci*sct%gam(iband2,info%ik,is))**2/(sct%gam(iband2,info%ik,is)* &
                        (enrgydiff(is)+ci*gamdiff(is))*(enrgydiff(is)-ci*gamplus(is))**2)

          M0_22=M0_22*(-pi/(4.q0*(enrgydiff(is)+ci*gamdiff(is))))

          M1_22=(enrgy(iband2,is)-ci*sct%gam(iband2,info%ik,is))**2/((enrgydiff(is)-ci*gamdiff(is))**2 &
                        *(enrgydiff(is)+ci*gamplus(is))**2*2.q0*sct%gam(iband2,info%ik,is)*ci)

          M2_22=1.q0/(ci*sct%gam(iband1,info%ik,is))-1.q0/(enrgydiff(is)-ci*gamdiff(is)) &
                        -1.q0/(enrgydiff(is)-ci*gamdiff(is))+2.q0/(enrgy(iband1,is)-ci*sct%gam(iband1,info%ik,is))

          M2_22=M2_22*(enrgy(iband1,is)-ci*sct%gam(iband1,info%ik,is))**2/(4.q0*sct%gam(iband1,info%ik,is)**2&
                         *(enrgydiff(is)-ci*gamdiff(is))*(enrgydiff(is)-ci*gamplus(is)))

          M3_22=(enrgy(iband1,is)-ci*sct%gam(iband1,info%ik,is))**2/(4.q0*sct%gam(iband1,info%ik,is)**2 &
                          *(enrgydiff(is)-ci*gamplus(is))*(enrgydiff(is)-ci*gamdiff(is)))

          

          if (abs(info%Temp) <= abs(tol)) then !zero temepature lowest order expansion of kernel function
            calc_xi = real(M1_22*DiGamLim(1,iband2,info%ik,is))+real(M2_22*DiGamLim(1,iband1,info%ik,is)) &
                            -aimag(M3_22*PolyGamLim(1,iband1,info%ik,is))
          else
            calc_xi = real(M1_22*DiGamma(1,iband2,info%ik,is))+real(M2_22*DiGamma(1,iband1,info%ik,is)) &
                            -info%beta2pQ*aimag(M3_22*PolyGamma(1,iband1,info%ik,is))
          endif
          
          calc_xi=2.q0*piQ*2.q0*sct%zqp(iband1,info%ik,is)**2 * sct%zqp(iband2,info%ik,is) * &
                            sct%gam(iband1,info%ik,is)**2*sct%gam(iband2,info%ik,is)/pi**3* &
                            (M0_11+calc_xi)


        if (algo%ldebug .and. (index(algo%dbgstr,"KernelsOnly") .ne. 0)) then
          cycle
        endif

        !only have three non-zero Berry curvature components to read in in the order

        ! |- 1 2 |
        ! |- - 3 |
        ! |- - - | 

        index1 = (/1,1,2/)
        index2 = (/2,3,3/)
        
        !(xy,xz,yz)

        ! multiply by Berry curvature 

        do idir = 1,3

          resp%s_full(index1(idir),index2(idir),iband1,is,info%ik) = &
          resp%s_full(index1(idir),index2(idir),iband1,is,info%ik) + calc_sigma &
              * edisp%BerryCurv(idir,iband2,iband1,is,info%ik)

          resp%a_full(index1(idir),index2(idir),iband1,is,info%ik) = &
          resp%a_full(index1(idir),index2(idir),iband1,is,info%ik) + calc_alpha &
              * edisp%BerryCurv(idir,iband2,iband1,is,info%ik)

          resp%x_full(index1(idir),index2(idir),iband1,is,info%ik) = &
          resp%x_full(index1(idir),index2(idir),iband1,is,info%ik) + calc_xi &
              * edisp%BerryCurv(idir,iband2,iband1,is,info%ik)
        enddo
      enddo
    enddo
  enddo


  do iband = edisp%nbopt_min, edisp%nbopt_max ! this is anti-symmetric transport so we need to anti-symmetrize this 

    resp%s_full(2,1,iband,:,info%ik) = -1.q0*resp%s_full(1,2,iband,:,info%ik)
    resp%s_full(3,1,iband,:,info%ik) = -1.q0*resp%s_full(1,3,iband,:,info%ik)
    resp%s_full(3,2,iband,:,info%ik) = -1.q0*resp%s_full(2,3,iband,:,info%ik)

    resp%a_full(2,1,iband,:,info%ik) = -1.q0*resp%a_full(1,2,iband,:,info%ik)
    resp%a_full(3,1,iband,:,info%ik) = -1.q0*resp%a_full(1,3,iband,:,info%ik)
    resp%a_full(3,2,iband,:,info%ik) = -1.q0*resp%a_full(2,3,iband,:,info%ik)

    resp%x_full(2,1,iband,:,info%ik) = -1.q0*resp%x_full(1,2,iband,:,info%ik)
    resp%x_full(3,1,iband,:,info%ik) = -1.q0*resp%x_full(1,3,iband,:,info%ik)
    resp%x_full(3,2,iband,:,info%ik) = -1.q0*resp%x_full(2,3,iband,:,info%ik)
  enddo

  deallocate(enrgy)
  deallocate(enrgydiff)
  deallocate(gamdiff)
  deallocate(gamplus)

end subroutine response_inter_anti_km_Q

end module Mantiresponse
