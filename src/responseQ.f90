module responseQ

contains

subroutine calc_responseQ(mu,it)


  use Mresponse
  implicit none
  real(8) :: mu
  integer it,ib
  real(16) :: beta2pQ,betaQ,facQ,facBQ

  betaQ=real(beta,16)
  beta2pQ=betaQ/(2.q0*piQ)

  s_respQ  = 0.q0
  sb_respQ = 0.q0
  a_respQ  = 0.q0
  ab_respQ = 0.q0

  !move somewhere more generic....
  nalpha=3
  if (icubic.eq.1) nalpha=1


  ! outer k-loop
  do ik=ikstart,ikend
     do ib=1,nband ! at this point everything is band-diagonal
        gammaQ=real(gam(it,ib),16)
        ! pre-compute all needed digamma functions.
        aqpQ=real(ek(ik,ib)-mu,16)
        zQ=real(zqp(ib),16)
        zargQ=0.5q0+beta2pQ*(ciQ*aqpQ+gammaQ)
        do ipg=1,3 ! do i need zeroth order? check peltier in B formula...    XXX
           ctmpQ=wpsipghp(zargQ,ipg)
           RePolyGammaQ(ipg)=real(ctmpQ,16)
           ImPolyGammaQ(ipg)=imag(ctmpQ)
        enddo

        ! compute transport kernels (omega-part)
        !

        tmpQ=zQ**2 / (4.q0*piQ**3) ! missing: beta/gamma (multiplied later, see below)

        s_kernelQ = tmpQ * (     RePolyGammaQ(1) - GammaQ*beta2pQ * RePolyGammaQ(2) )
        a_kernelQ = tmpQ * ( aqpQ * RePolyGammaQ(1) - aqpQ*GammaQ*beta2pQ*RePolyGammaQ(2) &
                                 - GammaQ**2.q0 * beta2pQ * ImPolyGammaQ(2) )


        tmpQ=tmpQ*3.q0*zQ/(4.q0*piQ) ! additionally missing: 1/gamma (multiplied later)

        sb_kernelQ = tmpQ * ( - RePolyGammaQ(1) + gammaQ * beta2pQ * RePolyGammaQ(2) &
                         + (beta2pQ*gammaQ)**2.q0/3.q0 * RePolyGammaQ(3)    )

        ab_kernelQ = tmpQ * ( aqpQ * RePolyGammaQ(1) - aqpQ * beta2pQ * gammaQ * RePolyGammaQ(2) &
                          + gammaQ**2.q0 / 3.q0 * beta2pQ * ImPolyGammaQ(2) &
             - aqpQ*gammaQ**2.q0 / 3.q0 * beta2pQ**2.q0 * ImPolyGammaQ(3) &
             + gammaQ**3.q0 / 3.q0 * beta2pQ**2.q0 * RePolyGammaQ(3) )


        ! multiply with matrix elements for all polarizations and add to final band-dependent response functions
        ! create an option for computing only polarization-diagonal elements in absence of B
        !

        ! B = 0
        do ialpha=1,nalpha
           do ibeta=ialpha,nalpha
              if ((icubic.eq.1).and.(ibeta.gt.ialpha)) exit
              if ((idiag.eq.1).AND.(ialpha.ne.ibeta)) exit

              tmpQ=real(vka(ik,ib,ialpha)*vka(ik,ib,ibeta),16)

              s_tmpQ(ib,ialpha,ibeta)=s_kernelQ * tmpQ
              a_tmpQ(ib,ialpha,ibeta)=a_kernelQ * tmpQ
           enddo !ibeta
        enddo ! ialpha


        ! B .ne. 0
        if (ivkab.eq.1) then
        do ialpha=1,nalpha
           do ibeta=ialpha+1,3
              if ((icubic.eq.1).and.(ibeta.gt.2)) exit

              tmpQ=real( vka(ik,ib,ialpha)*( vkab(ik,ib,ibeta,ialpha)*vka(ik,ib,ibeta) &
                          - vkab(ik,ib,ibeta,ibeta)*vka(ik,ib,ialpha)    ) ,16)

              sb_tmpQ(ib,ialpha,ibeta)=sb_kernelQ * tmpQ
              ab_tmpQ(ib,ialpha,ibeta)=ab_kernelQ * tmpQ

           enddo !ibeta
        enddo ! ialpha
        else
           sb_tmpQ=0.d0
           ab_tmpQ=0.d0
        endif
     enddo ! ib


     ! add to k--summed response functions (keep band-dependence here)
     !
     s_respQ  = s_respQ+s_tmpQ
     sb_respQ = sb_respQ+sb_tmpQ

     a_respQ  = a_respQ+a_tmpQ
     ab_respQ = ab_respQ+ab_tmpQ

  enddo ! outer k-loop

  ! reduce MPI_sum all contributions & output & sum over bands
  !


  if (nproc.gt.1) then

  do ib=1,nband
     do ialpha=1,3 ! reduce all here...
        do ibeta=1,3
           call MPI_REDUCE_QUAD(s_respQ(ib,ialpha,ibeta),s_tmpQ(ib,ialpha,ibeta))
           call MPI_REDUCE_QUAD(a_respQ(ib,ialpha,ibeta),a_tmpQ(ib,ialpha,ibeta))
           call MPI_REDUCE_QUAD(sb_respQ(ib,ialpha,ibeta),sb_tmpQ(ib,ialpha,ibeta))
           call MPI_REDUCE_QUAD(ab_respQ(ib,ialpha,ibeta),ab_tmpQ(ib,ialpha,ibeta))
        enddo
     enddo
  enddo

!  XXX
!  call MPI_REDUCE(s_resp,s_tmp,nband*3*3,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_WORLD,mpierr)
!  call MPI_REDUCE(sb_resp,sb_tmp,nband*3*3,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_WORLD,mpierr)
!  call MPI_BARRIER( MPI_COMM_WORLD, mpierr )
 !  call MPI_REDUCE(a_resp,a_tmp,nband*3*3,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_WORLD,mpierr)
!  call MPI_BARRIER( MPI_COMM_WORLD, mpierr )
!  call MPI_REDUCE(ab_resp,ab_tmp,nband*3*3,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_WORLD,mpierr)
 ! call MPI_BARRIER( MPI_COMM_WORLD, mpierr )

  if (myid.eq.master) then
     s_respQ  = s_tmpQ
     sb_respQ = sb_tmpQ
     a_respQ  = a_tmpQ
     ab_respQ = ab_tmpQ
  endif
  endif ! nproc>1



  ! Multiply with global constants here !!! XXX
  !
  ! TO DO

! J = sigma E - alpha NablaT
!
! Then S = alpha / sigma
!
! and I need to multiply alpha from above with (-kb*beta) wrt sigma




  if (myid.eq.master) then

     if (icubic.eq.1) then
        sb_respQ(:,1,3)=sb_respQ(:,1,2)
        sb_respQ(:,2,3)=sb_respQ(:,1,2)

        ab_respQ(:,1,3)=ab_respQ(:,1,2)
        ab_respQ(:,2,3)=ab_respQ(:,1,2)

        s_respQ(:,2,2)=s_respQ(:,1,1)
        s_respQ(:,3,3)=s_respQ(:,1,1)

        a_respQ(:,2,2)=a_respQ(:,1,1)
        a_respQ(:,3,3)=a_respQ(:,1,1)
    endif




     ! divide by number of k-points and multiply with global constants

     facQ  = 2.q0 * piQ * ( real(echarge,16) / real(vol*hbarevs,16)) * 1.q10
     facBQ = 2.q0 * piQ**2 * ( real(echarge,16) / real(vol*hbarevs,16)) *  (1.q-10 / real(hbareVs,16))


     s_respQ  = s_respQ/real(nk,16)  * facQ ! --> sigma in 1/(Ohm m)     [vk's are in eV*Angstroem]
     sb_respQ = sb_respQ/real(nk,16) * facBQ
     a_respQ  = a_respQ/real(nk,16)  * facQ * ( - betaQ * kbQ)  ! --> S=alpha/sigma in units V/K (below conversion to mV/K for output of S)
     ab_respQ = ab_respQ/real(nk,16) * facBQ * ( - betaQ * kbQ)

     s_resp_totQ  = 0.q0
     sb_resp_totQ = 0.q0
     a_resp_totQ  = 0.q0
     ab_resp_totQ = 0.q0


     ! multiply by missing gamma, beta factors and sum over bands
     do ib=1,nband
        gammaQ=real(gam(it,ib),16)

        s_respQ(ib,:,:)=s_respQ(ib,:,:)*betaQ/gammaQ  ! dont do this in Boltzmann mode... as gamma=zero
        a_respQ(ib,:,:)=a_respQ(ib,:,:)*betaQ/gammaQ
        s_resp_totQ=s_resp_totQ+s_respQ(ib,:,:)
        a_resp_totQ=a_resp_totQ+a_respQ(ib,:,:)


        sb_respQ(ib,:,:)=sb_respQ(ib,:,:)*betaQ/(gammaQ**2)
        ab_respQ(ib,:,:)=ab_respQ(ib,:,:)*betaQ/(gammaQ**2)
        sb_resp_totQ=sb_resp_totQ+sb_respQ(ib,:,:)
        ab_resp_totQ=ab_resp_totQ+ab_respQ(ib,:,:)
     enddo


     ! Output raw and combines (Seebeck/Nernst...) data
     !
     ! TO DO

! CHECK!
!     ab_resp_totQ = - betaQ * kbQ * ab_resp_totQ
!     ab_respQ     = - betaQ * kbQ * ab_respQ


     write(130,'(100E20.12)') T,(s_resp_totQ(ia,ia),ia=1,nalpha)
     write(131,'(100E20.12)') T,(s_respQ(:,ia,ia),ia=1,nalpha)

     write(140,'(100E20.12)') T,(a_resp_totQ(ia,ia),ia=1,nalpha)
     write(141,'(100E20.12)') T,(a_respQ(:,ia,ia),ia=1,nalpha)


! ACHTUNG XXX doesn't write out all possible offdiagonals...
     if (ivkab.eq.1) then
     write(150,'(100E20.12)') T,sb_resp_totQ(1,2)
     write(151,'(100E20.12)') T,sb_respQ(:,1,2)

     write(160,'(100E20.12)') T,ab_resp_totQ(1,2)
     write(161,'(100E20.12)') T,ab_respQ(:,1,2)

     endif

! In Seebeck: *1000. so as to yield [S]=mV/K

     do ialpha=1,3
        if (icubic.eq.0) then
           SeebeckQ(ialpha)=1000.q0*a_resp_totQ(ialpha,ialpha)/s_resp_totQ(ialpha,ialpha) ! CHECK CONSTANTS .... kbe ??? SIGN ? XXX
        else
           SeebeckQ(ialpha)=1000.q0*a_resp_totQ(1,1)/s_resp_totQ(1,1) ! CHECK CONSTANTS .... kbe ??? SIGN ? XXX
!           s_resp_totQ(ialpha,ialpha)=s_resp_totQ(1,1)
!           a_resp_totQ(ialpha,ialpha)=a_resp_totQ(1,1)
        endif
     enddo


!     1 = xy
!     2 = xz
!     3 = yz

     NernstQ(1) = ( ab_resp_totQ(1,2) * s_resp_totQ(1,1) - a_resp_totQ(1,1) * sb_resp_totQ(1,2) ) / ( s_resp_totQ(1,1)**2  )
     NernstQ(2) = ( ab_resp_totQ(1,3) * s_resp_totQ(1,1) - a_resp_totQ(1,1) * sb_resp_totQ(1,3) ) / ( s_resp_totQ(1,1)**2  )
     NernstQ(3) = ( ab_resp_totQ(2,3) * s_resp_totQ(2,2) - a_resp_totQ(2,2) * sb_resp_totQ(2,3) ) / ( s_resp_totQ(2,2)**2  )
     NernstQ = NernstQ * 1000.q0 ! V/K --> mV/K


     NernstpartQ(1) = (  ab_resp_totQ(1,2) * s_resp_totQ(1,1) ) / ( s_resp_totQ(1,1)**2  ) * 1000.q0
     NernstpartQ(2) = (- a_resp_totQ(1,1) * sb_resp_totQ(1,2) ) / ( s_resp_totQ(1,1)**2  ) * 1000.q0




     RHQ(1) = - sb_resp_totQ(1,2) / ( s_resp_totQ(1,1)*s_resp_totQ(2,2)    )
     RHQ(2) = - sb_resp_totQ(1,3) / ( s_resp_totQ(1,1)*s_resp_totQ(3,3)    )
     RHQ(3) = - sb_resp_totQ(2,3) / ( s_resp_totQ(3,3)*s_resp_totQ(3,3)    )

     RHQ = RHQ * 1.q+7 ! --> 10^-7 m^3/C


     write(170,'(100E20.12)') T,real(SeebeckQ(:),8) ! in mV/K

! ACHTUNG XXX the following assumes the conductivity tensor to be diagonal!
     write(175,'(100E20.12)') T,(real(1.q0/s_resp_totQ(ia,ia),16),ia=1,nalpha) ! resistivity in Ohm m

     if (ivkab.eq.1) then
     write(171,'(100E20.12)') T,real(NernstQ(:),8)  ! in mV/K
     write(172,'(100E20.12)') T,real(RHQ(:),8) ! in 10^-7 m^3/C

     write(173,'(100E20.12)') T, real( sb_resp_totQ(1,2) / s_resp_totQ(1,1) , 8) ! in 1/T
     write(174,'(100E20.12)') T, real( ab_resp_totQ(1,2) / a_resp_totQ(1,1) , 8) ! in 1/T


     write(180,'(100E20.12)') T,real(NernstpartQ(1),8)  ! in mV/K
     write(181,'(100E20.12)') T,real(NernstpartQ(2),8)  ! in mV/K

     endif

  endif ! master



end subroutine calc_responseQ

end module responseQ
