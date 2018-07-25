



module response



! eM: TODO: define all the variables inside newly defined types


use params
use mpi_org
use estruct !, only: delta,gam,gmax,iband_valence,gminall,emin,emax,ek,vka,vkab,zqp 
implicit none

integer ik,ipg,ialpha,ibeta,idiag,icubic,nalpha,ia

! DOUBLE PRECISION
real(8) :: RePolyGamma(0:3),ImPolyGamma(0:3),gamma,aqp,z,tmp
complex(8) :: ctmp,zarg,wpsipg
external wpsipg

!kernels
real(8) ::  s_kernel ! for conductivity
real(8) ::  sB_kernel ! nband,3,3 for conductivity in B-field
real(8) ::  a_kernel ! nband,3,3 for Peltier
real(8) ::  aB_kernel ! nband,3,3 for Peltier in B-field          
!response functions...
real(8), allocatable :: s_resp(:,:,:),s_tmp(:,:,:) ! nband,3,3 for conductivity
real(8), allocatable :: sB_resp(:,:,:),sb_tmp(:,:,:) ! nband,3,3 for conductivity in B-field
real(8), allocatable :: a_resp(:,:,:),a_tmp(:,:,:) ! nband,3,3 for Peltier
real(8), allocatable :: aB_resp(:,:,:),ab_tmp(:,:,:) ! nband,3,3 for Peltier in B-field           

real(8) :: s_resp_tot(3,3)
real(8) :: sb_resp_tot(3,3)
real(8) :: a_resp_tot(3,3)
real(8) :: ab_resp_tot(3,3)
real(8) :: Seebeck(3),Nernst(3),RH(3)

! QUAD PRECISION
real(16) :: RePolyGammaQ(0:3),ImPolyGammaQ(0:3),gammaQ,aqpQ,zQ,tmpQ
complex(16) :: ctmpQ,zargQ,wpsipghp
external wpsipghp

!kernels
real(16) ::  s_kernelQ ! for conductivity
real(16) ::  sB_kernelQ ! nband,3,3 for conductivity in B-field
real(16) ::  a_kernelQ ! nband,3,3 for Peltier
real(16) ::  aB_kernelQ ! nband,3,3 for Peltier in B-field          
!response functions...
real(16), allocatable :: s_respQ(:,:,:),s_tmpQ(:,:,:) ! nband,3,3 for conductivity
real(16), allocatable :: sB_respQ(:,:,:),sb_tmpQ(:,:,:) ! nband,3,3 for conductivity in B-field
real(16), allocatable :: a_respQ(:,:,:),a_tmpQ(:,:,:) ! nband,3,3 for Peltier
real(16), allocatable :: aB_respQ(:,:,:),ab_tmpQ(:,:,:) ! nband,3,3 for Peltier in B-field           

real(16) :: s_resp_totQ(3,3)
real(16) :: sb_resp_totQ(3,3)
real(16) :: a_resp_totQ(3,3)
real(16) :: ab_resp_totQ(3,3)
real(16) :: SeebeckQ(3),NernstQ(3),RHQ(3),NernstpartQ(2) ! nernstpart is axy sxx bzw axx sxy / sxx^2 


end module response
!contains

! eM: TODO: I guess that a good way to proceed here is to break up the calc_response routine into two 
! subroutines that compute the interband and the intraband contributions separately 


subroutine calc_response(mu,it)
  use response
  implicit none
  real(8) :: mu,fac,facB
  integer it,ib
    
  s_resp  = 0.d0
  sb_resp = 0.d0
  a_resp  = 0.d0
  ab_resp = 0.d0
  
  !put somewhere else...
  nalpha=3
  if (icubic.eq.1) nalpha=1
    
  ! outer k-loop   
  do ik=ikstart,ikend     
  ! eM: where have the values of ikstart and ikend been defined??
     do ib=1,nband ! at this point everything is band-diagonal 
     !(eM: this will not be true the moment one uses the optical matrix elements instead of the fermi velocities)
     ! and because of the tetrahedron method one has to trace over bands first and then k-points
        gamma=gam(it,ib)

        ! pre-compute all needed digamma functions.   
        aqp=ek(ik,ib)-mu
        z=zqp(ib)
        zarg=0.5d0+beta2p*(ci*aqp+gamma)
        do ipg=1,3 ! XXX need 0 for alphaxy ????
           ctmp=wpsipg(zarg,ipg)
           RePolyGamma(ipg)=real(ctmp,8)
           ImPolyGamma(ipg)=imag(ctmp)
        enddo
        
        ! compute transport kernels (omega-part)
        ! 
        
        tmp=z**2 / (4.d0*pi**3) ! missing: beta/gamma (multiplied later to keep number reasonable here)        
        s_kernel = tmp * (     RePolyGamma(1) - Gamma*beta2p * RePolyGamma(2) )
        a_kernel = tmp * ( aqp * RePolyGamma(1) - aqp*Gamma*beta2p*RePolyGamma(2) - Gamma**2.d0 * beta2p * ImPolyGamma(2) )
        
        
        tmp=tmp*3.d0*z/(4.d0*pi) ! additionally missing: 1/gamma (multiplied later) XXX
        
        sb_kernel = tmp * ( - RePolyGamma(1) - gamma * beta2p * RePolyGamma(2) - (beta2p*gamma)**2.d0/3.d0 * RePolyGamma(3)    )
        
        
        ab_kernel = tmp * ( aqp * RePolyGamma(1) - aqp * beta2p * gamma * RePolyGamma(2) + gamma**2.d0 / 3.d0  &
             * beta2p * ImPolyGamma(2) &
             - aqp*gamma**2.d0 / 3.d0 * beta2p**2.d0 * ImPolyGamma(3) &
             + gamma**3.d0 / 3.d0 * beta2p**2.d0 * RePolyGamma(3) )


        ! multiply with matrix elements for all polarizations and add to final band-dependent response functions
        ! create an option for computing only polarization-diagonal elements in absence of B 
        !                                                                                

        
        ! B = 0  
        do ialpha=1,nalpha
           do ibeta=ialpha,nalpha
              if ((icubic.eq.1).and.(ibeta.gt.ialpha)) exit
              if ((idiag.eq.1).AND.(ialpha.ne.ibeta)) exit

              tmp=vka(ik,ib,ialpha)*vka(ik,ib,ibeta)

              s_tmp(ib,ialpha,ibeta)=s_kernel * tmp
              a_tmp(ib,ialpha,ibeta)=a_kernel * tmp 
              
           enddo !ibeta  
        enddo ! ialpha           
        
        ! B .ne. 0 
        if (ivkab.eq.1) then ! file vkab exists
        do ialpha=1,nalpha
           do ibeta=ialpha+1,3
              if ((icubic.eq.1).and.(ibeta.gt.2)) exit
              
              tmp=vka(ik,ib,ialpha)*( vkab(ik,ib,ibeta,ialpha)*vka(ik,ib,ibeta) - vkab(ik,ib,ibeta,ibeta)*vka(ik,ib,ialpha)    )

              sb_tmp(ib,ialpha,ibeta)=sb_kernel * tmp
              ab_tmp(ib,ialpha,ibeta)=ab_kernel * tmp 
              
           enddo !ibeta  
        enddo ! ialpha
        else
           sb_tmp=0.d0
           ab_tmp=0.d0
        endif !ivkab

     enddo ! ib   

     ! add to k--summed response functions (keep band-dependence here)   
     ! 
     
     s_resp=s_resp+s_tmp
     sb_resp=sb_resp+sb_tmp

     a_resp=a_resp+a_tmp
     ab_resp=ab_resp+ab_tmp
     
  enddo ! outer k-loop  
  
  
  ! reduce MPI_sum all contributions & output & sum over bands   
  ! 


  if (nproc.gt.1) then
  s_tmp=0.d0
  sb_tmp=0.d0
  a_tmp=0.d0
  ab_tmp=0.d0
  
  call MPI_REDUCE(s_resp,s_tmp,nband*3*3,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_WORLD,mpierr)
!  call MPI_BARRIER( MPI_COMM_WORLD, mpierr )
  call MPI_REDUCE(sb_resp,sb_tmp,nband*3*3,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_WORLD,mpierr)
!  call MPI_BARRIER( MPI_COMM_WORLD, mpierr )
  call MPI_REDUCE(a_resp,a_tmp,nband*3*3,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_WORLD,mpierr)
!  call MPI_BARRIER( MPI_COMM_WORLD, mpierr )
  call MPI_REDUCE(ab_resp,ab_tmp,nband*3*3,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_WORLD,mpierr)
!  call MPI_BARRIER( MPI_COMM_WORLD, mpierr )
  
  if (myid.eq.master) then
     s_resp=s_tmp
     sb_resp=sb_tmp
     a_resp=a_tmp
     ab_resp=ab_tmp
  endif
  endif

  
  
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
        sb_resp(:,1,3)=sb_resp(:,1,2)
        sb_resp(:,2,3)=sb_resp(:,1,2)

        ab_resp(:,1,3)=ab_resp(:,1,2)
        ab_resp(:,2,3)=ab_resp(:,1,2)

        s_resp(:,2,2)=s_resp(:,1,1)
        s_resp(:,3,3)=s_resp(:,1,1)

        a_resp(:,2,2)=a_resp(:,1,1)
        a_resp(:,3,3)=a_resp(:,1,1)

     endif



 
     fac  = 2.d0 * pi * ( echarge / (vol*hbarevs)) * 1.d10 
     facB = 2.d0 * pi**2 * ( echarge / (vol*hbarevs) ) * (1.d-10 / hbareVs)


     s_resp=s_resp/real(nk,8)    * fac   ! --> sigma in 1/(Ohm m)     [vk's are in eV*Angstroem]
     sb_resp=sb_resp/real(nk,8)  * facB
     a_resp=a_resp/real(nk,8)    * fac * ( - beta * kb)  ! --> S=alpha/sigma in units V/K (below conversion to mV/K for output of S)
     ab_resp=ab_resp/real(nk,8)  * facB  * ( - beta * kb)
     
     
     s_resp_tot=0.d0
     sb_resp_tot=0.d0
     a_resp_tot=0.d0
     ab_resp_tot=0.d0
     
     do ib=1,nband
        gamma=gam(it,ib)

        s_resp(ib,:,:)=s_resp(ib,:,:)*beta/gamma  ! dont do this in Boltzmann mode... as gamma=zero 
        a_resp(ib,:,:)=a_resp(ib,:,:)*beta/gamma
        s_resp_tot=s_resp_tot+s_resp(ib,:,:)
        a_resp_tot=a_resp_tot+a_resp(ib,:,:)


        sb_resp(ib,:,:)=sb_resp(ib,:,:)*beta/gamma**2
        ab_resp(ib,:,:)=ab_resp(ib,:,:)*beta/gamma**2
        sb_resp_tot=sb_resp_tot+sb_resp(ib,:,:)
        ab_resp_tot=ab_resp_tot+ab_resp(ib,:,:)

     enddo
     
     
     ! Output raw and combines (Seebeck/Nernst...) data  
     ! 

! CHECK
!     ab_resp_tot = - beta * kb * ab_resp_tot
!     ab_resp     = - beta * kb * ab_resp


   !eM: all this writing will go in a subroutine with only one type being passed down
     
     write(30,'(100E20.12)') T,(s_resp_tot(ia,ia),ia=1,nalpha)
     write(31,'(100E20.12)') T,(s_resp(:,ia,ia),ia=1,nalpha)
     write(40,'(100E20.12)') T,(a_resp_tot(ia,ia),ia=1,nalpha)
     write(41,'(100E20.12)') T,(a_resp(:,ia,ia),ia=1,nalpha)
     
! XXX ACHTUNG not writing out all combinations...
     if (ivkab.eq.1) then
        write(50,'(100E20.12)') T,sb_resp_tot(1,2)
        write(51,'(100E20.12)') T,sb_resp(:,1,2)

        write(60,'(100E20.12)') T,ab_resp_tot(1,2)
        write(61,'(100E20.12)') T,ab_resp(:,1,2)
     endif

! In Seebeck: *1000. so as to yield [S]=mV/K        
     
     
     do ialpha=1,3
        if (icubic.eq.0) then
           Seebeck(ialpha)=1000.d0*a_resp_tot(ialpha,ialpha)/s_resp_tot(ialpha,ialpha) 
        else
           Seebeck(ialpha)=1000.d0*a_resp_tot(1,1)/s_resp_tot(1,1) 
!           s_resp_tot(ialpha,ialpha)=s_resp_tot(1,1)
!           a_resp_tot(ialpha,ialpha)=a_resp_tot(1,1)
        endif


     enddo

!     1 = xy
!     2 = xz
!     3 = yz

     Nernst(1) = ( ab_resp_tot(1,2) * s_resp_tot(1,1) - a_resp_tot(1,1) * sb_resp_tot(1,2) ) / ( s_resp_tot(1,1)**2.d0  ) 
     Nernst(2) = ( ab_resp_tot(1,3) * s_resp_tot(1,1) - a_resp_tot(1,1) * sb_resp_tot(1,3) ) / ( s_resp_tot(1,1)**2.d0  ) 
     Nernst(3) = ( ab_resp_tot(2,3) * s_resp_tot(2,2) - a_resp_tot(2,2) * sb_resp_tot(2,3) ) / ( s_resp_tot(2,2)**2.d0  ) 

     Nernst = Nernst * 1000.d0 ! V/K --> mV/K


     RH(1) = - sb_resp_tot(1,2) / ( s_resp_tot(1,1)*s_resp_tot(2,2)    )
     RH(2) = - sb_resp_tot(1,3) / ( s_resp_tot(1,1)*s_resp_tot(3,3)    )
     RH(3) = - sb_resp_tot(2,3) / ( s_resp_tot(3,3)*s_resp_tot(3,3)    )
 
     RH = RH * 1.d+7 ! --> 10^-7 m^3/C

    
     write(70,'(100E20.12)') T,Seebeck(:) ! in mV/K

     if (ivkab.eq.1) then
        write(71,'(100E20.12)') T,Nernst(:)  ! in mV/KT
        write(72,'(100E20.12)') T,RH(:) ! in 10^-7 m^3/C
        write(73,'(100E20.12)') T, sb_resp_tot(1,2) / s_resp_tot(1,1) ! in 1/T
        write(74,'(100E20.12)') T, ab_resp_tot(1,2) / a_resp_tot(1,1) ! in 1/T
     endif
! XXX ACHTUNG assumes diagonal conductivity tensor... !
     write(75,'(100E20.12)') T,(1.d0/s_resp_tot(ia,ia),ia=1,nalpha) ! resistivity in Ohm m
     
     
  endif ! master                                                                                                              
  
  
  
  
end subroutine calc_response



subroutine response_open_files()
use estruct, only: ivkab
implicit none

   open(30,file='sigma_tot.dat',status='unknown')
   open(31,file='sigma_band_xx.dat',status='unknown')
   open(40,file='peltier_tot.dat',status='unknown')
   open(41,file='peltier_band_xx.dat',status='unknown')
! 50 for sxy                                                                                
   if (ivkab.eq.1) then
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
   if (ivkab.eq.1) then
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

! QUAD PRECISION, BOLTZMANN MODE
   open(230,file='sigma_tot.datBQ',status='unknown')
   open(231,file='sigma_band_xx.datBQ',status='unknown')
   open(240,file='peltier_tot.datBQ',status='unknown')
   open(241,file='peltier_band_xx.datBQ',status='unknown')
! 50 for sxy                                                                                
   if (ivkab.eq.1) then
   open(250,file='sigmaB_tot.datBQ',status='unknown')
   open(251,file='sigmaB_band_xy.datBQ',status='unknown')
! 60 for axy                                                                                
   open(260,file='peltierB_tot.datBQ',status='unknown')
   open(261,file='peltierB_band_xy.datBQ',status='unknown')
   open(271,file='Nernst.datBQ',status='unknown')
   open(272,file='RH.datBQ',status='unknown')
   open(273,file='muH.datBQ',status='unknown') ! Hall mobility
   open(274,file='mut.datBQ',status='unknown') ! thermal counterpart of Hall mobility
   open(280,file='Nernst_part1.datBQ',status='unknown')
   open(281,file='Nernst_part2.datBQ',status='unknown')
   endif
! 70 for aux's: Seebeck, Hall, Nernst...                                                    
   open(270,file='Seebeck.datBQ',status='unknown')
   open(275,file='resistivity.datBQ',status='unknown')


end subroutine response_open_files

subroutine response_close_files()
use estruct, only : ivkab
implicit none

   close(30)
   close(31)
   close(40)
   close(41)
   if (ivkab.eq.1) then
   close(50)
   close(51)
   close(60)
   close(61)
   close(70)
   close(71)
   close(72)
   close(73)
   close(74)
   endif
   close(75)

   close(130)
   close(131)
   close(140)
   close(141)
   if (ivkab.eq.1) then
   close(150)
   close(151)
   close(160)
   close(161)
   close(170)
   close(171)
   close(172)
   close(173)
   close(174)
   endif
   close(175)

   close(180)
   close(181)

   close(230)
   close(231)
   close(240)
   close(241)
   if (ivkab.eq.1) then
   close(250)
   close(251)
   close(260)
   close(261)
   close(270)
   close(271)
   close(272)
   close(273)
   close(274)
   endif
   close(275)

   close(280)
   close(281)
end subroutine response_close_files


subroutine response_allocate_variables()
use response
implicit none

! allocate transport variables                                                              
allocate(s_tmp(nband,3,3))
allocate(a_tmp(nband,3,3))
allocate(sb_tmp(nband,3,3))
allocate(ab_tmp(nband,3,3))

allocate(s_resp(nband,3,3))
allocate(a_resp(nband,3,3))
allocate(sb_resp(nband,3,3))
allocate(ab_resp(nband,3,3))

allocate(s_tmpQ(nband,3,3))
allocate(a_tmpQ(nband,3,3))
allocate(sb_tmpQ(nband,3,3))
allocate(ab_tmpQ(nband,3,3))

allocate(s_respQ(nband,3,3))
allocate(a_respQ(nband,3,3))
allocate(sb_respQ(nband,3,3))
allocate(ab_respQ(nband,3,3))


end subroutine response_allocate_variables
