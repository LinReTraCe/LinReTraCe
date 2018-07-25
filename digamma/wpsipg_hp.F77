*
* $Id: wpsipg.F,v 1.2 2006/09/15 09:34:53 mclareni Exp $
*
* $Log: wpsipg.F,v $
* Revision 1.2  2006/09/15 09:34:53  mclareni
* Submitted mods for gcc4/gfortran and MacOSX, corrected to work also on slc4 with gcc3.4 and 4.1
*
* Revision 1.1.1.1  1996/04/01 15:02:01  mclareni
* Mathlib gen
*
*

c Jan Tomczak 2015. I have converted this to quad-precision and added higher orders
c                   in the Bernoulli coefficients. You need to chose the order kmax
c                   in the pre-compiling options.

#ifdef GFORTRAN
       module gfortran_fix
       integer dp,qp
       parameter(dp=16)
       parameter(qp=selected_real_kind (32))
       end module gfortran_fix
#endif

! kind=selected_real_kind (32)

       FUNCTION WPSIPGHP(Z,K)                                                      
#ifndef GFORTRAN
	use, intrinsic :: iso_fortran_env, dp =>real64, qp=>real128
#else
	use gfortran_fix
#endif

      IMPLICIT real(qp) (A-H,O-Z)                                       

c setting the expansion order of the digamma functions: 6 <= kmax <= 16
c
c if kmax is not defined via the precompiling command, e.g. "-D kmax=10"
#ifndef kmax
#define kmax 6
#else
c if expansion order was set to higher than is available -> set to maximum
#if kmax > 16
#undef kmax
#define kmax 16
#endif
c if expansion order was set to lower than 6, set it to 6 (lowest precision)
#if kmax < 6
#undef kmax
#define kmax 6
#endif
#endif

      complex(qp) WPSIPGHP,Z,U,V,H,R,P,GCMPLX                                      
      CHARACTER NAME*(*)                                                        
      CHARACTER*80 ERRTXT                                                       
      PARAMETER (NAME = 'CPSIPG/WPSIPG')                                        
      DIMENSION C(kmax,0:4),FCT(-1:4),SGN(0:4)                                     
                                                                                
!      PARAMETER (DELTA = 1Q-14)
      PARAMETER (DELTA = 1Q-20)
      PARAMETER (R1 = 1.q0, HF = R1/2.q0)                                             
      PARAMETER (PI=3.14159 26535 89793 23846 26433 83279 50288 41972Q0)                                 
      PARAMETER (C1=PI**2, C2=2.q0*PI**3, C3=2.q0*PI**4, C4=8.q0*PI**5)          
                                                                                
      DATA SGN /-1,+1,-1,+1,-1/
      DATA FCT /0,1,1,2,6,24/                             
                                                                                
      DATA C(1,0) / 8.33333 33333 33333 33333 33333 33333Q-2/                                  
      DATA C(2,0) /-8.33333 33333 33333 33333 33333 33333Q-3/                                  
      DATA C(3,0) / 3.96825 39682 53968 25396 82539 68254Q-3/                                  
      DATA C(4,0) /-4.16666 66666 66666 66666 66666 66667Q-3/                                  
      DATA C(5,0) / 7.57575 75757 57575 75757 57575 75758Q-3/                                  
      DATA C(6,0) /-2.10927 96092 79609 27960 92796 09280Q-2/                                  
#if kmax > 6
      DATA C(7,0) / 8.33333 33333 33333 33333 33333 33333Q-2/                                  
#if kmax > 7
      DATA C(8,0) /-4.43259 80392 15686 27450 98039 21569Q-1/
#if kmax > 8
      DATA C(9,0) / 3.05395 43302 70119 74380 39543 30270Q+0/
#if kmax > 9
      DATA C(10,0) /-2.64562 12121 21212 12121 21212 12121Q+1/
#if kmax > 10
      DATA C(11,0) / 2.81460 14492 75362 31884 05797 10145Q+2/
#if kmax > 11
      DATA C(12,0) /-3.60751 05463 98046 39804 63980 46398Q+3/
#if kmax > 12
      DATA C(13,0) / 5.48275 83333 33333 33333 33333 33333Q+3/
#if kmax > 13
      DATA C(14,0) /-9.74936 82385 05747 12643 67816 09195Q+5/
#if kmax > 14
      DATA C(15,0) / 2.00526 95796 68807 89461 43462 27249Q+7/
#if kmax > 15
      DATA C(16,0) /-4.72384 86772 16299 01960 78431 37255Q+8/
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
                                                                                
      DATA C(1,1) / 1.66666 66666 66666 66666 66666 66667Q-1/                                  
      DATA C(2,1) /-3.33333 33333 33333 33333 33333 33333Q-2/                                  
      DATA C(3,1) / 2.38095 23809 52380 95238 09523 80952Q-2/ 
      DATA C(4,1) /-3.33333 33333 33333 33333 33333 33333Q-2/                                  
      DATA C(5,1) / 7.57575 75757 57575 75757 57575 75758Q-2/                                  
      DATA C(6,1) /-2.53113 55311 35531 13553 11355 31136Q-1/
#if kmax > 6
      DATA C(7,1) / 1.16666 66666 66666 66666 66666 66667Q+0/
#if kmax > 7
      DATA C(8,1) /-7.09215 68627 45098 03921 56862 74510Q+0/
#if kmax > 8
      DATA C(9,1) / 5.49711 77944 86215 53884 71177 94486Q+1/
#if kmax > 9
      DATA C(10,1) /-5.29124 24242 42424 24242 42424 24242Q+2/
#if kmax > 10
      DATA C(11,1) / 6.19212 31884 05797 10144 92753 62319Q+3/
#if kmax > 11
      DATA C(12,1) /-8.65802 53113 55311 35531 13553 11355Q+4/
#if kmax > 12
      DATA C(13,1) / 1.42551 71666 66666 66666 66666 66667Q+6/
#if kmax > 13
      DATA C(14,1) /-2.72982 31067 81609 19540 22988 50575Q+7/
#if kmax > 14
      DATA C(15,1) / 6.01580 87390 06423 68384 30386 81748Q+8/
#if kmax > 15
      DATA C(16,1) /-1.51163 15767 09215 68627 45090 3922Q+10/
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
                                                                                
      DATA C(1,2) / 5.00000 00000 00000 00000 00000 00000Q-1/                                 
      DATA C(2,2) /-1.66666 66666 66666 66666 66666 66667Q-1/                                 
      DATA C(3,2) / 1.66666 66666 66666 66666 66666 66667Q-1/                                 
      DATA C(4,2) /-3.00000 00000 00000 00000 00000 00000Q-1/                                 
      DATA C(5,2) / 8.33333 33333 33333 33333 33333 33333Q-1/                                 
      DATA C(6,2) /-3.29047 61904 76190 47619 04761 90476Q+0/
#if kmax > 6
      DATA C(7,2) / 1.75000 00000 00000 00000 00000 00000Q+1/
#if kmax > 7
      DATA C(8,2) /-1.20566 66666 66666 66666 66666 66667Q+2/
#if kmax > 8
      DATA C(9,2) / 1.04445 23809 52380 95238 09523 80952Q+3/
#if kmax > 9
      DATA C(10,2) /-1.11116 09090 90909 09090 90909 09091Q+4/
#if kmax > 10
      DATA C(11,2) / 1.42418 83333 33333 33333 33333 33333Q+5/
#if kmax > 11
      DATA C(12,2) /-2.16450 63278 38827 83882 78388 27839Q+6/
#if kmax > 12
      DATA C(13,2) / 3.84889 63500 00000 00000 00000 00000Q+7/
#if kmax > 13
      DATA C(14,2) /-7.91648 70096 66666 66666 66666 66667Q+8/
#if kmax > 14
      DATA C(15,2) / 1.86490 07090 91991 34199 13419 91342Q+10/
#if kmax > 15
      DATA C(16,2) /-4.98838 42031 40411 76470 58823 52941Q+11/
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
               
                                                
      DATA C(1,3) / 2.00000 00000 00000 00000 00000 00000Q+0/                                  
      DATA C(2,3) /-1.00000 00000 00000 00000 00000 00000Q+0/                                  
      DATA C(3,3) / 1.33333 33333 33333 33333 33333 33333Q+0/                                  
      DATA C(4,3) /-3.00000 00000 00000 00000 00000 00000Q+0/                                  
      DATA C(5,3) / 1.00000 00000 00000 00000 00000 00000Q+1/                                  
      DATA C(6,3) /-4.60666 66666 66666 66666 66666 66667Q+1/                                  
#if kmax > 6
      DATA C(7,3) / 2.80000 00000 00000 00000 00000 00000Q+2/
#if kmax > 7
      DATA C(8,3) /-2.17020 00000 00000 00000 00000 00000Q+3/ 
#if kmax > 8
      DATA C(9,3) / 2.08890 47619 04761 90476 19047 61905Q+4/
#if kmax > 9
      DATA C(10,3) /-2.44455 40000 00000 00000 00000 00000Q+5/
#if kmax > 10
      DATA C(11,3) / 3.41805 20000 00000 00000 00000 00000Q+6/
#if kmax > 11
      DATA C(12,3) /-5.62771 64523 80952 38095 23809 52381Q+7/
#if kmax > 12
      DATA C(13,3) / 1.07769 09780 00000 00000 00000 00000Q+9/
#if kmax > 13
      DATA C(14,3) /-2.37494 61029 00000 00000 00000 00000Q+10/
#if kmax > 14
      DATA C(15,3) / 5.96768 22690 94372 29437 22943 72294Q+11/
#if kmax > 15
      DATA C(16,3) /-1.69605 06290 67740 00000 00000 00000Q+13/
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif

      DATA C(1,4) / 1.00000 00000 00000 00000 00000 00000Q+1/                                  
      DATA C(2,4) /-7.00000 00000 00000 00000 00000 00000Q+0/                                  
      DATA C(3,4) / 1.20000 00000 00000 00000 00000 00000Q+1/                                  
      DATA C(4,4) /-3.30000 00000 00000 00000 00000 00000Q+1/                                  
      DATA C(5,4) / 1.30000 00000 00000 00000 00000 00000Q+2/                                  
      DATA C(6,4) /-6.91000 00000 00000 00000 00000 00000Q+2/                                  
#if kmax > 6
      DATA C(7,4) / 4.76000 00000 00000 00000 00000 00000Q+3/
#if kmax > 7
      DATA C(8,4) /-4.12338 00000 00000 00000 00000 00000Q+4/ 
#if kmax > 8
      DATA C(9,4) / 4.38670 00000 00000 00000 00000 00000Q+5/
#if kmax > 9
      DATA C(10,4) /-5.62247 42000 00000 00000 00000 00000Q+6/
#if kmax > 10
      DATA C(11,4) / 8.54513 00000 00000 00000 00000 00000Q+7/
#if kmax > 11
      DATA C(12,4) /-1.51948 34421 42857 14285 71428 57143Q+9/
#if kmax > 12
      DATA C(13,4) / 3.12530 38362 00000 00000 00000 00000Q+10/
#if kmax > 13
      DATA C(14,4) /-7.36233 29189 90000 00000 00000 00000Q+11/
#if kmax > 14
      DATA C(15,4) / 1.96933 51488 01142 85714 28571 42857Q+13/
#if kmax > 15
      DATA C(16,4) /-5.93617 72017 37090 00000 00000 00000Q+14/
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
                                                                                
      GCMPLX(X,Y)=CMPLX(X,Y,qp)                                                   
                                                                                
      U=Z 
      X=real(U,qp)                                                                       
      A=ABS(X)                                                                  
      IF(K .LT. 0 .OR. K .GT. 4) THEN                                           
       H=0                                                                      
       WRITE(ERRTXT,101) K                                                      
       CALL MTLPRT(NAME,'C317.1',ERRTXT)                                        
      ELSEIF(ABS(AIMAG(U)) .LT. DELTA .AND. ABS(X+NINT(A)) .LT. DELTA)           
     1                                                        THEN              
       H=0                                                                      
       WRITE(ERRTXT,102) X                                                      
       CALL MTLPRT(NAME,'C317.2',ERRTXT)                                        
      ELSE                                                                      
       K1=K+1                                                                   
       IF(X .LT. 0.q0) U=-U                                                        
       V=U                                                                      
       H=0                                                                      
       IF(A .LT. 15) THEN                                                       
        H=1.q0/V**K1                                                               
        DO 1 I = 1,14-INT(A)                                                    
        V=V+1.q0                                                                   
    1   H=H+1.q0/V**K1                                                             
        V=V+1.q0                                                                   
       END IF                                                                   
       R=1.q0/V**2                                                                 
       P=R*C(kmax,K)                                                               
       DO 2 I = kmax-1,1,-1                                                          
    2  P=R*(C(I,K)+P)                                                           
       H=SGN(K)*(FCT(K)*H+(V*(FCT(K-1)+P)+HF*FCT(K))/V**K1)                     
       IF(K .EQ. 0) H=H+LOG(V)                                                  
       IF(X .LT. 0) THEN                                                        
        V=PI*U                                                                  
        X=real(V,qp)                                                                     
        Y=AIMAG(V)                                                               
        A=SIN(X)                                                                
        B=COS(X)                                                                
        T=TANH(Y)                                                               
        P=GCMPLX(B,-A*T)/GCMPLX(A,B*T)                                          
        IF(K .EQ. 0) THEN                                                       
         H=H+1.q0/U+PI*P                                                           
        ELSEIF(K .EQ. 1) THEN                                                   
         H=-H+1.q0/U**2+C1*(P**2+1.q0)                                                
        ELSEIF(K .EQ. 2) THEN                                                   
         H=H+2.q0/U**3+C2*P*(P**2+1.q0)                                               
        ELSEIF(K .EQ. 3) THEN                                                   
         R=P**2                                                                 
         H=-H+6.q0/U**4+C3*((3.q0*R+4.q0)*R+1.q0)                                           
        ELSEIF(K .EQ. 4) THEN                                                   
         R=P**2                                                                 
         H=H+24.q0/U**5+C4*P*((3.q0*R+5.q0)*R+2.q0)                                         
        ENDIF                                                                   
       ENDIF                                                                    
      ENDIF                                                                     
      WPSIPGHP=H                                                                  
      RETURN                                                                    
  101 FORMAT('K = ',I5,'  (< 0  OR  > 4)')                                      
  102 FORMAT('ARGUMENT EQUALS NON-POSITIVE INTEGER = ',F8.1)                    
      END                                                                       
