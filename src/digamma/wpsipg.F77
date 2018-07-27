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

c Jan Tomczak 2015. I have added higher orders
c                   in the Bernoulli coefficients. You need to chose the order kmax
c                   in the pre-compiling options.


      FUNCTION WPSIPG(Z,K)                                                      
      IMPLICIT real (kind=selected_real_kind (8)) (A-H,O-Z)                                       
      integer, parameter :: PRES = selected_real_kind(8)

c setting the expansion order of the digamma functions: 6 <= kmax <= 12
c
c if kmax is not defined via the precompiling command, e.g. "-D kmax=10"
#ifndef kmax
#define kmax 6
#else
c if expansion order was set to higher than is available -> set to maximum
#if kmax > 12
#undef kmax
#define kmax 12
#endif
c if expansion order was set to lower than 6, set it to 6 (lowest precision)
#if kmax < 6
#undef kmax
#define kmax 6
#endif
#endif

      complex (kind=PRES) WPSIPG,Z,U,V,H,R,P,GCMPLX                                      
      CHARACTER NAME*(*)                                                        
      CHARACTER*80 ERRTXT                                                       
      PARAMETER (NAME = 'CPSIPG/WPSIPG')                                        
      DIMENSION C(kmax,0:4),FCT(-1:4),SGN(0:4)                                     
                                                                                
      PARAMETER (DELTA = 1D-14)
      PARAMETER (R1 = 1.d0, HF = R1/2.d0)                                             
      PARAMETER (PI=3.14159 26535 89793 23846 26433 83279 50288 41972D0)                                 
      PARAMETER (C1=PI**2, C2=2.d0*PI**3, C3=2.d0*PI**4, C4=8.d0*PI**5)          
                                                                                
      DATA SGN /-1,+1,-1,+1,-1/
      DATA FCT /0,1,1,2,6,24/                             
                                                                                
      DATA C(1,0) / 8.33333 33333 33333 33333 33333 33333D-2/                                  
      DATA C(2,0) /-8.33333 33333 33333 33333 33333 33333D-3/                                  
      DATA C(3,0) / 3.96825 39682 53968 25396 82539 68254D-3/                                  
      DATA C(4,0) /-4.16666 66666 66666 66666 66666 66667D-3/                                  
      DATA C(5,0) / 7.57575 75757 57575 75757 57575 75758D-3/                                  
      DATA C(6,0) /-2.10927 96092 79609 27960 92796 09280D-2/                                  
#if kmax > 6
      DATA C(7,0) / 8.33333 33333 33333 33333 33333 33333D-2/                                  
#if kmax > 7
      DATA C(8,0) /-4.43259 80392 15686 27450 98039 21569D-1/
#if kmax > 8
      DATA C(9,0) / 3.05395 43302 70119 74380 39543 30270D+0/
#if kmax > 9
      DATA C(10,0) /-2.64562 12121 21212 12121 21212 12121D+1/
#if kmax > 10
      DATA C(11,0) / 2.81460 14492 75362 31884 05797 10145D+2/
#if kmax > 11
      DATA C(12,0) /-3.60751 05463 98046 39804 63980 46398D+3/
#endif
#endif
#endif
#endif
#endif
#endif
                                                                                
      DATA C(1,1) / 1.66666 66666 66666 66666 66666 66667D-1/                                  
      DATA C(2,1) /-3.33333 33333 33333 33333 33333 33333D-2/                                  
      DATA C(3,1) / 2.38095 23809 52380 95238 09523 80952D-2/ 
      DATA C(4,1) /-3.33333 33333 33333 33333 33333 33333D-2/                                  
      DATA C(5,1) / 7.57575 75757 57575 75757 57575 75758D-2/                                  
      DATA C(6,1) /-2.53113 55311 35531 13553 11355 31136D-1/
#if kmax > 6
      DATA C(7,1) / 1.16666 66666 66666 66666 66666 66667D+0/
#if kmax > 7
      DATA C(8,1) /-7.09215 68627 45098 03921 56862 74510D+0/
#if kmax > 8
      DATA C(9,1) / 5.49711 77944 86215 53884 71177 94486D+1/
#if kmax > 9
      DATA C(10,1) /-5.29124 24242 42424 24242 42424 24242D+2/
#if kmax > 10
      DATA C(11,1) / 6.19212 31884 05797 10144 92753 62319D+3/
#if kmax > 11
      DATA C(12,1) /-8.65802 53113 55311 35531 13553 11355D+4/
#endif
#endif
#endif
#endif
#endif
#endif
                                                                                
      DATA C(1,2) / 5.00000 00000 00000 00000 00000 00000D-1/                                 
      DATA C(2,2) /-1.66666 66666 66666 66666 66666 66667D-1/                                 
      DATA C(3,2) / 1.66666 66666 66666 66666 66666 66667D-1/                                 
      DATA C(4,2) /-3.00000 00000 00000 00000 00000 00000D-1/                                 
      DATA C(5,2) / 8.33333 33333 33333 33333 33333 33333D-1/                                 
      DATA C(6,2) /-3.29047 61904 76190 47619 04761 90476D+0/
#if kmax > 6
      DATA C(7,2) / 1.75000 00000 00000 00000 00000 00000D+1/
#if kmax > 7
      DATA C(8,2) /-1.20566 66666 66666 66666 66666 66667D+2/
#if kmax > 8
      DATA C(9,2) / 1.04445 23809 52380 95238 09523 80952D+3/
#if kmax > 9
      DATA C(10,2) /-1.11116 09090 90909 09090 90909 09091D+4/
#if kmax > 10
      DATA C(11,2) / 1.42418 83333 33333 33333 33333 33333D+5/
#if kmax > 11
      DATA C(12,2) /-2.16450 63278 38827 83882 78388 27839D+6/
#endif
#endif
#endif
#endif
#endif
#endif
               
                                                
      DATA C(1,3) / 2.00000 00000 00000 00000 00000 00000D+0/                                  
      DATA C(2,3) /-1.00000 00000 00000 00000 00000 00000D+0/                                  
      DATA C(3,3) / 1.33333 33333 33333 33333 33333 33333D+0/                                  
      DATA C(4,3) /-3.00000 00000 00000 00000 00000 00000D+0/                                  
      DATA C(5,3) / 1.00000 00000 00000 00000 00000 00000D+1/                                  
      DATA C(6,3) /-4.60666 66666 66666 66666 66666 66667D+1/                                  
#if kmax > 6
      DATA C(7,3) / 2.80000 00000 00000 00000 00000 00000D+2/
#if kmax > 7
      DATA C(8,3) /-2.17020 00000 00000 00000 00000 00000D+3/ 
#if kmax > 8
      DATA C(9,3) / 2.08890 47619 04761 90476 19047 61905D+4/
#if kmax > 9
      DATA C(10,3) /-2.44455 40000 00000 00000 00000 00000D+5/
#if kmax > 10
      DATA C(11,3) / 3.41805 20000 00000 00000 00000 00000D+6/
#if kmax > 11
      DATA C(12,3) /-5.62771 64523 80952 38095 23809 52381D+7/
#endif
#endif
#endif
#endif
#endif
#endif

      DATA C(1,4) / 1.00000 00000 00000 00000 00000 00000D+1/                                  
      DATA C(2,4) /-7.00000 00000 00000 00000 00000 00000D+0/                                  
      DATA C(3,4) / 1.20000 00000 00000 00000 00000 00000D+1/                                  
      DATA C(4,4) /-3.30000 00000 00000 00000 00000 00000D+1/                                  
      DATA C(5,4) / 1.30000 00000 00000 00000 00000 00000D+2/                                  
      DATA C(6,4) /-6.91000 00000 00000 00000 00000 00000D+2/                                  
#if kmax > 6
      DATA C(7,4) / 4.76000 00000 00000 00000 00000 00000D+3/
#if kmax > 7
      DATA C(8,4) /-4.12338 00000 00000 00000 00000 00000D+4/ 
#if kmax > 8
      DATA C(9,4) / 4.38670 00000 00000 00000 00000 00000D+5/
#if kmax > 9
      DATA C(10,4) /-5.62247 42000 00000 00000 00000 00000D+6/
#if kmax > 10
      DATA C(11,4) / 8.54513 00000 00000 00000 00000 00000D+7/
#if kmax > 11
      DATA C(12,4) /-1.51948 34421 42857 14285 71428 57143D+9/
#endif
#endif
#endif
#endif
#endif
#endif
                                                                                
      GCMPLX(X,Y)=CMPLX(X,Y,PRES)                                                   
                                                                                
      U=Z                                                                       
      X=real(U,PRES)                                                                       
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
       IF(X .LT. 0.d0) U=-U                                                        
       V=U                                                                      
       H=0                                                                      
       IF(A .LT. 15) THEN                                                       
        H=1.d0/V**K1                                                               
        DO 1 I = 1,14-INT(A)                                                    
        V=V+1.d0                                                                   
    1   H=H+1.d0/V**K1                                                             
        V=V+1.d0                                                                   
       END IF                                                                   
       R=1.d0/V**2                                                                 
       P=R*C(kmax,K)                                                               
       DO 2 I = kmax-1,1,-1                                                          
    2  P=R*(C(I,K)+P)                                                           
       H=SGN(K)*(FCT(K)*H+(V*(FCT(K-1)+P)+HF*FCT(K))/V**K1)                     
       IF(K .EQ. 0) H=H+LOG(V)                                                  
       IF(X .LT. 0) THEN                                                        
        V=PI*U                                                                  
        X=real(V,PRES)                                                                     
        Y=AIMAG(V)                                                               
        A=SIN(X)                                                                
        B=COS(X)                                                                
        T=TANH(Y)                                                               
        P=GCMPLX(B,-A*T)/GCMPLX(A,B*T)                                          
        IF(K .EQ. 0) THEN                                                       
         H=H+1.d0/U+PI*P                                                           
        ELSEIF(K .EQ. 1) THEN                                                   
         H=-H+1.d0/U**2+C1*(P**2+1.d0)                                                
        ELSEIF(K .EQ. 2) THEN                                                   
         H=H+2.d0/U**3+C2*P*(P**2+1.d0)                                               
        ELSEIF(K .EQ. 3) THEN                                                   
         R=P**2                                                                 
         H=-H+6.d0/U**4+C3*((3.d0*R+4.d0)*R+1.d0)                                           
        ELSEIF(K .EQ. 4) THEN                                                   
         R=P**2                                                                 
         H=H+24.d0/U**5+C4*P*((3.d0*R+5.d0)*R+2.d0)                                         
        ENDIF                                                                   
       ENDIF                                                                    
      ENDIF                                                                     
      WPSIPG=H                                                                  
      RETURN                                                                    
  101 FORMAT('K = ',I5,'  (< 0  OR  > 4)')                                      
  102 FORMAT('ARGUMENT EQUALS NON-POSITIVE INTEGER = ',F8.1)                    
      END                                                                       
