

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
c                   in the Bernoulli coefficients. You need to chose the order 7
c                   in the pre-compiling options.


      FUNCTION WPSIPG(Z,K)                                                      
      IMPLICIT real (kind=selected_real_kind (8)) (A-H,O-Z)                                       
      integer, parameter :: PRES = selected_real_kind(8)

c setting the expansion order of the digamma functions: 6 <= 7 <= 12
c
c if 7 is not defined via the precompiling command, e.g. "-D kmax=10"



c if expansion order was set to higher than is available -> set to maximum




c if expansion order was set to lower than 6, set it to 6 (lowest precision)






      complex (kind=PRES) WPSIPG,Z,U,V,H,R,P,GCMPLX                                      
      CHARACTER NAME*(*)                                                        
      CHARACTER*80 ERRTXT                                                       
      PARAMETER (NAME = 'CPSIPG/WPSIPG')                                        
      DIMENSION C(7,0:4),FCT(-1:4),SGN(0:4)                                     
                                                                                
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

      DATA C(7,0) / 8.33333 33333 33333 33333 33333 33333D-2/                                  
                                                                                
      DATA C(1,1) / 1.66666 66666 66666 66666 66666 66667D-1/                                  
      DATA C(2,1) /-3.33333 33333 33333 33333 33333 33333D-2/                                  
      DATA C(3,1) / 2.38095 23809 52380 95238 09523 80952D-2/ 
      DATA C(4,1) /-3.33333 33333 33333 33333 33333 33333D-2/                                  
      DATA C(5,1) / 7.57575 75757 57575 75757 57575 75758D-2/                                  
      DATA C(6,1) /-2.53113 55311 35531 13553 11355 31136D-1/
      DATA C(7,1) / 1.16666 66666 66666 66666 66666 66667D+0/
                                                                                
      DATA C(1,2) / 5.00000 00000 00000 00000 00000 00000D-1/                                 
      DATA C(2,2) /-1.66666 66666 66666 66666 66666 66667D-1/                                 
      DATA C(3,2) / 1.66666 66666 66666 66666 66666 66667D-1/                                 
      DATA C(4,2) /-3.00000 00000 00000 00000 00000 00000D-1/                                 
      DATA C(5,2) / 8.33333 33333 33333 33333 33333 33333D-1/                                 
      DATA C(6,2) /-3.29047 61904 76190 47619 04761 90476D+0/
      DATA C(7,2) / 1.75000 00000 00000 00000 00000 00000D+1/
               
                                                
      DATA C(1,3) / 2.00000 00000 00000 00000 00000 00000D+0/                                  
      DATA C(2,3) /-1.00000 00000 00000 00000 00000 00000D+0/                                  
      DATA C(3,3) / 1.33333 33333 33333 33333 33333 33333D+0/                                  
      DATA C(4,3) /-3.00000 00000 00000 00000 00000 00000D+0/                                  
      DATA C(5,3) / 1.00000 00000 00000 00000 00000 00000D+1/                                  
      DATA C(6,3) /-4.60666 66666 66666 66666 66666 66667D+1/                                  
      DATA C(7,3) / 2.80000 00000 00000 00000 00000 00000D+2/

      DATA C(1,4) / 1.00000 00000 00000 00000 00000 00000D+1/                                  
      DATA C(2,4) /-7.00000 00000 00000 00000 00000 00000D+0/                                  
      DATA C(3,4) / 1.20000 00000 00000 00000 00000 00000D+1/                                  
      DATA C(4,4) /-3.30000 00000 00000 00000 00000 00000D+1/                                  
      DATA C(5,4) / 1.30000 00000 00000 00000 00000 00000D+2/                                  
      DATA C(6,4) /-6.91000 00000 00000 00000 00000 00000D+2/                                  
      DATA C(7,4) / 4.76000 00000 00000 00000 00000 00000D+3/
                                                                                
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
       P=R*C(7,K)                                                               
       DO 2 I = 7-1,1,-1                                                          
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
