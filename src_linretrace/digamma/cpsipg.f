*
* $Id: cpsipg.F,v 1.1.1.1 1996/04/01 15:02:01 mclareni Exp $
*
* $Log: cpsipg.F,v $
* Revision 1.1.1.1  1996/04/01 15:02:01  mclareni
* Mathlib gen
*
*
      FUNCTION CPSIPG(Z,K)                                                      
      COMPLEX*16 CPSIPG,Z                                                          
      COMPLEX(kind=selected_real_kind (32)) WPSIPG,W                                                       
      real(kind=selected_real_kind (32)) D                                                        
                                                                                
      SROUND(D)=D+(D-SNGL(D))                                                   
      W=Z                                                                       
      W=WPSIPG(W,K)                                                             
      CPSIPG=CMPLX(SROUND(REAL(W)),SROUND(IMAG(W)))                           
      RETURN                                                                    
      END                                                                       
