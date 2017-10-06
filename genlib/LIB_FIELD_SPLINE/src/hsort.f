

      SUBROUTINE HSORT(K,N)                                          
C  HEAPSORT ALGORITHM FOR SORTING ON VECTOR OF KEYS K OF LENGTH N    
C  J F MONAHAN        TRANSCRIBED FROM KNUTH, VOL 2, PP 146-7.       
      REAL*8 K(1),KK                                                 
      INTEGER R                                                      
      IF(N.LE.1) RETURN                                              
      L=N/2+1                                                        
      R=N                                                            
  2   IF(L.GT.1) GO TO 1                                             
      KK=K(R)                                                        
      K(R)=K(1)                                                      
      R=R-1                                                          
      IF(R.EQ.1) GO TO 9                                             
      GO TO 3                                                        
  1   L=L-1                                                          
      KK=K(L)                                                        
  3   J=L                                                            
  4   I=J                                                            
      J=2*J                                                          
C     IF(J-R) 5,6,8                                                  
      if((J-R).gt.0) go to 8
      if((J-R).eq.0) go to 6
C
  5   IF(K(J).LT.K(J+1)) J=J+1                                       
  6   IF(KK.GT.K(J)) GO TO 8                                         
  7   K(I)=K(J)                                                      
      GO TO 4                                                        
  8   K(I)=KK                                                        
      GO TO 2                                                        
  9   K(1)=KK                                                        
      RETURN                                                         
      END                                                            
