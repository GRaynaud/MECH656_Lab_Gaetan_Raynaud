!     Last change:  SB   21 Jun 1999   10:36 am
      SUBROUTINE REALFT(dat,N,ISIGN)

      implicit none
      INTEGER, INTENT(IN):: N, isign
      REAL(8), INTENT(IN OUT):: DAT(2*N)
      REAL(8):: WR,WI,WPR,WPI,WTEMP,THETA
      REAL:: H1R, H2R, H1I, H2I, c1, c2, WRS, WIS
      INTEGER:: i,i1,i2,i3,i4,n2p3

      interface

      subroutine four1(DAT,L,i)
        implicit none
        INTEGER, INTENT(IN):: i,L
        REAL(8), INTENT(IN OUT):: DAT(2*L)
      end subroutine four1

      end interface

      THETA=6.28318530717959D0*0.5D0/DBLE(N)
      C1=0.5                                                                    
      IF (ISIGN == 1) THEN
        C2=-0.5                                                                 
        CALL FOUR1(DAT,N,+1)
      ELSE                                                                      
        C2=0.5                                                                  
        THETA=-THETA                                                            
      END IF
      WPR=-2.0D0*SIN(0.5D0*THETA)**2
      WPI=SIN(THETA)
      WR=1.0D0+WPR
      WI=WPI                                                                    
      N2P3=2*N+3                                                                
      DO I=2,N/2
        I1=2*I-1                                                                
        I2=I1+1                                                                 
        I3=N2P3-I2                                                              
        I4=I3+1                                                                 
        WRS=WR
        WIS=WI
        H1R=C1*(DAT(I1)+DAT(I3))
        H1I=C1*(DAT(I2)-DAT(I4))
        H2R=-C2*(DAT(I2)+DAT(I4))
        H2I=C2*(DAT(I1)-DAT(I3))
        DAT(I1)=H1R+WRS*H2R-WIS*H2I
        DAT(I2)=H1I+WRS*H2I+WIS*H2R
        DAT(I3)=H1R-WRS*H2R+WIS*H2I
        DAT(I4)=-H1I+WRS*H2I+WIS*H2R
        WTEMP=WR                                                                
        WR=WR*WPR-WI*WPI+WR                                                     
        WI=WI*WPR+WTEMP*WPI+WI                                                  
      END do
      IF (ISIGN == 1) THEN
        H1R=DAT(1)
        DAT(1)=H1R+DAT(2)
        DAT(2)=H1R-DAT(2)
      ELSE                                                                      
        H1R=DAT(1)
        DAT(1)=C1*(H1R+DAT(2))
        DAT(2)=C1*(H1R-DAT(2))
        CALL FOUR1(DAT,N,-1)
      END IF
      RETURN                                                                    
      END subroutine realft
