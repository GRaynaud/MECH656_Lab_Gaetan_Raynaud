!     Last change:  SB   21 Jun 1999   10:37 am
     SUBROUTINE FOUR1(DAT,NN,ISIGN)

      implicit none
      REAL(8):: WR,WI,WPR,WPI,WTEMP,THETA
      REAL:: TEMPR,TEMPI
      INTEGER, INTENT(IN):: NN, ISIGN
      REAL(8), INTENT(IN OUT):: DAT(2*NN)
      INTEGER:: N,M,J,I,MMAX,ISTEP
      N=2*NN                                                                    
      J=1                                                                       

      DO I=1,N,2
        IF(J>I)THEN
          TEMPR=DAT(J)
          TEMPI=DAT(J+1)
          DAT(J)=DAT(I)
          DAT(J+1)=DAT(I+1)         ! Bit reversal section
          DAT(I)=TEMPR              ! This is just a swap section...
          DAT(I+1)=TEMPI
        END IF
        M=nn
        do while((M>=2).AND.(J>M))
          J=J-M                                                                 
          M=M/2
        end do
        J=J+M                                                                   
      END do
      MMAX=2                                                                    
      DO WHILE (N>MMAX)
        ISTEP=2*MMAX                                                            
        THETA=6.28318530717959D0/(ISIGN*MMAX)                                   
        WPR=-2.0D0*SIN(0.5D0*THETA)**2
        WPI=SIN(THETA)
        WR=1.0D0
        WI=0.0D0
        DO M=1,MMAX,2
          DO I=M,N,ISTEP
            J=I+MMAX                                                            
            TEMPR=WR*DAT(J)-WI*DAT(J+1)
            TEMPI=WR*DAT(J+1)+WI*DAT(J)
            DAT(J)=DAT(I)-TEMPR
            DAT(J+1)=DAT(I+1)-TEMPI
            DAT(I)=DAT(I)+TEMPR
            DAT(I+1)=DAT(I+1)+TEMPI
          END do
          WTEMP=WR                                                              
          WR=WR*WPR-WI*WPI+WR                                                   
          WI=WI*WPR+WTEMP*WPI+WI                                                
        END do
        MMAX=ISTEP
      end do
      RETURN                                                                    
      END subroutine four1
