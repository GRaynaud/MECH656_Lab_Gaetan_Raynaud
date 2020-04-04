!     Last change:  SB   20 Mar 2003   10:27 am
    SUBROUTINE POWER_SPEC(sig)

    USE common_lab
    implicit none
    REAL(8), INTENT(OUT):: SIG(:)
    REAL(8):: PI, M
    INTEGER:: i,n_half

    interface

    subroutine realft(ans,lbuff,i)
        implicit none
        INTEGER, INTENT(IN):: lbuff,i
        REAL(8), INTENT(IN OUT):: ans(2*lbuff)
      end subroutine realft

    end interface

    PI = ACOS(-1.0)
    M = 0.1*lbuff
    DO i = 1,INT(M)
      rbuff(i)=0.5*rbuff(i)*(1.0-COS(i*PI/M))
      rbuff(lbuff-INT(M)+I)=0.5*rbuff(lbuff-INT(M)+I)*(1.0+COS(i*PI/M))
    END DO

    N_half = lbuff*0.5

    Call realft(rbuff,n_half,1)

    Sig(1)=rbuff(1)*rbuff(1)

    do i = 2,n_half
      sig(i)=rbuff(2*i-1)**2 + rbuff(2*i)**2
    end do

    sig(n_half+1) = rbuff(2)*rbuff(2)

    DO i = 1,n_half+1
      SIG(i)=2.0*SIG(i)/(0.875*(frequency*lbuff))
    END DO
    RETURN
    end subroutine power_spec
