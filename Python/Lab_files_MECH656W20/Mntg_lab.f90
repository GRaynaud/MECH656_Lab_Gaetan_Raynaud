!     Last change:  SB   20 Mar 2003   11:02 am
    SUBROUTINE MNTG()

! This subroutine calculates central moments of the blocks
! in the work file. The user must select which moments he/she wants
! (cross-moments or single moments).
    USE Common_lab
    implicit none
    REAL(8):: total
    REAL(8):: x(4),x1, xf(4)
    REAL(8):: rmsx,skewx,curtx
    INTEGER:: i,i2
    CHARACTER:: dot*1, blanc*1, file_name*40

! Variables initialization.
    do i = 1,4
      xf(i) = 0.0
    end do
    x1 = 0.0

! Inform the user that MNTG is being performed.
    WRITE(*,*) "          **********************************************"
    WRITE(*,*) "          *                                            *"
    WRITE(*,*) "          *             MNTG is computing              *"
    WRITE(*,*) "          *                  moments                   *"
    WRITE(*,*) "          *                                            *"
    WRITE(*,*) "          **********************************************"
    WRITE(99,*) "          **********************************************"
    WRITE(99,*) "          *                                            *"
    WRITE(99,*) "          *             MNTG is computing              *"
    WRITE(99,*) "          *                  moments                   *"
    WRITE(99,*) "          *                                            *"
    WRITE(99,*) "          **********************************************"

! Construct appropriate file name for writing results.
    dot = '.'
    blanc = ' '
    i = INDEX(fil_nam,dot)
    if (i==0) then
      i2 = INDEX(fil_nam,blanc)
      file_name = fil_nam(1:i2-1)//dot//today
    else
      file_name = fil_nam                       ! file type has been specified.
    END if

! Now, we can start computing the moments.
    OPEN(UNIT=lun(9),FILE=fil_nam_old,STATUS='old')
    total = lbuff*nblock

! Now, we calculate the moments around a zero mean series for greater accuracy.
    do i = 1,nblock
        x(2) = 0.0
        x(3) = 0.0
        x(4) = 0.0
        do i2 = 1,lbuff
          READ(lun(9),*) rbuff(i2)
          x1 = rbuff(i2)-average
          x(2) = x(2) + x1*x1
          x(3) = x(3) + x1**3     ! moments are computed for the first variable.
          x(4) = x(4) + x1**4
        END do
        xf(2) = xf(2) + x(2)
        xf(3) = xf(3) + x(3)
        xf(4) = xf(4) + x(4)
    end do
    xf(2) = xf(2)/total
    xf(3) = xf(3)/total
    xf(4) = xf(4)/total                  ! Here we compute the "mean" of the moments.
    rmsx = xf(2)**0.5
    skewx = xf(3)/(rmsx*xf(2))
    curtx = xf(4)/(xf(2)*xf(2))
    WRITE(lun(1),103) file_name, average, (xf(i2), i2=2,4), rmsx, skewx, curtx
103 FORMAT(a25,tr3,9(Es13.4,tr2))

    CLOSE(UNIT=lun(9))

    return
    END SUBROUTINE MNTG
