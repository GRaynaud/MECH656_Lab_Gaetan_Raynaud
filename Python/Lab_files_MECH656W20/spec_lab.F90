!     Last change:  SB   20 Mar 2003   11:03 am
    SUBROUTINE SPEC()

! 	This subroutine calculates the power spectrum using fast Fourrier transforms
! 	of a single variable in the work file.
    USE Common_lab
    implicit none
    REAL(8):: SIG1(lbuff),SIG2(lbuff)
    REAL(8):: ave, fq
    CHARACTER:: out_file*40, dot*1, blanc*1, fil_nam2*40
    INTEGER:: i,i2

    interface

    subroutine power_spec(sig)
      implicit none
      REAL(8), INTENT(OUT):: sig(:)
    end subroutine power_spec

    end interface

! Open the work file, read in variable number and save file.
    OPEN(UNIT= LUN(9), FILE= fil_nam_old, STATUS= 'old')
    READ(50,101) out_file
101 format (a25)

! Inform the user of what is processing...
    WRITE(*,*) "          **********************************************"
    WRITE(*,*) "          *                                            *"
    WRITE(*,*) "          *          Now computing the power           *"
    WRITE(*,*) "          *               spectrum of u                *"
    WRITE(*,*) "          *                                            *"
    WRITE(*,*) "          **********************************************"
    WRITE(99,*) "          **********************************************"
    WRITE(99,*) "          *                                            *"
    WRITE(99,*) "          *          Now computing the power           *"
    WRITE(99,*) "          *               spectrum of u                *"
    WRITE(99,*) "          *                                            *"
    WRITE(99,*) "          **********************************************"

! For each block, compute the average and the substract it from the values
! in the work file. Then we call power_spec.
! It is done block by block because the of the amount of data we send to power_spec.
! It needs to be a power of 2. We then, just have to take the average of all the
! spectra of the blocks.
    DO i = 1,NBLOCK
      ave=0.0
      DO i2= 1,LBUFF
        READ(LUN(9),*) rbuff(i2)
        ave = ave + rbuff(i2)
      END DO
      ave = ave/lbuff
      DO i2= 1,LBUFF ! center all the value around the average.
        rbuff(i2) = rbuff(i2) - ave
      END DO
! Now we can call power_spec
      call power_spec(sig1)

      IF(i==1) THEN
        DO i2=1,LBUFF/2+1
          SIG2(i2)=SIG1(i2)
        END DO
      else
        DO i2=1,LBUFF/2+1
          SIG2(i2)=SIG2(i2)+SIG1(i2)   ! Summation of all the spectra.
        END DO
      END IF
    end do ! End of the nblock do-loop.

    do i= 1,lbuff/2+1
      sig2(i) = sig2(i)/nblock   ! Divide by the number of block to obtain the average.
    end do

! CHECK IF FILE TYPE HAS BEEN SPECIFIED
    dot = '.'
    blanc = ' '
    i = INDEX(out_file,dot)
    if (i==0) then
      i2 = INDEX(out_file,blanc)
      fil_nam2 = out_file(1:i2-1)//dot//today(1:idate-1)//'.txt'
    else
      fil_nam2 = out_file//'.txt'        ! file type has been specified.
    END if
    WRITE(*,*) "The file name to which we will write data is:  ",fil_nam2
    WRITE(99,*) "The file name to which we will write data is:  ",fil_nam2

    OPEN(UNIT=lun(8), FILE='c:\results\'//fil_nam2, STATUS='replace')

! Now, we write data into the selected file.
    WRITE(LUN(8),*) "New buffer's length = ",lbuff*0.5
    DO i= 2,LBUFF/2+1
      FQ=(i-1.0)*FREQUENCY/LBUFF
      WRITE(LUN(8),103) FQ, SIG2(i)
103   FORMAT(tr5,es12.5,tr10,es12.5)
    END DO

! Program is almost finished, close unit and go back to main program.
    CLOSE(UNIT=lun(9))
    CLOSE(UNIT=lun(8))
    return
    end subroutine spec
