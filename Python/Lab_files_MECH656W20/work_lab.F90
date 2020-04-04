!     Last change:  LM   28 Mar 2003   12:24 pm
    SUBROUTINE work()
    USE Common_lab
    implicit none

! This subroutine translate the recorded voltages into real U velocities.

! Local data for the subroutine
    CHARACTER:: dot*1, blank1*1
    INTEGER:: Id, idbl
    INTEGER:: i,i2                    ! Dummy variables used in Do-loops.
    REAL(8):: ave, gain, ofset, vin
    CHARACTER:: FILE_name*50              ! File_name is the file to be open

! Now, we construct the appropriate name from where we read the data.
    dot = '.'
    blank1 = ' '
    READ(50,1001) fil_nam
1001 FORMAT (a15)
    Id = 0
    Id = INDEX(fil_nam,dot)
    IF(Id == 0) then
      idbl = INDEX(fil_nam, blank1)
      file_name = 'c:\data\'//fil_nam(1:idbl-1)//dot//today
    else                                    ! file type has been specified
      file_name = 'c:\data\'//fil_nam
    END if
    WRITE(*,*) "Here is the actual file name from where we will read data:"
    WRITE(*,*) "     ",file_name     ! inform the user what name has the data file.
    WRITE(99,*) "Here is the actual file name from where we will read data:"
    WRITE(99,*) "     ",file_name     ! inform the user what name has the data file.

! Now we are ready to read data. First, get specification parameters.
    OPEN(UNIT = lun(10), FILE = file_name, STATUS='old')
    READ(lun(10),*) nblock       ! number of block
    READ(lun(10),*) lbuff        ! buffer length
    READ(lun(10),*) n_var_sample
    READ(lun(10),*) n_chan
    READ(lun(10),*) frequency    ! sampling frequency

! inform the user on the specifications of his data file.
    WRITE(*,111) nblock
    WRITE(*,112) lbuff
    write(*,113) n_var_sample
    write(*,114) frequency
    WRITE(99,111) nblock
    WRITE(99,112) lbuff
    write(99,113) n_var_sample
    write(99,114) frequency
111 FORMAT(" Number of block in data file: ",i4)
112 FORMAT(" Buffer's length: ",i4)
113 FORMAT(" Number of different variables sampled: ",i1)
114 FORMAT(" The sampling frequency is: ",f8.1)


! Note: workcheck is there to erase previous work file if, and only if the work
!       subroutine is called more than once by one run of the main program.
    if (workcheck /= 0) then
      OPEN(UNIT = lun(9), FILE= fil_nam_old, STATUS = 'old')
      CLOSE(UNIT = lun(9), STATUS = 'delete')
    end if

! Now, we can construct the new work file.
    if (id == 0) then
      fil_nam_old = 'c:\results\U_series_'//fil_nam//'.txt'
    else
      fil_nam_old = 'c:\results\U_series_'//fil_nam(1:Id-1)//'.txt'
    END if
    OPEN(UNIT = lun(9), FILE= fil_nam_old, STATUS = 'new')

! Read the gain and offset for this data file:
    READ(50,*) gain
    READ(50,*) ofset
    IF(gain < 1.0) then
      WRITE(*,*) "Wrong gain entry. You need to enter a value larger than 1."
      WRITE(*,*) "Program stops!"
      stop
    END if
    IF(ofset < 0.0) then
      WRITE(*,*) "Wrong offset entry. You need to enter a value larger than 0."
      WRITE(*,*) "Program stops!"
      stop
    END if

! Let user know what is processing
      WRITE(*,*) "          **********************************************"
      WRITE(*,*) "          *                                            *"
      WRITE(*,*) "          *         Now translating voltages to        *"
      WRITE(*,*) "          *               real velocities              *"
      WRITE(*,*) "          *              Data written into:            *"
      WRITE(*,101) fil_nam_old
      WRITE(*,*) "          *                                            *"
      WRITE(*,*) "          **********************************************"
      WRITE(99,*) "          **********************************************"
      WRITE(99,*) "          *                                            *"
      WRITE(99,*) "          *         Now translating voltages to        *"
      WRITE(99,*) "          *               real velocities              *"
      WRITE(99,*) "          *              Data written into:            *"
      WRITE(99,101) fil_nam_old
      WRITE(99,*) "          *                                            *"
      WRITE(99,*) "          **********************************************"
101   FORMAT("           *  ",a40,"  *")

! We are ready to write data.

      do i = 1,nblock
        do i2 = 1,lbuff
          READ(lun(10),*) rbuff(i2)
          vin = rbuff(i2)/gain+ofset
          rbuff(i2) = ((vin*vin - cal(1))/cal(2))**(1.0/cal(3))
          WRITE(lun(9),*) rbuff(i2)
        end do
      end do

! Since we need the average of all variable in most of the subroutines,
! it is better to calculate it once and store it in the module.
    WRITE(*,*) "          **********************************************"
    WRITE(*,*) "          *                                            *"
    WRITE(*,*) "          *        Now computing the average for       *"
    WRITE(*,*) "          *         every variable in work file        *"
    WRITE(*,*) "          *                                            *"
    WRITE(*,*) "          **********************************************"
    WRITE(99,*) "          **********************************************"
    WRITE(99,*) "          *                                            *"
    WRITE(99,*) "          *        Now computing the average for       *"
    WRITE(99,*) "          *         every variable in work file        *"
    WRITE(99,*) "          *                                            *"
    WRITE(99,*) "          **********************************************"

! Compute the average now.
    REWIND(lun(9))
    average = 0.0
    do i = 1,nblock
      READ(lun(9),*) (rbuff(i2), i2 = 1,lbuff)
      ave = 0.0
      do i2 = 1,lbuff
          ave = ave + rbuff(i2)
      end do
      average = average + ave
    end do
    average = average/(nblock*lbuff)    ! Average is stored in Common_Variable now.
! Close units now.
    CLOSE(UNIT=lun(10))
    CLOSE(UNIT=lun(9))
    return
    END SUBROUTINE work
