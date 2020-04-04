!     Last change:  SB   20 Mar 2003   11:04 am
    PROGRAM ANALYSE_XT

! This program accepts commands from the user and calls subroutines to
! perform velocity or temperature data reduction.

! The list of commands so far are:

!    EXIT -- exits the program gracefully.
!    KONS -- stores an array of constants used to convert voltages to
!            physical variables.
!    MNTG -- calculates the central moments upto sixth order. (global mean)
!    SPEC -- calculates the spectrum of any/all variable(s) in the work file.
!    TALE -- shows the time series of one block of data.
!    WORK -- Changes the voltages in the data file into a work file containing
!            real values, ie. velocities and temperatures.

    USE Common_lab
    implicit none
    INTEGER:: Integer_com_parts(4)
    CHARACTER(LEN=4):: com
    CHARACTER(LEN=1):: Char_com_parts(4), A, C, blanc
    INTEGER:: i,ios,j                          ! dummy variable.

! Variable initialization:
    fil_nam = '               '
    fil_nam_old = '                    '
    do i=1,10
      lun(i) = i
    end do

! Open an output file so that we can check it after the program ran.
    OPEN(UNIT=99,FILE='c:\results\output.txt',STATUS='replace')

! Tell to the user the program starts well:
    WRITE(*,*) "     ******************************************************"
    WRITE(*,*) "     *                                                    *"
    WRITE(*,*) "     *       Program ANALYSE_U has been launched.         *"
    WRITE(*,*) "     *        Please, be sure that you're batch           *"
    WRITE(*,*) "     *      file contains appropriate informations        *"
    WRITE(*,*) "     *         otherwise the program will stop.           *"
    WRITE(*,*) "     *                                                    *"
    WRITE(*,*) "     ******************************************************"
! All the 99 writing are for the output file.
    WRITE(99,*) "     ******************************************************"
    WRITE(99,*) "     *                                                    *"
    WRITE(99,*) "     *       Program ANALYSE_U has been launched.         *"
    WRITE(99,*) "     *        Please, be sure that you're batch           *"
    WRITE(99,*) "     *      file contains appropriate informations        *"
    WRITE(99,*) "     *         otherwise the program will stop.           *"
    WRITE(99,*) "     *                                                    *"
    WRITE(99,*) "     ******************************************************"

! Open the batch file and start to read parameter into it:
    OPEN(UNIT=50, FILE='c:\data\batch.input', STATUS='old',IOSTAT=ios)
    IF (ios/=0) then
      WRITE(*,*) " Error in opening file. Program stops!"
      WRITE(99,*) " Error in opening file. Program stops!"
      stop
    END if
    READ(50,101) today
! Count how long is the date.
    blanc = ' '
    idate = INDEX(today,blanc)
101 FORMAT(A11)
    WRITE(*,*) "Threating data on experiment held on ",today
    WRITE(99,*) "Threating data on experiment held on ",today

! We have the date, so we can open saving files.
    OPEN(UNIT=lun(1), FILE='c:\results\u_mom.'//today(1:idate-1)//'.txt', STATUS='replace')

! Printing headings in these files.
      A = 'u'
      C = CHAR(ICHAR(A)-32)        ! C is for Capital letter for instateneous average.
      WRITE(lun(1),102) C,A,A,A,A,A,A
102   FORMAT("File name",tr22,"<",a1,">",tr12,"<",a1,"**2>",tr9,"<"&
             &,a1,"**3>",tr9,"<",a1,"**4>"&
             &,tr9,a1," rms",tr10,a1," skew",tr9,a1," kurt",/)

! Read calibration constant like wire temperature.
    do i=1,3
      READ(50,*) cal(i)
    end do

!  Read the commands now.
    workcheck = 0    ! To check if we perform the work subroutine more than once.

    i = 1
    do WHILE(i==1)
      READ(50,104) com
104   format (a4)

! Make sure that com is made of capital characters only.
      DO j = 1,4
        Char_com_parts(j) = com(j:j)
        if (IACHAR(Char_com_parts(j))>96) then
          integer_com_parts(j) = IACHAR(Char_com_parts(j))-32
          Char_com_parts(j)=CHAR(integer_com_parts(j))
        end if
      END DO
      com = Char_com_parts(1)//Char_com_parts(2)//Char_com_parts(3)//Char_com_parts(4)
! Now com is completly made of capitals.

! Select which subroutine to perform with a select case construct.
      select case (com)
          case ('    ')
          case ('MNTG')
            call MNTG()
          case ('SPEC')
            call spec()
!          case ('TIME')
!            call tale()
          case ('WORK')
            call work()
            workcheck = workcheck + 1
          case ('EXIT')
            i = 0

! Tell to the user the program run without errors:
    WRITE(*,*) "     ******************************************************"
    WRITE(*,*) "     *                                                    *"
    WRITE(*,*) "     *          Program ANALYSE_U is finish.              *"
    WRITE(*,*) "     *             It ran without errors!                 *"
    WRITE(*,*) "     *                                                    *"
    WRITE(*,*) "     ******************************************************"
    WRITE(99,*) "     ******************************************************"
    WRITE(99,*) "     *                                                    *"
    WRITE(99,*) "     *          Program ANALYSE_U is finish.              *"
    WRITE(99,*) "     *             It ran without errors!                 *"
    WRITE(99,*) "     *                                                    *"
    WRITE(99,*) "     ******************************************************"

          case default
            WRITE(*,*) "Invalid code entry. List of commands so far:"
            WRITE(*,*) "    EXIT        MNTG       WORK"
            WRITE(*,*) "    SPEC"
            WRITE(99,*) "Invalid code entry. List of commands so far:"
            WRITE(*,*) "    EXIT        MNTG       WORK"
            WRITE(*,*) "    SPEC"
            stop
      end select
    end do

! Program is almost finish. Only need to close unit.
    do i = 1,10
      CLOSE(UNIT = lun(i))
    end do
    CLOSE(UNIT=50)
    CLOSE(UNIT=99)
    stop
    END PROGRAM analyse_xt
