!!-----------------------------------------------------------------
!! Sets the logical unit numbers to use for output
!!-----------------------------------------------------------------

SUBROUTINE logun
 USE dimen_mhd
 USE loguns
 IMPLICIT NONE
      
! IF (ndim.EQ.1) THEN
    iprint = 6
! ELSE
!    iprint = 8  ! log file / screen
! ENDIF
 ievfile = 13
 idatfile = 12
 iread = 21  ! for reading from input file
 ireadf = 15 ! for reading filenames
      
 RETURN
END SUBROUTINE logun
