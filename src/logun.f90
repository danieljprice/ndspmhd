!!-----------------------------------------------------------------
!! Sets the logical unit numbers to use for output
!!-----------------------------------------------------------------

SUBROUTINE logun
 USE loguns
 IMPLICIT NONE
      
 iprint = 6		! log file / screen
 ievfile = 13
 idatfile = 12
 iread = 21		! for reading from input file
 ireadf = 15		! for reading filenames
      
 RETURN
END SUBROUTINE logun
