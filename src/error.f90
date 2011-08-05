!!--------------------------------------------------------------
!! Prints appropriate error messages and exits gracefully
!! (NOT YET IMPLEMENTED)
!!--------------------------------------------------------------

SUBROUTINE error(ierr,whereami)
 USE loguns
 IMPLICIT NONE
 INTEGER, INTENT(IN) :: ierr
 CHARACTER(len=10), INTENT(IN) :: whereami

 IF (whereami(1:4) = 'link') THEN
    IF (ierr.EQ.1) WRITE(iprint,10) whereami,' max h <=0'
 ENDIF

10 FORMAT(' Error in subroutine ',a10,/,1x,a)

 CALL quit
 
END SUBROUTINE error
