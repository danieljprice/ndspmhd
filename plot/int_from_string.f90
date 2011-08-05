!---------------------------------------------------------------------
!  function to convert a character string representing a number
!  to an integer value. Relies on ASCII character set.
!
!  returns zero on non-integer character
!  should do negative numbers
!
!  Daniel Price, Institute of Astronomy, Cambridge, May 2004
!---------------------------------------------------------------------

INTEGER FUNCTION int_from_string(string)
 IMPLICIT NONE
 INTEGER :: idigit,i,izero,ipower,maxdigits
 CHARACTER(LEN=*) :: string
 LOGICAL :: negative
 
 ipower = -1
 maxdigits = LEN(string)
 int_from_string = 0
 izero = IACHAR('0')    ! position of zero in ASCII character set
 negative = .false.

 DO i=maxdigits,1,-1    ! down through the characters
    IF (string(i:i).NE.' ' .AND. string(i:i).NE.'-') THEN
       ipower = ipower + 1
       idigit = IACHAR(string(i:i)) - izero
       IF ((idigit.LT.0).OR.(idigit.GT.10)) THEN
          !PRINT*, 'int_from_string: non-numeric character in input string'
          int_from_string = 0
          RETURN
       ELSE
          int_from_string = int_from_string + idigit*10**ipower
       ENDIF
    ELSEIF (string(i:i).EQ.'-') THEN
       negative = .true.
    ENDIF
 ENDDO   
 
 IF (negative) int_from_string = -int_from_string
 
END FUNCTION int_from_string
