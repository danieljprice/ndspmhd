!!-----------------------------------------------------------------
!! This subroutine returns the prsure and sound speed
!! from rho and/or u_therm via the equation of state
!!
!! size of array is specified on input, so can send in either whole
!! arrays or individual elements (just be consistent!)
!!-----------------------------------------------------------------

SUBROUTINE equation_of_state(pr,vsound,uu,rho,gamma,isize)
 USE dimen_mhd
 USE options
 USE loguns
 USE polyconst
!
!--define local variables
!
 IMPLICIT NONE
 INTEGER :: i
 INTEGER, INTENT(IN) :: isize
 REAL, INTENT(IN) :: gamma
 REAL, INTENT(IN), DIMENSION(isize) :: rho
 REAL, INTENT(OUT), DIMENSION(isize) :: pr
 REAL, INTENT(INOUT), DIMENSION(isize) :: uu,vsound
 REAL :: gamma1,B2

 gamma1 = gamma - 1.
!
!--exit gracefully if rho is negative
!
 IF (ANY(rho.LT.0.)) THEN
    WRITE(iprint,*) 'eos: rho -ve, exiting'
    DO i=1,isize
       IF (rho(i).LT.0.) WRITE(iprint,*) i,rho(i),uu(i)
    ENDDO
    CALL quit
 ELSEIF ((iener.NE.0).AND.ANY(uu.LT.0.)) THEN
    WRITE(iprint,*) 'eos: u_therm -ve, exiting'
    DO i=1,isize
       IF (uu(i).LT.0.) WRITE(iprint,*) i,rho(i),uu(i)
    ENDDO    
    CALL quit      
 ENDIF

 IF (iener.EQ.0) THEN	! polytropic (isothermal when gamma=1)
    WHERE (rho > 0.)
      pr = polyk*rho**gamma
      vsound = SQRT(gamma*pr/rho)
    END WHERE  
    IF (ABS(gamma1).GT.1.e-3) THEN
       WHERE (rho > 0.) 
       uu = pr/(gamma1*rho)    
       END WHERE
    ENDIF   
 ELSE		! adiabatic
    WHERE (rho > 0.)
      pr = gamma1*uu*rho
      vsound = SQRT(gamma*pr/rho)
    END WHERE  
 ENDIF
      
 RETURN
END SUBROUTINE equation_of_state
