!!-----------------------------------------------------------------
!! This subroutine returns the prsure and sound speed
!! from rho and/or u_therm via the equation of state
!!
!! This version does one variable at a time
!!-----------------------------------------------------------------

SUBROUTINE equation_of_state(pr,vsound,uu,rho,gamma)
 USE dimen_mhd
 USE options
 USE loguns
 USE polyconst
 IMPLICIT NONE
 REAL, INTENT(IN) :: gamma,rho
 REAL, INTENT(OUT) :: pr
 REAL, INTENT(INOUT) :: uu,vsound
 REAL :: gamma1

 gamma1 = gamma - 1.
!
!--check for errors
!
 IF (uu.LT.0. .or. rho.LE.0.) THEN
    WRITE(iprint,*) 'Error: eos: rho,uu = ',rho,uu
    CALL quit
 ENDIF
!
!--polytropic (isothermal when gamma=1)
! 
 IF (iener.EQ.0) THEN
    pr = polyk*rho**gamma
    vsound = SQRT(gamma*pr/rho)
    IF (ABS(gamma1).GT.1.e-3) THEN
       uu = pr/(gamma1*rho)    
    ENDIF   
 ELSE		
!
!--adiabatic
!
    pr = gamma1*uu*rho
    vsound = SQRT(gamma*pr/rho)
 ENDIF
      
 RETURN
END SUBROUTINE equation_of_state
