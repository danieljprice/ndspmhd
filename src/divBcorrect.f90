!!--------------------------------------------------------------------------
!! Interface subroutine for the div B correction by the projection method
!!
!!--------------------------------------------------------------------------

SUBROUTINE divBcorrect
 USE dimen_mhd
 USE debug
 USE loguns
 
 USE part
 USE derivB
 USE options
 IMPLICIT NONE
 REAL, PARAMETER :: pi = 3.1415926536
 REAL, PARAMETER :: fourpi1 = 1./(4.*pi)
 INTEGER :: ntot
 REAL, DIMENSION(ntot) :: source,phi
 REAL, DIMENSION(ndim,ntot) :: gradphi
!
!--allow for tracing flow
!      
 IF (trace) WRITE(iprint,*) ' Entering subroutine divBcorrect'

 SELECT CASE(idivBzero)
 
    CASE(1)
!
!--specify the source term for the Poisson equation
!    
       source(1:ntot) = pmass(1:ntot)*divBonrho(1:ntot)*fourpi1
!
!--calculate the correction to the magnetic field
!
       CALL direct_sum_poisson(x,source,phi,gradphi,ntot)
!
!--correct the magnetic field
!              
       Bfield(:,:) = Bfield(:,:) - gradphi(:,:)
       
    CASE(2)
!       CALL direct_sum_current
       
 END SELECT
 
 RETURN
END SUBROUTINE divBcorrect
