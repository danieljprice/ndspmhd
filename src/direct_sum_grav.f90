!!----------------------------------------------------------------------
!! Calculates the gravitational force by a direct summation over the 
!! particles. 
!!
!! Use this to check the accuracy of the tree code
!!
!!----------------------------------------------------------------------

SUBROUTINE direct_sum_grav
 USE dimen_mhd
 USE debug
 USE loguns
 
 USE part
 USE gravity
 IMPLICIT NONE
 INTEGER :: i,j
 REAL, DIMENSION(ndim) :: dx,term
 REAL :: rij,rij2,pmassi
!
!--allow for tracing flow
!      
 IF (trace) WRITE(iprint,*) ' Entering subroutine direct_sum_grav'
!
!--reset forces initially
!
 fgrav = 0.
!
!--calculate gravitational force by direct summation
!
 DO i=1,npart
    pmassi = pmass(i)
    
    DO j=i+1,npart
       
       dx = x(:,i) - x(:,j)
       rij2 = DOT_PRODUCT(dx,dx) + 0.01**2
       rij = SQRT(rij2)
       term(:) = dx(:)/(rij*rij2)
       
       fgrav(:,i) = fgrav(:,i) - pmass(j)*term(:) ! note that fgrav should be
       fgrav(:,j) = fgrav(:,j) + pmassi*term(:)   ! *added* to the force
       
    ENDDO
 ENDDO

 RETURN
END SUBROUTINE direct_sum_grav
