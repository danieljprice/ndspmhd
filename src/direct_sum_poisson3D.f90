!!----------------------------------------------------------------------------
!! Calculates the solution to any Poisson equation
!!
!! \nabla^2 \phi = \eta 
!!
!! by a direct summation over the particles. 
!!
!! Use this to check the accuracy of the tree code
!!
!! Input: 
!!
!!   x(ndim,ntot)  : co-ordinates of the particles
!!   source(ntot)  : quantity to be summed over 
!!                   for an arbitrary quantity \eta, this is given by
!!                   source = particle mass * eta / (4.*pi*density)
!!               ie. source = particle mass for gravity
!!
!! Output:
!!
!!   phitot          : total potential phi
!!   gradphi(ndim,ntot) : gradient of the potential (for gravity this = force)
!!----------------------------------------------------------------------------

SUBROUTINE direct_sum_poisson(x,source,phitot,gradphi,ntot)
 USE dimen_mhd
 USE debug
 USE loguns
 
 IMPLICIT NONE
 INTEGER, INTENT(IN) :: ntot
 REAL, DIMENSION(ndim,ntot), INTENT(IN) :: x
 REAL, DIMENSION(ntot), INTENT(IN) :: source
 REAL, DIMENSION(ntot) :: phi
 REAL, DIMENSION(ndim,ntot), INTENT(OUT) :: gradphi
 REAL, INTENT(OUT) :: phitot
 INTEGER :: i,j
 REAL, DIMENSION(ndim) :: dx,term
 REAL :: rij,rij2,sourcei
!
!--allow for tracing flow
!      
 IF (trace) WRITE(iprint,*) ' Entering subroutine direct_sum_poisson'
!
!--reset potential (BUT NOT FORCE) initially
!
 phi = 0.
!
!--calculate gravitational force by direct summation
!
 DO i=1,ntot
    sourcei = source(i)
    
    DO j=i+1,ntot
       dx = x(:,i) - x(:,j)
       rij2 = DOT_PRODUCT(dx,dx) + 0.1**2
       rij = SQRT(rij2)
       term(:) = dx(:)/(rij*rij2)
       
       phi(i) = phi(i) - source(j)/rij
       phi(j) = phi(j) - sourcei/rij
       
       gradphi(:,i) = gradphi(:,i) - source(j)*term(:)
       gradphi(:,j) = gradphi(:,j) + sourcei*term(:)
    ENDDO
 ENDDO

 phitot = 0.
 DO i=1,ntot
    phitot = phitot + 0.5*source(i)*phi(i)
 ENDDO
 print*,'phitot = ',phitot

 RETURN
END SUBROUTINE direct_sum_poisson
