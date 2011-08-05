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
!--reset forces initially
!
 gradphi = 0.
 phi = 0.
!
!--calculate gravitational force by direct summation
!
 DO i=1,ntot
    sourcei = source(i)
    
!!    DO j=i+1,ntot
    DO j=1,ntot       
       if (i.ne.j) then
       dx = x(:,i) - x(:,j)
       rij2 = DOT_PRODUCT(dx,dx) + 0.01**2
       rij = SQRT(rij2)
       term(:) = dx(:)/(rij*rij2)
       
!       phi(i) = phi(i) + 
!       phi(j) = phi(j) + 
       
       gradphi(:,i) = gradphi(:,i) + source(j)*term(:)
!!       gradphi(:,j) = gradphi(:,j) - sourcei*term(:)
       endif
    ENDDO
 ENDDO

 phitot = 0.

 RETURN
END SUBROUTINE direct_sum_poisson
