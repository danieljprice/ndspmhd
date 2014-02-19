!------------------------------------------------------------------------------!
! NDSPMHD: A Smoothed Particle (Magneto)Hydrodynamics code for (astrophysical) !
! fluid dynamics simulations in 1, 2 and 3 spatial dimensions.                 !
!                                                                              !
! (c) 2002-2014 Daniel Price                                                   !
!                                                                              !
! http://users.monash.edu.au/~dprice/ndspmhd                                   !
! daniel.price@monash.edu -or- dprice@cantab.net (forwards to current address) !
!                                                                              !
!  NDSPMHD comes with ABSOLUTELY NO WARRANTY.                                  !
!  This is free software; and you are welcome to redistribute                  !
!  it under the terms of the GNU General Public License                        !
!  (see LICENSE file for details) and the provision that                       !
!  this notice remains intact. If you modify this file, please                 !
!  note section 2a) of the GPLv2 states that:                                  !
!                                                                              !
!  a) You must cause the modified files to carry prominent notices             !
!     stating that you changed the files and the date of any change.           !
!                                                                              !
!  ChangeLog:                                                                  !
!------------------------------------------------------------------------------!

!!----------------------------------------------------------------------------
!! Calculates the 1D solution to any Poisson equation
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
 !!REAL, DIMENSION(ntot) :: phi
 REAL, DIMENSION(ndim,ntot), INTENT(OUT) :: gradphi
 REAL, INTENT(OUT) :: phitot
 INTEGER :: i,j
 REAL, DIMENSION(ndim) :: dx
 REAL :: sourcei
!
!--allow for tracing flow
!      
 IF (trace) WRITE(iprint,*) ' Entering subroutine direct_sum_poisson'
!
!--reset forces initially
!
 gradphi = 0.
 !!phi = 0.
!
!--calculate gravitational force by direct summation
!
 DO i=1,ntot
    sourcei = source(i)
    
!!    DO j=i+1,ntot
    DO j=1,ntot       
       if (i.ne.j) then
       dx = x(:,i) - x(:,j)
       
!       phi(i) = phi(i) + 
!       phi(j) = phi(j) + 
       
       gradphi(:,i) = gradphi(:,i) + source(j)*dx(:)
!!       gradphi(:,j) = gradphi(:,j) - sourcei*dx(:)
       endif
    ENDDO
 ENDDO

 phitot = 0.

 RETURN
END SUBROUTINE direct_sum_poisson
