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
!! Calculates the 2D solution to a vector Poisson equation
!!
!! \nabla^2 {\bf A} = {\bf \eta} 
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
!!
!! Output:
!!
!!   curlA(ndimV,ntot) : curl of the vector potential
!!----------------------------------------------------------------------------

SUBROUTINE direct_sum_poisson_vec(x,sourcevec,curlA,ntot)
 USE dimen_mhd
 USE debug
 USE loguns
 
 IMPLICIT NONE
 INTEGER, INTENT(IN) :: ntot
 REAL, DIMENSION(ndim,ntot), INTENT(IN) :: x
 REAL, DIMENSION(ndimV,ntot), INTENT(IN) :: sourcevec
 REAL, DIMENSION(ndimV,ntot), INTENT(OUT) :: curlA
 INTEGER :: i,j
 REAL, DIMENSION(ndimV) :: dx
 REAL :: rij2,term
!
!--allow for tracing flow
!      
 IF (trace) WRITE(iprint,*) ' Entering subroutine direct_sum_poisson_vec'
!
!--reset forces initially
!
 curlA = 0.
 if (ndimV.ne.3) then
    write(iprint,*) 'ERROR: vector correct won''t work for ndimV < 3'
 endif
!
!--calculate gravitational force by direct summation
!
 DO i=1,ntot
    
    DO j=i+1,ntot
       dx(1:ndim) = x(:,i) - x(:,j)
       rij2 = DOT_PRODUCT(dx,dx) !!+ 0.01**2
       term = 1./rij2
       if (ndimV.gt.ndim) dx(ndim+1:ndimV) = 0
       
       curlA(1,i) = curlA(1,i) - (sourcevec(2,j)*dx(3) - sourcevec(3,j)*dx(2))*term
       curlA(2,i) = curlA(2,i) - (sourcevec(3,j)*dx(1) - sourcevec(1,j)*dx(3))*term
       curlA(3,i) = curlA(3,i) - (sourcevec(1,j)*dx(2) - sourcevec(2,j)*dx(1))*term
       curlA(1,j) = curlA(1,j) + (sourcevec(2,i)*dx(3) - sourcevec(3,i)*dx(2))*term
       curlA(2,j) = curlA(2,j) + (sourcevec(3,i)*dx(1) - sourcevec(1,i)*dx(3))*term
       curlA(3,j) = curlA(3,j) + (sourcevec(1,i)*dx(2) - sourcevec(2,i)*dx(1))*term

    ENDDO
 ENDDO

 RETURN
END SUBROUTINE direct_sum_poisson_vec
