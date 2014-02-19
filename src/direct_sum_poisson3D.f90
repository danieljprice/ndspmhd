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
!!   psoft         : softening length for plummer softening
!!
!! Output:
!!
!!   phitot          : total potential phi
!!   gradphi(ndim,ntot) : gradient of the potential (for gravity this = force)
!!----------------------------------------------------------------------------

subroutine direct_sum_poisson(x,source,phitot,gradphi,psoft,ntot)
 use dimen_mhd, only:ndim
 use debug, only:trace
 use loguns, only:iprint
 
 implicit none
 integer, intent(in) :: ntot
 real, dimension(ndim,ntot), intent(in) :: x
 real, dimension(ntot), intent(in) :: source
 real, dimension(ntot) :: phi
 real, dimension(ndim,ntot), intent(out) :: gradphi
 real, intent(in) :: psoft
 real, intent(out) :: phitot
 integer :: i,j
 real, dimension(ndim) :: dx,term
 real :: rij,rij2,sourcei
!
!--allow for tracing flow
!      
 if (trace) write(iprint,*) ' Entering subroutine direct_sum_poisson'
!
!--reset potential (but not force) initially
!
 phi = 0.
!
!--calculate gravitational force by direct summation
!
 do i=1,ntot
    sourcei = source(i)
    
    do j=i+1,ntot
       dx = x(:,i) - x(:,j)
       rij2 = dot_product(dx,dx) + psoft**2
       rij = sqrt(rij2)
       term(:) = dx(:)/(rij*rij2)
       
       phi(i) = phi(i) - source(j)/rij
       phi(j) = phi(j) - sourcei/rij
       
       gradphi(:,i) = gradphi(:,i) - source(j)*term(:)
       gradphi(:,j) = gradphi(:,j) + sourcei*term(:)
    enddo
 enddo

 phitot = 0.
 do i=1,ntot
    phitot = phitot + 0.5*source(i)*phi(i)
 enddo

 return
end subroutine direct_sum_poisson
