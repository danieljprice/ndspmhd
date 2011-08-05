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
