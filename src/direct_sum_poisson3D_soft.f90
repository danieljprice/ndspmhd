!!----------------------------------------------------------------------------
!! Calculates the solution to the gravitational Poisson equation
!!
!! \nabla^2 \phi = 4 *pi*G \rho 
!!
!! With force softening from the kernel
!! by a direct summation over the particles. 
!!
!! Use this to check the accuracy of the tree code
!!
!! Input: 
!!
!!   x(ndim,ntot)  : co-ordinates of the particles
!!   pmass(ntot)   : particle masses
!!   hh(ntot)      : particle smoothing (softening) lengths
!!
!! Output:
!!
!!   phitot          : total potential phi
!!   fgrav(ndim,ntot) : gravitational force
!!----------------------------------------------------------------------------

subroutine direct_sum_poisson_soft(x,pmass,hh,phitot,fgrav,ntot)
 use dimen_mhd, only:ndim
 use debug, only:trace
 use loguns, only:iprint
 use kernels, only:interpolate_softening,radkern2,potensoft
 use hterms, only:gradh,gradhn,gradsoft
 
 implicit none
 integer, intent(in) :: ntot
 real, dimension(ndim,ntot), intent(in) :: x
 real, dimension(ntot), intent(in) :: pmass,hh
 real, dimension(ntot) :: phi
 real, dimension(ndim,ntot), intent(inout) :: fgrav
 real, intent(out) :: phitot
 integer :: i,j
 real, dimension(ndim) :: dx,dr
 real :: rij,rij1,rij2,rij21,pmassi,pmassj
 real :: gradhi,gradhni,gradsofti,grkerni,grkernj,dsofti,dsoftj
 real :: phii,phij,phiterm,fm,fmi,fmj
 real :: hi1,hj1,hi21,hj21,q2i,q2j
!
!--allow for tracing flow
!      
 if (trace) write(iprint,*) ' Entering subroutine direct_sum_poisson'
!
!--reset potential (but not force) initially
!
 phi = 0.
!
!--calculate gravitational force by direct summation on all particles
!
 do i=1,ntot
    pmassi = pmass(i)
    hi1 = 1./hh(i)
    hi21 = hi1*hi1
    gradhi = gradh(i)
    gradhni = gradhn(i)
    gradsofti = gradsoft(i)
    
    do j=i+1,ntot
       dx = x(:,i) - x(:,j)
       rij2 = dot_product(dx,dx)
       rij = sqrt(rij2)
       rij1 = 1./rij
       rij21 = rij1*rij1
       dr(:) = dx(:)*rij1
       hj1 = 1./hh(j)
       hj21 = hj1*hj1
       pmassj = pmass(j)
       
       q2i = rij2*hi21
       q2j = rij2*hj21
       if (q2i.lt.radkern2) then
          call interpolate_softening(q2i,phii,fmi,grkerni)
          phii = phii*hi1
          fmi = fmi*hi21
          grkerni = grkerni*hi1**(ndim+1)*gradhi !!(1. + gradhni*gradhi/pmassj)
          dsofti = 0.5*grkerni*gradsofti
          fgrav(:,i) = fgrav(:,i) - pmassj*dsofti*dr(:)
          fgrav(:,j) = fgrav(:,j) + pmassi*dsofti*dr(:)
       else
          phii = -rij1
          fmi = rij21
       endif
       if (q2j.lt.radkern2) then
          call interpolate_softening(q2j,phij,fmj,grkernj)
          phij = phij*hj1
          fmj = fmj*hj21
          grkernj = grkernj*hj1**(ndim+1)*gradh(j) !!(1. + gradhn(j)*gradh(j)/pmassi)
          dsoftj = 0.5*grkernj*gradsoft(j)
          fgrav(:,i) = fgrav(:,i) - pmassj*dsoftj*dr(:)
          fgrav(:,j) = fgrav(:,j) + pmassi*dsoftj*dr(:)
       else
          phij = -rij1
          fmj = rij21
       endif

       phiterm = 0.5*(phii + phij)
       phi(i) = phi(i) + pmassj*phiterm
       phi(j) = phi(j) + pmassi*phiterm

       fm = 0.5*(fmi + fmj)
       fgrav(:,i) = fgrav(:,i) - pmassj*dr(:)*fm
       fgrav(:,j) = fgrav(:,j) + pmassi*dr(:)*fm
    enddo
!
!--add self contribution to potential
!
    phi(i) = phi(i) + pmassi*potensoft(0)*hi1
 enddo

 phitot = 0.
 do i=1,ntot
    phitot = phitot + 0.5*pmass(i)*phi(i)
 enddo

 return
end subroutine direct_sum_poisson_soft
