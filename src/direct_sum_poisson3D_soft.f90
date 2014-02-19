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

subroutine direct_sum_poisson_soft(x,pmass,hh,phi,fgrav,phitot,ntot)
 use dimen_mhd, only:ndim
 use debug, only:trace
 use loguns, only:iprint
 use kernels, only:interpolate_softening,radkern2,potensoft
 use hterms, only:gradh,gradhn,gradsoft
 use options, only:igravity
 
 implicit none
 integer, intent(in) :: ntot
 real, dimension(ndim,ntot), intent(in) :: x
 real, dimension(ntot), intent(in) :: pmass,hh
 real, dimension(ntot), intent(out) :: phi
 real, dimension(ndim,ntot), intent(inout) :: fgrav
 real, intent(out) :: phitot
 integer :: i,j
 real, dimension(ndim) :: dx,dr,fgravi,xi
 real :: rij,rij1,rij2,rij21,pmassi,pmassj
 real :: gradhi,gradhni,gradsofti,grkerni,grkernj,dsofti,dsoftj
 real :: phii,phij,phiterm,fm,fmi,fmj,phitemp
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
    fgravi(:) = 0.
    xi(:) = x(:,i)
    phitemp = 0.
    
    do j=i+1,ntot
       dx(1) = xi(1) - x(1,j)
       dx(2) = xi(2) - x(2,j)
       dx(3) = xi(3) - x(3,j)
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
          if (igravity.ge.4) then
             grkerni = grkerni*hi1**(ndim+1)*gradhi !!(1. + gradhni*gradhi/pmassj)
             dsofti = 0.5*grkerni*gradsofti
             fgravi(:) = fgravi(:) - pmassj*dsofti*dr(:)
             fgrav(:,j) = fgrav(:,j) + pmassi*dsofti*dr(:)
          endif
       else
          phii = -rij1
          fmi = rij21
       endif
       if (q2j.lt.radkern2) then
          call interpolate_softening(q2j,phij,fmj,grkernj)
          phij = phij*hj1
          fmj = fmj*hj21
          if (igravity.ge.4) then
             grkernj = grkernj*hj1**(ndim+1)*gradh(j) !!(1. + gradhn(j)*gradh(j)/pmassi)
             dsoftj = 0.5*grkernj*gradsoft(j)
             fgravi(:) = fgravi(:) - pmassj*dsoftj*dr(:)
             fgrav(1,j) = fgrav(1,j) + pmassi*dsoftj*dr(1)
             fgrav(2,j) = fgrav(2,j) + pmassi*dsoftj*dr(2)
             fgrav(3,j) = fgrav(3,j) + pmassi*dsoftj*dr(3)
          endif
       else
          phij = -rij1
          fmj = rij21
       endif

       phiterm = 0.5*(phii + phij)
       phitemp = phitemp + pmassj*phiterm
       phi(j) = phi(j) + pmassi*phiterm

       fm = 0.5*(fmi + fmj)
       fgravi(1) = fgravi(1) - pmassj*dr(1)*fm
       fgravi(2) = fgravi(2) - pmassj*dr(2)*fm
       fgravi(3) = fgravi(3) - pmassj*dr(3)*fm
       fgrav(1,j) = fgrav(1,j) + pmassi*dr(1)*fm
       fgrav(2,j) = fgrav(2,j) + pmassi*dr(2)*fm
       fgrav(3,j) = fgrav(3,j) + pmassi*dr(3)*fm
    enddo
!
!--add self contribution to potential
!
    fgrav(:,i) = fgrav(:,i) + fgravi(:)
    phi(i) = phitemp + pmassi*potensoft(0)*hi1
 enddo

 phitot = 0.
 do i=1,ntot
    phitot = phitot + 0.5*pmass(i)*phi(i)
 enddo

 return
end subroutine direct_sum_poisson_soft
