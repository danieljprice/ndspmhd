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

!!------------------------------------------------------------------------!!
!!                                                                        !!
!!  Setup for the 2D current advection problem in MHD                     !!
!!   as in Gardiner & Stone (2005) J. Comp. Phys. 205, 509                !!
!!                                                                        !!
!!------------------------------------------------------------------------!!

subroutine setup
!
!--include relevant global variables
!
 use dimen_mhd, only:ndim,ndimV
 use debug, only:trace
 use loguns, only:iprint
 use bound
 use eos
 use options
 use part
 use setup_params
 
 use uniform_distributions
!
!--define local variables
!      
 implicit none
 integer :: ipart
 real :: denszero,przero,vzero,uuzero,Azero
 real :: totmass,gam1,massp,rr,rzero
!
!--check number of dimensions is right
!
 if (ndim.ne.2) stop ' ndim must be = 2 for current advection test'
!
!--set boundaries
!                        
 ibound = 3     ! periodic
 nbpts = 0      ! must use fixed particles if inflow/outflow at boundaries
 !costheta = 2./sqrt(5.)
 !sintheta = 1./sqrt(5.)
 xmin(1) = -1.0  ! x
 xmax(1) = -xmin(1)
 xmin(2) = -0.5  ! y
 xmax(2) = -xmin(2)
!
!--setup parameters
!
 vzero = sqrt(5.)
 przero = 1.0
 denszero = 1.0
 Azero = 1.e-3
 rzero = 0.3
 
 gam1 = gamma - 1.
 uuzero = przero/(gam1*denszero)

 write(iprint,*) 'Two dimensional current loop advection problem '
 write(iprint,10) 0.0,Azero,rzero,denszero,przero
10 format(/,' beta        = ',f6.3,', Azero = ',f6.3,/, &
            ' radius of loop = ',f6.3,', density = ',f6.3,', pressure = ',f6.3,/)
!
!--setup uniform density grid of particles (2D) with sinusoidal field/velocity
!  determines particle number and allocates memory
!
 call set_uniform_cartesian(2,psep,xmin,xmax)        ! 2 = close packed arrangement

 ntotal = npart
!
!--determine particle mass
!
 totmass = denszero*(xmax(2)-xmin(2))*(xmax(1)-xmin(1))
 massp = totmass/float(ntotal) ! average particle mass
!
!--now assign particle properties
! 
 do ipart=1,ntotal
    vel(1,ipart) = 2. !!vzero*costheta
    vel(2,ipart) = 1. !!vzero*sintheta
    if (ndimv.eq.3) vel(3,ipart) = 0.
    dens(ipart) = denszero
    pmass(ipart) = massp
    uu(ipart) = uuzero
    rr = sqrt(dot_product(x(:,ipart),x(:,ipart)))
    if (imhd.ge.1) then
       if (rr.lt.rzero) then
          Bfield(1,ipart) = Azero*(-x(2,ipart)/rr)
          Bfield(2,ipart) = Azero*(x(1,ipart)/rr)
       else
          Bfield(:,ipart) = 0.
       endif
       if (ndimv.eq.3) Bfield(3,ipart) = 0.0
    elseif (imhd.lt.0) then
!--vector potential setup
       Bevol(:,ipart) = 0.
       if (ndimV.lt.3) stop 'ndimV too small in setup_jadvect'
       Bevol(3,ipart) = MAX(Azero*(rzero-rr),0.)
    else
       Bfield(:,ipart) = 0.
    endif 
 enddo

!
!--allow for tracing flow
!
 if (trace) write(iprint,*) '  exiting subroutine setup'
            
 return
end subroutine setup

!
! use this routine to modify the dump upon code restart
!
subroutine modify_dump()
 implicit none

end subroutine modify_dump
