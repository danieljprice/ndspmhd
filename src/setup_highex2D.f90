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

!!-------------------------------------------------------------------------!!
!!                                                                         !!
!!  setup for 2D high explosives test                                      !!
!!                                                                         !!
!!  dense, high pressure blob at the origin,                               !!
!!  get lovely Richtmeyer-Meshkov instabilities                            !!
!!                                                                         !!
!!-------------------------------------------------------------------------!!

subroutine setup
!
!--include relevant global variables
!
 use dimen_mhd, only:ndim,ndimV
 use debug, only:trace
 use loguns, only:iprint
 use bound
 use eos
 use kernels, only:interpolate_kernel
 use options, only:damp,ibound
 use part
 use setup_params, only:psep,pi
 
 use uniform_distributions
!
!--define local variables
!      
 implicit none
 integer :: ipart,n1
 real :: denszero,przero,densc,uuc,uuzero,rc,prc
 real :: totmass,gam1,massp,psepc,radius
 real, dimension(ndim) :: xmaxtemp

!
!--set boundaries
!            	    
 ibound = 2     ! fixed ghosts
 nbpts = 0
 xmin(:) = 0.     ! same xmin in all dimensions
 xmax(:) = 2.5
 denszero = 1.293
 densc = 1630.
 rc = 0.0527
 !
 !--choose time unit in ms, so pressure unit is 10^6 Pa
 !
 przero = 0.101325   ! 1 atm = 101325 in pascal
 prc = 27603.94*przero
 
 write(iprint,10) ndim
 write(iprint,20) prc,rc,denszero,przero
10 format(/,1x,i1,'-dimensional high explosives test')
20 format(/,' central pressure  = ',f12.3,', blast radius = ',f6.3,/, &
            ' ext. density = ',f6.3,', ext. pressure = ',f9.3)

 gam1 = gamma - 1.
 if (abs(gam1).lt.1.e-3) stop 'eos cannot be isothermal for this setup'
 uuzero = przero/(gam1*denszero)
 uuc = prc/(gam1*densc)  !4.29e6 J ! should 
 write(iprint,30) uuc,uuzero
30 format(' central energy = ',f9.3,', ext. energy = ',f9.3,/)

!
!--setup uniform density grid of particles
!  (determines particle number and allocates memory)
!
 psepc = psep*(denszero/densc)**(1./ndim)
 xmaxtemp = rc
 call set_uniform_cartesian(1,psepc,xmin,xmaxtemp,rmax=rc) ! 2 = close packed arrangement
 n1 = npart
 call set_uniform_cartesian(1,psep,xmin,xmax,rmin=rc)
 ntotal = npart
!
!--determine particle mass
!
 totmass = denszero*(product(xmax(:)-xmin(:)) - pi*rc**2)   ! assumes cartesian boundaries
 massp = totmass/float(ntotal - n1)                  ! average particle mass
!
!--now assign particle properties
! 
 do ipart=1,ntotal
    vel(:,ipart) = 0.
!--uniform density and smoothing length
    dens(ipart) = denszero

    radius = sqrt(dot_product(x(:,ipart),x(:,ipart)))
    if (radius.le.rc) then
       dens(ipart) = densc
       uu(ipart) = uuc
    else
       dens(ipart) = denszero
       uu(ipart) = uuzero
    endif
    pmass(ipart) = massp
 enddo

!
!--allow for tracing flow
!
 if (trace) write(iprint,*) '  exiting subroutine setup'
            
 return
end

!
! use this routine to modify the dump upon code restart
!
subroutine modify_dump()
 implicit none

end subroutine modify_dump
