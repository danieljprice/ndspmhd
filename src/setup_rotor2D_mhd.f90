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
!!  Setup for the MHD Rotor problems in Toth (2000)                       !!
!!                                                                        !!
!!  dense, rotating disk of fluid at origin                               !!
!!  magnetic field is initially uniform in x-direction                    !!
!!  -> simplest way to set this up is to vary the particle mass           !!
!!                                                                        !!                                                                        !!
!!  Note for all MHD setups, only the magnetic field should be setup      !!
!!  Similarly the thermal energy is setup even if using total energy.     !!
!!                                                                        !!
!!------------------------------------------------------------------------!!

subroutine setup
!
!--include relevant global variables
!
 use dimen_mhd
 use debug
 use loguns
 use bound
 use eos
 use options
 use part
 use setup_params
 
 use uniform_distributions
 use cons2prim, only:primitive2conservative
!
!--define local variables
!      
 implicit none
 integer :: ipart
 real :: denszero,densdisk,przero,vzero,ftaper
 real :: pri,rdisk,rbuffer,radius
 real :: totmass,gam1,massp,const,psepdisk,totvol
 real, dimension(ndim) :: xorigin, dx, xmintemp, xmaxtemp
 real, dimension(ndimv) :: Bzero
 logical :: equalmass
!
!--check number of dimensions is right
!
 if (ndim.ne.2) stop ' ndim must be = 2 for mhd rotor problem'
 if (ndimv.ne.2) write(iprint,*) ' warning: best if ndimv=2 for this problem'
!
!--set boundaries
!                        
 ibound = 3        ! reflective ghosts (boundaries not important in this problem)
 nbpts = 0        ! no fixed particles
 xmin(:) = -0.5        ! unit square
 xmax(:) = 0.5
 const = sqrt(4.*pi) 
!
!--setup parameters for the problem
! 
 equalmass = .true.
 xorigin(:) = 0.0        ! co-ordinates of the centre of the initial blast
 rdisk = 0.1             ! radius of the initial disk
 rbuffer = 0. !!115      ! radius of the smoothed front
 vzero = 2.0             ! rotation speed of initial disk
 Bzero(:) = 0.
 if (imhd.ne.0) Bzero(1) = 5.0/const        ! uniform field in bx direction
 przero = 1.0              ! initial pressure
 denszero = 1.0            ! ambient density
 densdisk = 10.0           ! density of rotating disk
 
 gam1 = gamma - 1.

 write(iprint,*) 'two dimensional mhd rotor problem '
 write(iprint,10) densdisk,rdisk,bzero(1),vzero,przero
10 format(/,' central density  = ',f10.3,', disk radius = ',f6.3,/, &
            ' initial bx   = ',f6.3,', rotation = ',f6.3,', pressure = ',f6.3,/)
!
!--setup uniform density grid of particles (2d) 
!  (determines particle number and allocates memory)
!
 psepdisk = psep*(denszero/densdisk)**(1./ndim)
 write(iprint,*) 'psep in disk = ',psepdisk
 if (equalmass) then
    call set_uniform_cartesian(2,psep,xmin,xmax,rmin=rdisk)
    xmintemp = -rdisk
    xmaxtemp = rdisk
    call set_uniform_cartesian(2,psepdisk,xmintemp,xmaxtemp,rmax=rdisk)
 else
    call set_uniform_cartesian(2,psep,xmin,xmax)  ! 2 = close packed arrangement
 endif
 ntotal = npart
!
!--determine particle mass in ambient medium
!
 if (equalmass) then
    totvol = product(xmax(:)-xmin(:)) - pi*rdisk**2
    totmass = denszero*totvol + densdisk*pi*rdisk**2
    massp = totmass/float(ntotal)
 else
    totmass = denszero*product(xmax(:)-xmin(:))
    massp = totmass/float(ntotal) ! average particle mass
 endif
!
!--now assign particle properties
! 
 do ipart=1,ntotal
    dx(:) = x(:,ipart)-xorigin(:) 
    radius = sqrt(dot_product(dx,dx))
    if (radius.le.rdisk) then
       dens(ipart) = densdisk
       if (equalmass) then
          pmass(ipart) = massp
       else
          pmass(ipart) = massp*densdisk/denszero
       endif
       vel(1,ipart) = -vzero*(x(2,ipart)-xorigin(2))/rdisk
       vel(2,ipart) = vzero*(x(1,ipart)-xorigin(1))/rdisk
    elseif (radius.le.rbuffer) then        ! smooth edge with taper function (toth)
       ftaper = (rbuffer-radius)/(rbuffer - rdisk)
       dens(ipart) = denszero + (densdisk-denszero)*ftaper
       pmass(ipart) = massp*dens(ipart)/denszero
       vel(1,ipart) = -ftaper*vzero*(x(2,ipart)-xorigin(2))/radius
       vel(2,ipart) = ftaper*vzero*(x(1,ipart)-xorigin(1))/radius
    else
       pmass(ipart) = massp
       dens(ipart) = denszero
       vel(:,ipart) = 0.
    endif  
    pri = przero 
    Bfield(:,ipart) = Bzero(:)
    uu(ipart) = przero/(gam1*dens(ipart))
 enddo
!
!--make sure it is *really* in pressure equilibrium
!
 call primitive2conservative
 do ipart=1,ntotal   
    uu(ipart) = przero/(gam1*rho(ipart))
 enddo
 
 Bconst(:) = Bzero(:)
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
