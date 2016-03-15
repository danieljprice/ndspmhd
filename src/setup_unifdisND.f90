!------------------------------------------------------------------------------!
! NDSPMHD: A Smoothed Particle (Magneto)Hydrodynamics code for (astrophysical) !
! fluid dynamics simulations in 1, 2 and 3 spatial dimensions.                 !
!                                                                              !
! (c) 2002-2015 Daniel Price                                                   !
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

!----------------------------------------------------------------
!     Set up a uniform density cartesian grid of particles in ND
!----------------------------------------------------------------

subroutine setup
!
!--include relevant global variables
!
 use dimen_mhd
 use debug
 use loguns
 use bound
 use options
 use part
 use setup_params
 
 use uniform_distributions
!
!--define local variables
!            
 implicit none
 integer :: i
 real :: massp,volume,totmass
 real :: denszero
!
!--allow for tracing flow
!
 if (trace) write(iprint,*) ' entering subroutine setup(unifdis)'
!
!--set boundaries
! 	    
 ibound(:) = 3
 ibound(1) = 1     ! boundaries
 nbpts = 0      ! use ghosts not fixed
 xmin(:) = 0.   ! set position of boundaries
 xmax(:) = 1.
 
 call set_uniform_cartesian(1,psep,xmin,xmax,adjustbound=.true.)
 npart = ntotal
 print*,'npart =',npart
!
!--determine particle mass
!
 denszero = 1.0
 volume = product(xmax(:)-xmin(:))
 totmass = denszero*volume
 massp = totmass/float(ntotal) ! average particle mass
!
!--now assign particle properties
! 
 do i=1,ntotal
    vel(:,i) = 0.
    !vel(1,i) = vx(x(:,i))
    !vel(2,i) = vy(x(:,i))
    !vel(3,i) = vz(x(:,i))
    dens(i) = denszero
    pmass(i) = massp
    uu(i) = 1.0 ! isothermal
    if (idiffuse > 0) then
       !uu(i) = 1. + vx(x(:,i))
       if (x(1,i) < 0.5) then
          uu(i) = 1.
       else
          uu(i) = 2.
       endif
    endif
    Bfield(:,i) = 0.
 enddo 
!
!--allow for tracing flow
!
 if (trace) write(iprint,*) '  exiting subroutine setup'
  
 return
 
contains

!----------------------------------------------------------------
!+
!  functional form for v and its derivatives
!+
!----------------------------------------------------------------
real function vx(xyzhi)
 use bound
 use setup_params, only:pi
 real, intent(in) :: xyzhi(3)
 real :: dxbound
 
 dxbound = xmax(1) - xmin(1)

 vx = 0.5/pi*dxbound*sin(2.*pi*(xyzhi(1)-xmin(1))/dxbound)

end function vx

real function vy(xyzhi)
 use bound
 use setup_params, only:pi
 real, intent(in) :: xyzhi(3)
 real :: dxbound(3)
 
 dxbound(:) = xmax(:) - xmin(:)

 vy = 0.5/pi*dxbound(1)*sin(2.*pi*(xyzhi(1)-xmin(1))/dxbound(1)) &
     - 0.5/pi*dxbound(3)*sin(2.*pi*(xyzhi(3)-xmin(3))/dxbound(3))

end function vy

real function vz(xyzhi)
 use bound
 use setup_params, only:pi
 real, intent(in) :: xyzhi(3)
 real :: dxbound(3)
 
 dxbound(:) = xmax(:) - xmin(:)

 vz = 0.05/pi*dxbound(2)*cos(4.*pi*(xyzhi(2)-xmin(2))/dxbound(2))

end function vz

end subroutine setup

!
! use this routine to modify the dump upon code restart
!
subroutine modify_dump()
 implicit none

end subroutine modify_dump

