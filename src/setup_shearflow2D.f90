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
 use eos
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
 real :: denszero,uuzero,cs0,polyk0
!
!--allow for tracing flow
!
 if (trace) write(iprint,*) ' entering subroutine setup(unifdis)'
 write(iprint,*) '2d cartesian shear flow'
!
!--set boundaries
! 	    
 ibound = 3	! boundaries
 nbpts = 0	! use ghosts not fixed
 xmin(:) = 0.0  ! set position of boundaries
 xmax(:) = 1.0
!
!--set up the uniform density grid
!
 call set_uniform_cartesian(11,psep,xmin,xmax,.false.)
 ntotal = npart
!
!--determine particle mass
!
 denszero = 1.0
 volume = product(xmax(:)-xmin(:))
 totmass = denszero*volume
 massp = totmass/float(ntotal) ! average particle mass
!
!--sound speed
!
 cs0 = 0.05
 uuzero = cs0**2/(gamma*(gamma-1))
 polyk0 = cs0**2/(gamma*denszero**(gamma-1.))
 write(iprint,*) ' cs0 = ',cs0, ' u = ',uuzero
 write(iprint,*) ' polyk = ',polyk,' should be = ',polyk0
!
!--now assign particle properties
! 
 do i=1,ntotal
    vel(:,i) = 0.
    vel(2,i) = sin(2.*pi*x(1,i))
    dens(i) = denszero
    pmass(i) = massp
    uu(i) = uuzero ! isothermal
    bfield(:,i) = 0.
 enddo 
!
!--allow for tracing flow
!
 if (trace) write(iprint,*) '  exiting subroutine setup'
  
 return
end
