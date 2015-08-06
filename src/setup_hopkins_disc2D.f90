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
 use eos, only:gamma,polyk 
 use uniform_distributions
 use externf, only:eps2_soft
!
!--define local variables
!            
 implicit none
 integer :: i
 real :: massp,volume,totmass
 real :: denszero,r,phi,gamm1,omega
!
!--allow for tracing flow
!
 if (trace) write(iprint,*) ' entering subroutine setup(unifdis)'
!
!--set boundaries
! 	    
 ibound = 0     ! boundaries
 nbpts = 0      ! use ghosts not fixed
 xmin(:) = -2.   ! set position of boundaries
 xmax(:) = 2.
 iexternal_force = 2
 
 if (ndim /= 2) stop 'need ndim=2'
 
 call set_uniform_cartesian(1,psep,xmin,xmax,rmin=0.5,rmax=2.)
 npart = ntotal
 print*,'npart =',npart
!
!--determine particle mass
!
 denszero = 1.0
 volume = product(xmax(:)-xmin(:))
 totmass = denszero*volume
 massp = totmass/float(ntotal) ! average particle mass
 gamm1 = gamma - 1.
!
!--now assign particle properties
! 
 do i=1,ntotal
    vel(:,i) = 0.
    dens(i) = denszero
    Bfield(:,i) = 0.
    
    r = sqrt(dot_product(x(:,i),x(:,i)))
    phi = atan2(x(2,i),x(1,i))
    !
    !--keplerian velocity profile balancing pressure gradient
    !
    omega = sqrt(1./(r**2 + eps2_soft)**1.5)
    vel(1,i) = -r*sin(phi)*omega
    vel(2,i) = r*cos(phi)*omega
    if (gamm1 > 0.) then
       uu(i) = polyk/gamm1*dens(i)**gamm1
    else
       uu(i) = 1.e-5
    endif 
    pmass(i) = massp

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
