!------------------------------------------------------------------------------!
! NDSPMHD: A Smoothed Particle (Magneto)Hydrodynamics code for (astrophysical) !
! fluid dynamics simulations in 1, 2 and 3 spatial dimensions.                 !
!                                                                              !
! (c) 2002-2010 Daniel Price                                                   !
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
 use eos, only:gamma
 use cons2prim 
 use uniform_distributions
!
!--define local variables
!            
 implicit none
 integer :: i
 real :: massp,volume,totmass
 real :: denszero,rmin,rmax,spsound2
!
!--allow for tracing flow
!
 if (trace) write(iprint,*) ' entering subroutine setup(unifdis)'
!
!--set boundaries
! 	    
 ibound = 3     ! boundaries
 nbpts = 0      ! use ghosts not fixed
 xmin(:) = 0.   ! set position of boundaries
 xmax(:) = 1.
!
!--set up the uniform density grid
! 
 rmin = 0.
 rmax = 0.5

 call set_uniform_cartesian(2,psep,xmin,xmax,fill=.true.)
 npart = ntotal
 print*,'npart =',npart
!
!--determine particle mass
!
 denszero = 5.0
 volume = product(xmax(:)-xmin(:))
 totmass = denszero*volume
 massp = totmass/float(ntotal) ! average particle mass
 spsound2 = 1.0
!
!--now assign particle properties
! 
 do i=1,ntotal
    vel(:,i) = 0.
    dens(i) = denszero
    pmass(i) = massp
    if (gamma.lt.1.0001) then
       uu(i) = 1.0    ! isothermal
    else
       uu(i) = spsound2/(gamma*(gamma-1.))
    endif
    Bfield(:,i) = 0.
 enddo
 
 print*,'sound speed = ',sqrt(spsound2), ' sound crossing time = ',(xmax(1)-xmin(1))/sqrt(spsound2)
!
!--allow for tracing flow
!
 if (trace) write(iprint,*) '  exiting subroutine setup'
  
 return
end
