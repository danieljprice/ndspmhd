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
 use options
 use part
 use setup_params
 
 use uniform_distributions
!
!--define local variables
!            
 implicit none
 integer :: i,iprofile
 real :: massp,volume,totmass
 real :: denszero,dustfrac0,A,B,r,rmin
!
!--allow for tracing flow
!
 if (trace) write(iprint,*) ' entering subroutine setup(dustydiffusion)'
!
!--set boundaries
! 	    
 ibound = 3     ! boundaries
 nbpts = 0      ! use ghosts not fixed
 xmin(:) = -0.5   ! set position of boundaries
 xmax(:) = 0.5
 
 call set_uniform_cartesian(2,psep,xmin,xmax,adjustbound=.true.)
 npart = ntotal
 print*,'npart =',npart
 if (idust.eq.2) then
    itype(1:npart) = itypegas
    print*,' setting up dust particles...'
    !call set_uniform_cartesian(2,psep,xmin,xmax,adjustbound=.true.)
 endif
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
 A = 0.05
 B = 0.1
 dustfrac0 = 0.1
 iprofile = 2
 rmin = psep
 do i=1,ntotal
    vel(:,i) = 0.
    dens(i) = denszero
    pmass(i) = massp
    uu(i) = 1.0 ! isothermal
    Bfield(:,i) = 0.
    if (idust.eq.1 .or. idust.eq.3 .or. idust.eq.4) then
       r = sqrt(dot_product(x(:,i),x(:,i)))
       select case(iprofile)
       case(3)
          dustfrac(i) = 0.5*exp(-(r/0.1)**2)
       case(2)
          dustfrac(i) = max(A*(B)**(-0.6) - r**2/B,0.001)
       case default
          if (r < 0.5) then
             dustfrac(i) = min(sqrt(A**2/r + dustfrac0**2),0.9)      
          else
             dustfrac(i) = dustfrac0
          endif       
       end select
    else
       dustfrac(i) = 0.
    endif 
 enddo 
!
!--allow for tracing flow
!
 if (trace) write(iprint,*) '  exiting subroutine setup'
  
 return
end
