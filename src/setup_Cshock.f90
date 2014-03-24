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
! Set up C shock problem in MacLow et al. (1995), ApJ 442, 726
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
 use eos
 
 use uniform_distributions
!
!--define local variables
!            
 implicit none
 integer :: i,nparty
 real :: massp,volume,totmass,rhoi,shockl,domainl
 real :: denszero,Bzero,vamach,csmach,theta,cs,vs
!
!--allow for tracing flow
!
 if (trace) write(iprint,*) ' entering subroutine setup (C-shock)'
!
!--set boundaries
!
 ibound = 3  ! y-z boundaries are periodic 
 ibound(1) = 1  ! x-boundary is fixed
 nbpts = 0      ! use ghosts not fixed
 vamach = 5.
 csmach = 50.
 cs = 0.1
 vs = cs*csmach
 theta = 0.25*pi
 denszero = 1.0
 rhoi  = 1.e-5
 Bzero = 1.
 print "(/,1x,a)",'C-shock: '
 gamma = 1.
 shockl = Bzero/(gamma*rhoi*sqrt(denszero))
 print "(a,es10.3)",'      neutral density = ',denszero
 print "(a,es10.3)",'          ion density = ',rhoi
 print "(a,es10.3)",'                gamma = ',gamma
 print "(a,es10.3)",'                   B0 = ',Bzero
 print "(a,es10.3)",'         shock length = ',shockl
 print "(a,es10.3)",' Alfvenic mach number = ',vs/(Bzero**2/denszero)
 print "(a,es10.3)",'    sonic mach number = ',csmach
 print "(a,es10.3)",'          sound speed = ',cs
 print "(a,es10.3,/)",'               vshock = ',vs
 nparty = 8
 domainl = 30.*shockl
 xmin(:) = -domainl - nparty*psep   ! set position of boundaries
 xmax(:) = domainl + nparty*psep
 if (ndim.ge.2) then
    xmin(2:ndim) = 0.0
    xmax(2:ndim) = xmin(2:ndim) + nparty*psep
 endif
! 
!--set up the uniform density grid
!
 call set_uniform_cartesian(2,psep,xmin,xmax,adjustbound=.true.)
 npart = ntotal
 print*,'npart =',npart
!
!--determine particle mass
!
 volume = product(xmax-xmin)
 totmass = denszero*volume
 massp = totmass/float(ntotal) ! average particle mass
!
!--now assign particle properties
!
 do i=1,ntotal
    vel(:,i) = 0.
    if (x(1,i) < -domainl) then
       itype(i) = itypebnd
       nbpts = nbpts + 1
       vel(1,i) = vs
    elseif (x(1,i) > domainl) then
       nbpts = nbpts + 1
       itype(i) = itypebnd
       vel(1,i) = -vs
    elseif (x(1,i) < 0.) then
       itype(i) = itypegas
       vel(1,i) = vs
    else
       itype(i) = itypegas
       vel(1,i) = -vs    
    endif
    dens(i) = denszero
    pmass(i) = massp
    uu(i) = 1.0 ! isothermal
    Bfield(:,i) = 0.
    if (imhd > 0) then
       Bfield(1,i) = Bzero*cos(theta)
       Bfield(2,i) = Bzero*sin(theta)
    endif
 enddo 
!
!--allow for tracing flow
!
 if (trace) write(iprint,*) '  exiting subroutine setup'
  
 return
end
