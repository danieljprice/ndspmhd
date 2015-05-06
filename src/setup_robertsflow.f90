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

!!------------------------------------------------------------------------!!
!!                                                                        !!
!!  Setup for the Roberts flow with MHD                                   !!
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
 use mem_allocation, only:alloc
 
 use uniform_distributions
!
!--define local variables
!      
 implicit none
 integer :: ipart
 real :: betazero,denszero,przero,vzero,bzero,uuzero,machzero
 real :: totmass,gam1,massp,const,boxl,cszero
 real :: coskx,cosky,sinkx,sinky,kx,ky,fk,vz
!
!--check number of dimensions is right
!
 if (ndim.lt.2) stop ' ndim must be >= 2 for Roberts flow'
!
!--set boundaries
!                        
 ibound = 3     ! periodic
 nbpts = 0      ! must use fixed particles if inflow/outflow at boundaries
 boxl  = 1.
 xmin(:) = -0.5*boxl
 xmax(:) =  0.5*boxl
!
!--setup parameters
!
 gamma = 1.
 const = 4.*pi
 betazero = 100.
 machzero = 0.1
 denszero = 1.
 fk = 2.*pi/boxl
 vzero = 16.*pi
 cszero = vzero/machzero ! M = v/cs; cs = v/m
 przero = cszero**2*denszero
 
 gam1 = gamma - 1.
 if (gam1 <= 0.) then
    polyk = przero/denszero
    uuzero = 1.5*polyk
    print*,' polyk = ',polyk, ' mach = ',vzero/cszero
 else
    uuzero = przero/(gam1*denszero)
 endif
 print*,'end time = ',400./(vzero*fk)
 print*,' vzero = ',vzero,' uuzero = ',uuzero,cszero**2

 write(iprint,*) 'Three dimensional Roberts flow '
 if (ndim.ge.3) write(iprint,*) ' (in 3D...)'
 write(iprint,10) betazero,machzero,bzero,denszero,przero
10 format(/,' beta        = ',f6.3,', mach number = ',f6.3,/, &
            ' initial B   = ',f6.3,', density = ',f6.3,', pressure = ',f6.3,/)
!
!--setup uniform density grid of particles (2D) with sinusoidal field/velocity
!  determines particle number and allocates memory
!
 call set_uniform_cartesian(1,psep,xmin,xmax)        ! 2 = close packed arrangement
 ntotal = npart
!
!--determine particle mass
!
 totmass = denszero*product(xmax-xmin)
 massp = totmass/float(ntotal) ! average particle mass
!
!--now assign particle properties
! 

 do ipart=1,ntotal
    kx = fk*(x(1,ipart)-xmin(1))
    ky = fk*(x(2,ipart)-xmin(2))
    coskx = cos(kx)
    cosky = cos(ky)
    sinkx = sin(kx)
    sinky = sin(ky)
    vz = 0.  ! planar case
    !vz = coskx*cosky/sqrt(2.)
    vel(:,ipart) = vzero/fk*(/-sinky*coskx,sinkx*cosky,vz/)
    dens(ipart) = denszero
    pmass(ipart) = massp
    uu(ipart) = uuzero
    if (imhd.ge.1) then
       Bfield(:,ipart) = 0.
       Bfield(1,ipart) = Bzero
    elseif (imhd.lt.0) then
!--vector potential setup
       if (ndimV.lt.3) stop 'ndimV too small in setup(robertsflow)'
       Bevol(:,ipart) = 0.
       Bevol(3,ipart) = 0.
    else
       Bfield(:,ipart) = 0.
    endif 
 enddo
!
!--allow for tracing flow
!
 if (trace) write(iprint,*) '  exiting subroutine setup'

end subroutine setup

!
! use this routine to modify the dump upon code restart
!
subroutine modify_dump()
 implicit none

end subroutine modify_dump
