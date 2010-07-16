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
 if (trace) write(iprint,*) ' entering subroutine setup(2D shear)'
 write(iprint,*) '2D cartesian shear flow'
!
!--set boundaries
! 	    
 ibound = 3     ! boundaries
 nbpts = 0      ! use ghosts not fixed
 xmin(:) = 0.0  ! set position of boundaries
 xmax(:) = 1.0
!
!--set up the uniform density grid
!
 call set_uniform_cartesian(1,psep,xmin,xmax,fill=.true.)
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
 !cs0 = 0.05
 cs0 = 1.0
 uuzero = cs0**2/(gamma*(gamma-1))
 polyk0 = cs0**2/(gamma*denszero**(gamma-1.))
 write(iprint,*) ' cs0 = ',cs0, ' u = ',uuzero
 write(iprint,*) ' polyk = ',polyk,' should be = ',polyk0
!
!--now assign particle properties
! 
 do i=1,ntotal
    vel(:,i) = 0.
    vel(1,i) = 0.5*cos(2.*pi*x(2,i))
    !vel(2,i) = sin(2.*pi*x(1,i))
    dens(i) = denszero
    pmass(i) = massp
    uu(i) = uuzero
    bfield(:,i) = 0.
 enddo 
!
!--allow for tracing flow
!
 if (trace) write(iprint,*) '  exiting subroutine setup'
  
 return
end
