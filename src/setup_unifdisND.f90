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
 real :: denszero,rmin,rmax
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
 
! xmin(1) = -0.5
! xmax(1) = 0.5
! xmin(2) = -0.5*sqrt(3./4.)
! xmax(2) = 0.5*sqrt(3./4.)!
! xmin(3) = -0.5*sqrt(6.)/3.
! xmax(3) = 0.5*sqrt(6.)/3.!!
!--set up the uniform density grid
!
! npart = int((xmax(1)-xmin(1))/psep) !!int((1./psep)**3)
! call alloc(int(1.1*npart))
 
 rmin = 0.
 rmax = 0.5

!! call cp_distribute(rmin,rmax,psep,ntotal,x(1,1:npart),x(2,1:npart),x(3,1:npart),npart)
 call set_uniform_cartesian(2,psep,xmin,xmax,adjustbound=.true.)
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
    vel(1,i) = 0.01*sin(2.*pi*(x(1,i)-xmin(1)))
    dens(i) = denszero
    pmass(i) = massp
    uu(i) = 1.0 ! isothermal
    Bfield(:,i) = 0.
 enddo 
!
!--allow for tracing flow
!
 if (trace) write(iprint,*) '  exiting subroutine setup'
  
 return
end
