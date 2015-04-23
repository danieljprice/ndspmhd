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
 if (trace) write(iprint,*) ' entering subroutine setup(resistivity test)'
!
!--set boundaries
! 	    
 ibound = 3     ! boundaries
 nbpts = 0      ! use ghosts not fixed
 xmin(1) = 0.   ! set position of boundaries
 xmax(1) = 1.
if (ndim.ge.2) then
 xmin(2) = -0.5*sqrt(3./4.)
 xmax(2) = 0.5*sqrt(3./4.)
endif
if (ndim.ge.3) then
 xmin(3) = -0.5*sqrt(6.)/3.
 xmax(3) = 0.5*sqrt(6.)/3.
endif
!--set up the uniform density grid
!
! npart = int((xmax(1)-xmin(1))/psep) !!int((1./psep)**3)
! call alloc(int(1.1*npart))
 
 print*,' setting up resistivity test '

!! call cp_distribute(rmin,rmax,psep,ntotal,x(1,1:npart),x(2,1:npart),x(3,1:npart),npart)
 call set_uniform_cartesian(2,psep,xmin,xmax,adjustbound=.false.)
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
    dens(i) = denszero
    pmass(i) = massp
    uu(i) = 10.0 ! isothermal
    Bfield(2,i) = 0.01*sin(2.*pi*(x(1,i)-xmin(1)))
 enddo
 print*,' diffusion time = ',(xmax(1)-xmin(1))**2/etamhd
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
