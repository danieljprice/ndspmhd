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
 real :: massp,volume,totmass,dx
 real :: denszero,amp,wavekn
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
 dx = psep
if (ndim.ge.2) then
 xmin(2) = -3.*psep*sqrt(3./4.)
 xmax(2) = 3.*psep*sqrt(3./4.)
endif
if (ndim.ge.3) then
 xmin(3) =  -3.*psep*sqrt(6.)/3.
 xmax(3) = 3.*psep*sqrt(6.)/3.
endif
!--set up the uniform density grid
!
! npart = int((xmax(1)-xmin(1))/psep) !!int((1./psep)**3)
! call alloc(int(1.1*npart))
!
!--density perturbation
!
 amp = 1.0
 wavekn = 4. 

 print*,' setting up resistivity test '

!! call cp_distribute(rmin,rmax,psep,ntotal,x(1,1:npart),x(2,1:npart),x(3,1:npart),npart)
 call set_uniform_cartesian(2,psep,xmin,xmax,adjustbound=.true.,stretchfunc=rhofunc)
 npart = ntotal
 print*,'npart =',npart
!
!--determine particle mass
!
 denszero = 0.1
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
    uu(i) = 1.e-3 ! isothermal
    Bfield(1,i) = 0. !5
    Bfield(2,i) = 1.e-5*sin(2.*pi*(x(1,i)-xmin(1)))
 enddo
 print*,' diffusion time = ',(xmax(1)-xmin(1))**2/etamhd
!
!--allow for tracing flow
!
 if (trace) write(iprint,*) '  exiting subroutine setup'
  
 return

contains

 real function rhofunc(x)
  real, intent(in) :: x
  
  rhofunc = 1. + amp*sin(wavekn*2.*pi*x)
 
 end function rhofunc   

end

!
! use this routine to modify the dump upon code restart
!
subroutine modify_dump()
 implicit none

end subroutine modify_dump
