!!------------------------------------------------------------------------!!
!!                                                                        !!
!!  Generic setup for the 2D Orszag-Tang vortex test in MHD               !!
!!                                                                        !!
!!  Gives particles a variable separation so the mass per SPH particle    !!
!!  is constant. shock is smoothed slightly (simple smoothing).           !!
!!                                                                        !!
!!  Note for all MHD setups, only the magnetic field should be setup      !!
!!  similarly the thermal energy is setup even if using total energy.     !!
!!                                                                        !!
!!------------------------------------------------------------------------!!

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
 integer :: i,j,ntot,npartx,nparty,ipart
 real :: betazero,denszero,przero,vzero,bzero,uuzero,machzero
 real :: totmass,gam1,massp,const
!
!--check number of dimensions is right
!
 if (ndim.ne.2) stop ' ndim must be = 2 for orszag-tang vortex'
!
!--set boundaries
!            	    
 ibound = 3	! periodic
 nbpts = 0	! must use fixed particles if inflow/outflow at boundaries
 xmin(1) = 0.0		! x
 xmax(1) = 1.0
 xmin(2) = 0.0		! y
 xmax(2) = 1.0
!
!--setup parameters
!
 const = 4.*pi
 betazero = 10./3.
 machzero = 1.0
 vzero = 1.0
 bzero = 1.0/sqrt(const)
 przero = 0.5*bzero**2*betazero
 denszero = gamma*przero*machzero
! denszero = 1.0
 
 gam1 = gamma - 1.
! npartx = int((xmax(1)-xmin(1))/psep)
! nparty = int((xmax(2)-xmin(2))/psep)
 uuzero = przero/(gam1*denszero)

 write(iprint,*) 'two dimensional orszag-tang vortex problem '
 write(iprint,10) betazero,machzero,bzero,denszero,przero
10 format(/,' beta        = ',f6.3,', mach number = ',f6.3,/, &
            ' initial b   = ',f6.3,', density = ',f6.3,', pressure = ',f6.3,/)
!
!--allocate memory here
!
! ntot = npartx*nparty
! call alloc(ntot)
! npart = ntot
! ntotal = ntot
!
!--setup uniform density grid of particles (2d) with sinusoidal field/velocity
!  determines particle number and allocates memory
 call set_uniform_cartesian(2,psep,xmin,xmax,.false.)	! 2 = close packed arrangement

 ntotal = npart
!
!--determine particle mass
!
 totmass = denszero*(xmax(2)-xmin(2))*(xmax(1)-xmin(1))
 massp = totmass/float(ntotal) ! average particle mass
!
!--now assign particle properties
! 
 do ipart=1,ntotal
    vel(1,ipart) = -vzero*sin(2.*pi*x(2,ipart))
    vel(2,ipart) = vzero*sin(2.*pi*x(1,ipart))
    if (ndimv.eq.3) vel(3,ipart) = 0.
    dens(ipart) = denszero
    pmass(ipart) = massp
    uu(ipart) = uuzero
    if (imhd.ge.1) then 
       bfield(1,ipart) = -bzero*sin(2.*pi*x(2,ipart))
       bfield(2,ipart) = bzero*sin(4.*pi*x(1,ipart))
       if (ndimv.eq.3) bfield(3,ipart) = 0.0	
    else
       bfield(:,ipart) = 0.
    endif 
!       print*,ipart,x(:,ipart),dens(ipart),uu(ipart),pmass(ipart)
!       if (mod(i,1000).eq.0) read*
 enddo

!
!--allow for tracing flow
!
 if (trace) write(iprint,*) '  exiting subroutine setup'
            
 return
end
