!!------------------------------------------------------------------------!!
!!                                                                        !!
!!  Generic setup for the 2D Orszag-Tang vortex test in MHD               !!
!!                                                                        !!
!!  Note for all MHD setups, only the magnetic field should be setup      !!
!!  similarly the thermal energy is setup even if using total energy.     !!
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
 
 use uniform_distributions
!
!--define local variables
!      
 implicit none
 integer :: ipart
 real :: betazero,denszero,przero,vzero,bzero,uuzero,machzero
 real :: totmass,gam1,massp,const
!
!--check number of dimensions is right
!
 if (ndim.ne.2) stop ' ndim must be = 2 for Orszag-Tang vortex'
!
!--set boundaries
!                        
 ibound = 3     ! periodic
 nbpts = 0      ! must use fixed particles if inflow/outflow at boundaries
 xmin(1) = 0.0  ! x
 xmax(1) = 1.0
 xmin(2) = 0.0  ! y
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
 
 gam1 = gamma - 1.
 uuzero = przero/(gam1*denszero)

 write(iprint,*) 'Two dimensional Orszag-Tang vortex problem '
 write(iprint,10) betazero,machzero,bzero,denszero,przero
10 format(/,' beta        = ',f6.3,', mach number = ',f6.3,/, &
            ' initial B   = ',f6.3,', density = ',f6.3,', pressure = ',f6.3,/)
!
!--setup uniform density grid of particles (2D) with sinusoidal field/velocity
!  determines particle number and allocates memory
!
 call set_uniform_cartesian(22,psep,xmin,xmax)        ! 2 = close packed arrangement

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
       Bfield(1,ipart) = -Bzero*sin(2.*pi*x(2,ipart))
       Bfield(2,ipart) = Bzero*sin(4.*pi*x(1,ipart))
       if (ndimv.eq.3) Bfield(3,ipart) = 0.0
    elseif (imhd.lt.0) then
!--vector potential setup
       Bevol(:,ipart) = 0.
       Bevol(1,ipart) = 0.5/pi*Bzero*(cos(2.*pi*x(2,ipart)) + 0.5*cos(4.*pi*x(1,ipart)))
    else
       Bfield(:,ipart) = 0.
    endif 
 enddo

!
!--allow for tracing flow
!
 if (trace) write(iprint,*) '  exiting subroutine setup'
            
 return
end

subroutine modify_dump
 use loguns, only:iprint
 use part
 use options, only:imhd
 use timestep, only:time
 use setup_params, only:pi
 implicit none
 integer :: ipart
 real :: const,vzero,bzero
!
!--now assign particle properties
!
 write(iprint,*) 'modifying dump with Orzsag-Tang velocities/B field'
 vzero = 1.0
 const = 4.*pi
 bzero = 1.0/sqrt(const)
 do ipart=1,ntotal
    vel(1,ipart) = -vzero*sin(2.*pi*x(2,ipart))
    vel(2,ipart) = vzero*sin(2.*pi*x(1,ipart))
    if (ndimv.eq.3) vel(3,ipart) = 0.
    if (imhd.ge.1) then
       Bfield(1,ipart) = -Bzero*sin(2.*pi*x(2,ipart))
       Bfield(2,ipart) = Bzero*sin(4.*pi*x(1,ipart))
       if (ndimv.eq.3) Bfield(3,ipart) = 0.0
    elseif (imhd.lt.0) then
!--vector potential setup
       Bevol(:,ipart) = 0.
       Bevol(1,ipart) = 0.5/pi*Bzero*(cos(2.*pi*x(2,ipart)) + 0.5*cos(4.*pi*x(1,ipart)))
    else
       Bfield(:,ipart) = 0.
    endif 
 enddo
 time = 0.
 
end subroutine modify_dump
