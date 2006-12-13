!!------------------------------------------------------------------------!!
!!                                                                        !!
!!  Generic setup for the 2D current advection problem in MHD             !!
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
 real :: denszero,przero,vzero,uuzero,Azero
 real :: totmass,gam1,massp,rr,rzero,angle,anglev
!
!--check number of dimensions is right
!
 if (ndim.ne.2) stop ' ndim must be = 2 for current advection test'
!
!--set boundaries
!                        
 ibound = 3     ! periodic
 nbpts = 0      ! must use fixed particles if inflow/outflow at boundaries
 angle = 30.*pi/180.
 anglev = 60.*pi/180.
 xmin(1) = -1.0  ! x
 xmax(1) = -xmin(1)
 xmin(2) = -1.0/(2.*COS(angle))  ! y
 xmax(2) = -xmin(2)
!
!--setup parameters
!
 vzero = 1.0
 przero = 1.0
 denszero = 1.0
 Azero = 1.e-3
 rzero = 0.3
 
 gam1 = gamma - 1.
 uuzero = przero/(gam1*denszero)

 write(iprint,*) 'Two dimensional current loop advection problem '
 write(iprint,10) 0.0,Azero,rzero,denszero,przero
10 format(/,' beta        = ',f6.3,', Azero = ',f6.3,/, &
            ' radius of loop = ',f6.3,', density = ',f6.3,', pressure = ',f6.3,/)
!
!--setup uniform density grid of particles (2D) with sinusoidal field/velocity
!  determines particle number and allocates memory
!
 call set_uniform_cartesian(2,psep,xmin,xmax)        ! 2 = close packed arrangement

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
    vel(1,ipart) = SIN(anglev)
    vel(2,ipart) = COS(anglev)
    if (ndimv.eq.3) vel(3,ipart) = 0.
    dens(ipart) = denszero
    pmass(ipart) = massp
    uu(ipart) = uuzero
    if (imhd.ge.1) then
       Bfield(1,ipart) = 0.
       Bfield(2,ipart) = 0.
       if (ndimv.eq.3) Bfield(3,ipart) = 0.0
    elseif (imhd.lt.0) then
!--vector potential setup
       Bevol(:,ipart) = 0.
       if (ndimV.lt.3) stop 'ndimV too small in setup_jadvect'
       rr = sqrt(dot_product(x(:,ipart),x(:,ipart)))
       Bevol(3,ipart) = MAX(Azero*(rzero-rr),0.)
    else
       Bfield(:,ipart) = 0.
    endif 
 enddo

!
!--allow for tracing flow
!
 if (trace) write(iprint,*) '  exiting subroutine setup'
            
 return
end subroutine setup
