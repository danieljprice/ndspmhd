!!------------------------------------------------------------------------!!
!!                                                                        !!
!!  Setup for the MHD Rotor problems in Toth (2000)                       !!
!!                                                                        !!
!!  dense, rotating disk of fluid at origin                               !!
!!  magnetic field is initially uniform in x-direction                    !!
!!  -> simplest way to set this up is to vary the particle mass           !!
!!                                                                        !!                                                                        !!
!!  Note for all MHD setups, only the magnetic field should be setup      !!
!!  Similarly the thermal energy is setup even if using total energy.     !!
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
 real :: denszero,densdisk,przero,vzero,ftaper
 real :: pri,rdisk,rbuffer,radius
 real :: totmass,gam1,massp,const
 real, dimension(ndim) :: xorigin, dx
 real, dimension(ndimv) :: bzero
!
!--check number of dimensions is right
!
 if (ndim.ne.2) stop ' ndim must be = 2 for mhd rotor problem'
 if (ndimv.ne.2) write(iprint,*) ' warning: best if ndimv=2 for this problem'
!
!--set boundaries
!            	    
 ibound = 2	! reflective ghosts (boundaries not important in this problem)
 nbpts = 0	! no fixed particles
 xmin(:) = -0.5	! unit square
 xmax(:) = 0.5
 const = sqrt(4.*pi) 
!
!--setup parameters for the problem
! 
 xorigin(:) = 0.0	! co-ordinates of the centre of the initial blast
 rdisk = 0.1		! radius of the initial disk
 rbuffer = 0.115	! radius of the smoothed front
 vzero = 2.0		! rotation speed of initial disk
 bzero(:) = 0.
 if (imhd.ne.0) bzero(1) = 5.0/const	! uniform field in bx direction
 przero = 1.0		! initial pressure
 denszero = 1.0		! ambient density
 densdisk = 10.0		! density of rotating disk
 
 gam1 = gamma - 1.

 write(iprint,*) 'two dimensional mhd rotor problem '
 write(iprint,10) densdisk,rdisk,bzero(1),vzero,przero
10 format(/,' central density  = ',f10.3,', disk radius = ',f6.3,/, &
            ' initial bx   = ',f6.3,', rotation = ',f6.3,', pressure = ',f6.3,/)
!
!--setup uniform density grid of particles (2d) 
!  (determines particle number and allocates memory)
!
 call set_uniform_cartesian(1,psep,xmin,xmax,.false.)	! 2 = close packed arrangement

 ntotal = npart
!
!--determine particle mass in ambient medium
!
 totmass = denszero*product(xmax(:)-xmin(:))
 massp = totmass/float(ntotal) ! average particle mass
!
!--now assign particle properties
! 
 do ipart=1,ntotal
    dx(:) = x(:,ipart)-xorigin(:) 
    radius = sqrt(dot_product(dx,dx))
    if (radius.le.rdisk) then
       dens(ipart) = densdisk
       pmass(ipart) = massp*densdisk/denszero
       vel(1,ipart) = -vzero*(x(2,ipart)-xorigin(2))/rdisk
       vel(2,ipart) = vzero*(x(1,ipart)-xorigin(1))/rdisk
    elseif (radius.le.rbuffer) then	! smooth edge with taper function (toth)
       ftaper = (rbuffer-radius)/(rbuffer - rdisk)
       dens(ipart) = denszero + (densdisk-denszero)*ftaper
       pmass(ipart) = massp*dens(ipart)/denszero
       vel(1,ipart) = -ftaper*vzero*(x(2,ipart)-xorigin(2))/radius
       vel(2,ipart) = ftaper*vzero*(x(1,ipart)-xorigin(1))/radius
    else
       pmass(ipart) = massp
       dens(ipart) = denszero
       vel(:,ipart) = 0.
    endif  
    pri = przero 
    uu(ipart) = pri/(gam1*denszero)
    bfield(:,ipart) = bzero(:)
 enddo
!
!--allow for tracing flow
!
 if (trace) write(iprint,*) '  exiting subroutine setup'
            
 return
end
