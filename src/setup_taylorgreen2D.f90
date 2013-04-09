!!------------------------------------------------------------------------!!
!!                                                                        !!
!!  Generic setup for the 2D Taylor-Green vortex test                     !!
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
 integer :: ipart,i
 real :: betazero,denszero,przero,vzero,bzero,uuzero
 real :: totmass,gam1,massp
 real, dimension(ndim) :: xminregion, xmaxregion
 logical, parameter :: symmetric = .true.
!
!--check number of dimensions is right
!
 if (ndim.lt.2) stop ' ndim must be >= 2 for Taylor-Green vortex'
!
!--set boundaries
!                        
 ibound = 3     ! periodic
 nbpts = 0      ! must use fixed particles if inflow/outflow at boundaries
 if (symmetric) then
    xmin(1) = -0.5  ! x
    xmax(1) = -xmin(1)
    xmin(2) = -0.5  ! y
    xmax(2) = -xmin(2)
    if (ndim.ge.3) then
       xmin(3) = -0.0625
       xmax(3) = -xmin(3)
    endif
 else
    xmin(1) = 0.  ! x
    xmax(1) = 1.0
    xmin(2) = 0.0  ! y
    xmax(2) = 1.0
    if (ndim.ge.3) then
       xmin(3) = 0.
       xmax(3) = 0.125
    endif
 endif
!
!--setup parameters
!
 denszero = 1.
 przero = denszero/gamma
 vzero = 0.1
 bzero = 1.0/sqrt(4.*pi)
 if (gamma.gt.1.) then
    gam1 = gamma - 1.
    uuzero = przero/(gam1*denszero)
 else
    uuzero = 1.5*przero/denszero
 endif

 write(iprint,*) 'Two dimensional Taylor-Green vortex problem '
 if (ndim.ge.3) write(iprint,*) ' (in 3D...)'
 write(iprint,10) betazero,vzero,bzero,denszero,przero
10 format(/,' beta        = ',f6.3,', mach number = ',f6.3,/, &
            ' initial B   = ',f6.3,', density = ',f6.3,', pressure = ',f6.3,/)

 if (symmetric) then
    xminregion(:) = xmin(:)
    xmaxregion(:) = xmax(:)
    xmaxregion(1) = 0.5*(xmax(1)+xmin(1))
    xmaxregion(2) = 0.5*(xmax(2)+xmin(2))
    call set_uniform_cartesian(1,psep,xminregion,xmaxregion,fill=.true.)

    call modify_dump()
!
!--reallocate memory to new size of list
!
    call alloc(4*npart)
!
!--reflect particles across x=0 axis
!
    print*,'reflecting across x'
    ipart = npart
    do i=1,npart
       ipart = ipart + 1
       x(1,ipart) = -x(1,i)
       x(2,ipart) = x(2,i)
       if (ndim.eq.3) x(3,ipart) = x(3,i)
       vel(1,ipart) = vel(1,i)
       vel(2,ipart) = -vel(2,i)
       if (imhd.lt.0) then
          Bevol(3,ipart) = Bevol(3,i)
       endif
    enddo
    npart = ipart
!
!--reflect particles across y=0 axis
!
    print*,'reflecting across y'
    ipart = npart
    do i=1,npart
       ipart = ipart + 1
       x(1,ipart) = x(1,i)
       x(2,ipart) = -x(2,i)
       if (ndim.eq.3) x(3,ipart) = x(3,i)
       vel(1,ipart) = -vel(1,i)
       vel(2,ipart) = vel(2,i)
       if (imhd.lt.0) then
          Bevol(3,ipart) = Bevol(3,i)
       endif
    enddo
    npart = ipart
    ntotal = npart

 else
!
!--setup uniform density grid of particles (2D) with sinusoidal field/velocity
!  determines particle number and allocates memory
!
    call set_uniform_cartesian(2,psep,xmin,xmax)        ! 2 = close packed arrangement

 endif

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
    dens(ipart) = denszero
    pmass(ipart) = massp
    uu(ipart) = uuzero
 enddo
 call modify_dump()

!
!--allow for tracing flow
!
 if (trace) write(iprint,*) '  exiting subroutine setup'
            
 return
end subroutine setup

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
 write(iprint,*) 'modifying dump with Taylor-Green velocities/B field'
 vzero = 0.1
 const = 4.*pi
 bzero = 1.0/sqrt(const)
 do ipart=1,ntotal
    vel(1,ipart) = vzero*sin(2.*pi*x(1,ipart))*cos(2.*pi*x(2,ipart))
    vel(2,ipart) = -vzero*cos(2.*pi*x(1,ipart))*sin(2.*pi*x(2,ipart))
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
