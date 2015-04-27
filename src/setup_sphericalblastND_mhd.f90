!------------------------------------------------------------------------------!
! NDSPMHD: A Smoothed Particle (Magneto)Hydrodynamics code for (astrophysical) !
! fluid dynamics simulations in 1, 2 and 3 spatial dimensions.                 !
!                                                                              !
! (c) 2002-2014 Daniel Price                                                   !
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

!!-------------------------------------------------------------------------!!
!!                                                                         !!
!!  setup for spherical adiabatic (mhd) blast waves                        !!
!!  in 1,2 and 3 dimensions                                                !!
!!                                                                         !!
!!  density is set to unity all over, whilst the pressure (or equivalently !!
!!  the thermal energy) is set to some large quantity in a small circle    !!
!!  around the origin                                                      !!
!!                                                                         !!
!!-------------------------------------------------------------------------!!

subroutine setup
!
!--include relevant global variables
!
 use dimen_mhd, only:ndim,ndimV
 use debug, only:trace
 use loguns, only:iprint
 use bound
 use eos
 use kernels, only:interpolate_kernel
 use options, only:damp,ibound,geom
 use part
 use setup_params, only:psep
 
 use uniform_distributions
!
!--define local variables
!      
 implicit none
 integer :: ipart
 real :: denszero,przero
 real :: totmass,gam1,massp

!
!--set boundaries
!
 nbpts = 0
 ibound = 3     ! fixed ghosts
 xmin(:) = -0.5     ! same xmin in all dimensions
 xmax(:) = 0.5
 denszero = 1.0
 przero = 0.1
 
 gam1 = gamma - 1.
 if (abs(gam1).lt.1.e-3) stop 'eos cannot be isothermal for this setup'

!
!--setup uniform density grid of particles
!  (determines particle number and allocates memory)
!
 call set_uniform_cartesian(1,psep,xmin,xmax,offset=.true.) ! 2 = close packed arrangement

 ntotal = npart
!
!--determine particle mass
!
 totmass = denszero*product(xmax(:)-xmin(:))    ! assumes cartesian boundaries
 massp = totmass/float(ntotal)                  ! average particle mass
!
!--now assign particle properties
! 
 do ipart=1,ntotal
    vel(:,ipart) = 0.
!--uniform density and smoothing length
    dens(ipart) = denszero
    pmass(ipart) = massp
    uu(ipart) = przero/(gam1*denszero)
 enddo
!
!--add the blast wave if damping is off
!
 if (abs(damp).lt.tiny(damp)) call modify_dump

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
 use setup_params, only:pi,psep,hfact
 use eos, only:gamma
 use cons2prim, only:specialrelativity
 use kernels, only:interpolate_kernel
 implicit none
 integer :: ipart
 real :: prblast,pri,enblast,enzero,const !,uui
 real :: rblast,radius,przero,gam1,denszero
 real, dimension(ndim) :: xblast, dblast
 real, dimension(ndimV) :: Bzero
 real :: rbuffer, exx, hsmooth
 real :: q2, wab, grkern, uui
 logical, parameter :: dosedov = .true.
 
 write(iprint,*) 'modifying dump by adding blast'

 gam1 = gamma - 1.0
 if (abs(gam1).lt.1.e-3) stop 'eos cannot be isothermal for this setup'
!
!--setup parameters for the problem
! 
 xblast(:) = 0.0        ! co-ordinates of the centre of the initial blast
 bzero(:) = 0.0
 const = 1./sqrt(4.*pi) 
 if (imhd.ne.0) then
    bzero(1) = 3.0 !!sqrt(2.*pi)         !10.0*const	! uniform field in bx direction
!    bzero(2) = sqrt(2.*pi)
 endif
 if (specialrelativity) then
    rblast = 0.04
    przero = 1.0
    denszero = 1.
    prblast = 1000.
    enblast = 1.0
    prblast = gam1*enblast/(4./3.*pi*rblast**3)
    write(iprint,*) 'using setup for special relativistic problem'
 elseif (dosedov) then
    rblast = 2.*hfact*psep      ! radius of the initial blast
    przero = 0.0                ! initial pressure
    enblast = 1.0
    enzero = 0.
    prblast = gam1*enblast/(4./3.*pi*rblast**3)
    !enblast = enblast/massp   ! enblast is now the energy to put in a single particle
 else
    rblast = 0.1                ! radius of the initial blast
    przero = 0.1                ! initial pressure
    prblast = 10.0              ! initial pressure within rblast
 endif
 rbuffer = rblast       !+10.*psep      ! radius of the smoothed front
 denszero = 1.0
!
!--smoothing length for kernel smoothing
! 
 hsmooth = 2.*hfact*(pmass(1)/denszero)**dndim
! write(iprint,*) 'hsmooth = ',hsmooth

 write(iprint,10) ndim
 write(iprint,20) prblast,rblast,denszero,przero
10 format(/,1x,i1,'-dimensional adiabatic mhd blast wave problem')
20 format(/,' central pressure  = ',f10.3,', blast radius = ',f6.3,/, &
            ' density = ',f6.3,', pressure = ',f6.3,/)


 do ipart=1,ntotal
    vel(:,ipart) = 0.
!--uniform density and smoothing length
    dens(ipart) = denszero

    dblast(:) = x(:,ipart)-xblast(:) 
    radius = sqrt(dot_product(dblast,dblast))
!
!--smooth energy injection using the sph kernel
!    
    q2 = radius**2/hsmooth**2
    call interpolate_kernel(q2,wab,grkern)
    uui = enblast*wab/hsmooth**ndim
    if (radius.le.rblast) then
       pri = prblast
       !uui = enblast
    elseif (radius.lt.rbuffer) then ! smooth out front
       exx = exp((radius-rblast)/(psep))
       pri = (prblast + przero*exx)/(1.0+exx)
      !uui = (enblast + enzero*exx)/(1.0+exx)
    else
       pri = przero
       !uui = enzero
    endif   
    uu(ipart) = uui
!    uu(ipart) = pri/(gam1*denszero)
!
!--euler potentials setup (for use in other codes)
!    
    if (imhd.eq.-2 .and. ndim.eq.3) then
       if (abs(Bzero(1)).gt.tiny(Bzero)) then ! Bx
          Bevol(1,ipart) = -Bzero(1)*x(3,ipart)
          Bevol(2,ipart) = x(2,ipart)
       elseif (abs(Bzero(2)).gt.tiny(Bzero)) then ! By
          Bevol(1,ipart) = -Bzero(2)*x(1,ipart)
          Bevol(2,ipart) = x(3,ipart)
       elseif (abs(Bzero(3)).gt.tiny(Bzero)) then ! Bz
          Bevol(1,ipart) = -Bzero(3)*x(2,ipart)
          Bevol(2,ipart) = x(1,ipart)
       else
          Bevol(:,ipart) = 0.
       endif
    endif
    Bfield(:,ipart) = bzero(:)
 enddo
 
!
!--set the constant components of the mag field which can be subtracted
!
 bconst(:) = bzero(:)
 
 time = 0.
 
end subroutine modify_dump
