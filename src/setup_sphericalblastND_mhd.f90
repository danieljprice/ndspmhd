!!-------------------------------------------------------------------------!!
!!                                                                         !!
!!  setup for spherical adiabatic (mhd) blast waves                        !!
!!  in 1,2 and 3 dimensions                                                !!
!!                                                                         !!
!!  density is set to unity all over, whilst the pressure (or equivalently !!
!!  the thermal energy) is set to some large quantity in a small circle    !!
!!  around the origin                                                      !!
!!                                                                         !!
!!  magnetic field of strength 10g in the x-direction                      !!
!!                                                                         !!                                                                        !!
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
 use options
 use part
 use setup_params
 
 use uniform_distributions
!
!--define local variables
!      
 implicit none
 integer :: ipart
 real :: denszero,przero
 real :: prblast,pri,uui,rblast,radius,enblast,enzero
 real :: totmass,gam1,massp,const
 real, dimension(ndim) :: xblast, dblast
 real, dimension(ndimV) :: Bzero
 real :: rbuffer, exx, hsmooth
 real :: q2, wab, grkern
!
!--set boundaries
!            	    
 ibound = 3	! fixed ghosts
 nbpts = 0
 xmin(:) = -0.5		! same xmin in all dimensions
 xmax(:) = 0.5
 xmin(2) = -0.75
 xmax(2) = 0.75
 const = 1./sqrt(4.*pi) 
!
!--setup parameters for the problem
! 
 xblast(:) = 0.0	! co-ordinates of the centre of the initial blast
 rblast = 0.0     !!!05    !!!*psep		! radius of the initial blast
 rbuffer = rblast	!+10.*psep		! radius of the smoothed front
 bzero(:) = 0.0
 if (imhd.ne.0) then
    bzero(1) = sqrt(2.*pi)         !10.0*const	! uniform field in bx direction
    bzero(2) = sqrt(2.*pi)
 endif
 przero = 0.1		! initial pressure
 denszero = 1.0
 prblast = 10.0	! initial pressure within rblast
 enblast = 1.0
 enzero = 0.
 
 gam1 = gamma - 1.
 if (abs(gam1).lt.1.e-3) stop 'eos cannot be isothermal for this setup'

 write(iprint,10) ndim
 write(iprint,20) prblast,rblast,denszero,przero
 write(iprint,30) Bzero
10 format(/,1x,i1,'-dimensional adiabatic mhd blast wave problem')
20 format(/,' central pressure  = ',f10.3,', blast radius = ',f6.3,/, &
            ' density = ',f6.3,', pressure = ',f6.3,/)
30 format(' initial B   = ',3(f6.3,1x))
!
!--setup uniform density grid of particles
!  (determines particle number and allocates memory)
!
 call set_uniform_cartesian(1,psep,xmin,xmax,offset=.true.,perturb=0.5)	! 2 = close packed arrangement

 ntotal = npart
!
!--determine particle mass
!
 totmass = denszero*product(xmax(:)-xmin(:))	! assumes cartesian boundaries
 massp = totmass/float(ntotal) ! average particle mass
 !enblast = enblast/massp   ! enblast is now the energy to put in a single particle
 !enblast = enblast/(4./3.*pi*rblast**3)
!
!--smoothing length for kernel smoothing
! 
 hsmooth = hfact*(massp/denszero)**dndim
 write(iprint,*) 'hsmooth = ',hsmooth
!
!--now assign particle properties
! 
 do ipart=1,ntotal
    vel(:,ipart) = 0.
!--uniform density and smoothing length
    dens(ipart) = denszero
    pmass(ipart) = massp

    dblast(:) = x(:,ipart)-xblast(:) 
    radius = sqrt(dot_product(dblast,dblast))
!
!--smooth energy injection using the sph kernel
!    
!    q2 = radius**2/hsmooth**2
!    call interpolate_kernel(q2,wab,grkern)
!    uui = enblast*wab/hsmooth**ndim
    if (radius.lt.rblast) then
       pri = prblast
       !uui = enblast
    elseif (radius.lt.rbuffer) then	! smooth out front
       exx = exp((radius-rblast)/(psep))
       pri = (prblast + przero*exx)/(1.0+exx)
       !uui = (enblast + enzero*exx)/(1.0+exx)
    else
       pri = przero
       !uui = enzero
    endif   
    uu(ipart) = pri/(gam1*denszero)
    bfield(:,ipart) = bzero(:)
 enddo
!
!--set the constant components of the mag field which can be subtracted
!
 bconst(:) = bzero(:)

!
!--allow for tracing flow
!
 if (trace) write(iprint,*) '  exiting subroutine setup'
            
 return
end

subroutine modify_dump
 use loguns, only:iprint
 use part
 use timestep, only:time
 use setup_params, only:pi
 use eos, only:gamma
 implicit none
 integer :: ipart
 real :: rblast,radius,dblast(ndim),przero,prblast,pri,gam1,denszero
!
!--now assign particle properties
!
 write(iprint,*) 'modifying dump with velocities/B field'
!
!--now assign particle properties
! 
 rblast = 0.1
 prblast = 10.0
 przero = 0.1
 gam1 = gamma - 1.0
 denszero = 1.0
 do ipart=1,ntotal
    vel(:,ipart) = 0.

    dblast(:) = x(:,ipart)
    radius = sqrt(dot_product(dblast,dblast))
!
!--smooth energy injection using the sph kernel
!    
!    q2 = radius**2/hsmooth**2
!    call interpolate_kernel(q2,wab,grkern)
!    uui = enblast*wab/hsmooth**ndim
    if (radius.lt.rblast) then
       pri = prblast
       !uui = enblast
    else
       pri = przero
       !uui = enzero
    endif   
    uu(ipart) = pri/(gam1*denszero)
 enddo
 time = 0.
 
end subroutine modify_dump
