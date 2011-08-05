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
 use dimen_mhd
 use debug
 use loguns
 use bound
 use eos
 use kernel
 use options
 use part
 use setup_params
 
 use uniform_distributions
!
!--define local variables
!      
 implicit none
 integer :: i,j,ntot,npartx,nparty,ipart
 real :: denszero,przero,vzero
 real :: prblast,pri,uui,rblast,radius,enblast,enzero
 real :: totmass,gam1,massp,const
 real, dimension(ndim) :: xblast, dblast
 real, dimension(ndimb) :: bzero
 real :: rbuffer, exx, hsmooth
 real :: q2, wab, grkern
!
!--set boundaries
!            	    
 ibound = 3	! fixed ghosts
 nbpts = 0
 xmin(1) = -0.5		! same xmin in all dimensions
 xmax(1) = 0.5
 xmin(2) = -0.75
 xmax(2) = 0.75
 const = 1./sqrt(4.*pi) 
!
!--setup parameters for the problem
! 
 xblast(:) = 0.0	! co-ordinates of the centre of the initial blast
 rblast = 0.1     !!!05    !!!*psep		! radius of the initial blast
 rbuffer = rblast	!+10.*psep		! radius of the smoothed front
 bzero(:) = 0.0
 if (imhd.ne.0) then
    bzero(1) = sqrt(2.*pi)         !10.0*const	! uniform field in bx direction
    bzero(2) = sqrt(2.*pi)
 endif
 przero = 0.1		! initial pressure
 denszero = 1.0
 prblast = 100.0	! initial pressure within rblast
 enblast = 1.0
 enzero = 0.
 
 gam1 = gamma - 1.
 if (abs(gam1).lt.1.e-3) stop 'eos cannot be isothermal for this setup'

 write(iprint,10) ndim
 write(iprint,20) prblast,rblast,denszero,przero
 write(iprint,30) bzero
10 format(/,1x,i1,'-dimensional adiabatic mhd blast wave problem')
20 format(/,' central pressure  = ',f10.3,', blast radius = ',f6.3,/, &
            ' density = ',f6.3,', pressure = ',f6.3,/)
30 format(' initial b   = ',3(f6.3,1x))
!
!--setup uniform density grid of particles
!  (determines particle number and allocates memory)
!
 call set_uniform_cartesian(1,psep,xmin,xmax,.true.)	! 2 = close packed arrangement

 ntotal = npart
!
!--determine particle mass
!
 totmass = denszero*product(xmax(:)-xmin(:))	! assumes cartesian boundaries
 massp = totmass/float(ntotal) ! average particle mass
! enblast = enblast/massp   ! enblast is now the energy to put in a single particle
!
!--smoothing length for kernel smoothing
! 
 hsmooth = hfact*(massp/denszero)**dndim
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
