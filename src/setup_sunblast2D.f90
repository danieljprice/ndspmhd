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
 use options, only:damp,ibound
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
 ibound = 1	! fixed ghosts
 nbpts = 0
 xmin(:) = -0.1	
 xmax(:) = 0.1
! xmax(1) = 0.1
!! xmin(2) = -0.3
!! xmax(2) = 0.3	! same xmin in all dimensions
!! xmax(:) = 0.5
! xmin(2) = -0.75
! xmax(2) = 0.75
 denszero = 1.0
 przero = 0.1
 
 gam1 = gamma - 1.
 if (abs(gam1).lt.1.e-3) stop 'eos cannot be isothermal for this setup'

!
!--setup uniform density grid of particles
!  (determines particle number and allocates memory)
!
 call set_uniform_cartesian(1,psep,xmin,xmax,offset=.true.)	! 2 = close packed arrangement

 ntotal = npart
!
!--determine particle mass
!
 totmass = denszero*product(xmax(:)-xmin(:))	! assumes cartesian boundaries
 massp = totmass/float(ntotal) ! average particle mass
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
 use bound, only:xmin,xmax
 implicit none
 integer :: ipart
 real :: prblast,pri,uui,radius,enblast,enzero,const
 real :: rblast,przero,gam1,denszero
 real, dimension(ndim) :: xblast, dblast
 real, dimension(ndimV) :: Bzero
 real :: rbuffer, exx, hsmooth
 real :: q2, wab, grkern
 real :: Azdipole,Azradial,Bdipole,Brad
 real, parameter :: Rcen = 1.7, Rdipole = 0.3, eps = 0.3
 
 write(iprint,*) 'modifying dump by adding blast'

 gam1 = gamma - 1.0
 if (abs(gam1).lt.1.e-3) stop 'eos cannot be isothermal for this setup'
!
!--setup parameters for the problem
! 
 xblast(:) = 0.0	! co-ordinates of the centre of the initial blast
 rblast = 2.*hfact*psep		! radius of the initial blast
 rbuffer = rblast	!+10.*psep		! radius of the smoothed front
 bzero(:) = 0.0
 const = 1./sqrt(4.*pi) 
 if (imhd.ne.0) then
    bzero(1) = sqrt(2.*pi)         !10.0*const	! uniform field in bx direction
    bzero(2) = sqrt(2.*pi)
 endif
 przero = 0.0		! initial pressure
 denszero = 1.0
 prblast = 10.0	! initial pressure within rblast
 enblast = 1.0
 enzero = 0.
 prblast = gam1*enblast/(4./3.*pi*rblast**3)
 !enblast = enblast/massp   ! enblast is now the energy to put in a single particle
!
!--smoothing length for kernel smoothing
! 
 hsmooth = hfact*(pmass(1)/denszero)**dndim
! write(iprint,*) 'hsmooth = ',hsmooth

!--initial magnetic field strength
 Brad = 1.0
 Bdipole = 200. !!*Brad

 write(iprint,10) ndim
 write(iprint,20) prblast,rblast,denszero,przero
10 format(/,1x,i1,'-dimensional adiabatic mhd blast wave problem')
20 format(/,' central pressure  = ',f10.3,', blast radius = ',f6.3,/, &
            ' density = ',f6.3,', pressure = ',f6.3,/)

 nbpts = 0

 do ipart=1,ntotal
    vel(:,ipart) = 0.
!--uniform density and smoothing length
    dens(ipart) = denszero

    dblast(:) = x(:,ipart)-xblast(:) 
    radius = sqrt(dot_product(dblast,dblast))
!
!--smooth energy injection using the sph kernel
!    
!    q2 = radius**2/hsmooth**2
!    call interpolate_kernel(q2,wab,grkern)
!    uui = enblast*wab/hsmooth**ndim
    if (radius.le.rblast) then
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
    
    if (imhd.lt.0) then
       !!Azradial = Brad*(-(x(1,ipart)+Rcen)*x(2,ipart) + x(1,ipart)) !!sqrt((Rcen + x(1,ipart))**2 + x(2,ipart)**2)
       Azdipole = -Bdipole*Rdipole**3*(Rdipole + x(1,ipart))/ &
                   sqrt((Rdipole + x(1,ipart))**2 + x(2,ipart)**2 + (eps*Rdipole)**2)**3
       Bevol(3,ipart) = Azdipole
!       Bconst(1) = Brad
    elseif (imhd.gt.0) then
       stop 'not implemented'
       bfield(:,ipart) = bzero(:)
    endif
    if (ANY(abs(x(:,ipart)).GT.0.1-2.*hfact*psep)) then
       itype(ipart) = 1
       nbpts = nbpts + 1
    endif
 enddo
 
!
!--set the constant components of the mag field which can be subtracted
!
 bconst(1) = Brad
 
 time = 0.
 
end subroutine modify_dump
