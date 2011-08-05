!!-------------------------------------------------------------------------!!
!!                                                                         !!
!!  Setup for spherical adiabatic (MHD) blast waves                        !!
!!  in 1,2 and 3 dimensions                                                !!
!!                                                                         !!
!!  density is set to unity all over, whilst the pressure (or equivalently !!
!!  the thermal energy) is set to some large quantity in a small circle    !!
!!  around the origin                                                      !!
!!                                                                         !!
!!  Magnetic field of strength 10G in the x-direction                      !!
!!                                                                         !!                                                                        !!
!!-------------------------------------------------------------------------!!

SUBROUTINE setup
!
!--include relevant global variables
!
 USE dimen_mhd
 USE debug
 USE loguns
 USE bound
 USE eos
 USE kernel
 USE options
 USE part
 USE setup_params
!
!--define local variables
!      
 IMPLICIT NONE
 INTEGER :: i,j,ntot,npartx,nparty,ipart
 REAL :: denszero,przero,vzero
 REAL :: prblast,pri,uui,rblast,radius,enblast,enzero
 REAL :: totmass,gam1,massp,const
 REAL, DIMENSION(ndim) :: xblast, dblast
 REAL, DIMENSION(ndimB) :: Bzero
 REAL :: rbuffer, exx, hsmooth
 REAL :: q2, wab, grkern
!
!--set boundaries
!            	    
 ibound = 2	! fixed ghosts
 nbpts = 0
 xmin(:) = -0.5		! same xmin in all dimensions
 xmax(:) = 0.5
 const = 1./SQRT(4.*pi) 
!
!--setup parameters for the problem
! 
 xblast(:) = 0.0	! co-ordinates of the centre of the initial blast
 rblast = 0.5*psep		! radius of the initial blast
 rbuffer = rblast	!+10.*psep		! radius of the smoothed front
 Bzero(:) = 0.0
 IF (imhd.NE.0) THEN
    Bzero(1) = 10.0*const	! uniform field in Bx direction
 ENDIF
 przero = 1.0		! initial pressure
 denszero = 1.0
 prblast = 1000.0	! initial pressure within rblast
 enblast = 1.0
 enzero = 0.
 
 gam1 = gamma - 1.
 IF (abs(gam1).lt.1.e-3) STOP 'eos cannot be isothermal for this setup'

 WRITE(iprint,10) ndim
 WRITE(iprint,20) prblast,rblast,denszero,przero
 WRITE(iprint,30) Bzero
10 FORMAT(/,1x,i1,'-dimensional adiabatic MHD blast wave problem')
20 FORMAT(/,' Central pressure  = ',f10.3,', blast radius = ',f6.3,/, &
            ' density = ',f6.3,', pressure = ',f6.3,/)
30 FORMAT(' Initial B   = ',3(f6.3,1x))
!
!--setup uniform density grid of particles
!  (determines particle number and allocates memory)
!
 CALL set_uniform_cartesian(1,psep,xmin,xmax,.true.)	! 2 = close packed arrangement

 ntotal = npart
!
!--determine particle mass
!
 totmass = denszero*PRODUCT(xmax(:)-xmin(:))	! assumes cartesian boundaries
 massp = totmass/FLOAT(ntotal) ! average particle mass
! enblast = enblast/massp   ! enblast is now the energy to put in a single particle
!
!--smoothing length for kernel smoothing
! 
 hsmooth = hfact*(massp/denszero)**dndim
!
!--now assign particle properties
! 
 DO ipart=1,ntotal
    vel(:,ipart) = 0.
!--uniform density and smoothing length
    dens(ipart) = denszero
    pmass(ipart) = massp

    dblast(:) = x(:,ipart)-xblast(:) 
    radius = SQRT(DOT_PRODUCT(dblast,dblast))
!
!--smooth energy injection using the SPH kernel
!    
    q2 = radius**2/hsmooth**2
    CALL interpolate_kernel(q2,wab,grkern)
    uui = enblast*wab/hsmooth**ndim
!    IF (radius.LT.rblast) THEN
!       pri = prblast
!       uui = enblast
!    ELSEIF (radius.LT.rbuffer) THEN	! smooth out front
!       exx = exp((radius-rblast)/(psep))
!       pri = (prblast + przero*exx)/(1.0+exx)
!       uui = (enblast + enzero*exx)/(1.0+exx)
!    ELSE
!       pri = przero
!       uui = enzero
!    ENDIF   
    uu(ipart) = uui	!pri/(gam1*denszero)
    Bfield(:,ipart) = Bzero(:)
 ENDDO
!
!--allow for tracing flow
!
 IF (trace) WRITE(iprint,*) '  Exiting subroutine setup'
            
 RETURN
END
