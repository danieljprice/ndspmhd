!!------------------------------------------------------------------------!!
!!                                                                        !!
!!  Setup for spherical adiabatic (MHD) blast waves                       !!
!!  in 1,2 and 3 dimensions                                               !!
!!                                                                        !!
!!  density is set to unity all over, whilst the pressure is set          !!
!!  to 1000.0 in a small circle around the origin                         !!
!!                                                                        !!
!!  Magnetic field of strength 10G in the x-direction                     !!
!!                                                                        !!                                                                        !!
!!------------------------------------------------------------------------!!

SUBROUTINE setup
!
!--include relevant global variables
!
 USE dimen_mhd
 USE debug
 USE loguns
 USE bound
 USE eos
 USE options
 USE part
 USE part_in
 USE setup_params
!
!--define local variables
!      
 IMPLICIT NONE
 INTEGER :: i,j,ntot,npartx,nparty,ipart
 REAL :: rhozero,przero,vzero
 REAL :: prblast,pri,rblast,radius
 REAL :: totmass,gam1,massp,const
 REAL, DIMENSION(ndim) :: xblast, dblast
 REAL, DIMENSION(ndimB) :: Bzero
 REAL :: rbuffer, exx
!
!--set boundaries
!            	    
 ibound = 2	! reflective ghosts
 nbpts = 0	! must use fixed particles if inflow/outflow at boundaries
 xmin(:) = -0.5		! same xmin in all dimensions
 xmax(:) = 0.5
 const = 1./SQRT(4.*pi) 
!
!--setup parameters for the problem
! 
 xblast(:) = 0.0	! co-ordinates of the centre of the initial blast
 rblast = 0.05		! radius of the initial blast
 rbuffer = rblast+20.*psep		! radius of the smoothed front
 Bzero(:) = 0.0
 IF (imhd.NE.0) THEN
    Bzero(1) = 10.0*const	! uniform field in Bx direction
 ENDIF
 przero = 1.0		! initial pressure
 rhozero = 1.0
 prblast = 1000.0	! initial pressure within rblast
 
 gam1 = gamma - 1.
 IF (abs(gam1).lt.1.e-3) STOP 'eos cannot be isothermal for this setup'

 WRITE(iprint,*) ndim,' dimensional adiabatic MHD blast wave problem '
 WRITE(iprint,10) prblast,rblast,rhozero,przero
 WRITE(iprint,20) Bzero
10 FORMAT(/,' Central pressure  = ',f10.3,', blast radius = ',f6.3,/, &
            ,' density = ',f6.3,', pressure = ',f6.3,/)
20 FORMAT(' Initial B   = ',3(f6.3,1x))
!
!--setup uniform density grid of particles
!  (determines particle number and allocates memory)
!
 CALL set_uniform_cartesian(1,xmin,xmax,.false.)	! 2 = close packed arrangement

 ntotal = npart
!
!--determine particle mass
!
 totmass = rhozero*PRODUCT(xmax(:)-xmin(:))	! assumes cartesian boundaries
 massp = totmass/FLOAT(ntotal) ! average particle mass
!
!--now assign particle properties
! 
 DO ipart=1,ntotal
    velin(:,ipart) = 0.
    rhoin(ipart) = rhozero
    pmass(ipart) = massp
    dblast(:) = xin(:,ipart)-xblast(:) 
    radius = SQRT(DOT_PRODUCT(dblast,dblast))
    IF (radius.LT.rblast) THEN
       pri = prblast
    ELSEIF (radius.LT.rbuffer) THEN	! smooth out front
       exx = exp((radius-rblast)/(psep))
       pri = (prblast + przero*exx)/(1.0+exx)
    ELSE
       pri = przero
    ENDIF   
    uuin(ipart) = pri/(gam1*rhozero)
    hhin(ipart) = hfact*(massp/rhoin(ipart))**hpower	 ! ie constant everywhere
    Bin(:,ipart) = Bzero(:)
 ENDDO
!
!--allow for tracing flow
!
 IF (trace) WRITE(iprint,*) '  Exiting subroutine setup'
            
 RETURN
END
