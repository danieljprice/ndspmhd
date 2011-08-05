!!------------------------------------------------------------------------!!
!!                                                                        !!
!!  Setup for the 2D MHD blast waves in Balsara (1998)                    !!
!!                                                                        !!
!!  density is set to unity all over, whilst the pressure is set          !!
!!  to 1000.0 in a small circle around the origin                         !!
!!                                                                        !!
!!  Magnetic field of strength 10G in the x-direction                     !!
!!                                                                        !!                                                                        !!
!!  Note for all MHD setups, only the magnetic field should be setup      !!
!!  Similarly the thermal energy is setup even if using total energy.     !!
!!                                                                        !!
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
 REAL :: rhozero,przero,vzero,Bzero
 REAL :: prblast,pri,rblast,radius
 REAL :: totmass,gam1,massp,const
 REAL, DIMENSION(ndim) :: xblast, dblast
 REAL :: rbuffer, exx
!
!--check number of dimensions is right
!
 IF (ndim.NE.2) STOP ' ndim must be = 2 for MHD blast wave problem'
 IF (ndimV.NE.2) WRITE(iprint,*) ' WARNING: best if ndimV=2 for this problem'
!
!--set boundaries
!            	    
 ibound = 2	! reflective ghosts
 nbpts = 0	! must use fixed particles if inflow/outflow at boundaries
 xmin(1) = -0.5		! x
 xmax(1) = 0.5
 xmin(2) = -0.5		! y
 xmax(2) = 0.5
 const = 4.*pi 
!
!--setup parameters for the problem
! 
 xblast(1) = 0.0	! co-ordinates of the centre of the initial blast
 xblast(2) = 0.0
 rblast = 0.05		! radius of the initial blast
 rbuffer = rblast+20.*psep		! radius of the smoothed front
 vzero = 0.
 IF (imhd.NE.0) THEN
    Bzero = 10.0/SQRT(const)! uniform field in Bx direction
 ELSE
    Bzero = 0.0
 ENDIF
 przero = 1.0		! initial pressure
 rhozero = 1.0
 prblast = 1000.0	! initial pressure within rblast
 
 gam1 = gamma - 1.

 WRITE(iprint,*) 'Two dimensional MHD blast wave problem '
 WRITE(iprint,10) prblast,rblast,Bzero,rhozero,przero
10 FORMAT(/,' Central pressure  = ',f10.3,', blast radius = ',f6.3,/, &
            ' Initial B   = ',f6.3,', density = ',f6.3,', pressure = ',f6.3,/)
!
!--setup uniform density grid of particles (2D) 
!  (determines particle number and allocates memory)
!
 CALL set_uniform_cartesian(1,xmin,xmax,.false.)	! 2 = close packed arrangement

 ntotal = npart
!
!--determine particle mass
!
 totmass = rhozero*(xmax(2)-xmin(2))*(xmax(1)-xmin(1))
 massp = totmass/FLOAT(ntotal) ! average particle mass
!
!--now assign particle properties
! 
 DO ipart=1,ntotal
    velin(1,ipart) = 0.
    velin(2,ipart) = 0.
    IF (ndimV.EQ.3) velin(3,ipart) = 0.
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
    Bin(1,ipart) = Bzero
    Bin(2,ipart) = 0.
    IF (ndimV.EQ.3) Bin(3,ipart) = 0.0	
 ENDDO
!
!--allow for tracing flow
!
 IF (trace) WRITE(iprint,*) '  Exiting subroutine setup'
            
 RETURN
END
