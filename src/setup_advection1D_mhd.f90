!----------------------------------------------------------------
!     Set up for an advection test similar to that performed in
!     Stone, Hawley, Evans and Norman (1992) ApJ 388,415-437
!     (see also Evans & Hawley (1988) ApJ 332,659)
!
!     Optional smoothing of discontinuities (not necessary)
!
!     Note for all MHD setups, only the magnetic field should be setup
!----------------------------------------------------------------

SUBROUTINE setup
!
!--include relevant global variables
!
 USE dimen_mhd
 USE bound
 USE loguns
 USE eos
 USE options
 USE part
 USE part_in
 USE setup_params
!
!--define local variables
!
 IMPLICIT NONE
 INTEGER :: i,j,ntot
 REAL :: dist,xcentre,vxi,term,term2,rhozero,massp,gam1
 REAL :: xpulseleft,xpulseright,deltaleft,deltaright
 REAL :: exx,uuzero,Bypulse,velpulse,rhopulse,przero
 REAL :: ampl,spsoundi,velzero,dsmooth
 REAL, DIMENSION(ndimB) :: Bzero

 IF (ndim.NE.1) STOP ' advection test only setup for 1D : recompile '
 
 ibound = 3	! periodic
! nbpts = 10
 xmin(1) = -0.5
 xmax(1) = 0.5
 xcentre  = 0.0
 dsmooth = 0.0		! amount by which to smooth discontinuities ( < 20 )
 xpulseleft = xcentre - 25*psep	! 50 particles wide
 xpulseright = xcentre + 25*psep
 gam1 = gamma - 1.
 dist  = 0.5*psep
 rhozero =  10.0   
 massp = rhozero*psep                 ! the mass per sph particle
 npart = 0
 uuzero = 1.0/(gamma*gam1)
 ampl = 0.01		! amplitude of the density perturbation
 IF (imhd.NE.0) THEN
    Bypulse = 0.001
 ELSE
    WRITE(iprint,*) 'MHD not set, this test does nothing'
!   CALL quit
    Bypulse = 0.
 ENDIF       
 Bzero(1) = 0.
 Bzero(2) = Bypulse
 Bzero(3) = 0.

 velpulse = 5.0
 velzero = velpulse
 rhopulse = rhozero*(1.+ampl)
!
!--allocate variables
!
 ntot = NINT((xmax(1)-xmin(1))/psep) + 1
 CALL alloc(ntot,1)

! i = 1
 xin(1,1) = xmin(1) + 0.5*psep
 xin(1,2) = xmin(1) + psep + 0.5*psep
 rhoin(1:2) = rhozero
 uuin(1:2) = uuzero
 pmass(1:2) = massp
 velin(:,1:2) = 0.
 velin(1,1:2) = velzero
 
 Bin(:,1:2) = 0.
 hhin(1:2) = hfact*massp/rhozero

 i = 1
 DO WHILE (xin(1,i).lt.xmax(1))  
    i = i + 1
!   npart = npart + 1
    xin(1,i) = xmin(1) + (i-1)*psep + 0.5*psep
    deltaleft = (xin(1,i) - xpulseleft)/dist
    deltaright = (xin(1,i) - xpulseright)/dist
	
    Bin(:,i) = 0.		! initially zero, reset later
    velin(:,i) = 0.	 
    velin(1,i) = velzero
    rhoin(i) = rhozero

    IF ((deltaleft.gt.dsmooth).AND.(deltaright.le.-dsmooth)) THEN
       Bin(2,i) = Bypulse
       velin(1,i) = velpulse
!      rhoin(i) = rhopulse
    ELSEIF (deltaleft.le.-dsmooth) THEN
       Bin(2,i) = 0.
!      velin(1,i) = velzero
       rhoin(i) = rhozero
    ELSEIF ((deltaleft.gt.-dsmooth).AND.(deltaleft.le.dsmooth)) THEN
       exx = exp(deltaleft)
       Bin(2,i) = (Bypulse*exx)/(1 + exx)
!      velin(1,i) = (velpulse*exx)/(1+exx)
!      rhoin(i) = (rhozero + rhopulse*exx)/(1+exx)
    ELSEIF ((deltaright.gt.-dsmooth).AND.(deltaright.le.dsmooth)) THEN
       exx = exp(deltaright)
       Bin(2,i) = Bypulse - (Bypulse*exx)/(1 + exx)
!      velin(1,i) = velpulse - (velpulse*exx)/(1+exx)
!      rhoin(i) = rhopulse - ((rhopulse-rhozero)*exx)/(1+exx)
    ELSE
       Bin(2,i) = 0.  
       rhoin(i) = rhozero
!      velin(1,i) = velzero
    ENDIF

!   xin(1,i) = xin(1,i-2) + 2.*massp/rhoin(i-1)	 
    hhin(i) = hfact*massp/rhoin(i)
    uuin(i) = uuzero
    pmass(i) = massp
 ENDDO
 
 npart = i-1
 ntotal = npart
    
 RETURN
END
