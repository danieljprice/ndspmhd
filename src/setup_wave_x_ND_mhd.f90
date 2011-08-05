!----------------------------------------------------------------
!     Set up a travelling soundwave in the x direction
!     perturbs particles from a uniform density grid
!     should work in 1, 2 and 3 dimensions
!----------------------------------------------------------------

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
 INTEGER :: i
 INTEGER, PARAMETER :: itsmax = 100
! REAL, PARAMETER :: pi = 3.1415926536
 REAL, PARAMETER :: tol = 1.e-8
 INTEGER :: imax,its,iwave
 REAL, DIMENSION(ndimV) :: Bzero
 REAL :: xcentre,massp
 REAL :: ampl,wk,xlambda,dxmax,denom
 REAL :: dxi,dxprev,xmassfrac,func,fderiv
 REAL :: pri,spsoundi,valfven2i,vamplx,vamply,vamplz
 REAL :: vfast,vslow,vcrap,vwave,term,rhoin1
 REAL :: rhozero,uuzero,przero
!
!--allow for tracing flow
!
 IF (trace) WRITE(iprint,*) ' Entering subroutine setup'
 WRITE(iprint,*) '------------ 1D Wave setup ----------------'

 iwave = 1 		! preset wave parameters
!
!--setup parameters (could read in from a file)
!
 IF (iwave.EQ.1) THEN	! MHD slow wave
    ampl = 0.006
    rhozero = 1.0
    uuzero = 4.5
    Bzero(1:3) = SQRT(2.)
    WRITE(iprint,*) ' setting up for MHD slow wave...'
 ELSE			! MHD fast wave
    ampl = 0.0055
    rhozero = 1.0
    uuzero = 0.3
    Bzero(1:3) = 0.5  
    WRITE(iprint,*) ' setting up for MHD fast wave...'
 ENDIF
 xlambda = 1.0    
 wk = 2.0*pi/xlambda	! 	wave number
!
!--set boundaries
!
 ibound = 3	! periodic boundaries
 nbpts = 0		! use ghosts not fixed
 xmin(:) = 0.   ! set position of boundaries
 xmax(1) = 1.0 
 IF (ndim.GE.2) THEN
    xmax(2:ndim) = 0.25 ! would need to adjust this depending on grid setup
 ENDIF
!
!--initially set up a uniform density grid (also determines npart)
!  (the call to set_uniform_cartesian means this works in 1,2 and 3D)
!
 CALL set_uniform_cartesian(2,xmin,xmax,.false.)
 
! npart = INT((xmax(1)-xmin(1))/psep)
 massp = 1.0/FLOAT(npart)	! average particle mass
!
!--allocate memory here
!
! CALL alloc(npart,1)
!
!--setup uniform density grid of particles
! 
 DO i=1,npart
!    xin(1,i) = xmin(1) + (i-1)*psep  + 0.5*psep 
    velin(:,i) = 0.
    rhoin(i) = rhozero
    pmass(i) = massp
    uuin(i) = uuzero
    hhin(i) = hfact*(pmass(i)/rhoin(i))**hpower	 ! ie constant everywhere
    IF (imhd.GE.1) THEN 
       Bin(:,i) = Bzero
    ENDIF 
 ENDDO

 ntotal = npart
!
!--get sound speed from equation of state (want average sound speed, so
!  before the density is perturbed)
!
 CALL equation_of_state(przero,spsoundi,uuzero,rhozero,gamma,1)
!
!--work out MHD wave speeds
!
 rhoin1 = 1./rhozero
 
 valfven2i = DOT_PRODUCT(Bzero,Bzero)*rhoin1
 vfast = SQRT(0.5*(spsoundi**2 + valfven2i			&
                 + SQRT((spsoundi**2 + valfven2i)**2	&
                 - 4.*(spsoundi*Bzero(1))**2*rhoin1)))
 vslow = SQRT(0.5*(spsoundi**2 + valfven2i			&
                 - SQRT((spsoundi**2 + valfven2i)**2	&
                 - 4.*(spsoundi*Bzero(1))**2*rhoin1)))
 vcrap = SQRT(spsoundi**2 + valfven2i)

!----------------------------
! set the wave speed to use
!----------------------------
		    
 IF (iwave.EQ.1) THEN
    vwave = vslow
 ELSE
    vwave = vfast
 ENDIF
  
 dxmax = xmax(1) - xmin(1)
 denom = dxmax - ampl/wk*(COS(wk*dxmax)-1.0)
 WRITE (iprint,*) 'Wave number = ',wk
 WRITE (iprint,*) 'Amplitude = ',ampl
!
!--now perturb the particles appropriately
!
 DO i=1,npart
    
    dxi = xin(1,i)-xmin(1)
    dxprev = dxmax*2.
    xmassfrac = dxi/dxmax	! current mass fraction(for uniform density)
				! determines where particle should be
!
!--Use rootfinder on the integrated density perturbation
!  to find the new position of the particle
!    
    its = 0
	 
    DO WHILE ((ABS(dxi-dxprev).GT.tol).AND.(its.LT.itsmax))
       dxprev = dxi
       func = xmassfrac*denom - (dxi - ampl/wk*(COS(wk*dxi)-1.0))
       fderiv = -1.0 - ampl*SIN(wk*dxi)
       dxi = dxi - func/fderiv	! Newton-Raphson iteration
       its = its + 1 
!      PRINT*,'iteration',its,'dxi =',dxi 
    ENDDO
	 	 	 
    IF (its.GE.itsmax) THEN
       WRITE(iprint,*) 'Error: soundwave - too many iterations'
       CALL quit 
    ENDIF
	  
!    IF (idebug(1:5).EQ.'sound') THEN
!       WRITE(*,99002) i,its,dxi-dxprev,xin(1,i),xmin(1)+dxi
!99002  FORMAT('Particle',i5,' converged in ',i4,	&
!       ' iterations, error in x =',1(1pe10.2,1x),/,	&
!       'previous x = ',1(0pf8.5,1x),			&
!       'moved to x = ',1(0pf8.5,1x)) 
!    ENDIF
    xin(1,i) = xmin(1)+dxi

!
!--multiply by appropriate wave speed
!
    term = vwave**2 - Bin(1,i)**2/rhoin(i)
    vamplx = vwave*ampl
    vamply = -vamplx*Bin(1,i)*Bin(2,i)/(rhoin(i)*term)
    vamplz = -vamplx*Bin(1,i)*Bin(3,i)/(rhoin(i)*term)
    
    velin(1,i) = vamplx*SIN(wk*dxi)
    velin(2,i) = vamply*SIN(wk*dxi)
    velin(3,i) = vamplz*SIN(wk*dxi)
!
!--perturb internal energy if not using a polytropic equation of state 
!  (do this before density is perturbed)
!
   uuin(i) = uuin(i) + przero/rhoin(i)*ampl*SIN(wk*dxi)	! if not polytropic
!    
!--perturb density if not using summation
!
    Bin(2,i) = Bin(2,i) + vwave*Bin(2,i)*vamplx/term*SIN(wk*dxi)
    Bin(3,i) = Bin(3,i) + vwave*Bin(3,i)*vamplx/term*SIN(wk*dxi)
    rhoin(i) = rhoin(i)*(1.+ampl*SIN(wk*dxi))

    IF (ihvar.NE.0) THEN
       hhin(i) = hfact*(pmass(i)/rhoin(i))**hpower	! if variable smoothing length
    ENDIF

 ENDDO

 WRITE(iprint,*) ' Wave set: Amplitude = ',ampl,' Wavelength = ',xlambda,' k = ',wk
 WRITE(iprint,*) ' sound speed  = ',spsoundi
 WRITE(iprint,*) ' alfven speed = ',SQRT(valfven2i)
 WRITE(iprint,*) ' fast speed   = ',vfast
 WRITE(iprint,*) ' slow speed   = ',vslow
 WRITE(iprint,*) ' wave speed   = ',vwave,iwave
 
 RETURN
END
