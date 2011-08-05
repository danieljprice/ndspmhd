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
 REAL :: vfast,vslow,vcrap,vwave,term,dens1
 REAL :: denszero,uuzero,przero, Rzero
!
!--allow for tracing flow
!
 IF (trace) WRITE(iprint,*) ' Entering subroutine setup'
10 FORMAT(/,'-------------- ',a,' ----------------')

 iwave = 0 		! preset wave parameters
!
!--setup parameters (could read in from a file)
!
 Bzero = 0.
 IF (iwave.EQ.1) THEN        ! MHD slow wave
    ampl = 0.006
    denszero = 1.0
    uuzero = 4.5
    IF (imhd.ne.0) Bzero(1:3) = SQRT(2.)
    WRITE(iprint,10) 'MHD slow wave'
 ELSEIF (iwave.EQ.2) THEN    ! MHD fast wave
    ampl = 0.0055
    denszero = 1.0
    uuzero = 0.3
    IF (imhd.ne.0) Bzero(1:3) = 0.5  
    WRITE(iprint,10) 'MHD fast wave'
 ELSE                        ! Sound wave
    ampl = 0.005
    denszero = 1.0
    if (abs(gamma-1.).gt.1e-3) then
       uuzero = 1.0/((gamma-1.)*gamma)
    else
       uuzero = 1.0
    endif
    Rzero = -1.0  !  negative stress parameter
    OPEN(unit=20,ERR=30,file=trim(rootname)//'.rstress',status='old')
      READ(20,*) Rzero
    CLOSE(unit=20)
    PRINT*,'Read stress file: R parameter = ',Rzero
30  CONTINUE    
    Bzero(1) = sqrt(2.*(1.-Rzero))
    WRITE(iprint,10) 'sound wave'
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
    xmax(2:ndim) = 6.*psep ! would need to adjust this depending on grid setup
 ENDIF
!
!--initially set up a uniform density grid (also determines npart)
!  (the call to set_uniform_cartesian means this works in 1,2 and 3D)
!
 CALL set_uniform_cartesian(2,psep,xmin,xmax,.false.)
 
! npart = INT((xmax(1)-xmin(1))/psep)
 massp = 1.0/FLOAT(npart)	! average particle mass
!
!--allocate memory here
!
! CALL alloc(npart)
!
!--setup uniform density grid of particles
! 
 DO i=1,npart
!    x(1,i) = xmin(1) + (i-1)*psep  + 0.5*psep 
    vel(:,i) = 0.
    dens(i) = denszero
    pmass(i) = massp
    uu(i) = uuzero
    IF (imhd.GT.0) THEN 
       Bfield(:,i) = Bzero
    ELSE
       Bfield(:,i) = 0.
    ENDIF 
 ENDDO

 ntotal = npart
!
!--get sound speed from equation of state (want average sound speed, so
!  before the density is perturbed)
!
 CALL equation_of_state(przero,spsoundi,uuzero,denszero,gamma,polyk,1)
!
!--work out MHD wave speeds
!
 dens1 = 1./denszero
 
 valfven2i = DOT_PRODUCT(Bzero,Bzero)*dens1
 vfast = SQRT(0.5*(spsoundi**2 + valfven2i			&
                 + SQRT((spsoundi**2 + valfven2i)**2	&
                 - 4.*(spsoundi*Bzero(1))**2*dens1)))
 vslow = SQRT(0.5*(spsoundi**2 + valfven2i			&
                 - SQRT((spsoundi**2 + valfven2i)**2	&
                 - 4.*(spsoundi*Bzero(1))**2*dens1)))
 vcrap = SQRT(spsoundi**2 + valfven2i)

!----------------------------
! set the wave speed to use
!----------------------------
		    
 IF (iwave.EQ.1) THEN
    vwave = vslow
    IF (imhd.le.0) WRITE(iprint,*) 'Error: can''t have slow wave if no mhd'
 ELSEIF (iwave.EQ.2) THEN
    vwave = vfast
 ELSE
    vwave = spsoundi
 ENDIF
 IF (vwave.le.0.) THEN
    WRITE(iprint,*) 'Error in setup: vwave = ',vwave
    STOP
 ENDIF
  
 dxmax = xmax(1) - xmin(1)
 denom = dxmax - ampl/wk*(COS(wk*dxmax)-1.0)
 WRITE (iprint,*) 'Wave number = ',wk
 WRITE (iprint,*) 'Amplitude = ',ampl
!
!--now perturb the particles appropriately
!
 DO i=1,npart
    
    dxi = x(1,i)-xmin(1)
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
!       WRITE(*,99002) i,its,dxi-dxprev,x(1,i),xmin(1)+dxi
!99002  FORMAT('Particle',i5,' converged in ',i4,	&
!       ' iterations, error in x =',1(1pe10.2,1x),/,	&
!       'previous x = ',1(0pf8.5,1x),			&
!       'moved to x = ',1(0pf8.5,1x)) 
!    ENDIF
    x(1,i) = xmin(1)+dxi

!
!--multiply by appropriate wave speed
!
    term = vwave**2 - Bfield(1,i)**2/dens(i)
    vamplx = vwave*ampl
    vamply = -vamplx*Bfield(1,i)*Bfield(2,i)/(dens(i)*term)
    vamplz = -vamplx*Bfield(1,i)*Bfield(3,i)/(dens(i)*term)
    
    vel(1,i) = vamplx*SIN(wk*dxi)
!    vel(2,i) = vamply*SIN(wk*dxi)
!    vel(3,i) = vamplz*SIN(wk*dxi)
!
!--perturb internal energy if not using a polytropic equation of state 
!  (do this before density is perturbed)
!
   uu(i) = uu(i) + przero/dens(i)*ampl*SIN(wk*dxi)	! if not polytropic
!    
!--perturb density if not using summation
!
!    Bfield(2,i) = Bfield(2,i) + vwave*Bfield(2,i)*vamplx/term*SIN(wk*dxi)
!    Bfield(3,i) = Bfield(3,i) + vwave*Bfield(3,i)*vamplx/term*SIN(wk*dxi)
    dens(i) = dens(i)*(1.+ampl*SIN(wk*dxi))

 ENDDO

 WRITE(iprint,*) ' Wave set: Amplitude = ',ampl,' Wavelength = ',xlambda,' k = ',wk
 WRITE(iprint,*) ' sound speed  = ',spsoundi
 WRITE(iprint,*) ' alfven speed = ',SQRT(valfven2i)
 WRITE(iprint,*) ' fast speed   = ',vfast
 WRITE(iprint,*) ' slow speed   = ',vslow
 WRITE(iprint,*) ' wave speed   = ',vwave,iwave
 
 RETURN
END
