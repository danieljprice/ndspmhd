!----------------------------------------------------------------
!     Set up a wave in 2 or 3 dimensions
!     propagating in some arbitrary direction relative to the 
!     x, y and z axes
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
 INTEGER :: its
 REAL, DIMENSION(ndim) :: runit,dxi,dxmax
 REAL, DIMENSION(ndimV) :: velzero
 REAL, DIMENSION(ndimB) :: Bzero
 REAL :: massp,totmass,rhozero,gam1,uuzero,przero
 REAL :: anglexy,ampl,wk,xlambda,rmax,denom
 REAL :: ri,rprev,xmassfrac,func,fderiv
 REAL :: pri,spsoundi,valfven2i
 REAL :: vampl_par,vampl_perp,vamplz,vamply,vamplx
 REAL :: vparallel,vperp,vz
 REAL :: Bparallel,Bperp,Bz
 REAL :: Bamplx,Bamply,Bamplz,Bampl_perp
 REAL :: vfast,vslow,vcrap,vwave,term,rhoin1
!
!--allow for tracing flow
!
 IF (trace) WRITE(iprint,*) ' Entering subroutine setup'
 WRITE(iprint,*) '------------ Wave setup ----------------'
!
!--set direction of wave propagation (runit is unit vector in this direction)
!
 anglexy = 30	! angle in degrees x,y plane
! anglez = 45	! angle in degrees z plane
 anglexy = anglexy*pi/180.	! convert to radians

 IF (ndim.EQ.2) THEN
    runit(1) = COS(anglexy)
    runit(2) = SIN(anglexy)
 ELSE 
    STOP 'This wave setup only for 2D'       
!      runit(3) = 0.
 ENDIF
 WRITE(iprint,*) ' runit = ',runit
!
!--read/set wave parameters
! 
 ampl = 0.001
! WRITE (*,*) 'Enter amplitude of disturbance'
! READ (*,*) ampl
 
 xlambda = 1.0
! WRITE (*,*) 'Enter wavelength lambda'
! READ (*,*) xlambda
    
 wk = 2.0*pi/xlambda	! 	wave number
 
!
!--set boundaries
! 	    
 ibound = 3	! periodic boundaries
 nbpts = 0	! no fixed particles
 xmin(:) = 0.0	! set position of boundaries
 xmax(:) = 1.0/runit(:)
! PRINT*,'xmin,xmax = ',xmin,xmax
!
!--setup parameters
!
 vperp = 0.1
 vparallel = 0.1
 vz = 0.1
 rhozero = 1.0
 przero = 0.1
 Bparallel = 1.0
 Bperp = 0.1
 Bz = 0.1
 uuzero = 0.3
!
!--work out dependent parameters
!
 gam1 = gamma - 1.
! uuzero = przero/(gam1*rhozero)
 velzero(1) = vparallel*runit(1) - vperp*runit(2)
 velzero(2) = vparallel*runit(2) + vperp*runit(1)
 velzero(3) = vz
 Bzero(1) = Bparallel*runit(1) - Bperp*runit(2)
 Bzero(2) = Bparallel*runit(2) + Bperp*runit(1)
 Bzero(3) = Bz
!
!--initially set up a uniform density grid (also determines npart)
!
 CALL set_uniform_cartesian(2,xmin,xmax,.false.)	! 2 = close packed
!
!--determine particle mass
!
 dxmax(:) = xmax(:) - xmin(:)
 totmass = rhozero*PRODUCT(dxmax)
 massp = totmass/FLOAT(npart) ! average particle mass
 PRINT*,'npart,massp = ',npart,massp
 
 DO i=1,npart
    velin(:,i) = velzero
    rhoin(i) = rhozero
    pmass(i) = massp
    uuin(i) = uuzero
    hhin(i) = hfact*(pmass(i)/rhozero)**hpower	 ! ie constant everywhere
    IF (imhd.GE.1) THEN 
       Bin(:,i) = Bzero
    ENDIF 
 ENDDO

 ntotal = npart
  
 rmax = DOT_PRODUCT(dxmax,runit)
! PRINT*,'rmax = ',rmax
 denom = rmax - ampl/wk*(COS(wk*rmax)-1.0)
!
!--get sound speed from equation of state (want average sound speed, so
!  before the density is perturbed)
!
 CALL equation_of_state(przero,spsoundi,uuzero,rhozero,gamma,1)
!
!--work out MHD wave speeds
!
 rhoin1 = 1./rhozero
 
 valfven2i = Bparallel**2*rhoin1
 vfast = SQRT(0.5*(spsoundi**2 + valfven2i			&
                 + SQRT((spsoundi**2 + valfven2i)**2	&
                 - 4.*(spsoundi*Bparallel)**2*rhoin1)))
 vslow = SQRT(0.5*(spsoundi**2 + valfven2i			&
                 - SQRT((spsoundi**2 + valfven2i)**2	&
                 - 4.*(spsoundi*Bparallel)**2*rhoin1)))
 vcrap = SQRT(spsoundi**2 + valfven2i)

!----------------------------
! set the wave speed to use
!----------------------------
		    
 vwave = vfast

!----------------------------

!
!--now perturb the particles appropriately
!    
 DO i=1,npart
 
    dxi(:) = xin(:,i) - xmin(:)
    ri = DOT_PRODUCT(dxi,runit)

    rprev = 2.*rmax
    xmassfrac = ri/rmax	! current mass fraction(for uniform density)
				! determines where particle should be
!
!--Use rootfinder on the integrated density perturbation
!  to find the new position of the particle
!    
    its = 0
	 
    DO WHILE ((ABS(ri-rprev).GT.tol).AND.(its.LT.itsmax))
       rprev = ri
       func = xmassfrac*denom - (ri - ampl/wk*(COS(wk*ri)-1.0))
       fderiv = -1.0 - ampl*SIN(wk*ri)
       ri = ri - func/fderiv	! Newton-Raphson iteration
       its = its + 1 
!      PRINT*,'iteration',its,'ri =',ri 
    ENDDO
	 	 	 
    IF (its.GE.itsmax) THEN
       WRITE(iprint,*) 'Error: soundwave - too many iterations'
       CALL quit 
    ENDIF
	  
!    IF (idebug(1:5).EQ.'sound') THEN
!       WRITE(*,99002) i,its,ri-rprev,xin(1,i),xmin(1)+ri*runit(1)
99002  FORMAT('Particle',i5,' converged in ',i4,	&
       ' iterations, error in x =',1(1pe10.2,1x),/,	&
       'previous r = ',1(0pf8.5,1x),			&
       'moved to r = ',1(0pf8.5,1x)) 
!    ENDIF
    xin(:,i) = xmin(:) + ri*runit(:)
!
!--multiply by the appropriate amplitudes
!
    rhoin1 = 1./rhoin(i)
    term = 1./(vwave**2 - Bparallel**2*rhoin1)
    vampl_par = vwave*ampl
    vampl_perp = -vampl_par*Bparallel*Bperp*rhoin1*term
    vamplz = -vampl_par*Bparallel*Bin(3,i)*rhoin1*term
    
    vamplx = vampl_par*runit(1) - vampl_perp*runit(2)
    vamply = vampl_par*runit(2) + vampl_perp*runit(1)

    velin(1,i) = vamplx*SIN(wk*ri)
    velin(2,i) = vamply*SIN(wk*ri)
    velin(3,i) = vamplz*SIN(wk*ri)
!
!--perturb internal energy if not using a polytropic equation of state 
!  (do this before density is perturbed)
!
   uuin(i) = uuin(i) + pri/rhoin(i)*ampl*SIN(wk*ri)	! if not polytropic
!    
!--perturb density if not using summation
!
    Bampl_perp = vwave*Bperp*vampl_par*term*SIN(wk*ri)
    Bamplz = vwave*Bz*vampl_par*term*SIN(wk*ri)
    
    Bamplx = -Bampl_perp*runit(2)
    Bamply = Bampl_perp*runit(1)
    
    Bin(1,i) = Bzero(1) + Bamplx*SIN(wk*ri)
    Bin(2,i) = Bzero(2) + Bamply*SIN(wk*ri)
    Bin(3,i) = Bzero(3) + Bamplz*SIN(wk*ri)
    rhoin(i) = rhoin(i)*(1.+ampl*SIN(wk*ri))

    IF (ihvar.NE.0) THEN
       hhin(i) = hfact*(pmass(i)/rhoin(i))**hpower	! if variable smoothing length
    ENDIF

 ENDDO

 WRITE(iprint,*) ' Wave set: Amplitude = ',ampl,' Wavelength = ',xlambda,' k = ',wk
 WRITE(iprint,*) ' sound speed  = ',spsoundi
 WRITE(iprint,*) ' alfven speed = ',SQRT(valfven2i)
 WRITE(iprint,*) ' fast speed   = ',vfast
 WRITE(iprint,*) ' slow speed   = ',vslow
 WRITE(iprint,*) ' wave speed   = ',vwave
 
 RETURN
END
