!----------------------------------------------------------------
!     Set up a travelling soundwave in 1D
!     perturbs particles from a uniform density grid
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
 INTEGER :: imax,its
 REAL :: xcentre,massp
 REAL :: ampl,wk,xlambda,dxmax,denom
 REAL :: dxi,dxprev,xmassfrac,func,fderiv
 REAL :: pri,spsoundi,valfven2i,vsigx,vsigy,vsigz
 REAL :: vfast,vslow,vcrap,vwave
!
!--allow for tracing flow
!
 IF (trace) WRITE(iprint,*) ' Entering subroutine setup'
!
!--initially set up a uniform density grid
! 	    
 ibound = 3	! periodic boundaries
 nbpts = 0		! use ghosts not fixed
 xmin = 0.0	! set position of boundaries
 xmax = 1.0
 xcentre = (xmax + xmin)/2.0
 imax = INT((xmax-xmin)/psep)
 massp = 1.0/imax	! average particle mass
!
!--allocate memory here
!
 CALL alloc(imax)
 
 DO i=1,imax
    xin(i) = xmin + (i-1)*psep  + 0.5*psep 
    velin(:,i) = 0.
    rhoin(i) = 1.0
    pmass(i) = massp
    uuin(i) = 0.3
    hhin(i) = hfact*(massp/rhoin(i))**hpower	 ! ie constant everywhere
    IF (imhd.GE.1) THEN 
       Bin(1,i) = 0.5
       Bin(2,i) = 0.5
       Bin(3,i) = 0.5
    ENDIF 
!   print*,i,xin(i),rhoin(i),uuin(i),pmass(i),xin(i)-xin(i-1)
 ENDDO
  
 npart = imax
 ntotal = npart
  
 ampl = 0.0055
! WRITE (*,*) 'Enter amplitude of disturbance'
! READ (*,*) ampl
 
 xlambda = 1.0
! WRITE (*,*) 'Enter wavelength lambda'
! READ (*,*) xlambda
    
 wk = 2.0*pi/xlambda
 dxmax = xmax - xmin
 denom = dxmax - ampl/wk*(COS(wk*dxmax)-1.0)
 WRITE (iprint,*) 'Wave number = ',wk
 WRITE (iprint,*) 'Amplitude = ',ampl

 DO i=1,npart
    dxi = xin(i)-xmin
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
       WRITE(*,99002) i,its,dxi-dxprev,xin(i),xmin+dxi
99002  FORMAT('Particle',i5,' converged in ',i4,	&
       ' iterations, error in x =',1(1pe10.2,1x),/,	&
       'previous x = ',1(0pf8.5,1x),			&
       'moved to x = ',1(0pf8.5,1x)) 
!    ENDIF
    xin(i) = xmin+dxi
!
!--get sound speed from equation of state (want average sound speed, so
!  before the density is perturbed)
!
    CALL equation_of_state(pri,spsoundi,uuin(i),rhoin(i),gamma,1)
!
!--multiply by appropriate wave speed
!
!    PRINT*,' sound speed = ',spsoundi, ' alfven speed = ',SQRT(valfven2i)
    valfven2i = DOT_PRODUCT(Bin(:,i),Bin(:,i))/rhoin(i)
    vfast = SQRT(0.5*(spsoundi**2 + valfven2i			&
                      + SQRT((spsoundi**2 + valfven2i)**2	&
                      - 4.*(spsoundi*Bin(1,i))**2/rhoin(i))))
    vslow = SQRT(0.5*(spsoundi**2 + valfven2i			&
                      - SQRT((spsoundi**2 + valfven2i)**2	&
                      - 4.*(spsoundi*Bin(1,i))**2/rhoin(i))))
    vcrap = SQRT(spsoundi**2 + valfven2i)
!    PRINT*,' c_f = ',vfast, ' c_slow = ',vslow
!    PRINT*,'c_f_joe = ',0.5*(SQRT(spsoundi**2 + valfven2i 	&
!      		    - 2.*spsoundi*Bin(1,i)/SQRT(rhoin(i)) )	&
!                                +SQRT(spsoundi**2 + valfven2i 	&
!      		    + 2.*spsoundi*Bin(1,i)/SQRT(rhoin(i)) ))		      
!    PRINT*,'c_crap = ',vcrap,spsoundi + SQRT(valfven2i)
		    
    vwave = vfast
    vsigy = -Bin(1,i)*Bin(2,i)/rhoin(i)/(vwave**2 - Bin(1,i)**2/rhoin(i))
    vsigz = -Bin(1,i)*Bin(3,i)/rhoin(i)/(vwave**2 - Bin(1,i)**2/rhoin(i))
!    PRINT*,' vsig = ',vsigx,vsigy,vsigz
!    READ*
    
    velin(1,i) = vwave*ampl*SIN(wk*dxi)
    velin(2,i) = vsigy*ampl*SIN(wk*dxi)
    velin(3,i) = vsigz*ampl*SIN(wk*dxi)
!
!--perturb internal energy if not using a polytropic equation of state 
!  (do this before density is perturbed)
!
   uuin(i) = uuin(i) + pri/rhoin(i)*ampl*SIN(wk*dxi)	! if not polytropic
!    
!--perturb density if not using summation
!
    Bin(2,i) = Bin(2,i) + vwave*Bin(2,i)/(vwave**2.-Bin(1,i)**2/rhoin(i))*ampl*SIN(wk*dxi)
    Bin(3,i) = Bin(3,i) + vwave*Bin(3,i)/(vwave**2.-Bin(1,i)**2/rhoin(i))*ampl*SIN(wk*dxi)
    rhoin(i) = rhoin(i)*(1.+ampl*SIN(wk*dxi))

    IF (ihvar.NE.0) THEN
       hhin(i) = hfact*(pmass(i)/rhoin(i))**hpower	! if variable smoothing length
    ENDIF

 ENDDO

 
 RETURN
END
