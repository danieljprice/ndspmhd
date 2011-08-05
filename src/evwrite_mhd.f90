!!--------------------------------------------------------------------
!! Calculate conserved quantities etc and write to .ev file
!!--------------------------------------------------------------------
	 
SUBROUTINE evwrite(t,etot,momtot)
 USE dimen_mhd
 USE debug
 USE loguns
 
 USE derivB
 USE options
 USE part
 USE rates
 USE fmagarray
!
!--define local variables
!
 IMPLICIT NONE
 INTEGER :: i,j,ierr
 REAL, INTENT(IN) :: t
 REAL, INTENT(OUT) :: etot,momtot
 REAL :: ekin,etherm,emag,epot
 REAL :: pmassi,rhoi
 REAL, DIMENSION(ndimV) :: veli,mom
!
!--mhd
!
 REAL, DIMENSION(ndimB) :: Bi,Brhoi,fluxtot
 REAL :: B2i,Bmagi
 REAL :: fluxtotmag,crosshel
 REAL :: betamhdi,betamhdmin,betamhdmax,betamhdav
 REAL :: fdotBi,fdotBmax,fdotBav
 REAL :: forcemagi,force_erri,force_err_max,force_err_av
 REAL :: divBi,divBav,divBmax,divBtot
 REAL :: omegamhdi,omegamhdav,omegamhdmax
 REAL :: omegtol,fracdivBok
 REAL :: fmagabs
!
!--allow for tracing flow
!      
 IF (trace) WRITE(iprint,*) ' Entering subroutine evwrite'
    	  
 ekin = 0.0
 etherm = 0.0
 emag = 0.0
 etot = 0.0
 epot = 0.0
 mom(:) = 0.0
 momtot = 0.0
!
!--mhd parameters
!     
 IF (imhd.NE.0) THEN
    betamhdav = 0.
    betamhdmax = 0.
    betamhdmin = 1.E30
    divBmax = 0.
    divBav = 0.
    divBtot = 0.
    FdotBmax = 0.
    FdotBav = 0.
    force_err_max = 0.
    force_err_av = 0.
    omegamhdav = 0.
    omegamhdmax = 0.
    omegtol = 1.E-2
    fracdivBok = 0.
    fluxtot(:) = 0.
    fluxtotmag = 0.
    crosshel = 0.      
 ENDIF 

!
!--should really recalculate the thermal energy from the total energy here
!  (otherwise uu is from the half time step and same with Bfield)
! 
 CALL conservative2primitive
      
 DO i=1,npart

    pmassi = pmass(i)
    rhoi = rho(i)
    veli(:) = vel(:,i)	 
    mom(:) = mom(:) + pmassi*veli(:)
    ekin = ekin + 0.5*pmassi*DOT_PRODUCT(veli,veli)
    etherm = etherm + pmassi*uu(i)
!
!--potential energy from external forces
!    
    SELECT CASE(iexternal_force)
    CASE(1)	! toy star force (x^2 potential)
       epot = epot + 0.5*pmassi*DOT_PRODUCT(x(:,i),x(:,i))
    CASE(2)	! 1/r^2 force(1/r potential)
       epot = epot + pmassi/SQRT(DOT_PRODUCT(x(:,i),x(:,i)))
    CASE(3)	! potential from n point masses
       WRITE(iprint,*) 'potential not calculated for point masses'
    END SELECT
!
!--mhd parameters
!
    IF (imhd.NE.0) THEN
       IF (imhd.GE.11) THEN
	  Bi(:) = Bcons(:,i)
	  Brhoi(:) = Bi(:)/rhoi
       ELSE
          Brhoi(:) = Bcons(:,i)
	  Bi(:) = Brhoi(:)*rhoi
       ENDIF
       B2i = DOT_PRODUCT(Bi,Bi)
       Bmagi = SQRT(B2i)
       forcemagi = SQRT(DOT_PRODUCT(force(:,i),force(:,i)))
       divBi = abs(divB(i))
 
       emag = emag + 0.5*pmassi*B2i/rhoi
!
!--Plasma beta minimum/maximum/average
!	 
       IF (B2i.LT.1.e-5) THEN
	  betamhdi = 0.
       ELSE 
	  betamhdi = pr(i)/(0.5*B2i)	    
       ENDIF
       betamhdav = betamhdav + betamhdi
       IF (betamhdi.GT.betamhdmax) betamhdmax = betamhdi
       IF (betamhdi.LT.betamhdmin) betamhdmin = betamhdi
!
!--Maximum divergence of B
!	 
       IF (divBi.GT.divBmax) divBmax = divBi
       divBav = divBav + divBi
!
!--volume integral of div B (int B.dS)
!
       divBtot = divBtot + pmassi*divBi/rhoi
!
!--Max component of magnetic force in the direction of B (should be zero)
!
       fmagabs = SQRT(DOT_PRODUCT(fmag(:,i),fmag(:,i)))
       IF (fmagabs.GT.1.e-8) THEN
          fdotBi = ABS(DOT_PRODUCT(fmag(:,i),Bi(:)))/(fmagabs*Bmagi)	  
       ELSE
          FdotBi = 0.
       ENDIF
       fdotBav = fdotBav + fdotBi
       IF (fdotBi.GT.fdotBmax) fdotBmax = fdotBi  
!
!--Compute total error in the force due to the B(div B) term
!  only slight worry with this is that fmag is calculated in rates, whilst
!  B has been evolved a bit further since then. A possible solution is to
!  evaluate these quantities just after the call to rates.
!       
       IF (forcemagi.GT.1.e-8 .AND. Bmagi.GT.1e-8) THEN
          force_erri = ABS(DOT_PRODUCT(fmag(:,i),Bi(:)))/(forcemagi*Bmagi)
       ELSE
          force_erri = 0.
       ENDIF
       force_err_av = force_err_av + force_erri
       IF (force_erri.GT.force_err_max) force_err_max = force_erri
!
!--|div B| x smoothing length / |B| (see e.g. Cerqueira and Gouveia del Pino 1999) 
!  this quantity should be less than ~0.01.
!
       IF (Bmagi.EQ.0.) THEN
          omegamhdi = 0.
       ELSE
	  omegamhdi = divBi*hh(i)/Bmagi	    
       ENDIF	   
       IF (omegamhdi.LT.omegtol) fracdivBok = fracdivBok + 1.
       IF (omegamhdi.GT.omegamhdmax) omegamhdmax = omegamhdi
       omegamhdav = omegamhdav + omegamhdi	  
!
!--Conserved magnetic flux (int B dV)
!
       pmassi = pmass(i)
       fluxtot(:) = fluxtot(:) + pmassi*Brhoi(:)
!
!--Conserved Cross Helicity (int v.B dV)
!
       crosshel = crosshel + pmassi*DOT_PRODUCT(veli,Brhoi)

    ENDIF

 ENDDO
 
 etot = etherm + ekin + emag + epot
 momtot = SQRT(DOT_PRODUCT(mom,mom))

!
!--write line to .ev file
!     
 IF (imhd.NE.0) THEN      

    fluxtotmag = SQRT(DOT_PRODUCT(fluxtot,fluxtot))
    betamhdav = betamhdav/FLOAT(npart)
    fracdivBok = 100.*fracdivBok/FLOAT(npart)
    omegamhdav = omegamhdav/FLOAT(npart)
    divBav = divBav/FLOAT(npart)
    fdotBav = fdotBav/FLOAT(npart)
    force_err_av = force_err_av/FLOAT(npart)

!    print*,'t=',t,' emag =',emag,' etot = ',etot, 'ekin = ',ekin,' etherm = ',etherm

    WRITE(ievfile,30) t,ekin,etherm,emag,etot,momtot,fluxtotmag,	&
          crosshel,betamhdmin,betamhdav,betamhdmax,			&
	  divBav,divBmax,divBtot,					&
          fdotBav,FdotBmax,force_err_av,force_err_max,			&
	  omegamhdav,omegamhdmax,fracdivBok
30  FORMAT(20(1pe18.10,1x),1pe8.2)
      
 ELSE

    WRITE(ievfile,40) t,ekin,etherm,emag,etot,momtot
40  FORMAT(6(1pe20.13,1x))	       

 ENDIF

!
!--flush the buffer so that the line is written to the file immediately
!
! CALL flush(ievfile)
 
 RETURN
END SUBROUTINE evwrite
