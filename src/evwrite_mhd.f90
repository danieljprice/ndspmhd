!!--------------------------------------------------------------------
!! Calculate conserved quantities etc and write to .ev file
!!--------------------------------------------------------------------
	 
SUBROUTINE evwrite(t)
 USE dimen_mhd
 USE debug
 USE loguns
 
 USE derivB
 USE options
 USE part
 USE fmagarray
!
!--define local variables
!
 IMPLICIT NONE
 INTEGER :: i,j
 REAL, INTENT(IN) :: t
 REAL :: ekin,etherm,emag,etot,momtot
 REAL :: pmassi,rhoi
 REAL, DIMENSION(ndimV) :: veli,mom
!
!--mhd
!
 REAL, DIMENSION(ndimB) :: Bi,Brhoi,fluxtot
 REAL :: B2i,Bmagi,divBi,omegamhdi,betamhdi,FdotBi
 REAL :: fluxtotmag,crosshel
 REAL :: betamhdmin,betamhdmax,FdotBmax,betamhdav
 REAL :: divBav,divBmax,divBtot,omegtol,omegamhdav,omegamhdmax
 REAL :: fracdivBok,fmagabs
!
!--allow for tracing flow
!      
 IF (trace) WRITE(iprint,*) ' Entering subroutine evwrite'
    	  
 ekin = 0.0
 etherm = 0.0
 emag = 0.0
 etot = 0.0
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
    omegamhdav = 0.
    omegamhdmax = 0.
    omegtol = 1.E-2
    fracdivBok = 0.
    fluxtot(:) = 0.
    fluxtotmag = 0.
    crosshel = 0.      
 ENDIF 
      
 DO i=1,npart
    pmassi = pmass(i)
    rhoi = rho(i)
    veli(:) = vel(:,i)	 
    mom(:) = mom(:) + pmassi*veli(:)
    ekin = ekin + 0.5*DOT_PRODUCT(veli,veli)
    etherm = etherm + uu(i)
!
!--mhd parameters
!
    IF (imhd.NE.0) THEN
       IF (imhd.GE.11) THEN
	  Bi(:) = Bfield(:,i)
	  Brhoi(:) = Bi(:)/rhoi
       ELSE
          Brhoi(:) = Bfield(:,i)
	  Bi(:) = Brhoi(:)*rhoi
       ENDIF
       B2i = DOT_PRODUCT(Bi,Bi)
       Bmagi = SQRT(B2i)
       divBi = abs(divB(i))
 
       emag = emag + 0.5*B2i/rhoi
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
!--Max component of magnetic force in the direction of B
!
       fmagabs = 1.	!SQRT(DOT_PRODUCT(fmag(:,i),fmag(:,i)))
       IF (fmagabs.GT.1.e-6) THEN
!          FdotBi = DOT_PRODUCT(fmag(:,i),Bi(:))/(fmagabs*Bmagi)
       ELSE
          FdotBi = 0.
       ENDIF
       IF (FdotBi.GT.FdotBmax) THEN
          FdotBmax = FdotBi
!          PRINT*,' Fmag,B,fdotB        = ',fmag(:,i),Bi,	&
!	  DOT_PRODUCT(fmag(:,i),Bi(:))
!          PRINT*,' absF,absB,absF*absB = ',fmagabs,Bmagi,fmagabs*Bmagi
!	  PRINT*,' FdotBmax = ',FdotBi,FdotBmax	  
       ENDIF  
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
 
 etot = etherm + ekin + emag
 momtot = SQRT(DOT_PRODUCT(mom,mom))

!
!--write line to .ev file
!     
 IF (imhd.NE.0) THEN      

    fluxtotmag = SQRT(DOT_PRODUCT(fluxtot,fluxtot))
    betamhdav = betamhdav/FLOAT(npart)
    fracdivBok = 100.*fracdivBok/FLOAT(npart)
    omegamhdav = omegamhdav/FLOAT(npart)

    WRITE(ievfile,30) t,ekin,etherm,emag,etot,momtot,fluxtotmag,	&
          crosshel,betamhdmin,betamhdav,betamhdmax,			&
	  divBav,divBmax,divBtot,					&
          FdotBmax,omegamhdav,omegamhdmax,fracdivBok
30  FORMAT(17(1pe18.10,1x),1pe8.2)
      
 ELSE

    WRITE(ievfile,40) t,ekin,etherm,emag,etot,momtot
40  FORMAT(6(1pe20.13,1x))	       

 ENDIF

!
!--flush the buffer so that the line is written to the file immediately
!
 CALL flush(ievfile)
 
 RETURN
END SUBROUTINE evwrite
