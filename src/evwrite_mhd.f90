!!--------------------------------------------------------------------
!! Calculate conserved quantities etc and write to .ev file
!!--------------------------------------------------------------------
  
SUBROUTINE evwrite(t,etot,momtot)
 USE dimen_mhd, only:ndim,ndimV
 USE debug, only:trace
 USE loguns, only:iprint,ievfile
 
 USE derivB
 USE options
 USE part
 USE rates, only:force,potengrav
 USE fmagarray
 USE timestep, only:dt
!
!--define local variables
!
 IMPLICIT NONE
 INTEGER :: i
 REAL, INTENT(IN) :: t
 REAL, INTENT(OUT) :: etot,momtot
 REAL :: ekin,etherm,emag,epot
 REAL :: pmassi,rhoi,epoti,rr
 REAL, DIMENSION(ndim) :: rhat
 REAL, DIMENSION(ndimV) :: veli,mom,ang,angi
 REAL :: alphai,betai,alphatstarav,betatstarav
!
!--mhd
!
 REAL, DIMENSION(ndimB) :: Bi,Brhoi,fluxtot
 REAL :: B2i,Bmagi
 REAL :: fluxtotmag,crosshel
 REAL :: betamhdi,betamhdmin,betamhdav
 REAL :: fdotBi,fdotBmax,fdotBav
 REAL :: forcemagi,force_erri,force_err_max,force_err_av
 REAL :: divBi,divBav,divBmax,divBtot
 REAL :: omegamhdi,omegamhdav,omegamhdmax
 REAL :: fracdivBok
 REAL, PARAMETER :: omegtol = 1.E-2
 REAL :: fmagabs,rhomax,rhomean,rhomin
 REAL :: angtot,ekiny,emagp
!
!--allow for tracing flow
!      
 IF (trace) WRITE(iprint,*) ' Entering subroutine evwrite'
       
 ekin = 0.0
 etherm = 0.0
 emag = 0.0
 emagp = 0.
 etot = 0.0
 IF (igravity.NE.0) THEN
    epot = potengrav
 ELSE
    epot = 0.
 ENDIF
 mom(:) = 0.0
 momtot = 0.0
 ang(:) = 0.
 ekiny = 0.
! alphatstarav = 0.
! betatstarav = 0.
!
!--mhd parameters
!     
 IF (imhd.NE.0) THEN
    betamhdav = 0.
    betamhdmin = huge(betamhdmin)
    divBmax = 0.
    divBav = 0.
    divBtot = 0.
    FdotBmax = 0.
    FdotBav = 0.
    force_err_max = 0.
    force_err_av = 0.
    omegamhdav = 0.
    omegamhdmax = 0.
    fracdivBok = 0.
    fluxtot(:) = 0.
    fluxtotmag = 0.
    crosshel = 0.      
 ENDIF 

!
!--should really recalculate the thermal energy from the total energy here
!  (otherwise uu is from the half time step and same with Bfield)
! 
! CALL conservative2primitive
      
 DO i=1,npart

    pmassi = pmass(i)
    rhoi = rho(i)
    veli(:) = vel(:,i)  
    mom(:) = mom(:) + pmassi*veli(:)
    if (ndim.eq.3) then
       call cross_product3D(x(:,i),veli(:),angi(:))
       ang(:) = ang(:) + pmassi*angi(:)
    elseif (ndim.eq.2) then
       ang(3) = ang(3) + pmassi*(x(1,i)*veli(2) - x(2,i)*veli(1))
    endif
    ekin = ekin + 0.5*pmassi*DOT_PRODUCT(veli,veli)
    if (ndim.ge.2) ekiny = ekiny + 0.5*pmassi*vel(1,i)*vel(1,i)
    etherm = etherm + pmassi*uu(i)
!
!--potential energy from external forces
!    
    CALL external_potentials(iexternal_force,x(:,i),epoti,ndim)
    epot = epot + pmassi*epoti
    
!    rr = dot_product(x(1:ndim,i),x(1:ndim,i))
!    if (rr.gt.tiny(rr)) then
!       rhat(1:ndim) = x(1:ndim,i)/rr
!    else
!       rhat = 0.
!    endif
!    alphai = dot_product(vel(1:ndim,i),rhat(1:ndim))
!    betai = vel(2,i)*rhat(1) - vel(1,i)*rhat(2)
!    alphatstarav = alphatstarav + alphai
!    betatstarav = betatstarav + betai
!
!--mhd parameters
!
    IF (imhd.NE.0) THEN
       Bi(:) = Bfield(:,i)
       Brhoi(:) = Bi(:)/rhoi
       B2i = DOT_PRODUCT(Bi,Bi)
       Bmagi = SQRT(B2i)
       forcemagi = SQRT(DOT_PRODUCT(force(:,i),force(:,i)))
       divBi = abs(divB(i))
 
       emag = emag + 0.5*pmassi*B2i/rhoi
       emagp = emagp + 0.5*pmassi*dot_product(Bi(1:2),Bi(1:2))/rhoi
!
!--Plasma beta minimum/maximum/average
!  
       IF (B2i.LT.tiny(B2i)) THEN
          betamhdi = 0.
       ELSE 
          betamhdi = pr(i)/(0.5*B2i)     
       ENDIF
       betamhdav = betamhdav + betamhdi
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
       IF (fmagabs.GT.1.e-8 .and. Bmagi.gt.1.e-8) THEN
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
       IF (Bmagi.lt.1e-8) THEN
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
 
 etot = ekin + emag + epot
 IF (iprterm.ge.0 .or. iprterm.lt.-1) etot = etot + etherm
 momtot = SQRT(DOT_PRODUCT(mom,mom))
 CALL minmaxave(rho(1:npart),rhomin,rhomax,rhomean,npart)
 angtot = SQRT(DOT_PRODUCT(ang,ang))
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

!!    print*,'t=',t,' emag =',emag,' etot = ',etot, 'ekin = ',ekin,' etherm = ',etherm

    WRITE(ievfile,30) t,ekin,etherm,emag,epot,etot,momtot,angtot,rhomax,rhomean,dt, &
          emagp,crosshel,betamhdmin,betamhdav,  &
          divBav,divBmax,divBtot,     &
          fdotBav,FdotBmax,force_err_av,force_err_max,   &
          omegamhdav,omegamhdmax,fracdivBok
30  FORMAT(24(1pe18.10,1x),1pe8.2)
      
 ELSE
    !alphatstarav = alphatstarav/FLOAT(npart)
    !betatstarav = betatstarav/FLOAT(npart)
   !! print*,'t=',t,' emag =',emag,' etot = ',etot, 'ekin = ',ekin,' etherm = ',etherm

    WRITE(ievfile,40) t,ekin,etherm,emag,epot,etot,momtot,angtot,rhomax,rhomean,dt,ekiny
40  FORMAT(24(1pe18.10,1x))        

 ENDIF

!
!--flush the buffer so that the line is written to the file immediately
!
 CALL flush(ievfile)
 
 RETURN
END SUBROUTINE evwrite
