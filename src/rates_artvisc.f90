!!-------------------------------------------------------------------------
!! This subroutine computes all of the artificial dissipation terms
!! and adds them to the appropriate equations
!!
!! We should include this in the main rates file, rather than calling it as a
!! separate subroutine. This is for speed.
!!-------------------------------------------------------------------------

SUBROUTINE artvis_terms(i,j)
 USE dimen_mhd
 USE debug
 USE loguns
 
 USE artvi
 USE options
 USE rates
 USE part
!
!--variables shared with rates only
!
 USE rates_local
 USE rates_mhdlocal
!
!--local variables
!
 IMPLICIT NONE
 INTEGER, INTENT(IN) :: i,j
 REAL :: visc,alphaav
 REAL, DIMENSION(ndimB) :: Bvisc,dBdtvisc
 REAL :: v2i,v2j,B2i,B2j
 REAL :: eni,enj,ediff,ediffB,qdiff
 REAL :: vissv,vissB,vissu,envisc

!
!--artificial viscosity terms
!
 alphaav = 0.5*(alphai + alpha(j))		
!
!--viscosity
!
 visc = 0.5*alphaav*vsig*viss/rhoav*grkern	! viss defined in rates
!
!--add this to the momentum equation (viscous force)
! 
 IF (ndimV.GT.ndim) THEN
    force(:,i) = force(:,i) - pmassj*visc*vunit(:)
    force(:,j) = force(:,j) + pmassi*visc*vunit(:)
 ELSE
    force(:,i) = force(:,i) - pmassj*visc*dr(:)
    force(:,j) = force(:,j) + pmassi*visc*dr(:)		
 ENDIF
  
!
!--ohmic dissipation
!
 Bvisc(:) = dB(:) - dr(:)*projdB
 dBdtvisc(:) = 0.5*alphaav*vsig*Bvisc(:)/rhoav**2
!
!--add this to the induction equation
!
 IF (imhd.LE.10) THEN		! evolving B/rho
    dBfielddt(:,i) = dBfielddt(:,i) + rhoi*pmassj*dBdtvisc(:)*grkern		   
    dBfielddt(:,j) = dBfielddt(:,j) - rhoj*pmassi*dBdtvisc(:)*grkern
 ELSEIF (imhd.EQ.11) THEN	! evolving B
    dBfielddt(:,i) = dBfielddt(:,i) + rho2i*pmassj*dBdtvisc(:)*grkern		   
    dBfielddt(:,j) = dBfielddt(:,j) - rho2j*pmassi*dBdtvisc(:)*grkern
 ENDIF

!
!--add dissipation terms to the energy equation
!
 IF (iener.EQ.3) THEN
!
!--total energy equation
!
    IF (ndimV.GT.ndim) THEN
!      v2i = DOT_PRODUCT(veli,veli)	! total kinetic energy
!      v2j = DOT_PRODUCT(velj,velj)
       v2i = DOT_PRODUCT(veli,vunit)**2.	! along velocity line
       v2j = DOT_PRODUCT(velj,vunit)**2.	! (use in 1.5D)
    ELSE
       v2i = DOT_PRODUCT(veli,dr)**2	! energy along line
       v2j = DOT_PRODUCT(velj,dr)**2	! of sight		   
    ENDIF
    B2i = (DOT_PRODUCT(Bi,Bi) - DOT_PRODUCT(Bi,dr)**2)	!  "   " 
    B2j = (DOT_PRODUCT(Bj,Bj) - DOT_PRODUCT(Bj,dr)**2)
    enj = gconst*uu(j) + 0.5*v2j + 0.5*B2j/rhoav
    eni = gconst*uu(i) + 0.5*v2i + 0.5*B2i/rhoav
    ediffB = 0.5*(B2i-B2j)/rhoav	! needed if applying Bvisc everywhere
    ediff = eni - enj 
    qdiff = -0.5*alphaav*vsig*ediff*rij/(rhoav)*grkern
!
!--add to total energy equation
!
    dendt(i) = dendt(i) - pmassj*qdiff
    dendt(j) = dendt(j) + pmassi*qdiff     		

 ELSEIF (iener.GT.0) THEN
!
!--thermal energy equation
!    
!  (kinetic energy term)
    IF (ndimV.GT.ndim) THEN
       vissv = -0.5*(DOT_PRODUCT(veli,vunit) - DOT_PRODUCT(velj,vunit))**2
!      vissv = -0.5*DOT_PRODUCT(dvel,dvel)
    ELSE
       vissv = -0.5*(DOT_PRODUCT(veli,dr) - DOT_PRODUCT(velj,dr))**2		    
    ENDIF
    
!  (magnetic energy term)
    vissB = -0.5*(DOT_PRODUCT(dB,dB)-projdB**2)/rhoav
!  (thermal energy term - thermal conductivity set by gconst)
    vissu = gconst*(uu(i) - uu(j))
!  (total contribution to thermal energy equation)
    envisc = 0.5*alphaav*vsig*(vissv+vissB+vissu)/rhoav*grkern!*abs(rx)
!
!--add to thermal energy equation
! 
    dudt(i) = dudt(i) + pmassj*envisc
    dudt(j) = dudt(j) + pmassi*envisc
 
 ENDIF

 RETURN
END SUBROUTINE artvis_terms
