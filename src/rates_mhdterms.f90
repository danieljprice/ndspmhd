!!--------------------------------------------------------------------------
!! This subroutine computes the contribution to the magnetic field
!! terms from the interaction between two given particles i and j.
!! (called from rates)
!!--------------------------------------------------------------------------

SUBROUTINE mhd_terms(i,j)
 USE dimen_mhd
 USE debug
 USE loguns

 USE anticlumping
 USE derivB
 USE fmagarray
 USE options
 USE rates
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
 REAL, DIMENSION(ndimB) :: faniso,fmagi,fmagj,fcorr
 REAL, DIMENSION(ndimB) :: curlBi,divBvec
 REAL :: fiso,divBonrho,divBvec
 REAL :: BidotdB,BjdotdB
 REAL :: wabjoe
!
!--calculate curl B for current (only if field is 3D)
!
 IF (ndimB.EQ.3) THEN
    curlBi(1) = dB(2)*dr(3) - dB(3)*dr(2)
    curlBi(2) = dB(3)*dr(1) - dB(1)*dr(3)
    curlBi(3) = dB(1)*dr(2) - dB(2)*dr(1)
 ENDIF
!
!--Lorentz force (different formalisms)
!
 IF (imagforce.EQ.1) THEN	! vector form (dot products)
		  		  
    BidotdB = DOT_PRODUCT(Brhoi,dB)
    BjdotdB = DOT_PRODUCT(Brhoj,dB)
    fmagi(:) = grkern*(dr(:)*BidotdB - projBrhoi*dB(:))		  
    fmagj(:) = grkern*(dr(:)*BjdotdB - projBrhoj*dB(:))
		  		  			
 ELSEIF (imagforce.EQ.2) THEN	! tensor formalism
!
!--isotropic mag force (pressure gradient)
!		  
    fiso = 0.5*(Brho2i*grkerni/gradhi + Brho2j*grkernj/gradhj)
!
!--anisotropic force (including variable smoothing length terms)
!
    faniso(:) = Brhoi(:)*projBrhoi*grkerni/gradhi  &
	      + Brhoj(:)*projBrhoj*grkernj/gradhj		   
!
!--Joe's correction term (preserves momentum conservation)
!  *** in 1D only do this in the x-direction?
    IF (ianticlump.EQ.1) THEN
       wabjoe = wabjoei*hfacwab
       Rjoe = eps*(wab/wabjoe)**neps
!      IF (Rjoe.GT.0.1) PRINT*,'Rjoe = ',Rjoe,i,j,wab,wabjoe,x(i)-x(j)
!      faniso(:) = faniso(:) - 0.5*Rjoe*faniso(:)  
       faniso(1:ndim) = faniso(1:ndim) - 0.5*Rjoe*faniso(1:ndim)  
    ENDIF
!
!--add contributions to magnetic force
!	  
    fmagi(:) = faniso(:) - fiso*dr(:)
		  
 ELSEIF (imagforce.EQ.3) THEN	! D.Price version

    fiso = grkern*0.5*(Brho2i + Brho2j)
    faniso(:) = grkern*(Brhoi(:)*projBrhoj + Brhoj(:)*projBrhoi)     
    fmagi(:) = faniso(:) - fiso*dr(:)
		  		  
 ELSEIF (imagforce.EQ.4) THEN ! form that goes with alt. time deriv

    rhoij = rhoi*rhoj
    fiso = grkern*0.5*(Brho2i + Brho2j)
    faniso(:) = grkern*(Bi(:)*projBi + Bj(:)*projBj)/rhoij
    fmagi(:) = faniso(:) - fiso*dr(:)
		  	  
 ELSEIF (imagforce.EQ.5) THEN	! Morris' Hybrid form

    rhoij = rhoi*rhoj
    fiso = grkern*0.5*(Brho2i + Brho2j)
    faniso(:) = grkern*(Bj(:)*projBj - Bi(:)*projBi)/rhoij
    fmagi(:) = faniso(:) - fiso*dr(:)		  
		  
 ENDIF
!
!--divergence of B
!
 divB(i) = divB(i) - pmassj*DOT_PRODUCT(dB,dr)*grkern
 divB(j) = divB(j) - pmassi*DOT_PRODUCT(dB,dr)*grkern
!
!--compute current density J
!
 curlB(:,i) = curlB(:,i) - pmassj*curlBi(:)*grkern
 curlB(:,j) = curlB(:,j) - pmassi*curlBi(:)*grkern
!
!--stabilising correction term from Borve et al. (2001)
!
 IF (idivBzero.EQ.3) THEN	          
    divBvec(:) = Brhoi(:)/rhoi + Brhoj(:)/rhoj
    divBonrho = grkern*(DOT_PRODUCT(divBvec,dr))
    fcorr(:) = Bi(:)*divBonrho
    fmagi(:) = fmagi(:) - fcorr(:)
 ENDIF		   
!
!--add Lorentz force to total force
!        	
 IF (imagforce.EQ.1) THEN
    fmag(:,i) = fmag(:,i) + pmassj*fmagi(:)/rho2i
    fmag(:,j) = fmag(:,j) - pmassi*fmagj(:)/rho2j
    force(:,i) = force(:,i) + pmassj*fmagi(:)/rho2i
    force(:,j) = force(:,j) - pmassi*fmagj(:)/rho2j
 ELSEIF (imagforce.EQ.6) THEN	! Morris' Hybrid force
    fmag(:,i) = fmag(:,i) + pmassj*(faniso(:)-fiso*dr(:))
    fmag(:,j) = fmag(:,j) + pmassi*(faniso(:)+fiso*dr(:))
    force(:,i) = force(:,i) + pmassj*(faniso(:)-fiso*dr(:))
    force(:,j) = force(:,j) + pmassi*(faniso(:)+fiso*dr(:)) 			        
 ELSE 	! symmetric forces fmagxi = -fmagxj
    fmag(:,i) = fmag(:,i) + pmassj*fmagi(:)
    fmag(:,j) = fmag(:,j) - pmassi*fmagi(:)
    force(:,i) = force(:,i) + pmassj*fmagi(:)
    force(:,j) = force(:,j) - pmassi*fmagi(:) 
 ENDIF

!
!--time derivative of magnetic field (divide by rho later)
!
!
!--evolving B/rho as field variable
!
 IF (imhd.EQ.1) THEN   

    dBfielddt(:,i) = dBfielddt(:,i) - pmassj*(dvel(:)*projBrhoi)*grkerni/gradhi 
    dBfielddt(:,j) = dBfielddt(:,j) - pmassi*(dvel(:)*projBrhoj)*grkernj/gradhj

 ELSEIF (imhd.EQ.2) THEN	! good for rho discts

    dBfielddt(:,i) = dBfielddt(:,i) - pmassj*(dvel(:)*projBrhoi/rhoj)*grkern		   
    dBfielddt(:,j) = dBfielddt(:,j) - pmassi*(dvel(:)*projBrhoj/rhoi)*grkern
   
 ELSEIF (imhd.EQ.3) THEN	! add v*div B term

    dBfielddt(:,i) = dBfielddt(:,i) - pmassj*(dvel(:)*projBrhoi		&
                                      +veli(:)*DOT_PRODUCT(dB,dr))*grkern		   
    dBfielddt(:,j) = dBfielddt(:,j) - pmassi*(dvel(:)*projBrhoj		&
                                      +velj(:)*DOT_PRODUCT(dB,dr))*grkern
!
!--evolving B as field variable
!
 ELSEIF (imhd.EQ.11) THEN	! note divided by rho later		  

    dBfielddt(:,i) = dBfielddt(:,i)		&
                   + pmassj*(Bi(:)*dvdotr - dvel(:)*projBi)*grkerni/gradhi   
    dBfielddt(:,j) = dBfielddt(:,j) 		&
                   + pmassi*(Bj(:)*dvdotr - dvel(:)*projBj)*grkernj/gradhj
 ENDIF

RETURN		
END SUBROUTINE mhd_terms
