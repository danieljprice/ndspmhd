!!--------------------------------------------------------------------------
!! This subroutine computes the energy equation terms 
!! from the interaction between two given particles i and j.
!! (called from rates).
!!
!! enables different variables to be evolved for the energy
!!
!!--------------------------------------------------------------------------

SUBROUTINE energy_terms(i,j)
 USE dimen_mhd
 USE options
 USE rates
!
!--variables shared with rates only
!
 USE rates_local
 USE rates_mhdlocal
!
!--define local variables
!
 IMPLICIT NONE
 INTEGER, INTENT(IN) :: i,j
 REAL, DIMENSION(ndimV) :: prvterm
 REAL :: prvaniso,prvaniso2 ,projprv
 REAL :: Bidotvj,Bjdotvi


 IF (iener.EQ.2) THEN
!
!--1st law of thermodynamics (evolving thermal energy per unit mass)
!
    dudt(i) = dudt(i) + Prho2i*drhodti
    dudt(j) = dudt(j) + Prho2j*drhodtj
    
 ELSEIF (iener.GE.3) THEN	! total energy/particle
!
!--total energy equation
!    
    Bidotvj = DOT_PRODUCT(Brhoi,velj)
    Bjdotvi = DOT_PRODUCT(Brhoj,veli)
! (isotropic stress)
    prvterm(:) = (Prho2i+0.5*Brho2i)*velj(:)*grkerni/gradhi &
               + (Prho2j+0.5*Brho2j)*veli(:)*grkernj/gradhj
! (anisotropic stress)
    prvaniso =  - Bidotvj*projBrhoi*grkerni/gradhi	& 
                - Bjdotvi*projBrhoj*grkernj/gradhj
    projprv = DOT_PRODUCT(prvterm,dr)		   
!
! (Joe's anticlumping term)		   
!
    IF (ianticlump.EQ.1) THEN
!		      prvaniso = prvaniso - 0.5*Rjoe*prvaniso
! (if applied in x-direction only)
       Bidotvj = DOT_PRODUCT(Brhoi(1:ndim),velj(1:ndim))
       Bjdotvi = DOT_PRODUCT(Brhoj(1:ndim),veli(1:ndim))
       prvaniso = prvaniso - 0.5*Rjoe*(- Bidotvj*projBrhoi*grkerni/gradhi  &
                           - Bjdotvi*projBrhoj*grkernj/gradhj)
    ENDIF	   

    dendt(i) = dendt(i) - pmassj*(projprv+prvaniso)
    dendt(j) = dendt(j) + pmassi*(projprv+prvaniso)     		
 ELSE
!
!--default is 'symmetric' form of thermal energy equation in Monaghan (1992)
!
    dudt(i) = dudt(i) + 0.5*pmassj*prterm*dvdotr
    dudt(j) = dudt(j) + 0.5*pmassi*prterm*dvdotr 
 ENDIF

 RETURN
END SUBROUTINE energy_terms
