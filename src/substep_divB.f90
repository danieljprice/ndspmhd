!!------------------------------------------------------------------------
!! Computes one timestep
!! Change this subroutine to change the timestepping algorithm
!!
!! This subroutine is a cut down version of the 
!! predictor-corrector scheme to be used for sub-stepping
!! (operator splitting) of the hyperbolic divergence correction
!!
!!
!! NB: AT THE MOMENT, ONLY WORKS IF Bevol is B on input
!!
!!---------------------------------------------------------------------------

SUBROUTINE substep_divB(dtfull,nsubsteps,Bevol,psi,divB,gradpsi, &
                        x,hh,pmass,itype,npart,ntot)
 USE dimen_mhd
 USE debug
 USE loguns
 USE bound
 USE options
 USE timestep
!
!--define local variables
!
 IMPLICIT NONE
 INTEGER, INTENT(IN) :: nsubsteps,npart,ntot
 INTEGER, DIMENSION(ntot), INTENT(IN) :: itype
 REAL, INTENT(IN) :: dtfull
 REAL, DIMENSION(ndimV,ntot), INTENT(INOUT) :: Bevol,gradpsi
 REAL, DIMENSION(ntot), INTENT(INOUT) :: psi,divB
 REAL, DIMENSION(ndim,ntot), INTENT(IN) :: x
 REAL, DIMENSION(ntot), INTENT(IN) :: hh, pmass
 REAL, DIMENSION(ntot) :: psiin, rho, dpsidt
 REAL, DIMENSION(ndimV,ntot) :: Bevolin
 INTEGER :: i,j,istep,indexmax
 REAL :: dtsub,hdt,vsig2substep,vsigsubstep
!
!--allow for tracing flow
!      
 IF (trace) WRITE(iprint,*) ' Entering subroutine substep_divB'
!
!--take nsubsteps, moving the time forward by dtsub
!
 dtsub = dtfull/REAL(nsubsteps)
 vsigsubstep = nsubsteps*sqrt(vsig2max)
 vsig2substep = vsigsubstep**2
 
 print*,'Substepping: dtfull = ',dtfull,' dtsub = ',dtsub
 print*,' nsubsteps = ',nsubsteps, ' vsigmax = ',vsigsubstep,sqrt(vsig2max)
 print*,' maxdivB = ',maxval(divB(1:npart)), 'npart = ',npart,ntot
 
 psiin = psi
 Bevolin = Bevol 

 DO i=1,npart
    dpsidt(i) = -vsig2substep*divB(i) - psidecayfact*psi(i)*vsigsubstep/hh(i)
 ENDDO

 substeploop: DO istep = 1,nsubsteps
!
!--Mid-point Predictor step
!      
 hdt = 0.5*dtsub
   
 DO i=1,npart
    IF (itype(i).EQ.1) THEN    ! fixed particles
       Bevol(:,i) = Bevolin(:,i)     
       psi(i) = psiin(i)
    ELSE
       Bevol(:,i) = Bevolin(:,i) + hdt*gradpsi(:,i)
       psi(i) = psiin(i) + hdt*dpsidt(i)  
    ENDIF
 ENDDO
 
!
!--update psi and Bevol on ghosts
!
 IF (ANY(ibound.GE.2)) THEN
    DO i=npart+1,ntot
       j = ireal(i)
       psi(i) = psi(j)
       Bevol(:,i) = Bevol(:,j)
    ENDDO
 ENDIF
!
!--calculate grad psi and div B
!
 CALL get_divBgradpsi(divB,gradpsi,Bevol,psi,x,hh,pmass,rho,npart,ntot)
 print*,' max divB = ',maxval(divB(1:npart)), &
        rho(maxloc(divB(1:npart))),maxloc(divB(1:npart))

 DO i=1,npart
    dpsidt(i) = -vsig2substep*divB(i) - psidecayfact*psi(i)*vsigsubstep/hh(i)
 ENDDO
!
!--Mid-point Corrector step
!
 DO i=1,npart
    IF (itype(i).EQ.1) THEN
       Bevol(:,i) = Bevolin(:,i)
       psi(i) = psiin(i)
    ELSE
       Bevol(:,i) = Bevolin(:,i) + hdt*gradpsi(:,i)
       psi(i) = psiin(i) + hdt*dpsidt(i)          
    ENDIF      
 ENDDO        
           
 DO i=1,npart
    Bevolin(:,i) = 2.*Bevol(:,i) - Bevolin(:,i)
    Bevol(:,i) = Bevolin(:,i)
    psiin(i) = 2.*psi(i) - psiin(i)
    psi(i) = psiin(i)
 ENDDO

 ENDDO substeploop
 print*,'finished substeps'
 read*
!
!--return appropriate quantity (B or B/rho)
!
 
 IF (trace) WRITE (iprint,*) ' Exiting subroutine substep_divB'
      
 RETURN
END SUBROUTINE substep_divB
