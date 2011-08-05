!!------------------------------------------------------------------------
!! This subroutine is a cut down version of the 
!! predictor-corrector scheme to be used for sub-stepping
!! (operator splitting) of the hyperbolic divergence correction
!!
!! The predictor-corrector scheme is cut into 2*n+1 steps
!! such that the centred rates call is common to both the main (hydro)
!! stepping and the substepping
!! ie. for n=1 during main predictor we do:
!!      predictor->corrector->predictor->
!!     during main corrector step we do:
!!      ->corrector->predictor->corrector
!!
!! NB: AT THE MOMENT, ONLY WORKS IF Bevol is B on input
!!
!!---------------------------------------------------------------------------

SUBROUTINE substep_divB(icall,dtfull,nsubstepsin,Bevol,psi,divB,gradpsi, &
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
 INTEGER, INTENT(IN) :: icall,nsubstepsin,npart,ntot
 INTEGER, DIMENSION(ntot), INTENT(IN) :: itype
 REAL, INTENT(IN) :: dtfull
 REAL, DIMENSION(ndimV,ntot), INTENT(INOUT) :: Bevol,gradpsi
 REAL, DIMENSION(ntot), INTENT(INOUT) :: psi,divB
 REAL, DIMENSION(ndim,ntot), INTENT(IN) :: x
 REAL, DIMENSION(ntot), INTENT(IN) :: hh, pmass
 REAL, DIMENSION(ntot) :: psiin, rho, dpsidt
 REAL, DIMENSION(ndimV,ntot) :: Bevolin
 INTEGER :: i,j,istep,indexmax,nsubsteps
 REAL :: dtsub,hdt,vsig2substep,vsigsubstep
 LOGICAL :: nopredictor, docorrector
!
!--allow for tracing flow
!      
 IF (trace) WRITE(iprint,*) ' Entering subroutine substep_divB'
!
!--if nsubstepsin = 0 then just does a predictor step, followed by a corrector step
!
 nsubsteps = 2*nsubstepsin + 1
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

 hdt = 0.5*dtsub
 !
 !--work out whether or not to do predictor and corrector steps
 !
 nopredictor = (icall.eq.1) .or. (istep.ne.nsubsteps)
 docorrector = (icall.eq.2) .or. (istep.ne.nsubsteps)

 if (.not.nopredictor) then
    !
    !--Mid-point Predictor step
    !
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
 endif
 
 if (icall.eq.1 .and. istep.eq.nsubsteps) then
    !
    !--return after predictor step
    !
    return
 elseif (.not.(icall.eq.2 .and.istep.eq.1)) then
    !
    !--or else calculate grad psi and div B for corrector step
    !
    CALL get_divBgradpsi(divB,gradpsi,Bevol,psi,x,hh,pmass,rho,npart,ntot)
    print*,' max divB = ',maxval(divB(1:npart)), &
            rho(maxloc(divB(1:npart))),maxloc(divB(1:npart))

    DO i=1,npart
       dpsidt(i) = -vsig2substep*divB(i) - psidecayfact*psi(i)*vsigsubstep/hh(i)
    ENDDO
 endif
 
 if (.not.(icall.eq.1 .and. istep.eq.1)) then
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
 endif

 ENDDO substeploop
 print*,'finished substeps'
 read*
!
!--return appropriate quantity (B or B/rho)
!
 
 IF (trace) WRITE (iprint,*) ' Exiting subroutine substep_divB'
      
 RETURN
END SUBROUTINE substep_divB
