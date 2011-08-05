!!------------------------------------------------------------------------
!! This subroutine is a cut down version of the 
!! predictor-corrector scheme to be used for sub-stepping
!! (operator splitting) of the hyperbolic divergence correction
!!
!! The predictor-corrector scheme is cut into 2*n+1 steps
!! such that the centred rates call is common to both the main (hydro)
!! stepping and the substepping
!!
!! ie. n=0 does the following
!!     icall=1 : predictor->
!!     (return to main step to get total rates, including divB, grad psi)
!!     icall=2 : ->corrector
!!
!! n=1 does:
!!     icall=1 : predictor->(div B/grad psi)->corrector->predictor->
!!     (return to main step to get total rates)
!!     icall=2 : ->corrector->predictor->(div B/grad psi)->corrector
!!
!! n=2 does:
!!     icall=1 : pred->(divB/gradpsi)->corr->pred->(divBgradpsi)->corr->pred->
!!     (main rates)
!!     icall=2 : ->corr->pred->(divB/gradpsi)->corr->pred->(divBgradpsi)->corr
!!
!! and so on.
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
 REAL, DIMENSION(ntot) :: psiin, rho
 REAL, DIMENSION(ndimV,ntot) :: Bevolin, Bfield
 REAL :: dpsidti
 INTEGER :: i,j,istep,indexmax,nsubsteps,nloops
 REAL :: dtsub,hdt,vsig2substep,vsigsubstep
 LOGICAL :: dopredictor, docorrector
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
 
 print*,'Substepping: dtfull = ',dtfull,' dtsub = ',dtsub,' n = ',nsubstepsin
 print*,' nsubsteps = ',nsubsteps, ' vsigmax = ',vsigsubstep,sqrt(vsig2max)
 print*,' maxdivB = ',maxval(divB(1:npart)), 'npart = ',npart,ntot
 
 psiin = psi
 Bevolin = Bevol 

 nloops = nsubstepsin + 1

 substeploop: DO istep = 1,nloops

 hdt = 0.5*dtsub
 !
 !--work out whether or not to do predictor and corrector steps
 !
 dopredictor = .true.
 docorrector = .true.
 
 if ((icall.eq.2) .and. (istep.eq.1)) dopredictor = .false.
 if ((icall.eq.1) .and. (istep.eq.nloops)) docorrector = .false.

 if (dopredictor) then
    write(iprint,*) ' substep: predictor'
    !
    !--Mid-point Predictor step
    !
    DO i=1,npart
       dpsidti = -vsig2substep*divB(i) - psidecayfact*psi(i)*vsigsubstep/hh(i)
       IF (itype(i).EQ.1) THEN    ! fixed particles
          Bevol(:,i) = Bevolin(:,i)     
          psi(i) = psiin(i)
       ELSE
          Bevol(:,i) = Bevolin(:,i) + hdt*gradpsi(:,i)
          psi(i) = psiin(i) + hdt*dpsidti  
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
 
 if (.not.docorrector) then
    !
    !--return after predictor step, with predicted values of Bevol and psi
    !
    return
 elseif (dopredictor) then
    !
    !--or else calculate grad psi and div B for corrector step
    !
    DO i=1,ntot
       IF (imhd.ge.11) THEN
          Bfield(:,i) = Bevol(:,i) 
       ELSE
          Bfield(:,i) = Bevol(:,i)*rho(i)
       ENDIF
    ENDDO
    CALL get_divBgradpsi(divB,gradpsi,Bfield,psi,x,hh,pmass,rho,npart,ntot)
    print*,' max divB = ',maxval(divB(1:npart)), &
            rho(maxloc(divB(1:npart))),maxloc(divB(1:npart))

 endif
 
 if (docorrector) then
    write(iprint,*) 'substep: corrector'
    !
    !--Mid-point Corrector step
    !


    DO i=1,npart
       dpsidti = -vsig2substep*divB(i) - psidecayfact*psi(i)*vsigsubstep/hh(i)
       IF (itype(i).EQ.1) THEN
          Bevol(:,i) = Bevolin(:,i)
          psi(i) = psiin(i)
       ELSE
          Bevol(:,i) = Bevolin(:,i) + hdt*gradpsi(:,i)
          psi(i) = psiin(i) + hdt*dpsidti
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
 IF (trace) WRITE (iprint,*) ' Exiting subroutine substep_divB'
      
 RETURN
END SUBROUTINE substep_divB
