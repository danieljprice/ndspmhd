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
!!     (return to main step to get total rates, including divBi, grad psi)
!!     icall=2 : ->corrector
!!
!! n=1 does:
!!     icall=1 : predictor->(div B/grad psi)->corrector->predictor->
!!     (return to main step to get total rates)
!!     icall=2 : ->corrector->predictor->(div B/grad psi)->corrector
!!
!! n=2 does:
!!     icall=1 : pred->(divBi/gradpsi)->corr->pred->(divBigradpsi)->corr->pred->
!!     (main rates)
!!     icall=2 : ->corr->pred->(divBi/gradpsi)->corr->pred->(divBigradpsi)->corr
!!
!! and so on.
!!
!!---------------------------------------------------------------------------

subroutine substep_divb(icall,dtfull,nsubstepsin,bevol,psi,divbi,gradpsi, &
                        x,hh,pmass,rho,itype,npart,ntot)
 use dimen_mhd
 use debug
 use loguns
 use bound
 use options
 use timestep
 use derivb
!
!--define local variables
!
 implicit none
 integer, intent(in) :: icall,nsubstepsin,npart,ntot
 integer, dimension(ntot), intent(in) :: itype
 real, intent(in) :: dtfull
 real, dimension(ndimV,ntot), intent(inout) :: Bevol,gradpsi
 real, dimension(ntot), intent(inout) :: psi,divBi
 real, dimension(ndim,ntot), intent(in) :: x
 real, dimension(ntot), intent(in) :: hh, pmass, rho
 real, dimension(ntot) :: psiin
 real, dimension(ndimV,ntot) :: Bevolin, Bfield
 real :: dpsidti, alphasub
 integer :: i,j,istep,nsubsteps,nloops,isteptot
 integer, dimension(1) :: imaxpart
 real :: dtsub,hdt,vsig2substep,vsigsubstep, maxdivbi,crap1,crap2
 logical :: dopredictor, docorrector
 save isteptot
!
!--allow for tracing flow
!      
 IF (trace) WRITE(iprint,*) ' Entering subroutine substep_divBi'
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
 imaxpart = maxloc(abs(divBi(1:npart)))
 print*,' maxdivBi = ',divBi(imaxpart),' on part',imaxpart
 !!print*,'npart = ',npart,ntot
 
 print*,'min h = ',minval(hh(1:npart))
 print*,'dt should be = ',minval(hh(1:npart))/vsigsubstep
 
 psiin = psi
 Bevolin = Bevol 
! divBi(:) = 0.
! gradpsi(:,:) = 0.
 alphasub = 0.1

 nloops = nsubstepsin + 1

 substeploop: DO istep = 1,nloops


!! psidecayfact = psidecayfact - 0.01
 hdt = 0.5*dtsub
 !
 !--work out whether or not to do predictor and corrector steps
 !
 dopredictor = .true.
 docorrector = .true.
 
! if ((icall.eq.2) .and. (istep.eq.1)) dopredictor = .false.
! if ((icall.eq.1) .and. (istep.eq.nloops)) docorrector = .false.

 if (dopredictor) then
    write(iprint,*) ' substep: predictor ',Bevol(:,imaxpart) !!,Bevolin(:,2) !!,gradpsi(:,2),hdt
    !
    !--Mid-point Predictor step
    !
    DO i=1,npart
       dpsidti = -vsig2substep*divBi(i) - psidecayfact*psi(i)*vsigsubstep/hh(i)
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
    !IF (ANY(ibound.GE.2)) THEN
    !  DO i=npart+1,ntot
    !      j = ireal(i)
    !      psi(i) = psi(j)
    !      Bevol(:,i) = Bevol(:,j)
    !   ENDDO
    !ENDIF
 endif
 
 if (.not.docorrector) then
    !
    !--return after predictor step, with predicted values of Bevol and psi
    !
    print*,'substep: returning after predictor...'
    return
 elseif (dopredictor) then
    !
    !--or else calculate grad psi and div B for corrector step
    !
    DO i=1,npart
       IF (imhd.ge.11) THEN
          Bfield(:,i) = Bevol(:,i) 
       ELSE
          Bfield(:,i) = Bevol(:,i)*rho(i)
       ENDIF
    ENDDO
    IF (ANY(ibound.GE.2)) THEN
       DO i=npart+1,ntot
          j = ireal(i)
          Bevol(:,i) = Bevol(:,j)
          Bfield(:,i) = Bfield(:,j)
          psi(i) = psi(j)
       ENDDO
    ENDIF

    !maxdivBi = 0.
    !imaxpart = 0
    !do i=1,npart
    !   if (itype(i).eq.0) then
    !      if (divBi(i).gt.maxdivBi) then
    !         maxdivBi = divBi(i)
    !         imaxpart = i
    !      endif
    !   endif
    !enddo
    !print*,'before mini rates: max divBi = ',maxdivBi,imaxpart

    CALL get_divBgradpsi(divBi,gradpsi,Bfield,psi,x,hh,pmass,rho,npart,ntot,alphasub)
    
    maxdivBi = -1.0
    imaxpart = 1
    do i=1,npart
       if (itype(i).eq.0) then
          if (abs(divBi(i)).gt.maxdivBi) then
             maxdivBi = abs(divBi(i))
             imaxpart = i
          endif
       endif
       divB(i) = divBi(i)
    enddo
    do i=npart+1,ntot
       divB(i) = divB(ireal(i))
    enddo
    print*,'after mini rates: max divBi = ',maxdivBi,imaxpart,istep,real(isteptot)
    isteptot = isteptot + 1
    call evwrite(real(isteptot),crap1,crap2)
    call output(psidecayfact,1)
    !!print*,psi(2),Bevol(:,2),divBi(2),gradpsi(:,2),rho(2)
    !endif
 endif
 
 if (docorrector) then
    write(iprint,*) 'substep: corrector'
    !
    !--Mid-point Corrector step
    !
    DO i=1,npart
       dpsidti = -vsig2substep*divBi(i) - psidecayfact*psi(i)*vsigsubstep/hh(i)
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

! IF (ANY(ibound.GE.2)) THEN
!    do i=npart+1,ntot
!       j = ireal(i)
!       Bevol(:,i) = Bevol(:,j)
!       Bevolin(:,i) = Bevolin(:,j)
!    enddo
! ENDIF

 print*,'finished substeps',icall
 IF (trace) WRITE (iprint,*) ' Exiting subroutine substep_divBi'
      
 return
end subroutine substep_divb
