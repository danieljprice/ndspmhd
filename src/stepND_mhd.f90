!------------------------------------------------------------------------------!
! NDSPMHD: A Smoothed Particle (Magneto)Hydrodynamics code for (astrophysical) !
! fluid dynamics simulations in 1, 2 and 3 spatial dimensions.                 !
!                                                                              !
! (c) 2002-2014 Daniel Price                                                   !
!                                                                              !
! http://users.monash.edu.au/~dprice/ndspmhd                                   !
! daniel.price@monash.edu -or- dprice@cantab.net (forwards to current address) !
!                                                                              !
!  NDSPMHD comes with ABSOLUTELY NO WARRANTY.                                  !
!  This is free software; and you are welcome to redistribute                  !
!  it under the terms of the GNU General Public License                        !
!  (see LICENSE file for details) and the provision that                       !
!  this notice remains intact. If you modify this file, please                 !
!  note section 2a) of the GPLv2 states that:                                  !
!                                                                              !
!  a) You must cause the modified files to carry prominent notices             !
!     stating that you changed the files and the date of any change.           !
!                                                                              !
!  ChangeLog:                                                                  !
!------------------------------------------------------------------------------!

!!------------------------------------------------------------------------
!! Computes one timestep
!! Change this subroutine to change the timestepping algorithm
!!
!! This subroutine is a slightly modified version of the 
!! predictor-corrector scheme. Algorithm goes as follows:
!!
!! Predictor:
!!
!! v^{1/2}   = v^0 + dt/2*f^(-1/2)
!! x^{1/2}   = x^0 + dt/2*v^{1/2}   (note used updated v)
!! rho^{1/2} = rho^0 + dt/2*drhodt^{-1/2}
!!
!!  --> calculate f^{1/2}, drhodt^{1/2} using x^{1/2} and v^{1/2}
!!
!! Corrector:
!!
!! v^*  = v^0 + dt/2*f^{1/2}
!! x^*   = x^0 + dt/2*v^*          (note uses updated v)
!! rho^* = rho^0 + dt/2*drhodt^*
!!
!! v^1   = 2*v^* - v^0     = v^0 + dt*f^{1/2}
!! x^1   = 2*x^* - x^0     = x^0 + dt*v^*
!! rho^1 = 2*rho^* - rho^0 = rho^0 + dt*drhodt^{1/2}
!!
!! Energy, smoothing length, alpha and magnetic field follow density
!!
!! Gives good results (and good stability) on wave and shock-type problems
!! with a reasonable courant number (C_cour = 0.4 seems to be enough).
!! On these problems seems to do much better than leapfrog/symplectic,
!! mainly related to the derivatives which depend on velocity (for which the
!! other two methods are not reversible/symplectic).
!!
!! However, results for disks and orbits are slightly better with leapfrog
!!
!!---------------------------------------------------------------------------
         
SUBROUTINE step
 USE dimen_mhd
 USE debug
 USE loguns
 
 USE bound
 USE derivB
 USE eos
 USE hterms
 USE options
 USE part
 USE part_in
 USE rates
 USE timestep
 USE setup_params
 USE xsph
 USE kernels, only:eps
 USE utils, only:cross_product3D
!
!--define local variables
!
 IMPLICIT NONE
 INTEGER :: i
 REAL :: hdt
 REAL :: dtrhoi
 real, dimension(ndimV) :: vcrossB
!
!--allow for tracing flow
!      
 IF (trace) WRITE(iprint,*) ' Entering subroutine step'
!
!--Mid-point Predictor step
!      
 hdt = 0.5*dt

 DO i=1,npart
    xin(:,i) = x(:,i)
    velin(:,i) = vel(:,i)
    Bevolin(:,i) = Bevol(:,i)
    rhoin(i) = rho(i)
    hhin(i) = hh(i)
    enin(i) = en(i)
    alphain(:,i) = alpha(:,i)
    psiin(i) = psi(i)
 ENDDO         
   
 DO i=1,npart
    IF (itype(i).EQ.itypebnd .OR. itype(i).EQ.itypebnd2) THEN        ! fixed particles
       vel(:,i) = velin(:,i)
       IF (icty.GE.1) rho(i) = rhoin(i)
       if (imhd.lt.0) then
          call cross_product3D(vel(:,i),Bconst(:),vcrossB)
          Bevol(:,i) = Bevolin(:,i) + hdt*vcrossB(:)
       else
          Bevol(:,i) = Bevolin(:,i)             
       endif
       IF (iener.NE.0) en(i) = enin(i)
       hh(i) = hhin(i)            
       x(:,i) = xin(:,i) + hdt*(vel(1:ndim,i) + xsphfac*xsphterm(1:ndim,i))
       alpha(:,i) = alphain(:,i)
       psi(i) = psiin(i)
    ELSE
       vel(:,i) = (velin(:,i) + hdt*force(:,i))/(1.+damp)
       !--vertical damping
       if (ndim.ge.3) vel(3,i) = vel(3,i)/(1.+dampz)
       !--radial damping
       if (dampr.gt.tiny(dampr) .and.ndim.ge.2) call radialdamping(vel(1:2,i),x(1:2,i),dampr)    
       x(:,i) = xin(:,i) + hdt*(vel(1:ndim,i) + xsphfac*xsphterm(1:ndim,i))

       IF (imhd.NE.0) THEN
          IF (nsubsteps_divB.LT.0) THEN
             Bevol(:,i) = Bevolin(:,i) + hdt*gradpsi(:,i)
             psi(i) = (psiin(i) + hdt*dpsidt(i))/(1.+eps)
          ENDIF
          Bevol(:,i) = Bevol(:,i) + hdt*dBevoldt(:,i)
       ENDIF
       rho(i) = rhoin(i) + hdt*drhodt(i)
       IF (ihvar.EQ.1) THEN
!           hh(i) = hfact(pmass(i)/rho(i))**dndim        ! my version
          hh(i) = hhin(i)*(rhoin(i)/rho(i))**dndim    ! Joe's           
       ELSEIF (ihvar.EQ.2 .OR. ihvar.EQ.3) THEN
          hh(i) = hhin(i) + hdt*dhdt(i)
       ENDIF
       IF (iener.NE.0) en(i) = enin(i) + hdt*dendt(i)
       IF (ANY(iavlim.NE.0)) alpha(:,i) = alphain(:,i) + hdt*daldt(:,i)
    ENDIF

 ENDDO
 
!
!--calculate all derivatives
!
 call derivs

 dtrho = huge(dtrho)
 DO i=1,npart
    dtrhoi = abs(rho(i)/(drhodt(i) + 1.e-8))
    dtrho = min(dtrho,dtrhoi)
 ENDDO
 if (C_rho*dtrho/dtcourant .lt. C_cour) then
    write(iprint,*) 'dtrho equiv courant number = ',C_rho*dtrho/dtcourant
 endif
!
!--Mid-point Corrector step
!
 DO i=1,npart
    IF (itype(i).EQ.itypebnd .or. itype(i).EQ.itypebnd2) THEN
       vel(:,i) = velin(:,i)
       IF (icty.GE.1) rho(i) = rhoin(i)
       if (imhd.lt.0) then
          call cross_product3D(vel(:,i),Bconst(:),vcrossB)
          Bevol(:,i) = Bevolin(:,i) + dt*vcrossB(:)
       else
          Bevol(:,i) = Bevolin(:,i)
       endif
       IF (iener.NE.0) en(i) = enin(i)
       x(:,i) = xin(:,i) + dt*(vel(1:ndim,i) + xsphfac*xsphterm(1:ndim,i))
       alpha(:,i) = alphain(:,i)
       hh(i) = hhin(i)
       psi(i) = psiin(i)
    ELSE
       vel(:,i) = (velin(:,i) + dt*force(:,i))/(1.+damp)            
       !--vertical damping
       if (ndim.ge.3) vel(3,i) = vel(3,i)/(1.+dampz)
       !--radial damping
       if (dampr.gt.tiny(dampr) .and.ndim.ge.2) call radialdamping(vel(1:2,i),x(1:2,i),dampr)      
       x(:,i) = xin(:,i) + dt*(vel(1:ndim,i) + xsphfac*xsphterm(1:ndim,i))
       IF (ihvar.EQ.2 .OR. (ihvar.EQ.3 .and. itsdensity.eq.0)) THEN
          hh(i) = hhin(i) + dt*dhdt(i)
          IF (hh(i).LE.0.) THEN
             WRITE(iprint,*) 'step: hh -ve ',i,hh(i)
             CALL quit
          ENDIF
       ENDIF
       IF (icty.GE.1) THEN
          rho(i) = rhoin(i) + dt*drhodt(i)
          IF (rho(i).LE.0.) THEN
             WRITE(iprint,*) 'step: rho -ve ',i,rho(i)
             CALL quit
          ENDIF
       ENDIF
       IF (iener.NE.0) en(i) = enin(i) + dt*dendt(i)
       IF (ANY(iavlim.NE.0)) alpha(:,i) = alphain(:,i) + dt*daldt(:,i)           
       IF (imhd.NE.0) THEN
          IF (nsubsteps_divB.LT.0) THEN
             Bevol(:,i) = Bevolin(:,i) + dt*gradpsi(:,i)
             psi(i) = (psiin(i) + dt*dpsidt(i))/(1.+eps)
          ENDIF          
          Bevol(:,i) = Bevol(:,i) + dt*dBevoldt(:,i)
       ENDIF
    ENDIF 
              
 ENDDO

 IF (ANY(ibound.NE.0)) CALL boundary        ! inflow/outflow/periodic boundary conditions
!
!--set new timestep from courant/forces condition
!
 dt = min(C_force*dtforce,C_cour*dtcourant,C_force*dtdrag,C_force*dtvisc)

 IF (trace) WRITE (iprint,*) ' Exiting subroutine step'
      
 RETURN
END SUBROUTINE step

subroutine radialdamping(vxy,xy,dampr)
 implicit none
 real, intent(in), dimension(2) :: xy
 real, intent(inout), dimension(2) :: vxy
 real, intent(in) :: dampr
 real, dimension(2) :: rhat
 real :: vrad,vradnew
 
 rhat(:) = xy(:)/sqrt(dot_product(xy,xy))
 vrad = dot_product(vxy,rhat)
 vradnew = vrad/(1. + dampr)
 vxy(:) = vxy(:) + (vradnew - vrad)*rhat(:)
 
 return
end subroutine radialdamping
