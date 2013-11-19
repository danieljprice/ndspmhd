!!--------------------------------------------------------------------
!! Computes one timestep
!! Change this subroutine to change the timestepping algorithm
!!--------------------------------------------------------------------
         
subroutine step (integratorcheck)
 use dimen_mhd
 use debug
 use loguns
 
 use bound
 use eos
 use hterms
 use options
 use part
 use part_in
 use rates
 use timestep
 use setup_params
 use xsph

 use cons2prim, only:conservative2primitive,getv_from_pmom
!
!--define local variables
!
 implicit none
 integer :: i
 real :: hdt
 character (len=*), intent (inout) :: integratorcheck

 if (trim(integratorcheck).eq.'query') then
    integratorcheck = 'mhd_gr'
    return
 endif

!
!--allow for tracing flow
!      
 if (trace) write(iprint,*) ' Entering subroutine step'
!
!--set initial quantities
!
 hdt = 0.5*dt
 do i=1,npart
    xin(:,i) = x(:,i)
    pmomin(:,i) = pmom(:,i)
    Bevolin(:,i) = Bevol(:,i)
    rhoin(i) = rho(i)
    hhin(i) = hh(i)
    enin(i) = en(i)
    alphain(:,i) = alpha(:,i)
 enddo   
!
!--Mid-point Predictor step
!      
 do i=1,npart
    if (itype(i).EQ.itypebnd .or. itype(i).EQ.itypebnd2) then        ! fixed particles
       pmom(:,i) = pmomin(:,i)
       rho(i) = rhoin(i)
       Bevol(:,i) = Bevolin(:,i)     
       if (iener.NE.0) en(i) = enin(i)
       hh(i) = hhin(i)            
       x(:,i) = xin(:,i) + hdt*(vel(1:ndim,i) + xsphfac*xsphterm(1:ndim,i))
       alpha(:,i) = alphain(:,i)
    else
       pmom(:,i) = pmomin(:,i) + hdt*force(:,i)
       if (iener.NE.0) en(i) = enin(i) + hdt*dendt(i)
!--use updated v to change position
       !print*,i,'predictor'
       call getv_from_pmom(xin(:,i),pmom(:,i),vel(:,i),en(i),pr(i),rho(i),dens(i),uu(i))
       x(:,i) = xin(:,i) + hdt*(vel(1:ndim,i) + xsphfac*xsphterm(1:ndim,i))
!--do an 'iteration' of this to be on the safe side
       !call getv_from_pmom(x(:,i),pmom(:,i),vel(:,i),en(i),pr(i),rhoin(i),dens(i),uu(i))
       x(:,i) = xin(:,i) + hdt*(vel(1:ndim,i) + xsphfac*xsphterm(1:ndim,i))
       if (imhd.NE.0) Bevol(:,i) = Bevolin(:,i) + hdt*dBevoldt(:,i)
       if (ihvar.EQ.1) then
!           hh(i) = hfact(pmass(i)/rho(i))**dndim        ! my version
          hh(i) = hhin(i)*(rhoin(i)/rho(i))**dndim                ! Joe's           
       elseif (ihvar.EQ.2 .OR. ihvar.EQ.3) then
          hh(i) = hhin(i) + hdt*dhdt(i)
       endif
       if (icty.GE.1) rho(i) = rhoin(i) + hdt*drhodt(i)
       if (ANY(iavlim.NE.0)) alpha(:,i) = min(alphain(:,i) + hdt*daldt(:,i),1.0)           
    endif

 enddo
!
!--calculate all derivatives
!
 call derivs

 !print*,'corrector'
!
!--Mid-point Corrector step
!
 do i=1,npart
    if (itype(i).EQ.itypebnd .or. itype(i).eq.itypebnd2) then
       pmom(:,i) = pmomin(:,i)
       rho(i) = rhoin(i)
       Bevol(:,i) = Bevolin(:,i)
       if (iener.NE.0) en(i) = enin(i)
       x(:,i) = xin(:,i) + hdt*(vel(1:ndim,i) + xsphfac*xsphterm(1:ndim,i)) 
       alpha(:,i) = alphain(:,i)
       hh(i) = hhin(i)
    else
       pmom(:,i) = pmomin(:,i) + dt*force(:,i)
       if (iener.NE.0) en(i) = enin(i) + dt*dendt(i)
!--use updated v to change position
       call getv_from_pmom(xin(:,i),pmom(:,i),vel(:,i),en(i),pr(i),rho(i),dens(i),uu(i))
       x(:,i) = xin(:,i) + dt*(vel(1:ndim,i) + xsphfac*xsphterm(1:ndim,i))
!--do an 'iteration' of this (with updated x) to be on the safe side
       call getv_from_pmom(x(:,i),pmom(:,i),vel(:,i),en(i),pr(i),rho(i),dens(i),uu(i))
       x(:,i) = xin(:,i) + dt*(vel(1:ndim,i) + xsphfac*xsphterm(1:ndim,i))
       if (ihvar.EQ.2) then
          hh(i) = hhin(i) + dt*dhdt(i)
          if (hh(i).LE.0.) then
             write(iprint,*) 'step: hh -ve ',i,hh(i)
             call quit
          endif
       endif
       if (icty.GE.1) then
          rho(i) = rhoin(i) + dt*drhodt(i)
          if (rho(i).LE.0.) then
             write(iprint,*) 'step: rho -ve ',i,rho(i)
             call quit
          endif
       endif
       if (ANY(iavlim.NE.0)) alpha(:,i) = min(alphain(:,i) + dt*daldt(:,i),1.0)
       if (imhd.NE.0) Bevol(:,i) = Bevolin(:,i) + dt*dBevoldt(:,i)          
    endif
 enddo                 
!
!--update density using a full summation every so often
!         
 if (MOD(nsteps,ndirect).EQ.0) then
    call iterate_density
    do i=1,npart
       rhoin(i) = rho(i)
    enddo
 endif
!
!--if doing divergence correction then do correction to magnetic field
! 
! if (idivBzero.NE.0) call divBcorrect
!
 if (ANY(ibound.NE.0)) call boundary        ! inflow/outflow/periodic boundary conditions
!
!--set new timestep from courant/forces condition
!
 dt = min(C_force*dtforce,C_cour*dtcourant)
!
!--need to calculate velocity from particle momentum for next predictor step
!  (no need to call this from output if called here)
! 
 call conservative2primitive
  
 if (trace) write (iprint,*) ' Exiting subroutine step'
      
 return
end subroutine step
