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

!!--------------------------------------------------------------------
!! Calculate conserved quantities etc and write to .ev file
!!--------------------------------------------------------------------
  
subroutine evwrite(t,etot,momtot)
 use dimen_mhd, only:ndim,ndimV
 use debug, only:trace
 use loguns, only:iprint,ievfile
 
 use derivB
 use options
 use part
 use rates, only:force,potengrav
 use fmagarray
 use timestep, only:dt
 use utils,    only:cross_product3D,minmaxave
 use externf,  only:external_potentials
!
!--define local variables
!
 implicit none
 integer :: i
 real, intent(in) :: t
 real, intent(out) :: etot,momtot
 real :: ekin,etherm,emag,epot
 real :: pmassi,rhoi,epoti!,rr
! real, dimension(ndim) :: rhat
 real, dimension(ndimV) :: veli,mom,ang,angi
! real :: betai,alphai,alphatstarav,betatstarav
!
!--mhd
!
 real, dimension(ndimB) :: Bi,Brhoi,fluxtot
 real :: B2i,Bmagi
 real :: fluxtotmag,crosshel
 real :: betamhdi,betamhdmin,betamhdav
 real :: fdotBi,fdotBmax,fdotBav
 real :: forcemagi,force_erri,force_err_max,force_err_av
 real :: divBi,divBav,divBmax,divBtot
 real :: omegamhdi,omegamhdav,omegamhdmax
 real :: fracdivBok
 real, parameter :: omegtol = 1.E-2
 real :: fmagabs,rhomax,rhomean,rhomin
 real :: angtot,ekiny,emagp
!
!--one fluid dust
!
 real :: totmassgas,totmassdust,dtgi,dustfraci,dterm,ekindeltav

!
!--allow for tracing flow
!      
 if (trace) write(iprint,*) ' Entering subroutine evwrite'
       
 ekin = 0.0
 etherm = 0.0
 emag = 0.0
 emagp = 0.
 etot = 0.0
 if (igravity.ne.0) then
    epot = potengrav
 else
    epot = 0.
 endif
 mom(:) = 0.0
 momtot = 0.0
 ang(:) = 0.
 ekiny = 0.
 totmassdust = 0.
 totmassgas  = 0.
! alphatstarav = 0.
! betatstarav = 0.
!
!--mhd parameters
!     
 if (imhd.ne.0) then
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
 endif 

!
!--should really recalculate the thermal energy from the total energy here
!  (otherwise uu is from the half time step and same with Bfield)
! 
! CALL conservative2primitive
      
 do i=1,npart

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

    if (idust.eq.1 .or. idust.eq.3 .or. idust.eq.4) then
       dustfraci  = dustfrac(i)
       dtgi  = dustfraci/(1. - dustfraci)
       dterm = (1. - dustfraci)
       ekindeltav = 0.5*pmassi*dtgi*dterm**2*dot_product(deltav(:,i),deltav(:,i))
       ekin = ekin + ekindeltav
       ekiny = ekiny + ekindeltav
       etherm = etherm + pmassi*uu(i)*dterm
       totmassgas  = totmassgas  + pmassi*dterm
       totmassdust = totmassdust + pmassi*dtgi*dterm
    else
       if (ndim.ge.2) ekiny = ekiny + 0.5*pmassi*vel(1,i)*vel(1,i)
       etherm = etherm + pmassi*uu(i)
    endif
!
!--potential energy from external forces
!    
    call external_potentials(iexternal_force,x(:,i),epoti,ndim)
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
    if (imhd.ne.0) then
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
       if (B2i.LT.tiny(B2i)) then
          betamhdi = 0.
       else 
          betamhdi = pr(i)/(0.5*B2i)     
       endif
       betamhdav = betamhdav + betamhdi
       if (betamhdi.LT.betamhdmin) betamhdmin = betamhdi
!
!--Maximum divergence of B
!  
       if (divBi.GT.divBmax) divBmax = divBi
       divBav = divBav + divBi
!
!--volume integral of div B (int B.dS)
!
       divBtot = divBtot + pmassi*divBi/rhoi
!
!--Max component of magnetic force in the direction of B (should be zero)
!
       fmagabs = SQRT(DOT_PRODUCT(fmag(:,i),fmag(:,i)))
       if (fmagabs.GT.1.e-8 .and. Bmagi.gt.1.e-8) then
          fdotBi = ABS(DOT_PRODUCT(fmag(:,i),Bi(:)))/(fmagabs*Bmagi)   
       else
          FdotBi = 0.
       endif
       fdotBav = fdotBav + fdotBi
       if (fdotBi.GT.fdotBmax) fdotBmax = fdotBi  
!
!--Compute total error in the force due to the B(div B) term
!  only slight worry with this is that fmag is calculated in rates, whilst
!  B has been evolved a bit further since then. A possible solution is to
!  evaluate these quantities just after the call to rates.
!       
       if (forcemagi.GT.1.e-8 .AND. Bmagi.GT.1e-8) then
          force_erri = ABS(DOT_PRODUCT(fmag(:,i),Bi(:)))/(forcemagi*Bmagi)
       else
          force_erri = 0.
       endif
       force_err_av = force_err_av + force_erri
       if (force_erri.GT.force_err_max) force_err_max = force_erri
!
!--|div B| x smoothing length / |B| (see e.g. Cerqueira and Gouveia del Pino 1999) 
!  this quantity should be less than ~0.01.
!
       if (Bmagi.lt.1e-8) then
          omegamhdi = 0.
       else
          omegamhdi = divBi*hh(i)/Bmagi     
       endif    
       if (omegamhdi.LT.omegtol) fracdivBok = fracdivBok + 1.
       if (omegamhdi.GT.omegamhdmax) omegamhdmax = omegamhdi
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

    endif

 enddo
 
 etot = ekin + emag + epot
 if (iprterm.ge.0 .or. iprterm.lt.-1) etot = etot + etherm
 momtot = sqrt(dot_product(mom,mom))
 call minmaxave(rho(1:npart),rhomin,rhomax,rhomean,npart)
 angtot = sqrt(dot_product(ang,ang))
!
!--write line to .ev file
!     
 if (imhd.ne.0) then      

    fluxtotmag = SQRT(DOT_PRODUCT(fluxtot,fluxtot))
    betamhdav = betamhdav/FLOAT(npart)
    fracdivBok = 100.*fracdivBok/FLOAT(npart)
    omegamhdav = omegamhdav/FLOAT(npart)
    divBav = divBav/FLOAT(npart)
    fdotBav = fdotBav/FLOAT(npart)
    force_err_av = force_err_av/FLOAT(npart)

!!    print*,'t=',t,' emag =',emag,' etot = ',etot, 'ekin = ',ekin,' etherm = ',etherm

    write(ievfile,30) t,ekin,etherm,emag,epot,etot,momtot,angtot,rhomax,rhomean,dt, &
          emagp,crosshel,betamhdmin,betamhdav,  &
          divBav,divBmax,divBtot,     &
          fdotBav,FdotBmax,force_err_av,force_err_max,   &
          omegamhdav,omegamhdmax,fracdivBok
30  format(24(1pe18.10,1x),1pe8.2)
      
 else
    !alphatstarav = alphatstarav/float(npart)
    !betatstarav = betatstarav/float(npart)
   !! print*,'t=',t,' emag =',emag,' etot = ',etot, 'ekin = ',ekin,' etherm = ',etherm

    if (idust.eq.1) then
       write(ievfile,40) t,ekin,etherm,emag,epot,etot,momtot,angtot,rhomax,rhomean,dt,ekiny,&
                         totmassgas,totmassdust    
    else
       write(ievfile,40) t,ekin,etherm,emag,epot,etot,momtot,angtot,rhomax,rhomean,dt,ekiny
    endif
40  format(24(1pe18.10,1x))        

 endif

!
!--flush the buffer so that the line is written to the file immediately
!
 call flush(ievfile)
 
 return
end subroutine evwrite
