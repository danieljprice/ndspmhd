!!--------------------------------------------------------------------
!! calculate conserved quantities etc and write to .ev file
!!--------------------------------------------------------------------
  
subroutine evwrite(t,etot,momtot)
 use dimen_mhd
 use debug
 use loguns
 
 use derivb
 use options
 use part
 use rates
 use fmagarray
 
 use grutils
!
!--define local variables
!
 implicit none
 integer :: i
 real, intent(in) :: t
 real, intent(out) :: etot,momtot
 real :: ekin,etherm,emag,epot
 real :: pmassi,rhoi,epoti
 real, dimension(ndimV) :: veli,mom,gdiag
!
!--mhd
!
 real, dimension(ndimV) :: Bi,Brhoi,fluxtot
 real :: B2i,Bmagi
 real :: fluxtotmag,crosshel
 real :: betamhdi,betamhdmin,betamhdmax,betamhdav
 real :: FdotBi,Fdotbmax,Fdotbav
 real :: forcemagi,force_erri,force_err_max,force_err_av
 real :: divBi,divbav,divbmax,divbtot
 real :: omegamhdi,omegamhdav,omegamhdmax
 real :: omegtol,fracdivbok
 real :: fmagabs
!
!--allow for tracing flow
!      
 if (trace) write(iprint,*) ' entering subroutine evwrite'
       
 ekin = 0.0
 etherm = 0.0
 emag = 0.0
 etot = 0.0
 epot = 0.0
 mom(:) = 0.0
 momtot = 0.0
!
!--mhd parameters
!     
 if (imhd.ne.0) then
    betamhdav = 0.
    betamhdmax = 0.
    betamhdmin = 1.e30
    divbmax = 0.
    divbav = 0.
    divbtot = 0.
    fdotbmax = 0.
    fdotbav = 0.
    force_err_max = 0.
    force_err_av = 0.
    omegamhdav = 0.
    omegamhdmax = 0.
    omegtol = 1.e-2
    fracdivbok = 0.
    fluxtot(:) = 0.
    fluxtotmag = 0.
    crosshel = 0.      
 endif 

!
!--should really recalculate the thermal energy from the total energy here
!  (otherwise uu is from the half time step and same with bfield)
! 
 call conservative2primitive
      
 do i=1,npart
    call metric_diag(x(:,i),gdiag,sqrtg(i),ndim,ndimV,igeom)
    
    pmassi = pmass(i)
    rhoi = rho(i)
    veli(:) = vel(:,i)  
    mom(:) = mom(:) + pmassi*pmom(:,i)
    ekin = ekin + 0.5*pmassi*dot_product_gr(veli,veli,gdiag,ndimV)
    etherm = etherm + pmassi*uu(i)
!
!--potential energy from external forces
!    
    call external_potentials(iexternal_force,x(:,i),epoti,ndim)
    epot = epot + pmassi*epoti
!
!--mhd parameters
!
    if (imhd.ne.0) then
       if (imhd.ge.11) then
          Bi(:) = Bevol(:,i)
          Brhoi(:) = Bi(:)/rhoi
       else
          Brhoi(:) = Bevol(:,i)
          Bi(:) = Brhoi(:)*rhoi
       endif
       B2i = dot_product_gr(Bi,Bi,gdiag,ndimV)
       Bmagi = sqrt(B2i)
       forcemagi = sqrt(dot_product_gr(force(:,i),force(:,i),gdiag,ndimV))
       divBi = abs(divB(i))
 
       emag = emag + 0.5*pmassi*b2i/rhoi
!
!--plasma beta minimum/maximum/average
!  
       if (B2i.lt.1.e-5) then
          betamhdi = 0.
       else 
          betamhdi = pr(i)/(0.5*B2i)     
       endif
       betamhdav = betamhdav + betamhdi
       if (betamhdi.gt.betamhdmax) betamhdmax = betamhdi
       if (betamhdi.lt.betamhdmin) betamhdmin = betamhdi
!
!--maximum divergence of b
!  
       if (divBi.gt.divBmax) divBmax = divBi
       divbav = divbav + divBi
!
!--volume integral of div B (int B.ds)
!
       divbtot = divbtot + pmassi*divBi/rhoi
!
!--max component of magnetic force in the direction of B (should be zero)
!
       fmagabs = sqrt(dot_product(fmag(:,i),fmag(:,i)))
       if (fmagabs.gt.1.e-8 .and. Bmagi.gt.1.e-8) then
          fdotBi = abs(dot_product(fmag(:,i),Bi(:)))/(fmagabs*Bmagi)   
       else
          fdotBi = 0.
       endif
       fdotbav = fdotbav + fdotBi
       if (fdotBi.gt.fdotbmax) fdotbmax = fdotBi  
!
!--compute total error in the force due to the b(div B) term
!  only slight worry with this is that fmag is calculated in rates, whilst
!  B has been evolved a Bit further since then. a possible solution is to
!  evaluate these quantities just after the call to rates.
!       
       if (forcemagi.gt.1.e-8 .and. Bmagi.gt.1e-8) then
          force_erri = abs(dot_product_gr(fmag(:,i),Bi(:),gdiag,ndimV))/(forcemagi*Bmagi)
       else
          force_erri = 0.
       endif
       force_err_av = force_err_av + force_erri
       if (force_erri.gt.force_err_max) force_err_max = force_erri
!
!--|div B| x smoothing length / |b| (see e.g. cerqueira and gouveia del pino 1999) 
!  this quantity should be less than ~0.01.
!
       if (Bmagi.lt.1e-8) then
          omegamhdi = 0.
       else
          omegamhdi = divBi*hh(i)/Bmagi     
       endif    
       if (omegamhdi.lt.omegtol) fracdivbok = fracdivbok + 1.
       if (omegamhdi.gt.omegamhdmax) omegamhdmax = omegamhdi
       omegamhdav = omegamhdav + omegamhdi   
!
!--conserved magnetic flux (int B dv)
!
       pmassi = pmass(i)
       fluxtot(:) = fluxtot(:) + pmassi*Brhoi(:)
!
!--conserved cross helicity (int v.B dv)
!
       crosshel = crosshel + pmassi*dot_product_gr(veli,Brhoi,gdiag,ndimV)

    endif

 enddo
 
 etot = etherm + ekin + emag + epot
 momtot = sqrt(dot_product_gr(mom,mom,gdiag,ndimV))

!
!--write line to .ev file
!     
 if (imhd.ne.0) then      

    fluxtotmag = sqrt(dot_product(fluxtot,fluxtot))
    betamhdav = betamhdav/float(npart)
    fracdivbok = 100.*fracdivbok/float(npart)
    omegamhdav = omegamhdav/float(npart)
    divbav = divbav/float(npart)
    fdotbav = fdotbav/float(npart)
    force_err_av = force_err_av/float(npart)

!    print*,'t=',t,' emag =',emag,' etot = ',etot, 'ekin = ',ekin,' etherm = ',etherm

    write(ievfile,30) t,ekin,etherm,emag,etot,momtot,fluxtotmag, &
          crosshel,betamhdmin,betamhdav,betamhdmax,   &
          divbav,divbmax,divbtot,     &
          fdotbav,fdotbmax,force_err_av,force_err_max,   &
          omegamhdav,omegamhdmax,fracdivbok
30  format(20(1pe18.10,1x),1pe8.2)
      
 else

    write(ievfile,40) t,ekin,etherm,emag,etot,momtot
40  format(6(1pe18.10,1x))        

 endif

!
!--flush the buffer so that the line is written to the file immediately
!
! call flush(ievfile)
 
 return
end subroutine evwrite
