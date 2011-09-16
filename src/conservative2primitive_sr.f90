!---------------------------------------------------------------------
! these subroutines calculates the primitive variables from the
! conservative variables and vice versa.
!
! This version is for special relativistic hydro, using either
! the total energy or entropy variables and for
! an ideal gas equation of state
!
! This subroutine is more complicated in the relativistic case since the 
! Lorentz factor is involved everywhere.
!---------------------------------------------------------------------
module cons2prim
 implicit none
 logical, parameter :: specialrelativity = .true.

contains

subroutine conservative2primitive
  use dimen_mhd
  use debug
  use options
  use loguns  
  use bound
  use eos
  use part
  implicit none
  integer :: i,j,nerr,ierr,itsmax
  !real, dimension(ndimV) :: gdiag
  real :: lorentzfac

  if (trace) write(iprint,*) ' Entering subroutine conservative2primitive(sr)'

  nerr = 0
  geom = 'srel'
  itsmax = 0
!
!--this is the procedure for special relativity
!
  do i=1,npart
  
     if (iener.eq.0) then
        call solve_conservative2primitivei(pmom(:,i),polyk,rho(i),gamma, &
             vel(:,i),uu(i),dens(i),pr(i),lorentzfac,iener,ierr)
     elseif (iener.eq.1 .or. iener.eq.3) then
        call solve_conservative2primitivei(pmom(:,i),en(i),rho(i),gamma, &
             vel(:,i),uu(i),dens(i),pr(i),lorentzfac,iener,ierr)
     else
        stop 'energy equation not implemented for special relativity'
     endif
     if (ierr > 0) then
        print*,'conservative2primitive: could not solve rootfinding for particle ',i
        print*,' position = ',x(:,i)
        print*,' pmom=',pmom(:,i)
        print*,' rho*=',rho(i)
        if (iener.gt.0) print*,' en=',en(i)
        stop
     endif
     itsmax = max(abs(ierr),itsmax)
     !
     !--call equation of state to get pressure (needed for source terms)
     !
     call equation_of_state1(pr(i),spsound(i),uu(i),dens(i))
 
  enddo
  if (itsmax.gt.5) print*,' max iterations = ',itsmax
!
!--copy the primitive variables onto the ghost particles
! 
  if (any(ibound.gt.1)) then
     do i=npart+1,ntotal
        j = ireal(i)
        call copy_particle(i,j)
        where(ibound.eq.2)
           pmom(:,i) = -pmom(:,j)
           vel(:,i) = -vel(:,j)
        end where
     enddo
  endif
!
!--make fixed particles exact replicas of their closest particle
!
  if (any(ibound.eq.1)) then
     do i=1,npart
        if (itype(i).eq.itypebnd) then
           j = ireal(i)
           call copy_particle(i,j)
        endif
     enddo
  endif
 
  return
end subroutine conservative2primitive

!---------------------------------------------------------------------
! this subroutine is a cut-down version of conservative2primitive
! which just calculates v from pmom (used in timestepping)
!---------------------------------------------------------------------
subroutine getv_from_pmom(xi,pmomi,veli,eni,pri,rhoi,densi,uui)
 use dimen_mhd
 use options, only:geom,iener
 use eos, only:gamma,polyk
 real, dimension(ndim), intent(in) :: xi
 real, dimension(:), intent(in) :: pmomi
 real, intent(in) :: eni,pri,rhoi,uui
 real, intent(inout) :: densi
 real, dimension(:), intent(out) :: veli
 real :: prtemp,uutemp,lfac
 integer :: ierr
 
 uutemp = uui
 prtemp = pri
 select case(iener)
 case(0)
    call solve_conservative2primitivei(pmomi,polyk,rhoi,gamma,veli,uutemp,densi,prtemp,lfac,iener,ierr)
 case(1,3)
    call solve_conservative2primitivei(pmomi,eni,rhoi,gamma,veli,uutemp,densi,prtemp,lfac,iener,ierr) 
 case default
    stop 'invalid iener for special relativity'
 end select
 if (ierr > 0) then
   print*,'getv_from_pmom: could not solve root finding '
   print*,' position = ',xi(:)
   print*,' pmom=',pmomi(:)
   print*,' rho*=',rhoi
   if (iener.gt.0) print*,' en=',eni
   stop
endif


end subroutine getv_from_pmom

!---------------------------------------------------------------------
! this subroutine is called after setting up the initial conditions
! to set the initial values of the conserved variables
!---------------------------------------------------------------------
subroutine primitive2conservative
  use dimen_mhd
  use debug
  use loguns
  
  use bound, only:ireal
  use eos
  use options
  use part
  use setup_params
  use timestep
  use utils, only:minmaxave
  implicit none
  integer :: i,j,iktemp,ierr
  real :: lorentzfactor, v2i, hmin, hmax, hav, polyki, gam1

  if (trace) write(iprint,*) ' Entering subroutine primitive2conservative (sr)'

  call test_conservative2primitive
!
!--set initial h and rho (conservative) from dens (primitive) as given in setup/ data read
!
  do i=1,npart
  
     call solve_primitive2conservativei(pmom(:,i),en(i),rho(i),gamma, &
                                        vel(:,i),uu(i),dens(i),pr(i),spsound(i),iener,ierr)

     hh(i) = hfact*(pmass(i)/rho(i))**dndim
!
!--also work out what polyk should be if using iener = 0
!    
     if (iener.eq.0) then
        gam1 = gamma - 1.
        if (gamma.le.1.00001) then
           polyki = 2./3.*uu(i)
        else
           polyki = gam1*uu(i)/dens(i)**(gam1)
        endif
        if (abs(polyki-polyk)/polyk.gt.1.e-8) then
           write(iprint,*) 'NOTE: setting polyk = ',polyki,' (infile says ',polyk,')'
           polyk = polyki
        endif
     endif
  enddo
  if (ihvar.le.0) then
     call minmaxave(hh(1:npart),hmin,hmax,hav,npart)
     hh(1:ntotal) = hav
  endif
!
!--overwrite this with a direct summation
!  
  if (icty.eq.0 .or. ndirect.lt.100000) then
     write(iprint,*) 'Calculating initial density...' 
     if (ANY(ibound.GT.1)) call set_ghost_particles
     call set_linklist
     iktemp = ikernav
!        ikernav = 3            ! consistent with h for first density evaluation
     call iterate_density       ! evaluate density by direct summation
     ikernav = iktemp  
!!        hh(1:npart) = hfact*(pmass(1:npart)/rho(1:npart))**dndim
     if (ihvar.le.0) then
        call minmaxave(hh(1:npart),hmin,hmax,hav,npart)
        hh(1:npart) = hav
     endif
  endif
!
!--check for errors with memory allocation
!
  if (geom(1:4).ne.'cart' .and. .not.allocated(sourceterms)) then
     write(iprint,*) 'ERROR: non-cartesian geometry but source terms not allocated'
     stop    
  endif
!
!--this is the procedure for special relativity
!
  do i=1,npart
     !call metric_diag(x(:,i),gdiag,sqrtg(i),ndim,ndimV,geom)

     v2i = dot_product(vel(:,i),vel(:,i))
     if (v2i.gt.1.) then 
        print*,'here v2i > 1'
        stop
     endif
     lorentzfactor = 1./sqrt(1.-v2i)
     !dens(i) = rho(i)/lorentzfactor ! replace primitive variable after summation

     call solve_primitive2conservativei(pmom(:,i),en(i),rho(i),gamma, &
                                        vel(:,i),uu(i),dens(i),pr(i),spsound(i),iener,ierr)
  enddo
!
!--copy the conservative variables onto the ghost particles
!  
  if (any(ibound.gt.1)) then
     do i=npart+1,ntotal
        j = ireal(i)
        call copy_particle(i,j)
        where(ibound.eq.2)
           pmom(:,i) = -pmom(:,j)
           vel(:,i) = -vel(:,j)
        end where
     enddo
  endif
!
!--make fixed particles exact replicas of their closest particle
!
  if (any(ibound.eq.1)) then
     do i=1,npart
        if (itype(i).eq.itypebnd) then
           j = ireal(i)
           call copy_particle(i,j)
        endif
     enddo
  endif
!
!--call rates to get initial timesteps, div B etc
!
  call get_rates     

  return  
end subroutine primitive2conservative

!
! This is the internal workings for the above routine.
! Defines what the conservative variables are in terms of the primitive ones.
!
subroutine solve_primitive2conservativei(pmomi,eni,rhoi,gamma,veli,uui,densi,pri,spsoundi,iener,ierr)
  use eos, only:equation_of_state1
  implicit none
  real, intent(out), dimension(:) :: pmomi
  real, intent(out) :: eni,rhoi,pri,spsoundi
  real, intent(in) :: gamma
  real, intent(in) :: densi
  real, intent(inout) :: uui
  real, intent(in), dimension(:) :: veli
  integer, intent(in) :: iener ! choice of energy variable
  integer, intent(out) :: ierr
  real :: v2i,lorentzfactor,pmommag,vmag

  ierr = 0
  v2i = dot_product(veli(:),veli(:))
  if (v2i.ge.1.) then
     print*,' ERROR! v2 = ',v2i,' in special relativistic setup - means backwards in time!'
     print*,' vel = ',veli(:)
     ierr = 1
     stop
  endif
  lorentzfactor = 1./sqrt(1.-v2i)
  rhoi = densi*lorentzfactor

  !
  !--call equation of state to get pressure (needed to get momentum and energy variables)
  !  (strictly all this is unnecessary until after the density sum)
  !
  call equation_of_state1(pri,spsoundi,uui,densi)

  if (v2i.gt.tiny(v2i)) then
     vmag = sqrt(v2i)
     pmommag = lorentzfactor*vmag*(1. + uui + pri/densi)
     pmomi(:) = pmommag*(veli(:)/vmag)
  else
     pmomi(:) = 0.
  endif

  select case(iener)
  case(1,0)
     eni = pri/densi**gamma
  case(2)
     eni = uui
  case(3)
     eni = lorentzfactor*(1.+ uui + pri/densi) - pri/rhoi
  case default
     ierr = 2
     stop 'invalid iener in primitive2conservative(sr)'
  end select

  return
end subroutine solve_primitive2conservativei
!
!--this version solves for the enthalpy from the conserved momentum,
!  conserved density and either the total energy or an entropy variable
!  (only for ideal gas EOS)
!
subroutine solve_conservative2primitivei(pmomi,eni,rhoi,gamma,veli,uui,densi,pri,lorentzfactor,iener,ierr)
  implicit none
  real, intent(in), dimension(:) :: pmomi
  real, intent(in) :: eni,rhoi,gamma
  real, intent(inout) :: uui,densi,pri
  real, intent(out), dimension(size(pmomi)) :: veli
  integer, intent(in) :: iener ! choice of energy variable
  integer, intent(out) :: ierr
  integer :: its
  integer, parameter :: itsmax = 500
  integer, parameter :: itsmaxnr = 100
  real, parameter :: tol = 1.e-8
  logical :: converged
  real :: h,hnew,lorentzfactor
  real :: fh,dfh,gam1,pmommag,vmag,pmom2 !,grad
  real :: hmin,hmax,fhmin,dfhmin,hin
 
  pmom2 = dot_product(pmomi(:),pmomi(:))
  pmommag = sqrt(pmom2)
!
!--starting guess is the last stored values of the primitive variables  
!
  h = 1. + uui + pri/densi
  hin = h
  its = 0
  ierr = 0
  gam1 = gamma - 1. ! this is the adiabatic gamma
  converged = .false.
  hmin = 1.0 + tiny(hmin)
  hmax = 1.e8
  call eval_funch(hmin,fhmin,dfhmin)
  call eval_funch(hmax,fh,dfh)
  if (fhmin*fh .gt. 0) then
     print*,'fhmin,max = ',fhmin,fh
!     h = hmin
!     do its=1,500
!        h = h + 0.01*h
!        call eval_funch(h,fh,dfh)
!        write(8,*) h,fh
!     enddo
     print*,'en,pmom,etc =',eni,pmomi(1),rhoi,gamma,h
     print*,'error: root not bracketed for bisection: no solution'
     ierr = 1
     return
  endif
  
  iterate: do while (.not.converged .and. its < itsmax)
     its = its + 1

     call eval_funch(h,fh,dfh)
     
     !--Newton-Raphson iterations
     if (its.lt.itsmaxnr) then
        hnew = h - fh/dfh
        if (its.eq.itsmaxnr-1) print*,its,' hnew = ',h,' fh = ',fh,' dfh = ',dfh,' dh = ',fh/dfh

        !--watch out for large gradients...
        if (hnew.gt.1.2*h) then
           hnew = 1.2*h
        elseif (hnew.lt.0.8*h) then
          hnew = 0.8*h
        endif
     else
        if (its.eq.itsmaxnr) then
           h = 0.5*(hmax + hmin)
           call eval_funch(h,fh,dfh)
        endif
        if (fhmin*fh.lt.0.) then
           hmax = h
        else
           hmin = h
           call eval_funch(hmin,fhmin,dfhmin)
        endif
        hnew = 0.5*(hmax + hmin)
        print*,its,': bisection',h,hnew,fh
     endif
     
     converged = abs((hnew-h)/hin) < tol
     
  !   if (hnew.lt.1.) then
  !      print*,' hnew = ',hnew,' hold = ',h
  !      grad = fh/dfh
  !      do while (hnew < 1.)
  !         grad = grad*0.5
  !         hnew = h - grad
  !      enddo
     !   print*,' ERROR! enthalpy < 1 in conservative2primitive (sr)'
     !   print*,' fh = ',fh,' dfh = ',dfh,' dh = ',fh/dfh
     !   print*,' lorentzfactor = ',lorentzfactor
     !   print*,' enthalpy = ',h,hnew
     !   print*,' pmom = ',pmomi
     !   print*,' rho = ',rhoi
     !   print*,' en = ',eni
     !   !stop
    ! endif
     
     h = hnew
     
  enddo iterate
  
  !--if converged, set remaining primitive variables
  lorentzfactor = sqrt(pmom2/h**2 + 1.)
  densi = rhoi/lorentzfactor
  vmag = pmommag/(lorentzfactor*h)
  if (vmag.gt.1.) then
     print*,'EEK! v > 1 in conservative2primitive (sr), v = ',vmag
     vmag = 1.-tiny(vmag)
     ierr = 4
  endif
  if (pmommag.gt.tiny(pmommag)) then
     veli(:) = vmag*pmomi(:)/pmommag
     !print*,' vel = ',veli(:),vmag,vmag**2,pmommag,lorentzfactor*h
  else
     veli(:) = 0.
  endif

  uui = (h-1.)/gamma
  pri = (gamma-1.)*uui*densi
  
  if (its.ge.itsmax) then
     print*,' ERROR! Iterations not converged in conservative2primitive'
     print*,' lorentzfactor = ',lorentzfactor
     print*,' enthalpy = ',h
     print*,' pmom = ',pmomi, ' v = ',veli(:)
     print*,' rho = ',rhoi,' dens = ',densi
     print*,' en = ',eni,' uui = ',uui
     stop
     ierr = 1
  endif
  ierr = -its

  return 
  
contains

subroutine eval_funch(hi,funch,dfunch)
 implicit none
 real, intent(in) :: hi
 real, intent(out) :: funch,dfunch
 real :: lorz2,lorentzfactor,dlorentzfactor

 lorz2 = pmom2/hi**2 + 1.
 lorentzfactor = sqrt(lorz2)
 !print*,' lorentzfactor = ',lorentzfactor
 if (iener.eq.3) then ! total energy
    dlorentzfactor = -(pmom2)/(hi**3*lorentzfactor)
    funch = 1.- hi + (gamma*lorentzfactor*(lorentzfactor*hi - eni))/gam1
    dfunch = -1.+ (gamma*lorz2)/gam1 - (gamma*dlorentzfactor*(eni-2.*hi*lorentzfactor))/gam1
 else                 ! entropy
    funch = 1.- hi + (gamma*eni)/gam1*(rhoi**gam1)*(lorz2**(-0.5*gam1))
    dfunch = -1.+ gamma*eni*(rhoi**gam1)*(pmom2/hi**3)*lorz2**(-0.5*(gamma+1))
 endif
 
 return
end subroutine eval_funch

end subroutine solve_conservative2primitivei

subroutine test_conservative2primitive
 use dimen_mhd, only:ndimV
 implicit none
 real, dimension(ndimV) :: pmomi,veli,velout
 real :: gamma,eni,rhoi,uui,densi,pri,spsoundi
 real :: uuout,prout,densout,lorentzfac
 integer :: iener,ierr
 real, parameter :: tol = 1.e-5
 
 gamma = 1.66666666666667
 densi = 1.0
 uui = 1000.
 pri = (gamma-1.)*uui*densi
 veli(:) = 0.
 veli(1) = 0.9999999
 lorentzfac = 1./sqrt(1.-veli(1)*veli(1))
 print*,' lorentzfactor = ',lorentzfac

 uuout = uui
 prout = pri
 densout = densi
 pmomi(:) = 0.
 iener = 1
 
 print*,' testing conservative2primitive routines ... iener = ',iener
 
 call solve_primitive2conservativei(pmomi,eni,rhoi,gamma,veli,uui,densi,pri,spsoundi,iener,ierr)
 print*,' conservative variables: '
 print*,' en = ',eni,pri/densi**gamma
 print*,' pmomi = ',pmomi(1),lorentzfac*veli(1)*(1. + uui + pri/densi)
 print*,' spsoundi = ',spsoundi,sqrt(gamma*pri/densi)
 print*,' enthalpy = ',1. + uui + pri/densi
 print*,' -------------- '
 
!--throw it off the mark
 prout = 100000
 uuout = 0.001
 densout = 0.2
 call solve_conservative2primitivei(pmomi,eni,rhoi,gamma,velout,uuout,densout,prout,lorentzfac,iener,ierr)
 print*,'used ',abs(ierr),' iterations '
 print*,'lorentz factor = ',lorentzfac
 print*,'pr = ',pri,prout
 print*,'u = ',uui,uuout
 print*,'dens = ',densi,densout
 print*,'vx = ',veli(1),velout(1)
 if (abs(pri-prout).gt.tol .or. abs(uui-uuout).gt.tol .or. abs(densi-densout).gt.tol &
    .or.abs(veli(1)-velout(1)).gt.tol .or. ierr.gt.0) then
    stop 'error in conservative2primitive routines'
 endif

 iener = 3
 print*,' testing conservative2primitive routines ... iener = ',iener
 
 call solve_primitive2conservativei(pmomi,eni,rhoi,gamma,veli,uui,densi,pri,spsoundi,iener,ierr)
 print*,' conservative variables: '
 print*,' en = ',eni,pri/densi**gamma
 print*,' pmomi = ',pmomi(1),lorentzfac*veli(1)*(1. + uui + pri/densi)
 print*,' spsoundi = ',spsoundi,sqrt(gamma*pri/densi)
 print*,' enthalpy = ',1. + uui + pri/densi
 print*,' -------------- '

!--throw it off the mark
 prout = 0.
 uuout = 1.
 densout = 1.e-6

 call solve_conservative2primitivei(pmomi,eni,rhoi,gamma,velout,uuout,densout,prout,lorentzfac,iener,ierr)

 print*,'used ',abs(ierr),' iterations '
 print*,'lorentz factor = ',lorentzfac,sqrt(1./(1.-veli(1)*veli(1)))
 print*,'pr = ',pri,prout,abs(pri-prout)
 print*,'u = ',uui,uuout,abs(uui-uuout)
 print*,'dens = ',densi,densout,abs(densi-densout)
 print*,'vx = ',veli(1),velout(1),abs(veli(1)-velout(1))
 if (abs(pri-prout).gt.tol .or. abs(uui-uuout).gt.tol .or. abs(densi-densout).gt.tol &
     .or.abs(veli(1)-velout(1)).gt.tol .or. ierr.gt.0) then
    stop 'error in conservative2primitive routines (2)'
 endif

end subroutine test_conservative2primitive

end module cons2prim
