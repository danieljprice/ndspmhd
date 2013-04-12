module mem_allocation
 implicit none
 public :: alloc
 
 private

contains
!!----------------------------------------------------------------------
!! this subroutine allocates/reallocates memory for all arrays
!! given the new number of particles (should be initially called from setup).
!!----------------------------------------------------------------------
subroutine alloc(newsizein,sortlist)
!
!--include relevant global variables
!
 use dimen_mhd
 use debug
 use loguns
 
 use bound
 use derivb
 use eos
 use fmagarray
 use hterms
 use linklist
 use options
 use part
 use part_in
 use rates
 use xsph
 use matrixcorr
 use rkf
!
!--define local variables
!
 implicit none
 integer, intent(in) :: newsizein
 integer, dimension(newsizein), intent(in), optional :: sortlist
 integer, dimension(newsizein) :: iorder
 integer :: i,newsize,ioldsize,idumsize
 real, dimension(newsizein) :: dumpmass,dumrhoin,dumenin,dumpsiin
 real, dimension(newsizein) :: dumrho,dumdrhodt,dumuu,dumdudt,dumen,dumdendt
 real, dimension(3,newsizein) :: dumalpha,dumalphain,dumdaldt
 real, dimension(newsizein) :: dumpsi,dumrho0
 real, dimension(newsizein) :: dumhh,dumgradh,dumgradhn,dumgradsoft,dumgradgradh,dumpr,dumspsound
 real, dimension(newsizein) :: dumdivB,dumhhin,dumdhdt,dumdpsidt,dumpoten
 real, dimension(ndim,newsizein) :: dumxin,dumx
 real, dimension(ndimv,newsizein) :: dumvelin,dumvel,dumBevolin
 real, dimension(ndimv,newsizein) :: dumforce,dumBevol,dumdBevoldt,dumfmag
 real, dimension(ndimv,newsizein) :: dumBfield,dumcurlB,dumxsphterm,dumgradpsi
!--gr terms
 real, dimension(newsizein) :: dumsqrtg,dumdens
 real, dimension(ndimv,newsizein) :: dumsourceterms,dumpmom,dumpmomin
 integer, dimension(newsizein) :: idumireal,idumitype,idumnumneigh
 real, dimension(ndxdx,newsizein) :: dumdxdx
 real, dimension(ndim,newsizein) :: dumx0
!--dust
 real, dimension(newsizein)       :: dumdusttogas,dumdusttogasin,dumddusttogasdt
 real, dimension(ndimV,newsizein) :: dumdeltav,dumdeltavin,dumddeltavdt
!--integrator type check
 character (len=20) :: integratorcheck

 logical :: reallocate, isortparts
!
!--work out whether reallocating memory and what the current array size is
! 
 reallocate = allocated(x)
 if (allocated(rho)) then
    ioldsize = size(rho)  ! current array size
 else
    ioldsize = 0
 endif
 idumsize = 0
!
!--set default list order
!
 do i=1,size(iorder)
    iorder(i) = i
 enddo
 isortparts = .false.
!
!--check for errors
! 
 if (newsizein.le.0) then
    write(iprint,*) 'error: size<=0 in call to allocate : ',newsizein
    stop
 endif
!
!--set new size (add 5% if using ghosts so don't have to do this too often)
!
 if (any(ibound.ge.1) .and. newsizein.gt.ioldsize) then
    newsize = int(1.1*newsizein)
 elseif (newsizein.eq.ioldsize) then
    if (present(sortlist)) then
       iorder(:) = sortlist(:)
       isortparts = .true.
       reallocate = .false.
    else
       write(iprint,*) 'allocate: old size = new size: doing nothing'
       return
    endif
 else
    newsize = newsizein   
 endif

 if (reallocate) then
    write(iprint,*) 'reallocating memory for all arrays, old, new =',size(rho),newsize
 elseif (isortparts) then
    write(iprint,*) 'reshuffling all particles...'
 else
    write(iprint,*) 'allocating memory for all arrays, maxpart = ',newsize
 endif
 
 if (reallocate .or. isortparts) then
    idumsize = min(newsize,ioldsize)  ! dummy array size is minimum of old/new
!-----------------------------------------------------------------------------
!  if reallocating, copy everything into dummy arrays using the smallest size
!-----------------------------------------------------------------------------
    dumpmass(1:idumsize) = pmass(1:idumsize)
    dumrhoin(1:idumsize) = rhoin(1:idumsize)
    dumhhin(1:idumsize) = hhin(1:idumsize)
    dumenin(1:idumsize) = enin(1:idumsize)
    dumalphain(:,1:idumsize) = alphain(:,1:idumsize)
    dumpsiin(1:idumsize) = psiin(1:idumsize)
    
    dumxin(:,1:idumsize) = xin(:,1:idumsize)
    dumvelin(:,1:idumsize) = velin(:,1:idumsize)
    dumBevolin(:,1:idumsize) = Bevolin(:,1:idumsize)
    
    dumx(:,1:idumsize) = x(:,1:idumsize)
    dumvel(:,1:idumsize) = vel(:,1:idumsize)
    dumforce(:,1:idumsize) = force(:,1:idumsize)
    dumBevol(:,1:idumsize) = Bevol(:,1:idumsize)
    dumBfield(:,1:idumsize) = Bfield(:,1:idumsize)
    dumdBevoldt(:,1:idumsize) = dBevoldt(:,1:idumsize)
    dumfmag(:,1:idumsize) = fmag(:,1:idumsize)
    dumcurlB(:,1:idumsize) = curlB(:,1:idumsize)
    dumgradpsi(:,1:idumsize) = gradpsi(:,1:idumsize)
    dumxsphterm(:,1:idumsize) = xsphterm(:,1:idumsize)
        
    dumrho(1:idumsize) = rho(1:idumsize)
    dumdrhodt(1:idumsize) = drhodt(1:idumsize)
    dumuu(1:idumsize) = uu(1:idumsize)
    dumdudt(1:idumsize) = dudt(1:idumsize)
    dumen(1:idumsize) = en(1:idumsize)
    dumdendt(1:idumsize) = dendt(1:idumsize)
    dumalpha(:,1:idumsize) = alpha(:,1:idumsize)
    dumdaldt(:,1:idumsize) = daldt(:,1:idumsize)
    dumpsi(1:idumsize) = psi(1:idumsize)
    dumdpsidt(1:idumsize) = dpsidt(1:idumsize)    
    dumhh(1:idumsize) = hh(1:idumsize)
    dumdhdt(1:idumsize) = dhdt(1:idumsize)
    dumgradh(1:idumsize) = gradh(1:idumsize)
    dumgradhn(1:idumsize) = gradhn(1:idumsize)
    dumgradsoft(1:idumsize) = gradsoft(1:idumsize)
    dumgradgradh(1:idumsize) = gradgradh(1:idumsize)

    dumpr(1:idumsize) = pr(1:idumsize)
    dumspsound(1:idumsize) = spsound(1:idumsize)
!    dumll(1:idumsize) = ll(1:idumsize)
    dumdivB(1:idumsize) = divB(1:idumsize)
    idumireal(1:idumsize) = ireal(1:idumsize)
    idumitype(1:idumsize) = itype(1:idumsize)
    idumnumneigh(1:idumsize) = numneigh(1:idumsize)
    dumsqrtg(1:idumsize) = sqrtg(1:idumsize)
    dumdens(1:idumsize) = dens(1:idumsize)
    if (allocated(poten)) dumpoten(1:idumsize) = poten(1:idumsize)
    if (allocated(sourceterms)) dumsourceterms(:,1:idumsize) = sourceterms(:,1:idumsize)
    if (allocated(pmom)) dumpmom(:,1:idumsize) = pmom(:,1:idumsize)
    if (allocated(pmomin)) dumpmomin(:,1:idumsize) = pmomin(:,1:idumsize)
    
    if (allocated(dusttogas))    dumdusttogas(1:idumsize)    = dusttogas(1:idumsize)
    if (allocated(ddusttogasdt)) dumddusttogasdt(1:idumsize) = ddusttogasdt(1:idumsize)
    if (allocated(deltav))       dumdeltav(:,1:idumsize)     = deltav(:,1:idumsize)
    if (allocated(ddeltavdt))    dumddeltavdt(:,1:idumsize)  = ddeltavdt(:,1:idumsize)
    if (allocated(dusttogasin))    dumdusttogasin(1:idumsize) = dusttogasin(1:idumsize)
    if (allocated(deltavin))       dumdeltavin(:,1:idumsize)  = deltavin(:,1:idumsize)
    
    dumdxdx(:,1:idumsize) = dxdx(:,1:idumsize)
    dumx0(:,1:idumsize) = x0(:,1:idumsize)
    if (allocated(rho0)) dumrho0(1:idumsize) = rho0(1:idumsize)

!-----------------------------------------------------------------------------
!  deallocate the arrays
!-----------------------------------------------------------------------------

    if (reallocate) then
!
!--initial particle properties
!
    deallocate (pmass,xin,rhoin,hhin,enin,alphain,psiin)
    deallocate (velin)
    if (allocated(Bevolin)) deallocate (Bevolin)
!
!--Runge-Kutta-Fehlberg arrays
!
    if (allocated(v1)) deallocate(v1)
    if (allocated(v2)) deallocate(v2)
    if (allocated(v3)) deallocate(v3)
    if (allocated(v4)) deallocate(v4)
    if (allocated(v5)) deallocate(v5)
    if (allocated(v6)) deallocate(v6)
    if (allocated(f1)) deallocate(f1)
    if (allocated(f2)) deallocate(f2)
    if (allocated(f3)) deallocate(f3)
    if (allocated(f4)) deallocate(f4)
    if (allocated(f5)) deallocate(f5)
    if (allocated(f6)) deallocate(f6)
    if (allocated(dh1)) deallocate(dh1)
    if (allocated(dh2)) deallocate(dh2)
    if (allocated(dh3)) deallocate(dh3)
    if (allocated(dh4)) deallocate(dh4)
    if (allocated(dh5)) deallocate(dh5)
    if (allocated(dh6)) deallocate(dh6)
!
!--particle properties and derivatives
!
    deallocate(x,vel,force,rho,drhodt,uu,dudt,en,dendt)
    deallocate(alpha,daldt,psi,dpsidt,hh,dhdt,gradh,gradhn,gradsoft,pr)
    deallocate(gradgradh)
    if (allocated(zeta)) deallocate(zeta)
    if (allocated(poten)) deallocate(poten)
    if (allocated(Bfield)) deallocate(Bfield)
    if (allocated(Bevol)) deallocate(Bevol)
    if (allocated(dBevoldt)) deallocate(dBevoldt)
    if (allocated(gradpsi)) deallocate(gradpsi)
!
!--equation of state
!
    deallocate(spsound)
!
!--linklist and boundaries - note ifirstincell is reallocated in link
!  (depends on max # of link list cells)
!
    deallocate(ll,ireal,ifirstincell,iamincell)
!
!--itype is particle type (normal, fixed particle etc)
!
    if (allocated(itype)) deallocate(itype)
    if (allocated(numneigh)) deallocate(numneigh)
!
!--mhd quantities and derivatives
!
    if (allocated(fmag)) deallocate(fmag)
    if (allocated(divB)) deallocate(divB)
    if (allocated(curlB)) deallocate(curlB)
!
!--xsph
!
    if (allocated(xsphterm)) deallocate(xsphterm)
!
!--gr
!
    if (allocated(sqrtg)) deallocate(sqrtg)
    if (allocated(sourceterms)) deallocate(sourceterms)
    if (allocated(pmom)) deallocate(pmom)
    if (allocated(pmomin)) deallocate(pmomin)
    if (allocated(dens)) deallocate(dens)
    
    if (allocated(dxdx)) deallocate(dxdx)
    if (allocated(x0)) deallocate(x0)
    if (allocated(rho0)) deallocate(rho0)
!
!--dust
!
    if (allocated(dusttogas))    deallocate(dusttogas)
    if (allocated(deltav))       deallocate(deltav)
    if (allocated(ddusttogasdt)) deallocate(ddusttogasdt)
    if (allocated(ddeltavdt))    deallocate(ddeltavdt)
    if (allocated(dusttogasin)) deallocate(dusttogasin)
    if (allocated(deltavin))    deallocate(deltavin)
!
!--physical viscosity
!
    if (allocated(del2v)) deallocate(del2v)
    endif

 endif

 if (reallocate .or. .not.allocated(x)) then
!-----------------------------------------------------------------------------
!  allocate all arrays (for both first time and reallocation)
!-----------------------------------------------------------------------------
    idim = newsize
!
!--initial particle properties
!
    allocate (pmass(newsize))
    allocate (xin(ndim,newsize))
    allocate (rhoin(newsize),hhin(newsize))
    allocate (enin(newsize),psiin(newsize))   ! alpha done below
    allocate (velin(ndimv,newsize))
    allocate (Bevolin(ndimb,newsize))
!
!--particle properties and derivatives
!
    allocate(x(ndim,newsize))
    allocate(vel(ndimv,newsize),force(ndimv,newsize))
    allocate(rho(newsize),drhodt(newsize))
    allocate(uu(newsize),dudt(newsize),en(newsize),dendt(newsize))
    allocate(alpha(3,newsize),alphain(3,newsize),daldt(3,newsize))
    allocate(psi(newsize),dpsidt(newsize))
    allocate(hh(newsize))
    allocate(dhdt(newsize),gradh(newsize),gradhn(newsize),gradsoft(newsize))
    allocate(gradgradh(newsize))
    if (imhd.lt.0) allocate(zeta(newsize))
    allocate(pr(newsize))
    allocate(Bevol(ndimb,newsize))
    allocate(Bfield(ndimb,newsize),dBevoldt(ndimb,newsize))  ! mag field
    allocate(gradpsi(ndimb,newsize))
    if (igravity.ne.0) allocate(poten(newsize))
!
!--equation of state
!
    allocate(spsound(newsize))
!
!--linklist and boundaries - note ifirstincell is reallocated in link
!  (depends on max # of link list cells)
!
    allocate(ll(newsize),iamincell(newsize),ifirstincell(0),ireal(newsize))
!
!--particle type
!
    allocate(itype(newsize),numneigh(newsize))
    itype(:) = 0 ! give type 0 by default
!
!--mhd quantities and derivatives
!
    allocate(fmag(ndimb,newsize),divB(newsize),curlB(ndimb,newsize))
!
!--xsph
!
    allocate(xsphterm(ndimv,newsize))
!
!--gr
!
   allocate(sqrtg(newsize))
   allocate(dens(newsize))
   allocate(pmom(ndimv,newsize))
   allocate(pmomin(ndimv,newsize))
   if (geom(1:4).ne.'cart') then
      allocate(sourceterms(ndimv,newsize))
   endif
   allocate(dxdx(ndxdx,newsize))
   allocate(x0(ndim,newsize))
   allocate(rho0(newsize))
!
!--dust
!
   if (idust.eq.1) then
      allocate(dusttogas(newsize),dusttogasin(newsize))
      allocate(ddusttogasdt(newsize))
      allocate(deltav(ndimV,newsize),deltavin(ndimV,newsize))
      allocate(ddeltavdt(ndimV,newsize))
   endif
!
!--physical viscosity
!
   if (ivisc.gt.0) allocate(del2v(newsize))   
!
!--rkf extra storage
!
   integratorcheck = 'query'
   call step (integratorcheck)
   if (trim(integratorcheck).eq.'rkf') then
       allocate (v1(ndimv,newsize))
       allocate (v2(ndimv,newsize))
       allocate (v3(ndimv,newsize))
       allocate (v4(ndimv,newsize))
       allocate (v5(ndimv,newsize))
       allocate (v6(ndimv,newsize))
       allocate (f1(ndimv,newsize))
       allocate (f2(ndimv,newsize))
       allocate (f3(ndimv,newsize))
       allocate (f4(ndimv,newsize))
       allocate (f5(ndimv,newsize))
       allocate (f6(ndimv,newsize))
       allocate (dh1(newsize))
       allocate (dh2(newsize))
       allocate (dh3(newsize))
       allocate (dh4(newsize))
       allocate (dh5(newsize))
       allocate (dh6(newsize))
    endif
 endif
 
 if (reallocate .or. isortparts) then
!-----------------------------------------------------------------------------
!  copy properties back from old arrays to new arrays
!-----------------------------------------------------------------------------
    pmass(1:idumsize) = dumpmass(iorder(1:idumsize))
    rhoin(1:idumsize) = dumrhoin(iorder(1:idumsize))
    hhin(1:idumsize) = dumhhin(iorder(1:idumsize))
    enin(1:idumsize) = dumenin(iorder(1:idumsize))
    alphain(:,1:idumsize) = dumalphain(:,iorder(1:idumsize))
    psiin(1:idumsize) = dumpsiin(iorder(1:idumsize))
    
    xin(:,1:idumsize) = dumxin(:,iorder(1:idumsize))
    velin(:,1:idumsize) = dumvelin(:,iorder(1:idumsize))
    Bevolin(:,1:idumsize) = dumBevolin(:,iorder(1:idumsize))
    
    x(:,1:idumsize) = dumx(:,iorder(1:idumsize))
    vel(:,1:idumsize) = dumvel(:,iorder(1:idumsize))
    force(:,1:idumsize) = dumforce(:,iorder(1:idumsize))
    Bfield(:,1:idumsize) = dumBfield(:,iorder(1:idumsize))
    Bevol(:,1:idumsize) = dumBevol(:,iorder(1:idumsize))
    dBevoldt(:,1:idumsize) = dumdBevoldt(:,iorder(1:idumsize))
    gradpsi(:,1:idumsize) = dumgradpsi(:,iorder(1:idumsize))
    fmag(:,1:idumsize) = dumfmag(:,iorder(1:idumsize))
    curlB(:,1:idumsize) = dumcurlB(:,iorder(1:idumsize))
    xsphterm(:,1:idumsize) = dumxsphterm(:,iorder(1:idumsize))
        
    rho(1:idumsize) = dumrho(iorder(1:idumsize))
    drhodt(1:idumsize) = dumdrhodt(iorder(1:idumsize))
    uu(1:idumsize) = dumuu(iorder(1:idumsize))
    dudt(1:idumsize) = dumdudt(iorder(1:idumsize))
    en(1:idumsize) = dumen(iorder(1:idumsize))
    dendt(1:idumsize) = dumdendt(iorder(1:idumsize))
    alpha(:,1:idumsize) = dumalpha(:,iorder(1:idumsize))
    daldt(:,1:idumsize) = dumdaldt(:,iorder(1:idumsize))
    psi(1:idumsize) = dumpsi(iorder(1:idumsize))
    dpsidt(1:idumsize) = dumdpsidt(iorder(1:idumsize))
    hh(1:idumsize) = dumhh(iorder(1:idumsize))
    dhdt(1:idumsize) = dumdhdt(iorder(1:idumsize))
    gradh(1:idumsize) = dumgradh(iorder(1:idumsize))
    gradhn(1:idumsize) = dumgradhn(iorder(1:idumsize))
    gradsoft(1:idumsize) = dumgradsoft(iorder(1:idumsize))
    gradgradh(1:idumsize) = dumgradgradh(iorder(1:idumsize))
    if (allocated(poten)) poten(1:idumsize) = dumpoten(iorder(1:idumsize))

    pr(1:idumsize) = dumpr(iorder(1:idumsize))
    spsound(1:idumsize) = dumspsound(iorder(1:idumsize))
!    ll(1:idumsize) = dumll(1:idumsize)
    divB(1:idumsize) = dumdivB(iorder(1:idumsize))
    ireal(1:idumsize) = idumireal(iorder(1:idumsize))
    itype(1:idumsize) = idumitype(iorder(1:idumsize))
    numneigh(1:idumsize) = idumnumneigh(iorder(1:idumsize))
    sqrtg(1:idumsize) = dumsqrtg(iorder(1:idumsize))
    dens(1:idumsize) = dumdens(iorder(1:idumsize))
    pmom(:,1:idumsize) = dumpmom(:,iorder(1:idumsize))
    pmomin(:,1:idumsize) = dumpmomin(:,iorder(1:idumsize))
    if (allocated(sourceterms)) then
       sourceterms(:,1:idumsize) = dumsourceterms(:,iorder(1:idumsize))
    endif
    dxdx(:,1:idumsize) = dumdxdx(:,iorder(1:idumsize))
    x0(:,1:idumsize) = dumx0(:,iorder(1:idumsize))
    rho0(1:idumsize) = dumrho0(iorder(1:idumsize))
!
!--dust
!
    if (idust.eq.1) then
       dusttogas(1:idumsize)    = dumdusttogas(iorder(1:idumsize))
       ddusttogasdt(1:idumsize) = dumddusttogasdt(iorder(1:idumsize))
       deltav(:,1:idumsize)     = dumdeltav(:,iorder(1:idumsize))
       ddeltavdt(:,1:idumsize)  = dumddeltavdt(:,iorder(1:idumsize))
       dusttogasin(1:idumsize)  = dumdusttogasin(iorder(1:idumsize))
       deltavin(:,1:idumsize)   = dumdeltavin(:,iorder(1:idumsize))
    endif
    ! no need to copy physical viscosity stuff
 else
    itype(:) = 0 ! on first memory allocation, set all parts = normal
    numneigh(:) = 0
    dxdx(:,:) = 1.
 endif
 
!
!--debugging information
!
 if (trace) write(iprint,*) '  memory allocated, size(x,v)=',size(x),size(vel) 
 
 return

end subroutine alloc

end module mem_allocation
