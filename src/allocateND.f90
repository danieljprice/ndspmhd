!!----------------------------------------------------------------------
!! This subroutine allocates/reallocates memory for all arrays
!! given the new number of particles (should be initially called from setup).
!!----------------------------------------------------------------------
SUBROUTINE alloc(newsizein)
!
!--include relevant global variables
!
 USE dimen_mhd
 USE debug
 USE loguns
 
 USE bound
 USE derivB
 USE eos
 USE fmagarray
 USE gravity
 USE hterms
 USE linklist
 USE options
 USE part
 USE part_in
 USE rates
!! USE unityfunc
 USE xsph
!
!--define local variables
!
 IMPLICIT NONE
 INTEGER, INTENT(IN) :: newsizein
 INTEGER :: newsize,ioldsize,idumsize
 REAL, DIMENSION(newsizein) :: dumpmass,dumrhoin,dumenin,dumpsiin
 REAL, DIMENSION(newsizein) :: dumrho,dumdrhodt,dumuu,dumdudt,dumen,dumdendt
 REAL, DIMENSION(3,newsizein) :: dumalpha,dumalphain,dumdaldt
 REAL, DIMENSION(newsizein) :: dumpsi
 REAL, DIMENSION(newsizein) :: dumhh,dumgradh,dumpr,dumspsound
 REAL, DIMENSION(newsizein) :: dumdivB,dumhhin,dumdhdt,dumdpsidt
 REAL, DIMENSION(ndim,newsizein) :: dumxin,dumx,dumfgrav
 REAL, DIMENSION(ndimV,newsizein) :: dumvelin,dumvel,dumBconsin
 REAL, DIMENSION(ndimV,newsizein) :: dumforce,dumBcons,dumdBconsdt,dumfmag
 REAL, DIMENSION(ndimV,newsizein) :: dumBfield,dumcurlB,dumxsphterm
!--gr terms
 REAL, DIMENSION(newsizein) :: dumsqrtg,dumdens
 INTEGER, DIMENSION(newsizein) :: idumireal,idumitype
!--unity functions
!! REAL, DIMENSION(newsizein) :: dumunity
!! REAL, DIMENSION(ndim,newsizein) :: dumgradunity

 LOGICAL :: reallocate
!
!--set size of neighbour list
!
 nlistdim = newsizein/2	! ie. up to half of particles can be neighbours
!
!--work out whether reallocating memory and what the current array size is
! 
 reallocate = ALLOCATED(x)
 IF (allocated(rho)) THEN
    ioldsize = SIZE(rho)  ! current array size
 ELSE
    ioldsize = 0
 ENDIF
!
!--check for errors
! 
 IF (newsizein.LE.0) THEN
    WRITE(iprint,*) 'Error: size<=0 in call to allocate : ',newsizein
    CALL quit
 ENDIF
!
!--set new size (add 10% so don't have to do this too often)
!
 IF (ANY(ibound.GE.1) .AND. newsizein.GE.ioldsize) THEN
    newsize = INT(1.1*newsizein)
 ELSE
    newsize = newsizein   
 ENDIF

 IF (reallocate) THEN
    WRITE(iprint,*) 'Reallocating memory for all arrays, old, new =',SIZE(rho),newsize
 ELSEIF (trace) THEN
    WRITE(iprint,*) 'Allocating memory for all arrays',newsize
 ENDIF
 
 IF (reallocate) THEN
    idumsize = MIN(newsize,ioldsize)  ! dummy array size is minimum of old/new
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
    dumBconsin(:,1:idumsize) = Bconsin(:,1:idumsize)
    
    dumx(:,1:idumsize) = x(:,1:idumsize)
    dumvel(:,1:idumsize) = vel(:,1:idumsize)
    dumforce(:,1:idumsize) = force(:,1:idumsize)
    dumBcons(:,1:idumsize) = Bcons(:,1:idumsize)
    dumBfield(:,1:idumsize) = Bfield(:,1:idumsize)
    dumdBconsdt(:,1:idumsize) = dBconsdt(:,1:idumsize)
    dumfmag(:,1:idumsize) = fmag(:,1:idumsize)
    dumcurlB(:,1:idumsize) = curlB(:,1:idumsize)
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
    dumpr(1:idumsize) = pr(1:idumsize)
    dumspsound(1:idumsize) = spsound(1:idumsize)
!    dumll(1:idumsize) = ll(1:idumsize)
    dumdivB(1:idumsize) = divB(1:idumsize)
    idumireal(1:idumsize) = ireal(1:idumsize)
    idumitype(1:idumsize) = itype(1:idumsize)
    IF (ALLOCATED(fgrav)) dumfgrav(:,1:idumsize) = fgrav(:,1:idumsize)
    dumsqrtg(1:idumsize) = sqrtg(1:idumsize)
    dumdens(1:idumsize) = dens(1:idumsize)
    
!!    dumunity(1:idumsize) = unity(1:idumsize)
!!    dumgradunity(:,1:idumsize) = gradunity(:,1:idumsize)

!-----------------------------------------------------------------------------
!  deallocate the arrays
!-----------------------------------------------------------------------------
!
!--initial particle properties
!
    DEALLOCATE (pmass,xin,rhoin,hhin,enin,alphain,psiin)
    DEALLOCATE (velin)
    IF (ALLOCATED(Bconsin)) DEALLOCATE (Bconsin)
!
!--particle properties and derivatives
!
    DEALLOCATE(x,vel,force,rho,drhodt,uu,dudt,en,dendt)
    DEALLOCATE(alpha,daldt,psi,dpsidt,hh,dhdt,gradh,pr)
    IF (ALLOCATED(Bfield)) DEALLOCATE(Bfield)
    IF (ALLOCATED(Bcons)) DEALLOCATE(Bcons)
    IF (ALLOCATED(dBconsdt)) DEALLOCATE(dBconsdt)
!
!--equation of state
!
    DEALLOCATE(spsound)
!
!--linklist and boundaries - note ifirstincell is reallocated in link
!  (depends on max # of link list cells)
!
    DEALLOCATE(ll,ireal,ifirstincell,iamincell)
!
!--itype is particle type (normal, fixed particle etc)
!
    IF (ALLOCATED(itype)) DEALLOCATE(itype)
!
!--MHD quantities and derivatives
!
    IF (ALLOCATED(fmag)) DEALLOCATE(fmag)
    IF (ALLOCATED(divB)) DEALLOCATE(divB)
    IF (ALLOCATED(curlB)) DEALLOCATE(curlB)
!
!--XSPH
!
    IF (ALLOCATED(xsphterm)) DEALLOCATE(xsphterm)
!
!--gravity
!
    IF (ALLOCATED(fgrav)) DEALLOCATE(fgrav)
!
!--GR
!
    IF (ALLOCATED(sqrtg)) DEALLOCATE(sqrtg)
    IF (ALLOCATED(sourceterms)) DEALLOCATE(sourceterms)
    IF (ALLOCATED(pmom)) DEALLOCATE(pmom)
    IF (ALLOCATED(dens)) DEALLOCATE(dens)
!
!--kernel correction
!
!!    IF (ALLOCATED(unity)) DEALLOCATE(unity)
!!    IF (ALLOCATED(gradunity)) DEALLOCATE(gradunity)
 ENDIF

!-----------------------------------------------------------------------------
!  allocate all arrays (for both first time and reallocation)
!-----------------------------------------------------------------------------

!
!--initial particle properties
!
    ALLOCATE (pmass(newsize))
    ALLOCATE (xin(ndim,newsize))
    ALLOCATE (rhoin(newsize),hhin(newsize))
    ALLOCATE (enin(newsize),psiin(newsize))   ! alpha done below
    ALLOCATE (velin(ndimV,newsize))
    ALLOCATE (Bconsin(ndimB,newsize))
!
!--particle properties and derivatives
!
    ALLOCATE(x(ndim,newsize))
    ALLOCATE(vel(ndimV,newsize),force(ndimV,newsize))
    ALLOCATE(rho(newsize),drhodt(newsize))
    ALLOCATE(uu(newsize),dudt(newsize),en(newsize),dendt(newsize))
    ALLOCATE(alpha(3,newsize),alphain(3,newsize),daldt(3,newsize))
    ALLOCATE(psi(newsize),dpsidt(newsize))
    ALLOCATE(hh(newsize),dhdt(newsize),gradh(newsize),pr(newsize))
    ALLOCATE(Bcons(ndimB,newsize))
    ALLOCATE(Bfield(ndimB,newsize),dBconsdt(ndimB,newsize))	! mag field
!
!--equation of state
!
    ALLOCATE(spsound(newsize))
!
!--linklist and boundaries - note ifirstincell is reallocated in link
!  (depends on max # of link list cells)
!
    ALLOCATE(ll(newsize),iamincell(newsize),ifirstincell(0),ireal(newsize))
!
!--particle type
!
    ALLOCATE(itype(newsize))
!
!--MHD quantities and derivatives
!
    ALLOCATE(fmag(ndimB,newsize),divB(newsize),curlB(ndimB,newsize))
!
!--XSPH
!
    ALLOCATE(xsphterm(ndimV,newsize))
!
!--gravity
!   
    IF (igravity.NE.0) THEN
       ALLOCATE(fgrav(ndim,newsize))
    ENDIF
!
!--GR
!
   ALLOCATE(sqrtg(newsize))
   ALLOCATE(dens(newsize))
!
!--kernel correction
!   
!!   ALLOCATE(unity(newsize))
!!   ALLOCATE(gradunity(ndim,newsize))
   
 IF (reallocate) THEN
!-----------------------------------------------------------------------------
!  copy properties back from old arrays to new arrays
!-----------------------------------------------------------------------------
    pmass(1:idumsize) = dumpmass(1:idumsize)
    rhoin(1:idumsize) = dumrhoin(1:idumsize)
    hhin(1:idumsize) = dumhhin(1:idumsize)
    enin(1:idumsize) = dumenin(1:idumsize)
    alphain(:,1:idumsize) = dumalphain(:,1:idumsize)
    psiin(1:idumsize) = dumpsiin(1:idumsize)
    
    xin(:,1:idumsize) = dumxin(:,1:idumsize)
    velin(:,1:idumsize) = dumvelin(:,1:idumsize)
    Bconsin(:,1:idumsize) = dumBconsin(:,1:idumsize)
    
    x(:,1:idumsize) = dumx(:,1:idumsize)
    vel(:,1:idumsize) = dumvel(:,1:idumsize)
    force(:,1:idumsize) = dumforce(:,1:idumsize)
    Bfield(:,1:idumsize) = dumBfield(:,1:idumsize)
    Bcons(:,1:idumsize) = dumBcons(:,1:idumsize)
    dBconsdt(:,1:idumsize) = dumdBconsdt(:,1:idumsize)
    fmag(:,1:idumsize) = dumfmag(:,1:idumsize)
    curlB(:,1:idumsize) = dumcurlB(:,1:idumsize)
    xsphterm(:,1:idumsize) = dumxsphterm(:,1:idumsize)
        
    rho(1:idumsize) = dumrho(1:idumsize)
    drhodt(1:idumsize) = dumdrhodt(1:idumsize)
    uu(1:idumsize) = dumuu(1:idumsize)
    dudt(1:idumsize) = dumdudt(1:idumsize)
    en(1:idumsize) = dumen(1:idumsize)
    dendt(1:idumsize) = dumdendt(1:idumsize)
    alpha(:,1:idumsize) = dumalpha(:,1:idumsize)
    daldt(:,1:idumsize) = dumdaldt(:,1:idumsize)
    psi(1:idumsize) = dumpsi(1:idumsize)
    dpsidt(1:idumsize) = dumdpsidt(1:idumsize)
    hh(1:idumsize) = dumhh(1:idumsize)
    dhdt(1:idumsize) = dumdhdt(1:idumsize)
    gradh(1:idumsize) = dumgradh(1:idumsize)
    pr(1:idumsize) = dumpr(1:idumsize)
    spsound(1:idumsize) = dumspsound(1:idumsize)
!    ll(1:idumsize) = dumll(1:idumsize)
    divB(1:idumsize) = dumdivB(1:idumsize) 
    ireal(1:idumsize) = idumireal(1:idumsize)
    itype(1:idumsize) = idumitype(1:idumsize)
    IF (ALLOCATED(fgrav)) fgrav(:,1:idumsize) = dumfgrav(:,1:idumsize)
    sqrtg(1:idumsize) = dumsqrtg(1:idumsize)
    dens(1:idumsize) = dumdens(1:idumsize)    
    
!!    unity(1:idumsize) = dumunity(1:idumsize)
!!    gradunity(:,1:idumsize) = dumgradunity(:,1:idumsize)
    
 ELSE
    itype(:) = 0	! on first memory allocation, set all parts = normal
 ENDIF   
       
!
!--debugging information
!
 IF (trace) WRITE(iprint,*) '  Memory allocated, size(x,v)=',SIZE(x),SIZE(vel) 
 
 RETURN

END SUBROUTINE alloc
