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
 USE xsph
!
!--define local variables
!
 IMPLICIT NONE
 INTEGER, INTENT(IN) :: newsizein
 INTEGER :: newsize,ioldsize,idumsize
 REAL, DIMENSION(newsizein) :: dumpmass,dumrhoin,dumuuin,dumenin,dumalphain
 REAL, DIMENSION(newsizein) :: dumrho,dumdrhodt,dumuu,dumdudt,dumen,dumdendt
 REAL, DIMENSION(newsizein) :: dumalpha,dumdaldt,dumhh,dumgradh,dumpr,dumspsound
 REAL, DIMENSION(newsizein) :: dumdivB,dumhhin,dumdhdt
 REAL, DIMENSION(ndim,newsizein) :: dumxin,dumx,dumfgrav
 REAL, DIMENSION(ndimV,newsizein) :: dumvelin,dumvel,dumBfieldin,dumBin
 REAL, DIMENSION(ndimV,newsizein) :: dumforce,dumBfield,dumdBfielddt,dumfmag
 REAL, DIMENSION(ndimV,newsizein) :: dumcurlB,dumxsphterm
 INTEGER, DIMENSION(newsizein) :: idumireal,idumitype
 LOGICAL :: reallocate
!
!--set size of neighbour list
!
 nlistdim = newsize/2	! ie. up to half of particles can be neighbours
!
!--work out whether reallocating memory and what the current array size is
! 
 reallocate = ALLOCATED(x)
 ioldsize = SIZE(rho) 		! current array size
!
!--check for errors
! 
 IF (newsizein.EQ.0) THEN
    WRITE(iprint,*) 'Error: size=0 in call to allocate'
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
    dumuuin(1:idumsize) = uuin(1:idumsize)
    dumhhin(1:idumsize) = hhin(1:idumsize)
    dumenin(1:idumsize) = enin(1:idumsize)
    dumalphain(1:idumsize) = alphain(1:idumsize)
    
    dumxin(:,1:idumsize) = xin(:,1:idumsize)
    dumvelin(:,1:idumsize) = velin(:,1:idumsize)
    dumBfieldin(:,1:idumsize) = Bfieldin(:,1:idumsize)
    dumBin(:,1:idumsize) = Bin(:,1:idumsize)
    
    dumx(:,1:idumsize) = x(:,1:idumsize)
    dumvel(:,1:idumsize) = vel(:,1:idumsize)
    dumforce(:,1:idumsize) = force(:,1:idumsize)
    dumBfield(:,1:idumsize) = Bfield(:,1:idumsize)
    dumdBfielddt(:,1:idumsize) = dBfielddt(:,1:idumsize)
    dumfmag(:,1:idumsize) = fmag(:,1:idumsize)
    dumcurlB(:,1:idumsize) = curlB(:,1:idumsize)
    dumxsphterm(:,1:idumsize) = xsphterm(:,1:idumsize)
        
    dumrho(1:idumsize) = rho(1:idumsize)
    dumdrhodt(1:idumsize) = drhodt(1:idumsize)
    dumuu(1:idumsize) = uu(1:idumsize)
    dumdudt(1:idumsize) = dudt(1:idumsize)
    dumen(1:idumsize) = en(1:idumsize)
    dumdendt(1:idumsize) = dendt(1:idumsize)
    dumalpha(1:idumsize) = alpha(1:idumsize)
    dumdaldt(1:idumsize) = daldt(1:idumsize)
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
    
 
!-----------------------------------------------------------------------------
!  deallocate the arrays
!-----------------------------------------------------------------------------
!
!--initial particle properties
!
    DEALLOCATE (pmass,xin,rhoin,uuin,hhin,enin,alphain)
    DEALLOCATE (velin)
    IF (ALLOCATED(Bfieldin)) DEALLOCATE (Bfieldin)
    IF (ALLOCATED(Bin)) DEALLOCATE(Bin)
!
!--particle properties and derivatives
!
    DEALLOCATE(x,vel,force,rho,drhodt,uu,dudt,en,dendt)
    DEALLOCATE(alpha,daldt,hh,dhdt,gradh,pr)
    IF (ALLOCATED(Bfield)) DEALLOCATE(Bfield)
    IF (ALLOCATED(dBfielddt)) DEALLOCATE(dBfielddt)
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

 ENDIF

!-----------------------------------------------------------------------------
!  allocate all arrays (for both first time and reallocation)
!-----------------------------------------------------------------------------

!
!--initial particle properties
!
    ALLOCATE (pmass(newsize))
    ALLOCATE (xin(ndim,newsize))
    ALLOCATE (rhoin(newsize),uuin(newsize),hhin(newsize))
    ALLOCATE (enin(newsize),alphain(newsize))
    ALLOCATE (velin(ndimV,newsize))
    ALLOCATE (Bfieldin(ndimB,newsize),Bin(ndimB,newsize))
!
!--particle properties and derivatives
!
    ALLOCATE(x(ndim,newsize))
    ALLOCATE(vel(ndimV,newsize),force(ndimV,newsize))
    ALLOCATE(rho(newsize),drhodt(newsize))
    ALLOCATE(uu(newsize),dudt(newsize),en(newsize),dendt(newsize))
    ALLOCATE(alpha(newsize),daldt(newsize))
    ALLOCATE(hh(newsize),dhdt(newsize),gradh(newsize),pr(newsize))
    ALLOCATE(Bfield(ndimB,newsize),dBfielddt(ndimB,newsize))	! mag field
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
    
 IF (reallocate) THEN
!-----------------------------------------------------------------------------
!  copy properties back from old arrays to new arrays
!-----------------------------------------------------------------------------
    pmass(1:idumsize) = dumpmass(1:idumsize)
    rhoin(1:idumsize) = dumrhoin(1:idumsize)
    uuin(1:idumsize) = dumuuin(1:idumsize)
    hhin(1:idumsize) = dumhhin(1:idumsize)
    enin(1:idumsize) = dumenin(1:idumsize)
    alphain(1:idumsize) = dumalphain(1:idumsize)
    
    xin(:,1:idumsize) = dumxin(:,1:idumsize)
    velin(:,1:idumsize) = dumvelin(:,1:idumsize)
    Bfieldin(:,1:idumsize) = dumBfieldin(:,1:idumsize)
    Bin(:,1:idumsize) = dumBin(:,1:idumsize)
    
    x(:,1:idumsize) = dumx(:,1:idumsize)
    vel(:,1:idumsize) = dumvel(:,1:idumsize)
    force(:,1:idumsize) = dumforce(:,1:idumsize)
    Bfield(:,1:idumsize) = dumBfield(:,1:idumsize)
    dBfielddt(:,1:idumsize) = dumdBfielddt(:,1:idumsize)
    fmag(:,1:idumsize) = dumfmag(:,1:idumsize)
    curlB(:,1:idumsize) = dumcurlB(:,1:idumsize)
    xsphterm(:,1:idumsize) = dumxsphterm(:,1:idumsize)
        
    rho(1:idumsize) = dumrho(1:idumsize)
    drhodt(1:idumsize) = dumdrhodt(1:idumsize)
    uu(1:idumsize) = dumuu(1:idumsize)
    dudt(1:idumsize) = dumdudt(1:idumsize)
    en(1:idumsize) = dumen(1:idumsize)
    dendt(1:idumsize) = dumdendt(1:idumsize)
    alpha(1:idumsize) = dumalpha(1:idumsize)
    daldt(1:idumsize) = dumdaldt(1:idumsize)
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
 ELSE
    itype(:) = 0	! on first memory allocation, set all parts = normal
 ENDIF   
       
!
!--debugging information
!
 IF (trace) WRITE(iprint,*) '  Memory allocated, size(x,v)=',SIZE(x),SIZE(vel) 
 
 RETURN

END SUBROUTINE alloc
