!!----------------------------------------------------------------------
!! This subroutine allocates memory for all arrays once the number of
!! particles is known (should be called from setup).
!!
!! icall = 1 : allocate arrays for the first time
!! icall = 2 : reallocate arrays because particle # > array size
!!
!!----------------------------------------------------------------------
SUBROUTINE alloc(nsize,icall)
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
 INTEGER, INTENT(IN) :: nsize,icall
 INTEGER :: ntot
 REAL, DIMENSION(nsize) :: dumpmass,dumrhoin,dumuuin,dumenin,dumalphain
 REAL, DIMENSION(nsize) :: dumrho,dumdrhodt,dumuu,dumdudt,dumen,dumdendt
 REAL, DIMENSION(nsize) :: dumalpha,dumdaldt,dumhh,dumgradh,dumpr,dumspsound
 REAL, DIMENSION(nsize) :: dumdivB,dumhhin,dumdhdt
 REAL, DIMENSION(ndim,nsize) :: dumxin,dumx,dumfgrav
 REAL, DIMENSION(ndimV,nsize) :: dumvelin,dumvel,dumBfieldin,dumBin
 REAL, DIMENSION(ndimV,nsize) :: dumforce,dumBfield,dumdBfielddt,dumfmag
 REAL, DIMENSION(ndimV,nsize) :: dumcurlB,dumxsphterm
 INTEGER, DIMENSION(nsize) :: idumireal
!
!--debugging information
!
 IF (trace) WRITE(iprint,*) ' Allocating memory for all arrays'
!
!--set size of neighbour list
!
 nlistdim = nsize/2	! ie. up to half of particles can be neighbours
 IF (icall.EQ.1) THEN
!
!--leave extra space if ghost particles are used
!
    IF (ibound.GT.1) THEN
       ntot = nsize + 25 	! allow space for 25 ghosts
    ELSE
       ntot = nsize
    ENDIF
 
 ELSEIF (icall.EQ.2) THEN

!-----------------------------------------------------------------------------
!  copy everything into dummy arrays with the old size
!-----------------------------------------------------------------------------
    dumpmass = pmass(1:nsize)
    dumrhoin = rhoin(1:nsize)
    dumuuin = uuin(1:nsize)
    dumhhin = hhin(1:nsize)
    dumenin = enin(1:nsize)
    dumalphain = alphain(1:nsize)
    
    dumxin = xin(:,1:nsize)
    dumvelin = velin(:,1:nsize)
    dumBfieldin = Bfieldin(:,1:nsize)
    dumBin = Bin(:,1:nsize)
    
    dumx = x(:,1:nsize)
    dumvel = vel(:,1:nsize)
    dumforce = force(:,1:nsize)
    dumBfield = Bfield(:,1:nsize)
    dumdBfielddt = dBfielddt(:,1:nsize)
    dumfmag = fmag(:,1:nsize)
    dumcurlB = curlB(:,1:nsize)
    dumxsphterm = xsphterm(:,1:nsize)
        
    dumrho = rho(1:nsize)
    dumdrhodt = drhodt(1:nsize)
    dumuu = uu(1:nsize)
    dumdudt = dudt(1:nsize)
    dumen = en(1:nsize)
    dumdendt = dendt(1:nsize)
    dumalpha = alpha(1:nsize)
    dumdaldt = daldt(1:nsize)
    dumhh = hh(1:nsize)
    dumdhdt = dhdt(1:nsize)
    dumgradh = gradh(1:nsize)
    dumpr = pr(1:nsize)
    dumspsound = spsound(1:nsize)
!    dumll = ll(1:nsize)
    dumdivB = divB(1:nsize)
    idumireal = ireal(1:nsize)
    IF (ALLOCATED(fgrav)) dumfgrav(:,1:nsize) = fgrav(:,1:nsize)
    
 
!-----------------------------------------------------------------------------
!  deallocate the arrays
!-----------------------------------------------------------------------------
!
!--initial particle properties
!
    DEALLOCATE (pmass,xin,rhoin,uuin,hhin,enin,alphain)
    DEALLOCATE (velin,Bfieldin,Bin)
!
!--particle properties and derivatives
!
    DEALLOCATE(x,vel,force,rho,drhodt,uu,dudt,en,dendt)
    DEALLOCATE(alpha,daldt,hh,dhdt,gradh,pr)
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
!--MHD quantities and derivatives
!
    DEALLOCATE(Bfield,dBfielddt,fmag,divB,curlB)
!
!--XSPH
!
    DEALLOCATE(xsphterm)
!
!--gravity
!
    IF (ALLOCATED(fgrav)) DEALLOCATE(fgrav)


!-----------------------------------------------------------------------------
!  reallocate arrays to new size (add 10% so don't have to do this too often)
!-----------------------------------------------------------------------------
 
    IF (ibound.GT.1) THEN
       ntot = INT(1.1*nsize)
    ELSE 
       ntot = nsize
    ENDIF
    IF (ntot.GT.nsize) THEN   
       WRITE(iprint,*) ' old array size = ',nsize,' new array size = ',ntot
    ENDIF
 ENDIF

!-----------------------------------------------------------------------------
!  allocate all arrays (for both icall = 1 and icall = 2)
!-----------------------------------------------------------------------------

!
!--initial particle properties
!
    ALLOCATE (pmass(ntot))
    ALLOCATE (xin(ndim,ntot))
    ALLOCATE (rhoin(ntot),uuin(ntot),hhin(ntot))
    ALLOCATE (enin(ntot),alphain(ntot))
    ALLOCATE (velin(ndimV,ntot),Bfieldin(ndimB,ntot),Bin(ndimB,ntot))
!
!--particle properties and derivatives
!
    ALLOCATE(x(ndim,ntot))
    ALLOCATE(vel(ndimV,ntot),force(ndimV,ntot))
    ALLOCATE(rho(ntot),drhodt(ntot))
    ALLOCATE(uu(ntot),dudt(ntot),en(ntot),dendt(ntot))
    ALLOCATE(alpha(ntot),daldt(ntot))
    ALLOCATE(hh(ntot),dhdt(ntot),gradh(ntot),pr(ntot))
!
!--equation of state
!
    ALLOCATE(spsound(ntot))
!
!--linklist and boundaries - note ifirstincell is reallocated in link
!  (depends on max # of link list cells)
!
    ALLOCATE(ll(ntot),iamincell(ntot),ifirstincell(0),ireal(ntot))
!
!--MHD quantities and derivatives
!
    ALLOCATE(Bfield(ndimB,ntot),dBfielddt(ndimB,ntot),fmag(ndimB,ntot))	! mag field
    ALLOCATE(divB(ntot),curlB(ndimB,ntot))
!
!--XSPH
!
    ALLOCATE(xsphterm(ndimV,ntot))
    
    IF (igravity.NE.0) THEN
       ALLOCATE(fgrav(ndim,ntot))
    ENDIF
    
!-----------------------------------------------------------------------------
!  copy properties back to new, bigger arrays
!-----------------------------------------------------------------------------

 IF (icall.EQ.2) THEN
    pmass(1:nsize) = dumpmass(1:nsize)
    rhoin(1:nsize) = dumrhoin(1:nsize)
    uuin(1:nsize) = dumuuin(1:nsize)
    hhin(1:nsize) = dumhhin(1:nsize)
    enin(1:nsize) = dumenin(1:nsize)
    alphain(1:nsize) = dumalphain(1:nsize)
    
    xin(:,1:nsize) = dumxin(:,1:nsize)
    velin(:,1:nsize) = dumvelin(:,1:nsize)
    Bfieldin(:,1:nsize) = dumBfieldin(:,1:nsize)
    Bin(:,1:nsize) = dumBin(:,1:nsize)
    
    x(:,1:nsize) = dumx(:,1:nsize)
    vel(:,1:nsize) = dumvel(:,1:nsize)
    force(:,1:nsize) = dumforce(:,1:nsize)
    Bfield(:,1:nsize) = dumBfield(:,1:nsize)
    dBfielddt(:,1:nsize) = dumdBfielddt(:,1:nsize)
    fmag(:,1:nsize) = dumfmag(:,1:nsize)
    curlB(:,1:nsize) = dumcurlB(:,1:nsize)
    xsphterm(:,1:nsize) = dumxsphterm(:,1:nsize)
        
    rho(1:nsize) = dumrho(1:nsize)
    drhodt(1:nsize) = dumdrhodt(1:nsize)
    uu(1:nsize) = dumuu(1:nsize)
    dudt(1:nsize) = dumdudt(1:nsize)
    en(1:nsize) = dumen(1:nsize)
    dendt(1:nsize) = dumdendt(1:nsize)
    alpha(1:nsize) = dumalpha(1:nsize)
    daldt(1:nsize) = dumdaldt(1:nsize)
    hh(1:nsize) = dumhh(1:nsize)
    dhdt(1:nsize) = dumdhdt(1:nsize)
    gradh(1:nsize) = dumgradh(1:nsize)
    pr(1:nsize) = dumpr(1:nsize)
    spsound(1:nsize) = dumspsound(1:nsize)
!    ll(1:nsize) = dumll(1:nsize)
    divB(1:nsize) = dumdivB(1:nsize) 
    ireal(1:nsize) = idumireal(1:nsize)
    IF (ALLOCATED(fgrav)) fgrav(:,1:nsize) = dumfgrav(:,1:nsize)
 ENDIF   
       
!
!--debugging information
!
 IF (trace) WRITE(iprint,*) '  Memory allocated, size(x,v)=',SIZE(x),SIZE(vel) 
 
 RETURN

END SUBROUTINE alloc
