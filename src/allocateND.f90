!!----------------------------------------------------------------------
!! This subroutine allocates memory for all arrays once the number of
!! particles is known (should be called from setup).
!!
!! oldsize = 0 : allocate arrays for the first time (oldsize irrelevant)
!! oldsize .ne. 0 : reallocate arrays because particle # > array size
!!
!!----------------------------------------------------------------------
SUBROUTINE alloc(oldsize)
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
 INTEGER, INTENT(IN) :: oldsize
 INTEGER :: newsize
 REAL, DIMENSION(oldsize) :: dumpmass,dumrhoin,dumuuin,dumenin,dumalphain
 REAL, DIMENSION(oldsize) :: dumrho,dumdrhodt,dumuu,dumdudt,dumen,dumdendt
 REAL, DIMENSION(oldsize) :: dumalpha,dumdaldt,dumhh,dumgradh,dumpr,dumspsound
 REAL, DIMENSION(oldsize) :: dumdivB,dumhhin,dumdhdt
 REAL, DIMENSION(ndim,oldsize) :: dumxin,dumx,dumfgrav
 REAL, DIMENSION(ndimV,oldsize) :: dumvelin,dumvel,dumBfieldin,dumBin
 REAL, DIMENSION(ndimV,oldsize) :: dumforce,dumBfield,dumdBfielddt,dumfmag
 REAL, DIMENSION(ndimV,oldsize) :: dumcurlB,dumxsphterm
 INTEGER, DIMENSION(oldsize) :: idumireal,idumitype
 LOGICAL :: reallocate, increasingsize
!
!--set size of neighbour list
!
 nlistdim = oldsize/2	! ie. up to half of particles can be neighbours
 reallocate = .false.
 IF (ALLOCATED(x)) reallocate = .true.
 IF (trace .AND. reallocate) THEN
    WRITE(iprint,*) 'Reallocating memory for all arrays'
 ELSEIF (trace) THEN
    WRITE(iprint,*) 'Allocating memory for all arrays'
 ENDIF
 increasingsize = .true.
 IF (oldsize.LT.SIZE(rho)) THEN
    increasingsize = .false.
    WRITE(iprint,*) 'decreasing array size'
 ENDIF
 IF (oldsize.EQ.0) THEN
    WRITE(iprint,*) 'Error: size=0 in call to allocate'
    CALL quit
 ENDIF
 
 IF (.not.reallocate) THEN
!
!--leave extra space if ghost particles are used
!
    IF (ibound.GT.1) THEN
       newsize = oldsize + 25 	! allow space for 25 ghosts
    ELSE
       newsize = oldsize		! else just use size given
    ENDIF
 
 ELSE

!-----------------------------------------------------------------------------
!  copy everything into dummy arrays with the old size
!-----------------------------------------------------------------------------
    dumpmass = pmass(1:oldsize)
    dumrhoin = rhoin(1:oldsize)
    dumuuin = uuin(1:oldsize)
    dumhhin = hhin(1:oldsize)
    dumenin = enin(1:oldsize)
    dumalphain = alphain(1:oldsize)
    
    dumxin = xin(:,1:oldsize)
    dumvelin = velin(:,1:oldsize)
    dumBfieldin = Bfieldin(:,1:oldsize)
    dumBin = Bin(:,1:oldsize)
    
    dumx = x(:,1:oldsize)
    dumvel = vel(:,1:oldsize)
    dumforce = force(:,1:oldsize)
    dumBfield = Bfield(:,1:oldsize)
    dumdBfielddt = dBfielddt(:,1:oldsize)
    dumfmag = fmag(:,1:oldsize)
    dumcurlB = curlB(:,1:oldsize)
    dumxsphterm = xsphterm(:,1:oldsize)
        
    dumrho = rho(1:oldsize)
    dumdrhodt = drhodt(1:oldsize)
    dumuu = uu(1:oldsize)
    dumdudt = dudt(1:oldsize)
    dumen = en(1:oldsize)
    dumdendt = dendt(1:oldsize)
    dumalpha = alpha(1:oldsize)
    dumdaldt = daldt(1:oldsize)
    dumhh = hh(1:oldsize)
    dumdhdt = dhdt(1:oldsize)
    dumgradh = gradh(1:oldsize)
    dumpr = pr(1:oldsize)
    dumspsound = spsound(1:oldsize)
!    dumll = ll(1:oldsize)
    dumdivB = divB(1:oldsize)
    idumireal = ireal(1:oldsize)
    idumitype = itype(1:oldsize)
    IF (ALLOCATED(fgrav)) dumfgrav(:,1:oldsize) = fgrav(:,1:oldsize)
    
 
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
!--itype is particle type (normal, fixed particle etc)
!
    IF (ALLOCATED(itype)) DEALLOCATE(itype)
!
!--MHD quantities and derivatives
!
    DEALLOCATE(Bfield,dBfielddt,fmag,divB,curlB)
!
!--XSPH
!
    IF (ALLOCATED(xsphterm)) DEALLOCATE(xsphterm)
!
!--gravity
!
    IF (ALLOCATED(fgrav)) DEALLOCATE(fgrav)


!-----------------------------------------------------------------------------
!  reallocate arrays to new size (add 10% so don't have to do this too often)
!-----------------------------------------------------------------------------
 
    IF (ibound.GE.1 .AND. increasingsize) THEN
       newsize = INT(1.1*oldsize)
    ELSE 
       newsize = oldsize	! decreasing size
    ENDIF
    IF (newsize.NE.oldsize) THEN   
       WRITE(iprint,*) ' old array size = ',oldsize,' new array size = ',newsize
    ENDIF
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
    ALLOCATE (velin(ndimV,newsize),Bfieldin(ndimB,newsize),Bin(ndimB,newsize))
!
!--particle properties and derivatives
!
    ALLOCATE(x(ndim,newsize))
    ALLOCATE(vel(ndimV,newsize),force(ndimV,newsize))
    ALLOCATE(rho(newsize),drhodt(newsize))
    ALLOCATE(uu(newsize),dudt(newsize),en(newsize),dendt(newsize))
    ALLOCATE(alpha(newsize),daldt(newsize))
    ALLOCATE(hh(newsize),dhdt(newsize),gradh(newsize),pr(newsize))
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
    ALLOCATE(Bfield(ndimB,newsize),dBfielddt(ndimB,newsize),fmag(ndimB,newsize))	! mag field
    ALLOCATE(divB(newsize),curlB(ndimB,newsize))
!
!--XSPH
!
    ALLOCATE(xsphterm(ndimV,newsize))
    
    IF (igravity.NE.0) THEN
       ALLOCATE(fgrav(ndim,newsize))
    ENDIF
    
!-----------------------------------------------------------------------------
!  copy properties back to new, bigger arrays
!-----------------------------------------------------------------------------

 IF (reallocate) THEN
    pmass(1:oldsize) = dumpmass(1:oldsize)
    rhoin(1:oldsize) = dumrhoin(1:oldsize)
    uuin(1:oldsize) = dumuuin(1:oldsize)
    hhin(1:oldsize) = dumhhin(1:oldsize)
    enin(1:oldsize) = dumenin(1:oldsize)
    alphain(1:oldsize) = dumalphain(1:oldsize)
    
    xin(:,1:oldsize) = dumxin(:,1:oldsize)
    velin(:,1:oldsize) = dumvelin(:,1:oldsize)
    Bfieldin(:,1:oldsize) = dumBfieldin(:,1:oldsize)
    Bin(:,1:oldsize) = dumBin(:,1:oldsize)
    
    x(:,1:oldsize) = dumx(:,1:oldsize)
    vel(:,1:oldsize) = dumvel(:,1:oldsize)
    force(:,1:oldsize) = dumforce(:,1:oldsize)
    Bfield(:,1:oldsize) = dumBfield(:,1:oldsize)
    dBfielddt(:,1:oldsize) = dumdBfielddt(:,1:oldsize)
    fmag(:,1:oldsize) = dumfmag(:,1:oldsize)
    curlB(:,1:oldsize) = dumcurlB(:,1:oldsize)
    xsphterm(:,1:oldsize) = dumxsphterm(:,1:oldsize)
        
    rho(1:oldsize) = dumrho(1:oldsize)
    drhodt(1:oldsize) = dumdrhodt(1:oldsize)
    uu(1:oldsize) = dumuu(1:oldsize)
    dudt(1:oldsize) = dumdudt(1:oldsize)
    en(1:oldsize) = dumen(1:oldsize)
    dendt(1:oldsize) = dumdendt(1:oldsize)
    alpha(1:oldsize) = dumalpha(1:oldsize)
    daldt(1:oldsize) = dumdaldt(1:oldsize)
    hh(1:oldsize) = dumhh(1:oldsize)
    dhdt(1:oldsize) = dumdhdt(1:oldsize)
    gradh(1:oldsize) = dumgradh(1:oldsize)
    pr(1:oldsize) = dumpr(1:oldsize)
    spsound(1:oldsize) = dumspsound(1:oldsize)
!    ll(1:oldsize) = dumll(1:oldsize)
    divB(1:oldsize) = dumdivB(1:oldsize) 
    ireal(1:oldsize) = idumireal(1:oldsize)
    itype(1:oldsize) = idumitype(1:oldsize)
    IF (ALLOCATED(fgrav)) fgrav(:,1:oldsize) = dumfgrav(:,1:oldsize)
 ELSE
    itype(:) = 0	! on first memory allocation, set all parts = normal
 ENDIF   
       
!
!--debugging information
!
 IF (trace) WRITE(iprint,*) '  Memory allocated, size(x,v)=',SIZE(x),SIZE(vel) 
 
 RETURN

END SUBROUTINE alloc
