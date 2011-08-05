!!------------------------------------------------------------------------!!
!!                                                                        !!
!!  Generic setup for a uniform density distribution of particles         !!
!!   in cartesian geometry for 1, 2 and 3 dimensions                      !!
!!                                                                        !!
!!  In 3D, lattice can be close packed, body centred or cubic             !!
!!     2D, lattice can be close packed or cubic                           !!
!!     1D, all setups are the same                                        !!
!!                                                                        !!
!!------------------------------------------------------------------------!!

SUBROUTINE set_uniform_cartesian(idistin,psep,xmin,xmax,offset)
!
!--include relevant global variables
!
 USE dimen_mhd
 USE debug
 USE loguns

 USE options
 USE part
!
!--define local variables
!  (note we read boundaries of region as input, so that more than one region
!  can be defined in the particle setup)
! 
 IMPLICIT NONE
 INTEGER, INTENT(IN) :: idistin
 REAL, DIMENSION(ndim), INTENT(INOUT) :: xmin, xmax
 REAL, INTENT(IN) :: psep
 LOGICAL, INTENT(IN) :: offset
 
 INTEGER :: i,j,k,ntot,npartin,npartx,nparty,npartz,ipart,iseed
 INTEGER :: idist
 REAL :: xstart,ystart,deltax,deltay
 REAL :: psepx,psepy
 REAL :: ran1,ampl
 REAL, DIMENSION(ndim) :: xran
!
!--allow for tracing flow
!
 IF (trace) WRITE(iprint,*) ' Entering subroutine uniform_cartesian ',idistin
 WRITE(iprint,*) 'Uniform cartesian distribution '

 IF (ndim.EQ.1) THEN 
    idist = 1   ! in 1D use default version
 ELSE
    idist = idistin
 ENDIF
 
 npartin = npart    ! this is how many particles have already been set up
 
 SELECT CASE(idist)
!
!--hexagonal close packed arrangement 
!  (in 3D, 2D this is the same as body centred)
!
 CASE(2,12)
 
    IF (ndim.EQ.3) STOP 'close packed not implemented in 3D'
!
!--determine number of particles
! 
    deltax = psep
    deltay = 0.5*SQRT(3.)*psep
    npartx = INT((xmax(1)-xmin(1))/deltax)
    nparty = INT((xmax(2)-xmin(2))/deltay)
    npartz = 1
    
!--for periodic boundaries, ymax needs to be divisible by 2
    if (ibound(2).eq.3) then
       nparty = 2*INT(nparty/2)
       print*,' periodic boundaries: adjusting nparty = ',nparty
    endif
!
!--adjust psep so that particles fill the volume
!
    PRINT*,' npartx,y = ',npartx,nparty  !!,deltax,deltay
!    deltax = (xmax(1)-xmin(1))/(FLOAT(npartx))
!    deltay = (xmax(2)-xmin(2))/(FLOAT(nparty))
!    PRINT*,' adjusted ',deltax,deltay
!
!--or adjust the boundaries appropriately
!
    xmax(2) = xmin(2) + nparty*deltay
    PRINT*,' adjusted y boundary : ymax  = ',xmax(2)
!
!--allocate memory here
!
    ntot = npartx*nparty + npartin
    
    WRITE(iprint,*) ' hexagonal close packed distribution, npart = ',ntot
    
    CALL alloc(ntot)
    npart = ntot

    ystart = 0.5*deltay      !psepy
    ipart = npartin
    DO k=1,npartz
       DO j=1,nparty
          xstart = 0.25*psep
          IF (MOD(j,2).EQ.0) xstart = xstart + 0.5*psep
          DO i = 1,npartx
             ipart = ipart + 1      !(k-1)*nparty + (j-1)*npartx + i
             x(1,ipart) = xmin(1) + (i-1)*deltax + xstart
             IF (ndim.GE.2) x(2,ipart) = xmin(2) + (j-1)*deltay + ystart
          ENDDO
       ENDDO
    ENDDO 
!
!--body centred lattice
!
 CASE(3,13)
    IF (ndim.EQ.3) STOP 'body centred not implemented in 3D'
!
!--determine number of particles
! 
    npartx = INT((xmax(1)-xmin(1))/(SQRT(2.)*psep))
    nparty = INT((xmax(2)-xmin(2))/(SQRT(2.)*psep))
    npartz = 1
!
!--adjust psep so that particles fill the volume
!
    psepx = (xmax(1)-xmin(1))/(FLOAT(npartx)*SQRT(2.))
    psepy = (xmax(2)-xmin(2))/(FLOAT(nparty)*SQRT(2.))
!    PRINT*,'psep = ',psepx,psepy,psep
    deltax = SQRT(2.)*psepx
    deltay = 0.5*SQRT(2.)*psepy    
!
!--allocate memory here
!
    ntot = 2*npartx*nparty
    WRITE(iprint,*) 'Body centred distribution, npart = ',ntot
    
    CALL alloc(ntot)
    npart = ntot

    ystart = 0.5*deltay      !psepy
    ipart = 0
    DO k=1,npartz
       DO j=1,2*nparty
          xstart = 0.25*deltax
          IF (MOD(j,2).EQ.0) xstart = xstart + 0.5*deltax
          DO i = 1,npartx
             ipart = ipart + 1      !(k-1)*nparty + (j-1)*npartx + i
             x(1,ipart) = xmin(1) + (i-1)*deltax + xstart
             IF (ndim.GE.2) x(2,ipart) = xmin(2) + (j-1)*deltay + ystart
          ENDDO
       ENDDO
    ENDDO 
!
!--random particle distribution
!  (uses random number generator ran1 in module random)
! 
 CASE(4,14)
    ntot = INT(PRODUCT((xmax(:)-xmin(:))/psep))
    npart = ntot
    WRITE(iprint,*) 'Random particle distribution, npart = ',ntot 
    CALL alloc(ntot)
!
!--initialise random number generator
!    
    iseed = -87682
    xran(1) = ran1(iseed)
    
    DO i=1,ntot
       DO j=1,ndim
          xran(j) = ran1(iseed)
       ENDDO
       x(:,i) = xmin(:) + xran(:)*(xmax(:)-xmin(:))
    ENDDO

 CASE DEFAULT
!----------------------
!  cubic lattice
!----------------------

    IF (offset) THEN
       WRITE(iprint,*) 'Offset lattice'
       xmin(:) = xmin(:) + 0.5*psep
       xmax(:) = xmax(:) + 0.5*psep
    ENDIF
!
!--determine number of particles
! 
    npartx = INT((xmax(1)-xmin(1))/psep)    
    IF (ndim.GE.2) THEN
       nparty = INT((xmax(2)-xmin(2))/psep)    
    ELSE
       nparty = 1
    ENDIF
    IF (ndim.GE.3) THEN
       npartz = INT((xmax(3)-xmin(3))/psep)    
    ELSE
       npartz = 1    
    ENDIF
!
!--allocate memory here
!
    ntot = npartx*nparty*npartz
    ipart = npartin
    npart = npartin + ntot   ! add to particles already setup
    WRITE(iprint,*) 'Cubic lattice, adding ',ntot,', npart = ',npart,npartx,nparty,npartz
    WRITE(iprint,*) 'xmin = ',xmin, ' xmax = ',xmax
    CALL alloc(npart)

    DO k=1,npartz
       DO j=1,nparty
          DO i = 1,npartx
             ipart = ipart + 1   !(k-1)*nparty + (j-1)*npartx + i
             x(1,ipart) = xmin(1) + (i-1)*psep + 0.5*psep
             IF (ndim.GE.2) THEN
                x(2,ipart) = xmin(2) + (j-1)*psep + 0.5*psep
                IF (ndim.GE.3) THEN
                   x(3,ipart) = xmin(3) + (k-1)*psep + 0.5*psep
                ENDIF
             ENDIF
!           print*,'new particle ',ipart,'x =', x(:,ipart)
         ENDDO
       ENDDO
    ENDDO
 END SELECT
 
!
!--if idist > 10 apply a small random perturbation to the particle positions
! 
 IF (idist.GE.10) THEN
!
!--initialise random number generator
!    
    iseed = -268
    xran(1) = ran1(iseed)
    ampl = 0.1*psep   ! perturbation amplitude
    
    DO i=npartin+1,npart  ! apply to new particles only
       DO j=1,ndim
          xran(j) = ran1(iseed)
       ENDDO
       x(:,i) = x(:,i) + ampl*(xran(:)-0.5)
    ENDDO    
 ENDIF
 
!
!--set total number of particles = npart
!
 ntotal = npart
!
!--allow for tracing flow
!
 IF (trace) WRITE(iprint,*) '  Exiting subroutine uniform_cartesian'
            
 RETURN
END SUBROUTINE
