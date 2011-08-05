!!------------------------------------------------------------------------
!! Iterates the density calculation so that rho and h are
!! calculated self-consistently. 
!!
!! rho_a is a function of h_a and vice versa
!! via the summation rho_a = sum_b m_b W_ab(h_a)
!! we assume h \propto 1/rho^(1/ndim)
!!
!! we determine the change in h and recalculate the density only 
!! for those particles where h changes significantly after the density
!! is calculated.
!!
!! For details see Price and Monaghan 2003b
!!------------------------------------------------------------------------

SUBROUTINE iterate_density
 USE dimen_mhd
 USE debug
 USE loguns

 USE bound
 USE hterms
 USE options
 USE part
 USE part_in	! for rhoin on fixed particles
 USE setup_params
!
!--define local variables
!
 IMPLICIT NONE
 INTEGER :: i,j,itsdensity,itsdensitymax
 INTEGER :: ncalc,ncalcprev,isize
 INTEGER, DIMENSION(:), ALLOCATABLE :: redolist
 REAL :: tol,hnew
 REAL, DIMENSION(:), ALLOCATABLE :: rho_old,gradh_old
 LOGICAL :: converged,redolink
!
!--allow for tracing flow
!      
 IF (trace) WRITE(iprint,*) ' Entering subroutine iterate_density'  
!
!--allocate memory for local array copies
!
 isize = SIZE(rho)
 ALLOCATE( rho_old(isize), gradh_old(isize) )
 ALLOCATE( redolist(ntotal) )
!
!--set maximum number of iterations to perform
! 
 IF ((ikernav.EQ.3).AND.(ihvar.NE.0)) THEN
    itsdensitymax = 20	! perform 1 fixed point iteration
 ELSE
    itsdensitymax = 0	! no iterations
 ENDIF
!
!--Loop to find rho and h self-consistently (if using Springel/Hernquist)
!
 itsdensity = 0
 tol = 1.e-2
 ncalc = npart	! number of particles to calculate density on
 redolink = .false.
 ncalcprev = 0

 iterate: DO WHILE ((ncalc.GT.0).AND.(itsdensity.LE.itsdensitymax)    &
                    .AND. ncalc.NE.ncalcprev)

    itsdensity = itsdensity + 1
    ncalcprev = ncalc
    IF (isize.NE.SIZE(rho)) THEN
       DEALLOCATE(rho_old,gradh_old)
       isize = SIZE(rho)
       ALLOCATE(rho_old(isize),gradh_old(isize))
    ENDIF
    rho_old = rho
    gradh_old = gradh
    IF (redolink) THEN
       IF (ANY(ibound.GT.1)) CALL set_ghost_particles
       IF (idebug(1:3).EQ.'den') WRITE(iprint,*) 'relinking...'
       CALL link
    ENDIF
    
    IF (ncalc.EQ.npart) THEN	! calculate density on all particles
       CALL density		! do symmetrically
    ELSE	! calculate density on partial list of particles              
       CALL density		!_partial(ncalc,redolist)
!       DO j=1,ncalc
!          i = redolist(j)
!          PRINT*,' normal  rho(',i,' =',rho(i),gradh(i)
!       ENDDO
!       gradh = gradh_old      
!       print*,'rho(717) 1 = ',rho(717)
!       rho = rho_old
!       print*,'rho(717) 2 = ',rho(717)
!       ncalc = npart
!       IF (ncalc.EQ.npart) THEN
!          DO j=1,npart     
!             redolist(j) = j
!          ENDDO
!       ENDIF
!       CALL density_partial(ncalc,redolist)
!       DO j=1,ncalc
!          i = redolist(j)
!          PRINT*,' partial rho(',i,' =',rho(i),gradh(i)
!       ENDDO             
!       READ*
    ENDIF

    ncalc = 0
    redolink = .false.

    IF (ihvar.NE.0) THEN
       DO i=1,npart
          IF (itype(i).NE.1) THEN
             IF (rho(i).LE.1.e-3) THEN
	        IF (rho(i).LE.0.) THEN
	           WRITE(iprint,*) 'rho: rho -ve, dying rho(',i,') = ',rho(i)
	           CALL quit
	        ELSE
	           WRITE(iprint,*) 'Warning : rho < 1e-5 '
	        ENDIF
             ENDIF

             gradh(i) = gradh(i)/rho_old(i)    ! now that rho is known
             IF (abs(1.-gradh(i)).LT.1.e-5) PRINT*,'eek'
!
!--perform Newton-Raphson iteration on rho
!      
!             rho(i) = rho_old(i) - (rho_old(i) - rho(i))/(1.-gradh(i))
!
!--work out new smoothing length h
!
             hnew = hfact*(pmass(i)/rho(i))**hpower	! ie h proportional to 1/rho^dimen
!
!--if this particle is not converged, add to list of particles to recalculate
!
!            PRINT*,'hnew - hh(i) = ',abs(hnew-hh(i))/hh(i)

             converged = abs((hnew-hh(i))/hh(i)) < tol	  
	     IF (.not.converged) THEN
	        ncalc = ncalc + 1
	        redolist(ncalc) = i
!	        PRINT*,'not converged',i,abs(hnew-hh(i))/hh(i),rho(i),	&
!     	        ncalc,redolist(ncalc)
!
!--update smoothing length only if taking another iteration
!
 	        IF (itsdensity.LT.itsdensitymax) hh(i) = hnew
	        IF (hnew.GT.hhmax) THEN
		   redolink = .true.
	        ENDIF	
	     ENDIF
          ENDIF	! itype .NE. 1
       ENDDO

       IF ((idebug(1:3).EQ.'den').AND.(ncalc.GT.0)) THEN
	  WRITE(iprint,*) ' density, iteration ',itsdensity,' ncalc = ',ncalc
!         DO i=1,ncalc
!	     PRINT*,' i = ',redolist(i)
!	  ENDDO
       ENDIF  

    ENDIF

    IF (ANY(ibound.GT.1)) THEN
       DO i=npart+1,ntotal		! update ghosts
          j = ireal(i)
          rho(i) = rho(j)
          hh(i) = hh(j)
          gradh(i) = gradh(j)       
       ENDDO
    ELSEIF (ANY(ibound.EQ.1)) THEN		! update fixed particles
       WHERE (itype(:).EQ.1)
          rho(:) = rhoin(:)
       END WHERE
    ENDIF
            
 ENDDO iterate

 IF (itsdensity.GT.1) WRITE(iprint,*) ' Finished density, iterations = ',itsdensity
 
 RETURN
END SUBROUTINE iterate_density
