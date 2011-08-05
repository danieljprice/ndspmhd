!!--------------------------------------------------------------------
!!  Subroutine to write output to data file
!!  (change this to change format of output)
!!-------------------------------------------------------------------

SUBROUTINE output(t,nstep)
 USE dimen_mhd
 USE debug
 USE loguns
 
 USE eos
 USE options
 USE part
 USE setup_params
 USE derivB
 USE rates
!
!--define local variables
!
 IMPLICIT NONE
 REAL, INTENT(IN) :: t
! REAL, PARAMETER :: pi=3.1415926536
 INTEGER, INTENT(IN) :: nstep
 INTEGER :: nprint,ndata
 INTEGER :: i,ierr(3)
 INTEGER, EXTERNAL :: flush
!
!--allow for tracing flow
!      
 IF (trace) WRITE(iprint,*) ' Entering subroutine output'

 IF (idumpghost.EQ.1) THEN
    nprint = ntotal
 ELSE
    nprint = npart
 ENDIF

!
!--write timestep to log file
!      
 WRITE (iprint,10) t,abs(nstep),npart,ntotal-npart
10 FORMAT('| time = ',f7.3,' | timesteps = ',i8,' | npart = ',i7,	&
             ' | nghost = ',i5,' |')
!
!--write header line to this data block in data file
!
 IF (imhd.NE.0) THEN
    ndata = ndim + 10 + 2*ndimB + ndimV	+ndimV ! number of columns apart from co-ords
 ELSE
    ndata = ndim + 7 + ndimV + ndimV
 ENDIF
 WRITE(idatfile,20) t,npart,nprint,gamma,hfact,ndim,ndimV,ndata
20 FORMAT(e12.5,1x,i8,1x,i8,1x,f14.12,1x,f6.2,1x,i1,1x,i1,1x,i3)

!--calculate primitive variables from conservatives for output
!  nstep <0 means do not do this as we are on a quit dump
 IF (nstep.ge.0) CALL conservative2primitive  ! also calls equation of state
!
!--write data for this time to data file
!      
! scale = MAXVAL(Bfield(2,:))
 
 DO i=1,nprint
!
!--write the data (primitive variables) to the .dat file
!
    IF (imhd.NE.0) THEN	! MHD

       WRITE(idatfile,30) x(:,i),vel(:,i),dens(i),pr(i),uu(i),hh(i),   &
        pmass(i),alpha(:,i),Bfield(:,i),divB(i),curlB(:,i),psi(i),force(:,i)

    ELSE   ! non-MHD

       WRITE(idatfile,30) x(:,i),vel(:,i),dens(i),pr(i),uu(i),hh(i),   &                        
        pmass(i),alpha(1:2,i),force(:,i)

    ENDIF

 ENDDO
!
!--flush the buffer so that whole timestep is printed even if program crashes
! 
 ierr(1) = flush(idatfile)
 ierr(2) = flush(ievfile)
 ierr(3) = flush(iprint)
 IF ( ANY(ierr.NE.0) ) WRITE(*,*) 'Error flushing files, see Makefile ',ierr 
  			  
30    FORMAT (25(1pe21.14,1x),:)	! make sure the format statement has >/=	
				! max number of columns in the write statement

 RETURN
END SUBROUTINE output
      
