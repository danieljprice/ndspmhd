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
 REAL :: vpar,Bpar,vperp,Bperp,angle
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
 WRITE (iprint,99002) t,nstep,npart,ntotal-npart
99002 FORMAT('| time = ',f7.3,' | timesteps = ',i8,' | npart = ',i7,	&
             ' | nghost = ',i5,' |')
!
!--write header line to this data block in data file
!
 IF (imhd.NE.0) THEN
    ndata = ndim + 11 + 2*ndimB + ndimV	! number of columns apart from co-ords
 ELSE
    ndata = ndim + 8 + ndimV
 ENDIF
 WRITE(idatfile,*) t,npart,nprint,gamma,hfact,ndim,ndimV,ndata
!
!--write data for this time to data file
!      
! scale = MAXVAL(Bfield(2,:)) 
 DO i=1,nprint
    angle = 30.*pi/180.

    IF (imhd.GE.11) THEN	! if Bfield is B

       vpar = vel(2,i)*SIN(angle) + vel(1,i)*COS(angle)
       vperp = vel(2,i)*COS(angle) - vel(1,i)*SIN(angle)
       Bpar = Bfield(2,i)*SIN(angle) + Bfield(1,i)*COS(angle)
       Bperp = Bfield(2,i)*COS(angle) - Bfield(1,i)*SIN(angle)

       WRITE(idatfile,30) x(:,i),vel(:,i),rho(i),pr(i),uu(i),hh(i),   &
        pmass(i),alpha(i),Bfield(:,i),divB(i),curlB(:,i),	      &
	vpar,vperp,Bpar,Bperp

    ELSEIF (imhd.GT.0) THEN	! if Bfield is B/rho

       WRITE(idatfile,30) x(:,i),vel(:,i),rho(i),pr(i),uu(i),hh(i),   &                        
        pmass(i),alpha(i),Bfield(:,i)*rho(i),divB(i),curlB(:,i),      &
	vpar,vperp,Bpar,Bperp

    ELSE   ! non-MHD

       WRITE(idatfile,30) x(:,i),vel(:,i),rho(i),pr(i),uu(i),hh(i),   &                        
        pmass(i),alpha(i),vpar,vperp

    ENDIF

 ENDDO
!
!--flush the buffer so that whole timestep is printed even if program crashes
! 
 ierr(1) = flush(idatfile)
 ierr(2) = flush(ievfile)
 ierr(3) = flush(iprint)
 IF ( ANY(ierr.NE.0) ) WRITE(*,*) 'Error flushing files, see Makefile ',ierr 
  			  
30    FORMAT (23(1pe14.6,1x),:)	! make sure the format statement has >/=	
				! max number of columns in the write statement

 RETURN
END SUBROUTINE output
      
