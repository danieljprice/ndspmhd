!!------------------------------------------------------------------------!!
!!                                                                        !!
!!  Setup by reading a dump from a file                                   !!
!!  Should be compatible with the output format given in output           !!
!!                                                                        !!
!!------------------------------------------------------------------------!!

SUBROUTINE setup
!
!--include relevant global variables
!
 USE dimen_mhd
 USE debug
 USE loguns
 
 USE bound
 USE derivB
 USE eos
 USE options
 USE part
 USE setup_params
!
!--define local variables
!      
 IMPLICIT NONE
 INTEGER :: i,j,ndimfile,ndimVfile,npartfile,nprintfile,ncolumns
 REAL :: tfile,gammafile,hfactfile
 CHARACTER(LEN=30) :: setupfile
 LOGICAL :: iexist, mhdfile

!
!--set name of setup file
!
 setupfile = 'tstar2D_static.dat'

 WRITE(iprint,*) 'Reading setup from file: ',TRIM(setupfile)
 INQUIRE (file=setupfile,exist=iexist)
 IF (.not.iexist) THEN
    WRITE(iprint,*) 'Can''t find setup file: ',TRIM(setupfile)    
    STOP
 ENDIF
!
!--open setup file
!
 OPEN(UNIT=ireadf,FILE=setupfile,ERR=666,STATUS='old',FORM='formatted')
!
!--read header line
!
 READ(ireadf,*,ERR=667) tfile,npartfile,nprintfile,gammafile,hfactfile, &
      ndimfile,ndimVfile,ncolumns
 WRITE(iprint,*) 'time = ',tfile,' in setup file '
!
!--check for compatibility with current settings
!
 IF (ndimfile.NE.ndim) STOP 'x dimensions not equal between setup file and code'
 IF (ndimVfile.NE.ndimV) STOP 'v dimensions not equal between setup file and code'
 IF (abs(gammafile-gamma).gt.1.e-3) WRITE(iprint,10) 'gamma',gammafile,gamma
 IF (abs(hfactfile-hfact).gt.1.e-3) WRITE(iprint,10) 'hfact',hfactfile,hfact
10 FORMAT(/,'WARNING: ',a,' changed from original setup: old = ',f9.6,' new = ',f9.6,/)
 
 IF ((nprintfile.NE.npartfile).AND.(ALL(ibound.LE.1))) THEN
    WRITE(iprint,*) 'WARNING: setup file contains ghosts, but none set'
 ENDIF
 IF (ncolumns.EQ.(ndim+6+ndimV)) THEN
    mhdfile = .false.
    WRITE(iprint,*) '(non-MHD input file)'
    IF (imhd.GT.0) WRITE(iprint,*) 'WARNING: reading non-MHD file, but MHD is on'
 ELSEIF (ncolumns.EQ.(ndim+7+3*ndimV)) THEN
    mhdfile = .true.
    WRITE(iprint,*) '(MHD input file)'
    IF (imhd.LE.0) WRITE(iprint,*) 'WARNING: MHD input file, but MHD is off'       
 ELSE
    STOP 'error: unknown identity of columns in input file'
 ENDIF
!
!--allocate memory for the number of particles
!  do not read any ghost particles
!
 npart = npartfile
 ntotal = npartfile
 
 CALL alloc(ntotal)
!
!--read one step (only) from file
!
 DO i=1,npart
    IF (mhdfile) THEN
       READ(ireadf,*,ERR=668) x(:,i),vel(:,i),rho(i),pr(i),uu(i),hh(i),   &
            pmass(i),alpha(:,i),Bfield(:,i),divB(i),curlB(:,i)
       IF (imhd.EQ.0) THEN
          Bfield(:,i) = 0.
          divB(i) = 0.
          curlB(:,i) = 0.
       ENDIF
    ELSE
       READ(ireadf,*,ERR=668) x(:,i),vel(:,i),rho(i),pr(i),uu(i),hh(i),   &                        
            pmass(i),alpha(:,i)
       IF (imhd.NE.0) THEN
          Bfield(:,i) = 0.    ! could change this
          divB(i) = 0.
          curlB(:,i) = 0.
       ENDIF
    ENDIF
 ENDDO
!
!--close the file
!
 CLOSE(UNIT=ireadf)
 WRITE(iprint,*) 'Finished reading setup file: everything is AOK'
!
!--set bounds of setup
!                   
 xmax = 0.   ! irrelevant
 xmin = 0.
 ibound = 0	! no boundaries 
 iexternal_force = 1	! use toy star force
!
!--now change things according to the specific setup required
!
 vel(:,:) = 0.
 WHERE (rho > 0.)
    hh = hfact*(pmass(i)/rho(:))**hpower
 END WHERE
 
 RETURN

666 STOP 'error opening setup file'
667 STOP 'error reading header line from setup file'
668 STOP 'error reading timestep from setup file'

END SUBROUTINE setup
