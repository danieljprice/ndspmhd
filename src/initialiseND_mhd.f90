!!-----------------------------------------------------------------
!! This subroutine initialises everything needed for the run
!!-----------------------------------------------------------------

SUBROUTINE initialise
!
!--include relevant global variables
!
! USE dimen_mhd
 USE debug
 USE loguns
 
 USE artvi
 USE bound
 USE derivB
 USE eos
 USE fmagarray
 USE hterms
 USE rates
 USE timestep
 USE options
 USE part
 USE part_in
 USE setup_params
!
!--define local variables
!      
 IMPLICIT NONE
 CHARACTER(LEN=20) :: infile,evfile
 CHARACTER(LEN=20) :: datfile,logfile   ! rootname is global in loguns
 INTEGER :: i,j
!
!--add extensions to set filenames (remove trailing blanks)
!
 infile = TRIM(rootname)//'.in'
 evfile = TRIM(rootname)//'.ev'
 datfile = TRIM(rootname)//'.dat'
 logfile = TRIM(rootname)//'.log'      
!
!--Initialise log file (if used)
!
 IF (iprint.NE.6) THEN
    OPEN(UNIT=iprint,FILE=logfile,ERR=669,   &
         STATUS='replace',FORM='formatted')
 ELSE
    logfile = 'screen   '
 ENDIF
!
!--start tracing flow into logfile
!
 IF (trace) WRITE(iprint,*) ' Entering subroutine initialise'
!
!--set some global parameters based on ndim
!
 IF (ndim.GT.3 .OR. ndim.LT.0 .OR. ndimV.GT.3 .OR. ndimV.LT.0) THEN
    STOP 'Error ndim <0 or >3: We leave string theory for later'
 ELSE
    dndim = 1./REAL(ndim)
    hpower = dndim 
 ENDIF
 time = 0.
 nsteps = 0
 xmin = 0.
 xmax = 0.
 nbpts = 0
!
!--read parameters from the infile
!
 CALL set_default_options   ! set the default options
 CALL read_infile(infile)
!
!--Open data/ev files
!
 OPEN(UNIT=idatfile,ERR=667,FILE=datfile,STATUS='replace',FORM='unformatted')
 OPEN(UNIT=ievfile,ERR=668,FILE=evfile,STATUS='replace',FORM='formatted') 
!
!--work out multiplication factor for source term in Morris and Monaghan scheme
!  (just from gamma)

 IF (abs(gamma-1.).GT.1.e-3) THEN   ! adiabatic
    avfact = LOG(4.)/(LOG((gamma+1.)/(gamma-1.)))
 ELSE
    avfact = 1.0   ! isothermal 
 ENDIF
!
!--write first header to logfile/screen
!
 CALL write_header(1,infile,datfile,evfile,logfile)    
 
 CALL setkern      ! setup kernel tables
 npart = 0
 CALL setup          ! setup particles, allocation of memory is called
          ! could replace this with a call to read_dump
 CALL check_setup       ! check for errors in the particle setup
!
!--setup additional quantities that are not done in setup
!
 alpha(1,:) = alphamin
 alpha(2,:) = alphaumin
 alpha(3,:) = alphaBmin
 gradh = 0.
 divB = 0.
 curlB = 0.
 fmag = 0.
 psi = 0.
 sqrtg = 1.
 if (imhd.eq.0) Bfield = 0.  ! zero mag field if turned off
!
!--if using fixed particle boundaries, set them up
!
 IF (ANY(ibound.EQ.1)) CALL set_fixedbound
!
!--calculate the conservative quantities (rho, en, B/rho)
!
 call primitive2conservative
!
!--Set derivatives to zero until calculated
!      
 drhodt = 0.
 dudt = 0.
 dendt = 0.
 force = 0.
 dhdt = 0.
 daldt = 0.    
 dBevoldt = 0.
 dpsidt = 0.
 gradpsi = 0.
 DO i=1,npart
    xin(:,i) = x(:,i)
    velin(:,i) = vel(:,i)
    Bevolin(:,i) = Bevol(:,i)
    rhoin(i) = rho(i)
    hhin(i) = hh(i)
    enin(i) = en(i)
    alphain(:,i) = alpha(:,i)
    psiin(i) = psi(i)
 ENDDO         
!
!--make sure ghost particle quantities are the same
!  (copy both conservative and primitive here)
!
 IF (ANY(ibound.GT.1)) THEN
    DO i=npart+1,ntotal   
       j = ireal(i)
       call copy_particle(i,j)
    ENDDO
 ENDIF
!
!--write second header to logfile/screen
!
 CALL write_header(2,infile,datfile,evfile,logfile)    
      
 RETURN
!
!--Error control
!
666   WRITE(iprint,*) 'Initialise: Error reading filenames, exiting...'
      CALL quit
667   WRITE(iprint,*) 'Initialise: Error opening data file, exiting...'
      CALL quit 
668   WRITE(iprint,*) 'Initialise: Error opening ev file, exiting...'
      CALL quit 
669   WRITE(iprint,*) 'Initialise: Error opening log file, exiting...'
      CALL quit       
      
END SUBROUTINE initialise
