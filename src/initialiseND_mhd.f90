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
 USE eos
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
 CHARACTER(LEN=20) :: datfile,logfile	! rootname is global in loguns
 INTEGER :: ntot
 INTEGER :: i,j,iktemp
 REAL :: B2i,v2i
!
!--try reading rootname from the command line
!
 CALL getarg(1,rootname)
! 
!--If nothing on command exit and print usage
!
 IF (rootname(1:1).EQ.' ') THEN
    WRITE(6,*) 'Enter name of run:'
    READ*,rootname
!    STOP 'Usage: spmhd runname '
!    OPEN(UNIT=ireadf,ERR=666,FILE='runname',STATUS='old',FORM='formatted')
!    READ(ireadf,'(A)') rootname
!    CLOSE(UNIT=ireadf)
 ENDIF
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
    OPEN(UNIT=iprint,FILE=logfile,ERR=669,	&
         STATUS='replace',FORM='formatted')
 ELSE
    logfile = 'screen   '
 ENDIF
!
!--start tracing flow into logfile
!
 IF (trace) WRITE(iprint,*) ' Entering subroutine initialise'
!
!--prevent ridiculous compilations
!
 IF (ndim.GT.3 .OR. ndim.LT.0 .OR. ndimV.GT.3 .OR. ndimV.LT.0) THEN
    STOP 'Error ndim <0 or >3: We leave string theory for later'
 ENDIF   
!
!--Open data/ev files
!
 OPEN(UNIT=idatfile,ERR=667,FILE=datfile,STATUS='replace',FORM='formatted')
 OPEN(UNIT=ievfile,ERR=668,FILE=evfile,STATUS='replace',FORM='formatted') 
!
!--read parameters from the infile
!
 CALL read_infile(infile)
!
!--work out multiplication factor for source term in Morris and Monaghan scheme
!  (just from gamma)

 IF (abs(gamma-1.).GT.1.e-3) THEN	! adiabatic
    avfact = 2.0*LOG(4.)/(LOG((gamma+1.)/(gamma-1.)))
 ELSE
    avfact = 1.0	! isothermal 
 ENDIF
 itoystar = 0	! unless toy star is setup
!
!--write first header to logfile/screen
!
 CALL write_header(1,infile,datfile,evfile,logfile)    
 
 ikernel = 0		! set type of kernel (this could be read as input options)
 CALL setkern		! setup kernel tables
 CALL setup    		! setup particles, allocation of memory is called
!
!--set boundaries to lie 1/2 a particle spacing from the end
!   (this should only be done if the grid is uniform near the boundary!!)
! IF (ibound.GT.1) THEN
!    CALL adjust_boundaries(xmin,xmax,xin(1:npart),npart)
! ENDIF
!
!--set particle properties equal to their initial values
!    
 DO i=1,npart
    x(:,i) = xin(:,i)
    rho(i) = rhoin(i)
    uu(i) = uuin(i)
    en(i) = enin(i)
    vel(:,i) = velin(:,i)
!
!--set smoothing length (h) according to the initial value of the density
!    
    hh(i) = hfact*(pmass(i)/rhoin(i))**hpower
!    hh(i) = hhin(i)	! uncomment if h set in setup
    hhin(i) = hh(i)
    IF (hh(i).LE.0) WRITE(iprint,*) 'Error in setup: h <= 0 :',i,hh(i)
    IF (imhd.NE.0) THEN
       Bfield(:,i) = Bin(:,i)
    ELSE
       Bfield(:,i) = 0.
    ENDIF
    IF (iavlim.NE.1) THEN
       alphain(i) = alphamin
    ELSE
       alphain(i) = alphamin       	! remove if set in setup
    ENDIF
    alpha(i) = alphain(i)
    gradh(i) = 0.
 ENDDO
!
!--set ghost particles initially
!
 IF (ibound.GT.1) CALL set_ghost_particles
!
!--If using direct sum for density, compute initial density from direct sum
!  also do this if using cty equation but doing a direct sum every so often
!      
 IF (icty.EQ.0 .OR. ndirect.LT.100000) THEN
    CALL link
    iktemp = ikernav
    ikernav = 3		! consistent with h for first density evaluation
    WRITE(iprint,*) 'Calculating initial density...'
    CALL iterate_density	! evaluate density by direct summation
    ikernav = iktemp
 ENDIF
!! CALL smooth_initial_conditions(uuin(:),uu(:),SIZE(uuin))

!
!--set initial pressure, vsound from eos (also u if polytropic)
!      
 DO i=1,npart	! not ghosts
    IF (i.LE.nbpts) THEN	! set fixed particles to value of density,h
       rho(i) = rho(nbpts+1)	! calculated just outside fixed zone
       hh(i) = hh(nbpts+1)
       gradh(i) = gradh(nbpts+1)
    ELSEIF (i.GT.(npart-nbpts)) THEN
       rho(i) = rho(npart-nbpts)
       hh(i) = hh(npart-nbpts)
       gradh(i) = gradh(npart-nbpts)
    ENDIF 
    
    rhoin(i) = rho(i)		! set initial rho and h equal to values calculated  
    hhin(i) = hh(i)
    CALL equation_of_state(pr(i),spsound(i),uu(i),rho(i),gamma,1)
 ENDDO
!
!--calculate the total energy per SPH particle, and
!  if using B/rho as the magnetic field variable, divide by rho
!      
 DO i=1,npart
    IF (imhd.GE.11) THEN		! B 
       B2i = DOT_PRODUCT(Bfield(:,i),Bfield(:,i))/rho(i)
    ELSEIF (imhd.NE.0) THEN		! B/rho
       Bfield(:,i) = Bfield(:,i)/rho(i)
       B2i = DOT_PRODUCT(Bfield(:,i),Bfield(:,i))*rho(i)
    ENDIF
    v2i = DOT_PRODUCT(vel(:,i),vel(:,i))
    en(i) = 0.5*v2i + uu(i) + 0.5*B2i
 ENDDO	! (note energy and B/rho are not set for ghosts yet)
!
!--make sure ghost particle quantities are the same
!
 IF (ibound.GT.1) THEN
    DO i=npart+1,ntotal	
       j = ireal(i)
       rho(i) = rho(j)
       hh(i) = hh(j)
       gradh(i) = gradh(j)
       pr(i) = pr(j)
       en(i) = en(j)
       uu(i) = uu(j)
       Bfield(:,i) = Bfield(:,j)
    ENDDO
 ENDIF                      
!
!--Set derivatives to zero until calculated
!      
 DO i=1,ntotal
    drhodt(i) = 0.
    dudt(i) = 0.
    dendt(i) = 0.
    force(:,i) = 0.
    daldt(i) = 0.
    dhdt(i) = 0.
    dBfielddt(:,i) = 0.
    Bfieldin(:,i) = Bfield(:,i)	! if using midpoint
    enin(i) = en(i)		! " " and not set in setup
 ENDDO
!
!--write second header to logfile/screen
!
 CALL write_header(2,infile,datfile,evfile,logfile)    
!
!--Write header line to data file
!
! WRITE(idatfile,*) npart,gamma,hfact 
      
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
