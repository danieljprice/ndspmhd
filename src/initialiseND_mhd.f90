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
 CHARACTER(LEN=20) :: datfile,logfile	! rootname is global in loguns
 INTEGER :: ntot,ierr
 INTEGER :: i,j,iktemp,npart1
 REAL :: hmin,hmax,have
 LOGICAL :: set_fixed_particles
!
!--try reading rootname from the command line
!
 CALL getarg(1,rootname)
! 
!--If nothing on command exit and print usage
!
 IF (rootname(1:1).EQ.' ') THEN
10  WRITE(6,*) 'Enter name of run:'
    READ(*,*,ERR=10) rootname
!    STOP 'Usage: spmhd runname '
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
!--set some global parameters based on ndim
!
 dndim = 1./REAL(ndim)
 hpower = dndim
!
!--read parameters from the infile
!
 CALL set_default_options   ! set the default options
 CALL read_infile(infile)
!
!--Open data/ev files
!
 OPEN(UNIT=idatfile,ERR=667,FILE=datfile,STATUS='replace',FORM='formatted')
 OPEN(UNIT=ievfile,ERR=668,FILE=evfile,STATUS='replace',FORM='formatted') 
!
!--work out multiplication factor for source term in Morris and Monaghan scheme
!  (just from gamma)

 IF (abs(gamma-1.).GT.1.e-3) THEN	! adiabatic
    avfact = 2.0*LOG(4.)/(LOG((gamma+1.)/(gamma-1.)))
 ELSE
    avfact = 1.0	! isothermal 
 ENDIF
!
!--write first header to logfile/screen
!
 CALL write_header(1,infile,datfile,evfile,logfile)    
 
 CALL setkern		! setup kernel tables
 npart = 0
 CALL setup    		! setup particles, allocation of memory is called
 			! could replace this with a call to read_dump
!
!--calculate conservative quantities (needs to be here for rho -> h)
! 
 CALL primitive2conservative
!
!--set particle properties that are not set in setup 
!  (initial smoothing length, dissipation parameter, gradh terms)
!
 hmin = 1.e6
 hmax = 0.
 have = 0.
 DO i=1,npart
    hh(i) = hfact*(pmass(i)/rho(i))**hpower   ! set using density
    hmin = min(hmin,hh(i))
    hmax = max(hmax,hh(i))
    have = have + hh(i)
    hhin(i) = hh(i)
    IF (hh(i).LE.0) WRITE(iprint,*) 'Error in setup: h <= 0 :',i,hh(i)
    alpha(i) = alphamin
    gradh(i) = 0.
    sqrtg(i) = 1.
 ENDDO
 have = have/REAL(npart)
 WRITE(iprint,20) hmin,hmax,have
20 FORMAT (/,' Smoothing lengths: min = ',1pe8.2,' max = ',1pe8.2,' ave = ',1pe8.2,/)
 divB = 0.
 curlB = 0.
 fmag = 0.
!
!--work out whether or not to set up fixed boundary particles
!
 set_fixed_particles = (ANY(ibound.EQ.1) .AND. ALL(itype.EQ.0)  &
                        .AND..NOT.(ndim.EQ.1 .AND. nbpts.GT.0) )
 IF (set_fixed_particles) WRITE(iprint,*) 'setting up fixed particles'
!
!--set ghost particles initially for the density calculation
!  (also do this if using fixed particle boundaries 
!   and no fixed particles have been set)
!
 IF (ANY(ibound.GT.1) .OR. set_fixed_particles) CALL set_ghost_particles
!
!  in 1D setup can just specify the number of particles to hold fixed
!  this is for backwards compatibility of the code
!
 IF (ALL(ibound.EQ.1) .AND. ndim.EQ.1 .AND. nbpts.GT.0) THEN
    itype(1:nbpts) = 1
    itype((npart-nbpts+1):ntotal) = 1
 ENDIF
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
    !! CALL smooth_initial_conditions(uuin(:),uu(:),SIZE(uuin))

    DO i=1,npart	! not ghosts
    !
    !--for 1D fixed particle boundaries, set values of rho, h for the boundary
    !  particles
    !
       IF (itype(i).EQ.1 .AND. ndim.EQ.1) THEN	! only for 1D fixed particles
          IF (i.LE.nbpts) THEN	! set fixed particles to value of density,h
             rho(i) = rho(nbpts+1)	! calculated just outside fixed zone
             hh(i) = hh(nbpts+1)
             gradh(i) = gradh(nbpts+1)
          ELSEIF (i.GT.(npart-nbpts)) THEN
             rho(i) = rho(npart-nbpts)
             hh(i) = hh(npart-nbpts)
             gradh(i) = gradh(npart-nbpts)
          ENDIF
       ENDIF
       rhoin(i) = rho(i)		! set initial rho and h equal to values calculated  
       hhin(i) = hh(i)
       !
       !--set initial pressure, vsound from eos (also u if polytropic)
       !      
       CALL equation_of_state(pr(i),spsound(i),uu(i),rho(i),gamma)
       
    ENDDO	! (note energy and B/rho are not set for ghosts yet)
    !
    !--recalculate the conserved variables since rho has changed
    !
    CALL primitive2conservative
 ENDIF
!
!--make sure ghost particle quantities are the same
!  (copy both conservative and primitive here)
!
 IF (ANY(ibound.GT.1) .OR. set_fixed_particles) THEN
    DO i=npart+1,ntotal	
       j = ireal(i)
       rho(i) = rho(j)
       hh(i) = hh(j)
       gradh(i) = gradh(j)
       pr(i) = pr(j)
       en(i) = en(j)
       uu(i) = uu(j)
       Bfield(:,i) = Bfield(:,j)
       Bcons(:,i) = Bcons(:,j)
    ENDDO
    IF (set_fixed_particles) THEN
       ! fix the ghost particles that have been set
       npart1 = npart + 1
       itype(npart1:ntotal) = 1	! set all these particles to be fixed
       npart = ntotal		! no ghosts
    ENDIF
 ENDIF                      
!
!--Set derivatives to zero until calculated
!      
 drhodt = 0.
 dudt = 0.
 dendt = 0.
 force = 0.
 dhdt = 0.
 daldt = 0.    
 dBconsdt = 0.
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
