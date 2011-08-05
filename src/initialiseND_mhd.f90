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
 
 USE infiles, only:read_infile
 USE dumpfiles, only:read_dump
 USE convert, only:convert_setup
!
!--define local variables
!      
 IMPLICIT NONE
 CHARACTER(LEN=len(rootname)+3) :: infile,evfile,dumpfile
 CHARACTER(LEN=len(rootname)+6) :: logfile   ! rootname is global in loguns
 INTEGER :: i,j,idash,idot,ierr
 LOGICAL :: iexist
!
!--if filename is of the form crap_xxxxx.dat restart from a given dump
!
 idash = index(rootname,'_')
 idot = index(rootname,'.dat')
 if (idash.ne.0) then
    if (idot.eq.0) then
       idot = len_trim(rootname)+1
       dumpfile = trim(rootname)//'.dat'
       !!write(dumpfile,"(a,i5.5,'.dat')") trim(rootname),ifile
    else
       dumpfile = trim(rootname)
    endif
    read(rootname(idash+1:idot-1),*,iostat=ierr) ifile   
    if (ierr /=0) then
       print*,'***could not extract dump number ***'
       print*,'***starting new run using dumpfile***'
       ifile = 0
    endif
    rootname = rootname(1:idash-1)
    print*,'ifile = ',ifile,'rootname = ',trim(rootname)
    !!ifile = ifile - 1
    !stop
 else
    ifile = -1
 endif
!
!--add extensions to set filenames (remove trailing blanks)
!
 infile = TRIM(rootname)//'.in'
 evfile = TRIM(rootname)//'.ev'
!! datfile = TRIM(rootname)//'.dat'
 logfile = TRIM(rootname)//'.log'      
!
!--Initialise log file (if used)
!
 if (iprint.ne.6) then
    if (ifile.le.0) then
       open(unit=iprint,file=logfile,err=669,status='replace',form='formatted')
    else
       inquire(file=logfile,exist=iexist)
       if (iexist) then
          i = 0
          do while(iexist)
             i = i + 1
             write(logfile,"(a,i2.2,'.log')") trim(rootname),i
             inquire(file=logfile,exist=iexist)
          enddo
       endif
       open(unit=iprint,file=logfile,err=669,status='new',form='formatted')
    endif
 else
    logfile = 'screen   '
 endif
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
 ENDIF
 time = 0.
 nsteps = 0
 xmin = 0.
 xmax = 0.
 nbpts = 0
 Bconst(:) = 0.
 vsig2max = 0.
!
!--read parameters from the infile
!
 CALL set_default_options   ! set the default options
 CALL read_infile(infile)
 !!igeom = 2
!
!--Open data/ev files
!
!! OPEN(UNIT=idatfile,ERR=667,FILE=datfile,STATUS='replace',FORM='unformatted')
 if (ifile.gt.0) then
    inquire(file=evfile,exist=iexist)
    if (iexist) then
       i = 0
       do while(iexist)
          i = i + 1
          write(evfile,"(a,i2.2,'.ev')") trim(rootname),i
          inquire(file=evfile,exist=iexist)
       enddo
    endif
 endif
 open(unit=ievfile,err=668,file=evfile,status='replace',form='formatted')
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
 CALL write_header(1,infile,evfile,logfile)    
 
 CALL setkern      ! setup kernel tables
 npart = 0
 
 if (ifile.lt.0) then
    call setup          ! setup particles, allocation of memory is called
 else
    call read_dump(trim(dumpfile),time)      ! or read from dumpfile
 endif
!
!--change coordinate systems if necessary
!
 if (ifile.ge.0) call modify_dump
 print*,'igeom = ',igeomsetup,igeom
 if (igeomsetup.ne.igeom) call convert_setup(igeomsetup,igeom)
 
 call check_setup       ! check for errors in the particle setup
!
!--setup additional quantities that are not done in setup
!
 alpha(1,:) = alphamin
 alpha(2,:) = alphaumin
 alpha(3,:) = alphaBmin
 gradh = 1.
 divB = 0.
 curlB = 0.
 fmag = 0.
 psi = 0.
 sqrtg = 1.
 if (imhd.eq.0) Bfield = 0.  ! zero mag field if turned off
!
!--set minimum density if using variable particle masses
!  and density iterations
!
 if (any(pmass(1:npart).ne.pmass(1)) .and. icty.eq.0 .and. ikernav.eq.3) then
    rhomin = 0.05
 else
    rhomin = 0.
 endif

!
!--if using fixed particle boundaries, set them up
!
 IF (ANY(ibound.EQ.1)) CALL set_fixedbound
!
!--calculate the conservative quantities (rho, en, B/rho)
!  this also sets the smoothing length
!
 call primitive2conservative

! call check_neighbourlist
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
 CALL write_header(2,infile,evfile,logfile)    
      
 RETURN
!
!--Error control
!
668   WRITE(iprint,*) 'Initialise: Error opening ev file, exiting...'
      CALL quit 
669   WRITE(iprint,*) 'Initialise: Error opening log file, exiting...'
      CALL quit       
      
END SUBROUTINE initialise
