!!-----------------------------------------------------------------
!! This subroutine writes header to logfile / screen
!! with the main parameters used for the run
!!
!! icall = 1 (before particle have been set up)
!! icall = 2 (after particle setup)
!!-----------------------------------------------------------------

SUBROUTINE write_header(icall,infile,evfile,logfile)
 USE dimen_mhd
 USE debug
 USE loguns
 
 USE anticlumping
 USE artvi
 USE bound
 USE eos
 USE options
 USE setup_params
 USE part
 USE timestep
 USE versn
 USE xsph
!
!--define local variables
!
 IMPLICIT NONE
 INTEGER :: i
 INTEGER, INTENT(IN) :: icall
 CHARACTER(LEN=*), INTENT(IN) :: infile,evfile,logfile
 CHARACTER(LEN=19) :: boundtype
 CHARACTER(LEN=1), DIMENSION(3) :: coord
 CHARACTER(LEN=10) :: startdate, starttime
 REAL :: have,hmin,hmax
!
!--allow for tracing flow
!      
 IF (trace) WRITE(iprint,*) ' Entering subroutine write_header'

!-----------------------------------------------------------------------
! 1st header after options have been read, but before particle setup
!-----------------------------------------------------------------------

 IF (icall.EQ.1) THEN
!
!--write version to file (for making backups of this code version)
!
    OPEN(UNIT=ireadf,file='version',status='replace')
       WRITE(ireadf,"(a)") TRIM(version)
    CLOSE(UNIT=ireadf)
!
!--Write parameters to screen/logfile
!     
    WRITE (iprint,10) version
    IF (iprint.NE.6) WRITE (*,10) version
     
10 FORMAT(                                                                &
      "                                                  _         _ ",/,  &
      "    ___ _   _ _ __   ___ _ __ ___ _ __  _ __ ___ | |__   __| |",/,  &
      "   / __| | | | '_ \ / _ \ '__/ __| '_ \| '_ ` _ \| '_ \ / _` |",/,  &
      "   \__ \ |_| | |_) |  __/ |  \__ \ |_) | | | | | | | | | (_| |",/,  &
      "   |___/\__,_| .__/ \___|_|  |___/ .__/|_| |_| |_|_| |_|\__,_|",/,  &
      "             |_|                 |_|                          ",/,  &
      '      _   _     _   _   _   _   _   _     _   _   _   _   _   ',/,  &
      '     / \ / \   / \ / \ / \ / \ / \ / \   / \ / \ / \ / \ / \  ',/,  &
      '    ( B | y ) ( D | a | n | i | e | l ) ( P | r | i | c | e ) ',/,  &
      '     \_/ \_/   \_/ \_/ \_/ \_/ \_/ \_/   \_/ \_/ \_/ \_/ \_/  ',/,  &
      '                                                                     ',/,  &
      /,' Version: ',a30,/)
      
    CALL date_and_time(startdate,starttime)
    startdate = startdate(7:8)//'/'//startdate(5:6)//'/'//startdate(1:4)
    starttime = starttime(1:2)//':'//starttime(3:4)//':'//starttime(5:)
    WRITE(iprint,"(' Run started on ',a,' at ',a)") startdate,starttime
    
    WRITE(iprint, 20) infile,evfile,logfile
    IF (iprint.NE.6) WRITE(*, 20) infile,evfile,logfile
20 FORMAT(/,                                &
      ' Read input from   : ',a,/,        &
      ' Writing energy to : ',a,/,        &
      ' Writing log to    : ',a,/)
!
!--time stepping options
!
    WRITE (iprint, 30) tmax, nmax, tout, nout
30 FORMAT(' Maximum time   = ',f7.3, ' or ',i7,' timesteps',/,        &
             ' Output every t = ',f7.3, ' or ',i7,' timesteps',/)
!
!--number of dimensions
!
    WRITE (iprint, 40) ndim,ndimV
40 FORMAT(' Number of spatial dimensions = ',i1,                &
             ', velocity varies in ',i1,' dimensions.',/)
!
!--general options on form of equations
!
    WRITE (iprint,50) iener, icty, iprterm, iav, imhd
50 FORMAT(' Options: ',/,                                                &
      6x,' Energy equation  : ',i2,5x,' Continuity Equation  :',i2,/,        &
      6x,' Pressure term    : ',i2,5x,' Artificial viscosity :',i2,/,        &
      6x,' Magnetic fields  : ',i2,/)

    WRITE (iprint,60) ihvar, ikernav, hfact, ndim
60 FORMAT(' Variable smoothing length: ',/,                                &
      6x,' h varied using method : ',i2,4x,' Kernel averaging :',i2,/,  &
      6x,' h = ',f4.2,'*m/rho**(1/',i1,')',/)

    IF (imhd.NE.0) THEN
       WRITE (iprint,70) imagforce,idivBzero,ianticlump,eps,neps
70  FORMAT(' MHD options: ',/,                                        &
         6x,' Lorentz force    : ',i2,5x,' div B correction  :',i2,/,        &
         6x,' Artificial stress: ',i2,5x,' with eps = ',f4.2,', neps = ',i2,/)
 ENDIF
!
!--artificial viscosity options
!     
    IF (iav.NE.0) THEN
       WRITE (iprint,80) alphamin,alphaumin,alphaBmin,beta, &
                         iavlim(:),avdecayconst,avfact
80  FORMAT(' Artificial dissipative terms: ',/,                        &
         6x,' alpha (min) = ',f6.2,f6.2,f6.2,' beta = ',f6.2,/,                &
         6x,' viscosity limiter   : ',i2,/, &
         6x,' conduction limiter  : ',i2,/, &
         6x,' resistivity limiter : ',i2,/, &
         6x,' decay constant = ',f6.2,', av source term x',f10.6, /)
    ENDIF
!
!--general constants
!     
    WRITE (iprint, 90) gamma
90 FORMAT(' Equation of state: ',/,6x,' gamma = ',f10.6,/)

!
!--timestepping
!     
    WRITE (iprint, 100) C_cour,C_force
100 FORMAT(' Timestepping conditions : ',/, &
        6x,' C_courant = ',f4.2,5x,' C_force = ',f4.2,/)


    WRITE (iprint, 110)
110 FORMAT(/,' Particle setup: ',/)

!----------------------------------------------
! 2nd header after particles have been setup
!----------------------------------------------

 ELSEIF (icall.EQ.2) THEN
 
!
!--number of particles
!     
    WRITE (iprint, 200) npart     
200 FORMAT(' Number of particles = ',i6,/)

!
!--boundary options
!
    coord(1) = 'x'        ! this should change with different co-ordinate systems
    coord(2) = 'y'
    coord(3) = 'z'
    DO i=1,ndim
       IF (ibound(i).EQ.0) boundtype = 'None'
       IF (ibound(i).EQ.1) boundtype = 'Fixed particles'
       IF (ibound(i).EQ.2) boundtype = 'Reflective ghosts'
       IF (ibound(i).EQ.3) boundtype = 'Periodic (ghosts)'     
       WRITE (iprint, 210) coord(i),TRIM(boundtype)             
       WRITE (iprint, 220) coord(i),xmin(i),coord(i),xmax(i)
    ENDDO
    WRITE (iprint, 230) nbpts
210 FORMAT(1x,a1,' boundary:',1x,a)
220 FORMAT(13x,a1,'min = ',f10.6,1x,a1,'max = ',f10.6)
230 FORMAT(/,1x,'Number of fixed particles = ',i4,/)      
!
!--print out smoothing length information
!
    CALL minmaxave(hh(1:npart),hmin,hmax,have,npart)
    WRITE(iprint,240) hmin,hmax,have
240 FORMAT (/,' Smoothing lengths: min = ',1pe8.2,' max = ',1pe8.2,' ave = ',1pe8.2,/)

 ENDIF
      
RETURN
END SUBROUTINE write_header
