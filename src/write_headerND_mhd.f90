!!-----------------------------------------------------------------
!! This subroutine writes header to logfile / screen
!! with the main parameters used for the run
!!
!! icall = 1 (before particle have been set up)
!! icall = 2 (after particle setup)
!!-----------------------------------------------------------------

SUBROUTINE write_header(icall,infile,datfile,evfile,logfile)
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
 CHARACTER(LEN=*), INTENT(IN) :: infile,evfile,datfile,logfile
 CHARACTER(LEN=19) :: boundtype
 CHARACTER(LEN=1), DIMENSION(3) :: coord
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
       WRITE(ireadf,99001) TRIM(version)
    CLOSE(UNIT=ireadf)
99001 FORMAT(a)
!
!--Write parameters to screen/logfile
!     
    WRITE (iprint,99002) version
    IF (iprint.NE.6) WRITE (*,99002) version
     
99002 FORMAT(      							  &
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
      '								     ',/,  &
      /,' Version: ',a30,/)
    WRITE(iprint, 99003) infile,datfile,evfile,logfile
    IF (iprint.NE.6) WRITE(*, 99003) infile,datfile,evfile,logfile
99003 FORMAT(/,				&
      ' Read input from   : ',a,/,	&
      ' Writing output to : ',a,/,	&
      ' Writing energy to : ',a,/,	&
      ' Writing log to    : ',a,/)
!
!--time stepping options
!
    WRITE (iprint, 99004) tmax, nmax, tout, nout
99004 FORMAT(' Maximum time   = ',f7.3, ' or ',i7,' timesteps',/,	&
             ' Output every t = ',f7.3, ' or ',i7,' timesteps',/)
!
!--number of dimensions
!
    WRITE (iprint, 99104) ndim,ndimV
99104 FORMAT(' Number of spatial dimensions = ',i1,		&
             ', velocity varies in ',i1,' dimensions.',/)
!
!--general options on form of equations
!
    WRITE (iprint,99009) iener, icty, ialtform, iav, imhd
99009 FORMAT(' Options: ',/,						&
      6x,' Energy equation  : ',i2,5x,' Continuity Equation  :',i2,/,	&
      6x,' Pressure term    : ',i2,5x,' Artificial viscosity :',i2,/,	&
      6x,' Magnetic fields  : ',i2,/)

    WRITE (iprint,99010) ihvar, ikernav, hfact, ndim
99010 FORMAT(' Variable smoothing length: ',/,				&
      6x,' h varied using method : ',i2,4x,' Kernel averaging :',i2,/,  &
      6x,' h = ',f4.2,'*m/rho**(1/',i1,')',/)

    IF (imhd.NE.0) THEN
       WRITE (iprint,99011) imagforce,idivBzero,ianticlump,eps,neps
99011  FORMAT(' MHD options: ',/,					&
         6x,' Lorentz force    : ',i2,5x,' div B correction  :',i2,/,	&
         6x,' Artificial stress: ',i2,5x,' with eps = ',f4.2,', neps = ',i2,/)
 ENDIF
!
!--artificial viscosity options
!     
    IF (iav.NE.0) THEN
       WRITE (iprint,99012) alphamin,beta,iavlim,avconst
99012  FORMAT(' Artificial viscosity: ',/,			&
         6x,' alpha (min) = ',f10.6,' beta = ',f10.6,/,		&
         6x,' Morris & Monaghan AV limiter : ',i2,		&
         '   with constant = ',f10.6,/)
    ENDIF
!
!--general constants
!     
    WRITE (iprint, 99013) gamma
99013 FORMAT(' Equation of state: ',/,6x,' gamma = ',f10.6,/)

!
!--timestepping
!     
    WRITE (iprint, 99014) C_cour,C_force
99014 FORMAT(' Timestepping conditions : ',/, &
        6x,' C_courant = ',f4.2,5x,' C_force = ',f4.2,/)


    WRITE (iprint, 99015)
99015 FORMAT(/,' Particle setup: ',/)

!----------------------------------------------
! 2nd header after particles have been setup
!----------------------------------------------

 ELSEIF (icall.EQ.2) THEN
 
!
!--number of particles
!     
    WRITE (iprint, 99005) npart     
99005 FORMAT(' Number of particles = ',i6,/)

!
!--boundary options
!
    coord(1) = 'x'	! this should change with different co-ordinate systems
    coord(2) = 'y'
    coord(3) = 'z'
    DO i=1,ndim
       IF (ibound(i).EQ.0) boundtype = 'None'
       IF (ibound(i).EQ.1) boundtype = 'Fixed particles'
       IF (ibound(i).EQ.2) boundtype = 'Reflective ghosts'
       IF (ibound(i).EQ.3) boundtype = 'Periodic (ghosts)'     
       WRITE (iprint, 99006) coord(i),TRIM(boundtype)             
       WRITE (iprint, 99007) coord(i),xmin(i),coord(i),xmax(i)
    ENDDO
    WRITE (iprint, 99008) nbpts
99006 FORMAT(1x,a1,' boundary:',1x,a)
99007 FORMAT(13x,a1,'min = ',f10.6,1x,a1,'max = ',f10.6)
99008 FORMAT(/,1x,'Number of fixed particles = ',i4,/)      
!
!--print out smoothing length information
!
    CALL minmaxave(hh(1:npart),hmin,hmax,have,npart)
    WRITE(iprint,99099) hmin,hmax,have
99099 FORMAT (/,' Smoothing lengths: min = ',1pe8.2,' max = ',1pe8.2,' ave = ',1pe8.2,/)
      
!
!--write header for timestep table
!      
    WRITE (iprint, 99020)
99020 FORMAT (76('_'))       

 ENDIF
      
RETURN
END SUBROUTINE write_header
