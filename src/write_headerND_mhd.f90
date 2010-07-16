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
 
 USE kernels, only:ianticlump,eps,neps,radkern
 USE artvi
 USE bound
 USE hterms, only:rhomin
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
 CHARACTER(LEN=21) :: boundtype
 CHARACTER(LEN=5), DIMENSION(3) :: coord
 CHARACTER(LEN=10) :: startdate, starttime
 REAL :: have,hmin,hmax,v2i,B2i,fNneigh
 REAL, DIMENSION(npart) :: dummy
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
    
    WRITE(iprint, 20) trim(infile),trim(evfile),trim(logfile)
    IF (iprint.NE.6) WRITE(*, 20) trim(infile),trim(evfile),trim(logfile)
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
    WRITE (iprint,50) iener, icty, iprterm, iav, imhd, iexternal_force
50 FORMAT(' Options: ',/,                                                &
      6x,' Energy equation  : ',i2,5x,' Continuity Equation  :',i2,/,        &
      6x,' Pressure term    : ',i2,5x,' Artificial viscosity :',i2,/,        &
      6x,' Magnetic fields  : ',i2,5x,' External forces      :',i2/)
      
    select case(igravity)
    case(0)
       write(iprint,60) 'OFF'    
    case(1)
       write(iprint,65) 'ON, fixed plummer softening',hsoft    
    case(2)
       write(iprint,65) 'ON, fixed kernel softening ',hsoft
    case(3)
       write(iprint,60) 'ON, adaptive softening'    
    case(4)
       write(iprint,60) 'ON, adaptive softening with energy conservation'
    case default
       write(iprint,60) 'ON, unknown igravity setting'    
    end select
60  format(' Self gravity is ',a)
65  format(' Self gravity is ',a,': h =',1pe9.2,/)

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
90 FORMAT(' Equation of state: gamma = ',f10.6,/)

!
!--timestepping
!     
    WRITE (iprint, 100) C_cour,C_force
100 FORMAT(' Timestepping conditions: C_cour = ',f4.2,', C_force = ',f4.2,/)

    WRITE (iprint,"(' Particle setup: ',/)")

!----------------------------------------------
! 2nd header after particles have been setup
!----------------------------------------------

 ELSEIF (icall.EQ.2) THEN
 
!
!--number of particles
!     
    WRITE (iprint,"(/,' Number of particles = ',i9,/)") npart
!
!--boundary options
!
    select case(geom(1:6))
    case('cylrpz')
       coord(1) = 'r'
       coord(2) = 'phi'
       coord(3) = 'z'
    case('cylrzp')
       coord(1) = 'r'
       coord(2) = 'z'
       coord(3) = 'phi'
    case('sphrpt')
       coord(1) = 'r'
       coord(2) = 'phi'
       coord(3) = 'theta'    
    case('sphlog')
       coord(1) = 'log r'
       coord(2) = 'phi'
       coord(3) = 'theta'    
    case default
       coord(1) = 'x'
       coord(2) = 'y'
       coord(3) = 'z'
    end select
    
    DO i=1,ndim
       boundtype = 'unknown'
       IF (ibound(i).EQ.0) boundtype = 'None'
       IF (ibound(i).EQ.1) boundtype = 'Fixed particles'
       IF (ibound(i).EQ.2) boundtype = 'Reflective ghosts'
       IF (ibound(i).EQ.3) boundtype = 'Periodic (ghosts)'     
       IF (ibound(i).EQ.5) boundtype = 'Shearing box (ghosts)'
       WRITE (iprint, 210) coord(i),TRIM(boundtype)             
       IF (ibound(i).NE.0) WRITE (iprint, 220) coord(i),xmin(i),coord(i),xmax(i)
    ENDDO
    WRITE (iprint, 230) nbpts
210 FORMAT(1x,a1,' boundary:',1x,a)
220 FORMAT(13x,a1,'min = ',f10.6,1x,a1,'max = ',f10.6)
230 FORMAT(/,1x,'Number of fixed particles = ',i4,/)      
!
!--print out smoothing length information
!
    select case(ndim)
    case(1)
       fNneigh = 2.*radkern*hfact 
    case(2)
       fNneigh = pi*(radkern*hfact)**2
    case(3)
       fNneigh = 4./3.*pi*(radkern*hfact)**3
    end select
    
    WRITE (iprint,240) ihvar, ikernav, hfact, rhomin, ndim, tolh, fNneigh
240 FORMAT(' Variable smoothing length: ',/,                                &
      6x,' h varied using method : ',i2,4x,' Kernel averaging :',i2,/,  &
      6x,' h = ',f4.2,'*[m/(rho + ',f7.5,')]^(1/',i1,'); htol = ',es8.2,/ &
      6x,' Number of neighbours = ',f8.2)
!
!--print out diagnostics of run
!
    WRITE(iprint,*)
    CALL minmaxave(hh(1:npart),hmin,hmax,have,npart)
    WRITE(iprint,250) 'h',hmin,hmax,have
    CALL minmaxave(dens(1:npart),hmin,hmax,have,npart)
    WRITE(iprint,250) 'dens',hmin,hmax,have
    CALL minmaxave(uu(1:npart),hmin,hmax,have,npart)
    WRITE(iprint,250) 'u',hmin,hmax,have
    DO i=1,ndimV
       CALL minmaxave(vel(i,1:npart),hmin,hmax,have,npart)
       WRITE(iprint,250) 'v'//trim(coord(i)),hmin,hmax,have
    ENDDO
!--mach number
    IF (iprterm.GE.0) THEN
       DO i=1,npart
          v2i = dot_product(vel(:,i),vel(:,i))
          dummy(i) = sqrt(v2i)/spsound(i)
       ENDDO
       CALL minmaxave(dummy(1:npart),hmin,hmax,have,npart)
       WRITE(iprint,250) 'Mach #',hmin,hmax,have
    ENDIF    
    IF (imhd.NE.0) THEN
       DO i=1,ndimV
          CALL minmaxave(Bfield(i,1:npart),hmin,hmax,have,npart)
          WRITE(iprint,250) 'B'//trim(coord(i)),hmin,hmax,have
       ENDDO
!--plasma beta
       DO i=1,npart
          B2i = dot_product(Bfield(:,i),Bfield(:,i))
          dummy(i) = pr(i)/(0.5*B2i)
       ENDDO
       CALL minmaxave(dummy(1:npart),hmin,hmax,have,npart)
       WRITE(iprint,250) 'beta',hmin,hmax,have
    ENDIF

250 FORMAT (1x,a6,' min = ',1pe9.2,' max = ',1pe9.2,' av = ',1pe9.2)
    WRITE(iprint,*)

 ENDIF
      
RETURN
END SUBROUTINE write_header
