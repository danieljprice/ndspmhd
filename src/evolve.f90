!!---------------------------------------------------------------------!!
!! evolves the simulation through all timesteps			       !!
!! this subroutine contains the main timestepping loop and calls       !! 
!!  the output routines at the appropriate times		       !!
!!---------------------------------------------------------------------!!

SUBROUTINE evolve
 USE dimen_mhd
 USE debug
 USE loguns
 USE options
 USE timestep
!
!--define local variables
!
 IMPLICIT NONE
 REAL :: tprint
 INTEGER :: noutput,nevwrite
!
!--allow for tracing flow
!      
 IF (trace) WRITE(iprint,*) ' Entering subroutine evolve'
!
!--Set initial timestep
!
 dt = 0.
 tprint = 0.0
 time = 0.0
 nsteps = 0
 nevwrite = 1	! frequency of writing to .ev file (could be read as parameter)
!
!--write the initial conditions to the output and evolution files
!
 CALL output(time,nsteps)
 CALL evwrite(time)
 noutput = 1
 tprint = tout
! CALL quit
!
! --------------------- Main loop ----------------------------------------
!
 dostep: DO WHILE ((time.LT.tmax).AND.(nsteps.LT.nmax))

    hdt = 0.5*dt
    time = time + dt
    nsteps = nsteps + 1

    CALL step 	 		!  Evolve data for one timestep

    dt = min(dtforce,0.25*dtcourant) ! new timestep from force/Courant condition
!
!--write log every step in 2D/3D
!
    IF (ndim.GE.2) THEN
       WRITE(iprint,10) time,dtforce,dtcourant
10     FORMAT(' t = ',f9.4,' dtforce = ',1pe10.3,' dtcourant = ',1pe10.3)
    ENDIF
    
    IF (dt.LT.1e-8) THEN
       WRITE(iprint,*) 'main loop: timestep too small, dt = ',dt
       CALL quit
    ENDIF
!
!--Write to data file if time is right
! 
    IF (   (time.GE.tprint) 				&
       .OR.(time.GE.tmax)				&
       .OR.((MOD(nsteps,nout).EQ.0).AND.(nout.GT.0))	&
       .OR.(nsteps.GE.nmax)  ) THEN
 
!--if making movies and need ghosts to look right uncomment the line below, 
       IF (idumpghost.EQ.1 .AND. ANY(ibound.GE.2)) CALL set_ghost_particles
       CALL output(time,nsteps)
       noutput = noutput + 1
       tprint = noutput*tout
    ENDIF
!    
!--calculate total energy etc and write to ev file    
!
    IF (MOD(nsteps,nevwrite).EQ.0) CALL evwrite(time)

    IF (dt.GE.(tprint-time)) dt = tprint-time	! reach tprint exactly
	 	 
 ENDDO dostep

!------------------------------------------------------------------------
 RETURN 
END SUBROUTINE evolve
