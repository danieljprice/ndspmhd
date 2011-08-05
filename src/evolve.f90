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
 REAL :: t_start,t_end,t_used
 REAL :: etot, momtot, etotin, momtotin, detot, dmomtot
 CHARACTER(LEN=10) :: finishdate, finishtime
!
!--allow for tracing flow
!      
 IF (trace) WRITE(iprint,*) ' Entering subroutine evolve'
!
!--Set initial timestep
!
 dt = 0.
 dt0 = 0.
 tprint = 0.
 t_start = 0.
 t_end = 0.
 time = 0.
 nsteps = 0
 nevwrite = 1	! frequency of writing to .ev file (could be read as parameter)
 detot = 0.
 dmomtot = 0.
!
!--calculate initial values of conserved quantities
!
 CALL evwrite(time,etotin,momtotin)
 write(iprint,5) etotin, momtotin
5 format(' Total energy = ',1pe8.2,2x,'Linear momentum = ',1pe8.2) 
!
!--write header for timestep table
!      
 write (iprint,6)
6 format (76('_'))    
!
!--write initial conditions to output file
!  
 CALL output(time,nsteps)

 noutput = 1
 tprint = tout
! CALL quit
!
!--get starting CPU time
!
 CALL cpu_time(t_start)
!
! --------------------- Main loop ----------------------------------------
!
 timestepping: DO WHILE ((time.LT.tmax).AND.(nsteps.LT.nmax))

    time = time + dt
    nsteps = nsteps + 1

    CALL step 	 		!  Evolve data for one timestep
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
    IF (MOD(nsteps,nevwrite).EQ.0) THEN
       CALL evwrite(time,etot,momtot)
       detot = max(detot,abs(etot-etotin))
       dmomtot = max(dmomtot,abs(momtot-momtotin))
    ENDIF
!
!--reach tprint exactly. Must take this out for integrator to be symplectic
!
!    IF (dt.GE.(tprint-time)) dt = tprint-time	! reach tprint exactly
	 	 
 ENDDO timestepping

!------------------------------------------------------------------------

 write(iprint,6)
!
!--get ending CPU time
!
 CALL cpu_time(t_end)
!
!--print out total energy and momentum conservation
!
 if (abs(momtotin).lt.1.e-6)  momtotin = 1.0
 WRITE(iprint,20) detot/etotin,dmomtot/momtotin
20 FORMAT(/,' Max energy error   : ',1pe8.2,3x,' Max momentum error : ',1pe8.2)
!
!--Now print out total code timings:
!
 t_used = t_end - t_start
 IF (nsteps.gt.0) WRITE(iprint,30) t_used,t_used/nsteps
30  FORMAT(' Total CPU time : ',f8.4,'s',2x,' Average time per step : ',f8.4,'s',/)
!
!--print time and date of finishing
!
 CALL date_and_time(finishdate,finishtime)
 finishdate = finishdate(7:8)//'/'//finishdate(5:6)//'/'//finishdate(1:4)
 finishtime = finishtime(1:2)//':'//finishtime(3:4)//':'//finishtime(5:)
 WRITE(iprint,"(' Run finished on ',a,' at ',a,/)") finishdate,finishtime

 RETURN 
END SUBROUTINE evolve
