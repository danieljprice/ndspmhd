!!---------------------------------------------------------------------!!
!! evolves the simulation through all timesteps                        !!
!! this subroutine contains the main timestepping loop and calls       !! 
!!  the output routines at the appropriate times                       !!
!!---------------------------------------------------------------------!!

subroutine evolve
 use dimen_mhd, only:ndim
 use debug, only:trace
 use loguns, only:iprint,ifile
 use options, only:idumpghost,ibound,iexternal_force
 use timestep
 use part, only:x,uu,pmass,npart
!
!--define local variables
!
 implicit none
 real :: tprint
 integer :: i,noutput,nevwrite
 real :: t_start,t_end,t_used,tzero
 real :: etot, momtot, etotin, momtotin, detot, dmomtot, epoti
 character(len=10) :: finishdate, finishtime
!
!--allow for tracing flow
!      
 if (trace) write(iprint,*) ' Entering subroutine evolve'
!
!--set initial timestep
!  (nb: time is set in initialise and should not be changed here
!   as can be non-zero from reading a dumpfile)
!
 dt = 0.
!
!--get total potential
!
 Omega0 = 0.
 do i=1,npart
    call external_potentials(iexternal_force,x(:,i),epoti,ndim)
    Omega0 = Omega0 + pmass(i)*(uu(i) + epoti)
 enddo
 dt0 = min(C_cour*dtcourant,C_force*dtforce)
 w0 = 1./Omega0
 dtscale = 1.0
 dtrho = huge(dtrho)
 tprint = 0.
 t_start = 0.
 t_end = 0.
 nsteps = 0
 nevwrite = 1   ! frequency of writing to .ev file (could be read as parameter)
 detot = 0.
 dmomtot = 0.
 tzero = time
!
!--calculate initial values of conserved quantities
!
 call evwrite(time,etotin,momtotin)
 write(iprint,5) etotin, momtotin
5 format(' Total energy = ',1pe10.2,2x,'Linear momentum = ',1pe10.2) 
!
!--write header for timestep table
!      
 write (iprint,6)
6 format (76('_'))    
!
!--write initial conditions to output file
!  
 if (ifile.eq.0) ifile = ifile - 1 ! for runs from init dump of different name
 if (ifile.lt.0) call output(time,nsteps)
 noutput = 1
 tprint = tzero + tout
! call quit
!
!--get starting cpu time
!
 call cpu_time(t_start)
!
! --------------------- main loop ----------------------------------------
!
 timestepping: do while ((time.lt.tmax).and.(nsteps.lt.nmax))

    time = time + dt
    nsteps = nsteps + 1

    call step           !  evolve data for one timestep
!
!--write log every step in 2D/3D
!
    if (ndim.ge.2) then
       if (C_force*dtforce.lt.C_cour*dtcourant) then
          write(iprint,10) time,C_force*dtforce     
       else
          write(iprint,15) time,C_cour*dtcourant
       endif
10     format(' t = ',f9.4,' dtforce = ',1pe10.3)
15     format(' t = ',f9.4,' dtcourant = ',1pe10.3)
    endif
    
    if (dt.lt.1e-8) then
       write(iprint,*) 'main loop: timestep too small, dt = ',dt
       call quit
    endif
!
!--write to data file if time is right
! 
    if (   (time.ge.tprint)             &
       .or.(time.ge.tmax)            &
       .or.((mod(nsteps,nout).eq.0).and.(nout.gt.0))   &
       .or.(nsteps.ge.nmax)  ) then
 
!--if making movies and need ghosts to look right uncomment the line below, 
       if (idumpghost.eq.1 .and. any(ibound.ge.2)) call set_ghost_particles
       call output(time,nsteps)
       noutput = noutput + 1
       tprint = tzero + noutput*tout
    endif
!    
!--calculate total energy etc and write to ev file    
!
    if (mod(nsteps,nevwrite).eq.0) then
       call evwrite(time,etot,momtot)
       detot = max(detot,abs(etot-etotin))
       dmomtot = max(dmomtot,abs(momtot-momtotin))
    endif
!
!--reach tprint exactly. must take this out for integrator to be symplectic
!
!    if (dt.ge.(tprint-time)) dt = tprint-time   ! reach tprint exactly
        
 enddo timestepping

!------------------------------------------------------------------------

!
!--get ending cpu time
!
 call cpu_time(t_end)
!
!--print out total energy and momentum conservation
!
 if (abs(momtotin).lt.tiny(momtotin))  momtotin = 1.0
 if (abs(etotin).lt.tiny(etotin)) etotin = 1.0
 write(iprint,20) detot/abs(etotin),dmomtot/abs(momtotin)
20 format(/,' Max energy error   : ',1pe10.2,3x,' Max momentum error : ',1pe10.2)
!
!--now print out total code timings:
!
 t_used = t_end - t_start
 if (nsteps.gt.0) write(iprint,30) t_used,t_used/nsteps
30  format(' Total cpu time : ',1pe10.4,'s',2x,' Average time per step : ',f10.4,'s',/)
!
!--print time and date of finishing
!
 call date_and_time(finishdate,finishtime)
 finishdate = finishdate(7:8)//'/'//finishdate(5:6)//'/'//finishdate(1:4)
 finishtime = finishtime(1:2)//':'//finishtime(3:4)//':'//finishtime(5:)
 write(iprint,"(' Run finished on ',a,' at ',a,/)") finishdate,finishtime

 return 
end subroutine evolve
