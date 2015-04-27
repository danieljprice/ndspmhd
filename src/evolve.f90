!------------------------------------------------------------------------------!
! NDSPMHD: A Smoothed Particle (Magneto)Hydrodynamics code for (astrophysical) !
! fluid dynamics simulations in 1, 2 and 3 spatial dimensions.                 !
!                                                                              !
! (c) 2002-2014 Daniel Price                                                   !
!                                                                              !
! http://users.monash.edu.au/~dprice/ndspmhd                                   !
! daniel.price@monash.edu -or- dprice@cantab.net (forwards to current address) !
!                                                                              !
!  NDSPMHD comes with ABSOLUTELY NO WARRANTY.                                  !
!  This is free software; and you are welcome to redistribute                  !
!  it under the terms of the GNU General Public License                        !
!  (see LICENSE file for details) and the provision that                       !
!  this notice remains intact. If you modify this file, please                 !
!  note section 2a) of the GPLv2 states that:                                  !
!                                                                              !
!  a) You must cause the modified files to carry prominent notices             !
!     stating that you changed the files and the date of any change.           !
!                                                                              !
!  ChangeLog:                                                                  !
!------------------------------------------------------------------------------!

!!---------------------------------------------------------------------!!
!! evolves the simulation through all timesteps                        !!
!! this subroutine contains the main timestepping loop and calls       !! 
!!  the output routines at the appropriate times                       !!
!!---------------------------------------------------------------------!!

subroutine evolve
 use dimen_mhd, only:ndim
 use debug, only:trace
 use loguns, only:iprint,ifile
 use options, only:idumpghost,ibound
 use timestep
!
!--define local variables
!
 implicit none
 real :: tprint
 integer :: noutput,nevwrite,nsort
 real :: t_start,t_end,t_used,tzero
 real :: etot, momtot, etotin, momtotin, detot, dmomtot
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
 dt0 = min(C_cour*dtcourant,C_force*dtforce)
 if (dtfixed) dt = dt0
 dtscale = 1.0
 dtrho = huge(dtrho)
 tprint = 0.
 t_start = 0.
 t_end = 0.
 nsteps = 0
 nevwrite = 1   ! frequency of writing to .ev file (could be read as parameter)
 nsort = 1 ! frequency of sorting (again, should be parameter)
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
 timestepping: do while ((time.lt.tmax).and.(nsteps.lt.nmax).and.(time.ge.0.))

    time = time + dt
    nsteps = nsteps + 1
!    
!--sort particles at regular intervals  
!
!    if (nsort.gt.0 .and. mod(nsteps,nsort).eq.0) then
!       call sort_particles()
!    endif
    
    call step ()      !  evolve data for one timestep
!
!--write log every step in 2D/3D
!
    if (ndim.ge.1) then
       if (abs(dt-C_force*dtforce).lt.epsilon(0.)) then
          write(iprint,10) time,C_force*dtforce
       elseif (abs(dt-C_cour*dtcourant).lt.epsilon(0.)) then
          write(iprint,15) time,C_cour*dtcourant
       elseif (abs(dt-C_force*dtdrag).lt.epsilon(0.)) then
          write(iprint,16) time,C_force*dtdrag
       elseif (abs(dt-C_force*dtvisc).lt.epsilon(0.)) then
          write(iprint,17) time,C_force*dtvisc
       else
          write(iprint,18) time,dt
       endif
10     format(' t = ',f12.4,' dtforce = ',es10.3)
15     format(' t = ',f12.4,' dtcourant = ',es10.3)
16     format(' t = ',f12.4,' dtdrag = ',es10.3)
17     format(' t = ',f12.4,' dtvisc = ',es10.3)
18     format(' t = ',f12.4,' dt(unknown) = ',es10.3)
    endif
    
!    if (abs(dt).lt.1e-8) then
!       write(iprint,*) 'main loop: timestep too small, dt = ',dt
!       call quit
!    endif
!
!--write to data file if time is right
! 
    if (dt.gt.0.) then
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
    else
      if (   (time.le.tprint)             &
       .or.(time.le.0.)            &
       .or.((mod(nsteps,nout).eq.0).and.(nout.gt.0))   &
       .or.(nsteps.ge.nmax)  ) then
 
!--if making movies and need ghosts to look right uncomment the line below, 
       if (idumpghost.eq.1 .and. any(ibound.ge.2)) call set_ghost_particles
       call output(time,nsteps)
       noutput = noutput + 1
       tprint = tprint - tout
      endif    
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
    if (.not.dtfixed .and. dt.ge.(tprint-time)) dt = tprint-time   ! reach tprint exactly
        
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
