!!---------------------------------------------------------------------!!
!! evolves the simulation through all timesteps                        !!
!! this subroutine contains the main timestepping loop and calls       !! 
!!  the output routines at the appropriate times                       !!
!!---------------------------------------------------------------------!!

subroutine evolve
 use dimen_mhd, only:ndim
 use debug, only:trace
 use loguns, only:iprint,ifile
 use timestep, only:nsteps,nmax,time,nout,iseedMC
 use part !!, only:x,uu,pmass,npart
!
!--define local variables
!
 implicit none
 integer :: noutput,nevwrite
 real :: t_start,t_end,t_used,ran1
 real :: etot, momtot, etotin, momtotin, detot, dmomtot
 character(len=10) :: finishdate, finishtime
!
!--allow for tracing flow
!      
 if (trace) write(iprint,*) ' Entering subroutine evolve (Monte Carlo)'
!
!--get total potential
!
 t_start = 0.
 t_end = 0.
 nsteps = 0
 time = 0.
 nevwrite = 1
 detot = 0.
 dmomtot = 0.
!
!--calculate initial values of conserved quantities
!
 call evwrite(0.0,etotin,momtotin)
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
!
!--get starting cpu time
!
 call cpu_time(t_start)
!
! --------------------- main loop ----------------------------------------
!
 timestepping: do while (nsteps.lt.nmax)

    nsteps = nsteps + 1

!
!--get iseed from ran1
!
    iseedMC = -abs(int(1e8*ran1(iseedMC)))
    call setup           !  evolve data for one timestep
!
!--calculate the conservative quantities (rho, en, B/rho)
!  this also sets the smoothing length and calculates forces
!
    call primitive2conservative
!
!--write to data file if time is right
! 
    if (((mod(nsteps,nout).eq.0).and.(nout.gt.0))   &
       .or.(nsteps.ge.nmax)  ) then

       call output(time,nsteps)
       noutput = noutput + 1
    endif
!    
!--calculate total energy etc and write to ev file    
!
    if (mod(nsteps,nevwrite).eq.0) then
       call evwrite(real(nsteps),etot,momtot)
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
