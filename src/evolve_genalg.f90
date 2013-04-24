!!---------------------------------------------------------------------!!
!! evolves the simulation through all timesteps                        !!
!! this subroutine contains the main timestepping loop and calls       !! 
!!  the output routines at the appropriate times                       !!
!!---------------------------------------------------------------------!!

subroutine evolve
 use dimen_mhd, only:ndim
 use debug, only:trace
 use loguns, only:iprint,ifile,rootname
 use timestep, only:nsteps,nmax,time,nout,iseedMC
 use part !!, only:x,uu,pmass,npart
 use errors, only:calculate_errors
 use kernels
 use cons2prim, only:conservative2primitive
 use options,   only:ikernel,ikernelalt
!
!--define local variables
!
 implicit none
 integer :: i,ierr,noutput,nevwrite,ierr1,ierr2
 real :: t_start,t_end,t_used,ran1
 character(len=10) :: finishdate, finishtime
 integer, parameter :: nexact = 10001
 real :: err,t1,t2
!
!--allow for tracing flow
!      
 if (trace) write(iprint,*) ' Entering subroutine evolve (Genetic algorithm)'
!
!--get total potential
!
 t_start = 0.
 t_end = 0.
 nsteps = 0
 time = 0.
 nevwrite = 1
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
 open(unit=423,file='kernels.error',status='replace',form='formatted')
 open(unit=424,file='kernels.rank',status='replace',form='formatted')
 timestepping: do while (nsteps.lt.nmax)

    nsteps = nsteps + 1
!
!--get iseed from ran1
!
    iseedMC = -abs(int(1e8*ran1(iseedMC)))
    !call setup           !  evolve data for one timestep

    ikernel = nsteps
    call setkernels(ikernel,ikernel,ndim,ierr1,ierr2)
    write(iprint,"(/,' Smoothing kernel ',i3,' = ',a)") ikernel,trim(kernelname)
    select case(ikernel)
    case(8,19,35,53,67,69)
       print*,' SKIPPING...'
       write(423,*) ikernel,0.
       cycle timestepping
    end select
    
    if (ierr1.ne.0 .or. ierr2.ne.0) then
       print*,'ERROR: skipping'
       write(423,*) ikernel,0.
       cycle timestepping
    endif
!
!--calculate the conservative quantities (rho, en, B/rho)
!  this also sets the smoothing length and calculates forces
!
    call cpu_time(t1)
    call derivs
    call cpu_time(t2)
    print*,'cpu time used = ',t2-t1
!
!--calculate errors
!
    err = 0.
    do i=1,npart
       err = err + (del2v(i) - 1.)**2
       !print*,i,del2v(i)
    enddo
    err = sqrt(err/real(npart))
    print*,' kernel = ',ikernel,' error = ',err,' R = ',radkern
    write(423,*) ikernel,err
    write(424,"(f12.8,1x,i2,1x,a,f6.3)") err,ikernel,trim(kernelname),radkern
    call flush(423)
    call flush(424)
    !read*
!
!--write to data file if time is right
! 
    if (((mod(nsteps,nout).eq.0).and.(nout.gt.0))   &
       .or.(nsteps.ge.nmax)  ) then

       call output(time,nsteps)
       noutput = noutput + 1
    endif
        
 enddo timestepping
 close(423)
 close(424)
!------------------------------------------------------------------------

!
!--get ending cpu time
!
 call cpu_time(t_end)
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
