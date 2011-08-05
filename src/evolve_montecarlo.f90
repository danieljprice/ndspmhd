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
 use options, only:hsoft,igravity
 use setup_params, only:hfact,pi
 use densityprofiles, only:exact_densityprofiles
 use errors, only:calculate_errors
 use rates, only:force,poten
 use kernels, only:radkern
 use plummer_setup, only:npart1,npart2,idist
 use cons2prim, only:primitive2conservative
!
!--define local variables
!
 implicit none
 integer :: i,ierr,noutput,nevwrite
 real :: t_start,t_end,t_used,ran1
 real :: etot, momtot, etotin, momtotin, detot, dmomtot
 character(len=10) :: finishdate, finishtime
 integer, parameter :: nexact = 10001
 real, dimension(nexact) :: xexact, yexact
 real, dimension(npart) :: rr,fmag,residual
 real, dimension(2) :: msphere,rsoft
 real :: errL1,errL2,errLinf,toterrL2,toterrLinf,xmin,xmax,dxexact,t1,t2
 real :: errL2_1,errL2_2,toterrL2_1,toterrL2_2
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
 toterrL2 = 0.
 toterrLinf = 0.
 toterrL2_1 = 0.
 toterrL2_2 = 0.
 open(unit=423,file=trim(rootname)//'.err',status='replace',form='formatted')
 print*,'opening ',trim(rootname)//'.err',' for output'
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
    call cpu_time(t1)
    call primitive2conservative
    call cpu_time(t2)
    print*,'cpu time used = ',t2-t1
!
!--compute exact solution for comparison
!
    do i=1,npart
       rr(i) = sqrt(dot_product(x(:,i),x(:,i)))
       fmag(i) = sqrt(dot_product(force(:,i),force(:,i)))
    enddo
    xmin = 0.
    xmax = 200.
    if (maxval(rr(1:npart)).gt.xmax) stop 'error rrmax > xmax'
    dxexact = (xmax - xmin)/real(nexact-1)
    do i=1,nexact
       xexact(i) = xmin + (i-1)*dxexact
    enddo

    msphere(1) = 1.0
    msphere(2) = 0.0
    rsoft(1) = 1.0
    rsoft(2) = 0.1
    call exact_densityprofiles(3,idist,msphere,rsoft,xexact,yexact,ierr)
    call calculate_errors(xexact,yexact,rr,fmag,residual, &
         errL1,errL2,errLinf)
!    call force_error_densityprofiles(1,npart,x(1:ndim,1:npart),force(1:ndim,1:npart), &
!                                    poten(1:npart),errL2,errL1,ierr)
!    errLinf = 0.
    if (ierr /= 0) stop 'error in exact solution calculation'
    toterrL2 = toterrL2 + errL2
    toterrLinf = toterrLinf + errLinf
    print*,nsteps,' L2 error = ',errL2,' mean = ',toterrL2/real(nsteps), &
        ' Linf = ',errLinf,' mean = ',toterrLinf/real(nsteps)
    write(423,*) nsteps,errL2,toterrL2/real(nsteps),toterrLinf/real(nsteps)
!
!--print errors in each component
!
    if (npart2.gt.0) then
       call force_error_densityprofiles(1,npart1,x(1:ndim,1:npart1),force(1:ndim,1:npart1), &
                                       poten(1:npart1),errL2_1,errL1,ierr)
!--this is contribution of 1st component to total MASE, *NOT* MASE for 1st component
       errL2_1 = errL2_1*npart1/real(npart)
       toterrL2_1 = toterrL2_1 + errL2_2
       print*,nsteps,' L2 err (1) = ',errL2_1,' mean (1) = ',toterrL2_1/real(nsteps), &
           ' Linf (1) = ',errLinf
       call force_error_densityprofiles(1,npart2,x(1:ndim,npart1+1:npart2),force(1:ndim,npart1+1:npart2), &
                                       poten(npart1+1:npart2),errL2_2,errL1,ierr)
!--same here, but for second component
       errL2_2 = errL2_2*npart2/real(npart)
       toterrL2_2 = toterrL2_2 + errL2_2
       print*,nsteps,' L2 err (2) = ',errL2_2,' mean (2) = ',toterrL2_2/real(nsteps), &
           ' Linf (2) = ',errLinf
    endif
!
!--calculate errors
!
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
 close(unit=423)
!
!--write final error values to file
!
 open(unit=99,file=trim(rootname)//'.results',status='replace',action='write', &
     form='formatted')
     print*,' --- final results ---'
     if (igravity.ge.3) then
        write(99,*) hfact,4./3.*pi*(radkern*hfact)**3,toterrL2/real(nsteps),toterrLinf/real(nsteps),nsteps     
        print*,hfact,4./3.*pi*(radkern*hfact)**3,toterrL2/real(nsteps),toterrLinf/real(nsteps),nsteps     
     else
        write(99,*) hsoft,hsoft,toterrL2/real(nsteps),toterrLinf/real(nsteps),nsteps
        print*,hsoft,hsoft,toterrL2/real(nsteps),toterrLinf/real(nsteps),nsteps
     endif
 close(unit=99)

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
