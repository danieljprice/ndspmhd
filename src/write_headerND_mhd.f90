!!-----------------------------------------------------------------
!! This subroutine writes header to logfile / screen
!! with the main parameters used for the run
!!
!! icall = 1 (before particle have been set up)
!! icall = 2 (after particle setup)
!!-----------------------------------------------------------------

subroutine write_header(icall,infile,evfile,logfile)
 use dimen_mhd
 use debug
 use loguns
 
 use kernels, only:ianticlump,eps,neps,radkern
 use artvi
 use bound
 use hterms, only:rhomin
 use eos
 use options
 use setup_params
 use part
 use timestep
 use versn
 use xsph
 use utils, only:minmaxave
!
!--define local variables
!
 implicit none
 integer :: i
 integer, intent(in) :: icall
 character(len=*), intent(in) :: infile,evfile,logfile
 character(len=21) :: boundtype
 character(len=5), dimension(3) :: coord
 character(len=10) :: startdate, starttime
 real :: have,hmin,hmax,v2i,B2i,fnneigh
 real, dimension(npart) :: dummy
!
!--allow for tracing flow
!      
 if (trace) write(iprint,*) ' entering subroutine write_header'

!-----------------------------------------------------------------------
! 1st header after options have been read, but before particle setup
!-----------------------------------------------------------------------

 if (icall.eq.1) then
!
!--write version to file (for making backups of this code version)
!
    open(unit=ireadf,file='version',status='replace')
       write(ireadf,"(a)") trim(version)
    close(unit=ireadf)
!
!--write parameters to screen/logfile
!     
    write (iprint,10) version
    if (iprint.ne.6) write (*,10) version
     
10 format(                                                                &
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
      
    call date_and_time(startdate,starttime)
    startdate = startdate(7:8)//'/'//startdate(5:6)//'/'//startdate(1:4)
    starttime = starttime(1:2)//':'//starttime(3:4)//':'//starttime(5:)
    write(iprint,"(' Run started on ',a,' at ',a)") startdate,starttime
    
    write(iprint, 20) trim(infile),trim(evfile),trim(logfile)
    if (iprint.ne.6) write(*, 20) trim(infile),trim(evfile),trim(logfile)
20 format(/,                                &
      ' Read input from   : ',a,/,        &
      ' Writing energy to : ',a,/,        &
      ' Writing log to    : ',a,/)
!
!--time stepping options
!
    write (iprint, 30) tmax, nmax, tout, nout
30 format(' Maximum time   = ',f7.3, ' or ',i7,' timesteps',/,        &
             ' Output every t = ',f7.3, ' or ',i7,' timesteps',/)
!
!--number of dimensions
!
    write (iprint, 40) ndim,ndimV
40 format(' Number of spatial dimensions = ',i1,                &
             ', velocity varies in ',i1,' dimensions.',/)
!
!--general options on form of equations
!
    write (iprint,50) iener, icty, iprterm, iav, imhd, iexternal_force
50 format(' Options: ',/,                                                &
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

    if (imhd.ne.0) then
       write (iprint,70) imagforce,idivBzero,ianticlump,eps,neps
70  format(' MHD options: ',/,                                        &
         6x,' Lorentz force    : ',i2,5x,' div B correction  :',i2,/,        &
         6x,' Artificial stress: ',i2,5x,' with eps = ',f4.2,', neps = ',i2,/)
 endif
!
!--artificial viscosity options
!     
    if (iav.ne.0) then
       write (iprint,80) alphamin,alphaumin,alphaBmin,beta, &
                         iavlim(:),avdecayconst,avfact
80  format(' Artificial dissipative terms: ',/,                        &
         6x,' alpha (min) = ',f6.2,f6.2,f6.2,' beta = ',f6.2,/,                &
         6x,' viscosity limiter   : ',i2,/, &
         6x,' conduction limiter  : ',i2,/, &
         6x,' resistivity limiter : ',i2,/, &
         6x,' decay constant = ',f6.2,', av source term x',f10.6, /)
    endif
!
!--general constants
!     
    write (iprint, 90) gamma
90 format(' Equation of state: gamma = ',f10.6,/)

!
!--timestepping
!     
    write (iprint, 100) C_cour,C_force
100 format(' Timestepping conditions: C_cour = ',f4.2,', C_force = ',f4.2,/)

    write (iprint,"(' Particle setup: ',/)")

!----------------------------------------------
! 2nd header after particles have been setup
!----------------------------------------------

 elseif (icall.eq.2) then
 
!
!--number of particles
!     
    write (iprint,"(/,' Number of particles = ',i9,/)") npart
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
    
    do i=1,ndim
       boundtype = 'Unknown'
       if (ibound(i).eq.0) boundtype = 'None'
       if (ibound(i).eq.1) boundtype = 'Fixed particles'
       if (ibound(i).eq.2) boundtype = 'Reflective ghosts'
       if (ibound(i).eq.3) boundtype = 'Periodic (ghosts)'     
       if (ibound(i).eq.5) boundtype = 'Shearing box (ghosts)'
       write (iprint, 210) coord(i),trim(boundtype)             
       if (ibound(i).ne.0) write (iprint, 220) coord(i),xmin(i),coord(i),xmax(i)
    enddo
    write (iprint, 230) nbpts
210 format(1x,a1,' boundary:',1x,a)
220 format(13x,a1,'min = ',f10.6,1x,a1,'max = ',f10.6)
230 format(/,1x,'Number of fixed particles = ',i4,/)      
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
    
    write (iprint,240) ihvar, ikernav, hfact, rhomin, ndim, tolh, fNneigh
240 format(' Variable smoothing length: ',/,                                &
      6x,' h varied using method : ',i2,4x,' Kernel averaging :',i2,/,  &
      6x,' h = ',f4.2,'*[m/(rho + ',f7.5,')]^(1/',i1,'); htol = ',es8.2,/ &
      6x,' Number of neighbours = ',f8.2)
!
!--print out diagnostics of run
!
    write(iprint,*)
    call minmaxave(hh(1:npart),hmin,hmax,have,npart)
    write(iprint,250) 'h',hmin,hmax,have
    call minmaxave(dens(1:npart),hmin,hmax,have,npart)
    write(iprint,250) 'dens',hmin,hmax,have
    call minmaxave(uu(1:npart),hmin,hmax,have,npart)
    write(iprint,250) 'u',hmin,hmax,have
    do i=1,ndimv
       call minmaxave(vel(i,1:npart),hmin,hmax,have,npart)
       write(iprint,250) 'v'//trim(coord(i)),hmin,hmax,have
    enddo
!--mach number
    if (iprterm.ge.0) then
       do i=1,npart
          v2i = dot_product(vel(:,i),vel(:,i))
          dummy(i) = sqrt(v2i)/spsound(i)
       enddo
       call minmaxave(dummy(1:npart),hmin,hmax,have,npart)
       write(iprint,250) 'Mach #',hmin,hmax,have
    endif    
    if (imhd.ne.0) then
       do i=1,ndimv
          call minmaxave(bfield(i,1:npart),hmin,hmax,have,npart)
          write(iprint,250) 'B'//trim(coord(i)),hmin,hmax,have
       enddo
!--plasma beta
       do i=1,npart
          b2i = dot_product(bfield(:,i),bfield(:,i))
          dummy(i) = pr(i)/(0.5*b2i)
       enddo
       call minmaxave(dummy(1:npart),hmin,hmax,have,npart)
       write(iprint,250) 'beta',hmin,hmax,have
    endif

250 format (1x,a6,' min = ',1pe9.2,' max = ',1pe9.2,' av = ',1pe9.2)
    write(iprint,*)

 endif
      
return
end subroutine write_header
