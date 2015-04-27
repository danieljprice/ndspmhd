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

!!-----------------------------------------------------------------
!! This subroutine initialises everything needed for the run
!!-----------------------------------------------------------------

subroutine initialise
!
!--include relevant global variables
!
! use dimen_mhd
 use debug, only:trace
 use loguns, only:iprint,ievfile,ifile,rootname
 
 use artvi
 use bound
 use derivb
 use eos
 use fmagarray
 use kernels,only:setkern,setkernels,kernelname,kernelnamealt,radkern,&
                  setkerndrag,kernelnamedrag
 use hterms
 use rates
 use timestep
 use options
 use part
 use part_in
 use setup_params
 
 use infiles,   only:read_infile
 use dumpfiles, only:read_dump
 use cons2prim, only:primitive2conservative
 use dust,      only:init_drag
!
!--define local variables
!      
 implicit none
 character(len=len(rootname)+3) :: infile,evfile,dumpfile
 character(len=len(rootname)+6) :: logfile   ! rootname is global in loguns
 integer :: i,j,idash,idot,ierr,idotin,ierr1,ierr2
 logical :: iexist
 real :: radkernold
!
!--if filename is of the form crap_xxxxx.dat restart from a given dump
!
 idash = index(rootname,'_')
 idot = index(rootname,'.dat')
 if (idash.ne.0) then
    if (idot.eq.0) then
       idot = len_trim(rootname)+1
       dumpfile = trim(rootname)//'.dat'
       !!write(dumpfile,"(a,i5.5,'.dat')") trim(rootname),ifile
    else
       dumpfile = trim(rootname)
    endif
    read(rootname(idash+1:idot-1),*,iostat=ierr) ifile   
    if (ierr /=0) then
       print*,'***could not extract dump number ***'
       print*,'***starting new run using dumpfile***'
       ifile = 0
    endif
    rootname = rootname(1:idash-1)    
    print*,'ifile = ',ifile,'rootname = ',trim(rootname)
    !!ifile = ifile - 1
    !stop
 else
    !
    !--trim the .in if present
    !
    idotin = index(rootname,'.in')
    if (idotin.gt.0) rootname = rootname(1:idotin-1)
    ifile = -1
 endif
!
!--add extensions to set filenames (remove trailing blanks)
!
 infile = TRIM(rootname)//'.in'
 evfile = TRIM(rootname)//'.ev'
!! datfile = TRIM(rootname)//'.dat'
 logfile = TRIM(rootname)//'.log'      
!
!--Initialise log file (if used)
!
 if (iprint.ne.6) then
    if (ifile.le.0) then
       open(unit=iprint,file=logfile,err=669,status='replace',form='formatted')
    else
       inquire(file=logfile,exist=iexist)
       if (iexist) then
          i = 0
          do while(iexist)
             i = i + 1
             write(logfile,"(a,i2.2,'.log')") trim(rootname),i
             inquire(file=logfile,exist=iexist)
          enddo
       endif
       open(unit=iprint,file=logfile,err=669,status='new',form='formatted')
    endif
 else
    logfile = 'screen   '
 endif
!
!--start tracing flow into logfile
!
 if (trace) write(iprint,*) ' entering subroutine initialise'
!
!--set some global parameters based on ndim
!
 if (ndim.gt.3 .or. ndim.le.0 .or. ndimv.gt.3 .or. ndimv.le.0) then
    stop 'ERROR ndim <=0 or >3: We leave string theory for later'
 else
    dndim = 1./real(ndim)
 endif
 time = 0.
 nsteps = 0
 xmin = 0.
 xmax = 0.
 nbpts = 0
 bconst(:) = 0.
 vsig2max = 0.
!
!--read parameters from the infile
!
 call set_default_options   ! set the default options
 call read_infile(infile)
! geom = 'cylrpz' 
!
!--Open data/ev files
!
 if (ifile.gt.0) then
    inquire(file=evfile,exist=iexist)
    if (iexist) then
       i = 0
       do while(iexist)
          i = i + 1
          write(evfile,"(a,i2.2,'.ev')") trim(rootname),i
          inquire(file=evfile,exist=iexist)
       enddo
    endif
 endif
 open(unit=ievfile,err=668,file=evfile,status='replace',form='formatted')
!
!--work out multiplication factor for source term in morris and monaghan scheme
!  (just from gamma)

 if (abs(gamma-1.).gt.1.e-3) then   ! adiabatic
    avfact = log(4.)/(log((gamma+1.)/(gamma-1.)))
 else
    avfact = 1.0   ! isothermal 
 endif
!
!--write first header to logfile/screen
!
 call write_header(1,infile,evfile,logfile)    
!
!--setup kernel tables
!
 ikernelalt = ikernel
 call setkernels(ikernel,ikernelalt,ndim,ierr1,ierr2)
 write(iprint,"(/,' Smoothing kernel = ',a)") trim(kernelname)
 radkernold = radkern
!
!--setup drag kernel ONLY if idrag_nature > 0
!  also radkern MUST match between kernels
!
 ierr = 0
 if (idust /= 0) then
    if (ikernel.eq.0) then
       call setkerndrag(42,ndim,ierr)
    elseif (ikernel.eq.3) then
       call setkerndrag(43,ndim,ierr)
    elseif (ikernel.eq.99) then
        call setkerndrag(99,ndim,ierr)
      print*,'CAREFUL, HEXIC/HEPTIC DOUBLE HUMP NOT IMPLEMENTED YET'    
    elseif (ikernel.eq.100) then
       call setkerndrag(100,ndim,ierr)
       print*,'CAREFUL, HEXIC/HEPTIC DOUBLE HUMP NOT IMPLEMENTED YET'
    else
       call setkerndrag(42,ndim,ierr)
       if (idrag_nature.gt.0) stop 'cannot get matching drag kernel'
    endif
    if (radkern.ne.radkernold) stop 'drag kernel has different support radius'
    write(iprint,"(' Drag kernel = ',a)") trim(kernelnamedrag)
 endif
 if (ikernelalt.ne.ikernel) write(iprint,"(' Number density kernel = ',a)") trim(kernelnamealt)
 if (ierr1.ne.0 .or. ierr2.ne.0 .or. ierr.ne.0) stop 'error with kernel setup' 
 npart = 0
!
!--initialise drag terms
!
 if (idust /= 0) call init_drag(ierr,gamma)

 if (ifile.lt.0) then
    write(iprint,"(1x,80('-'))")
    call setup          ! setup particles, allocation of memory is called
    write(iprint,"(1x,80('-'))")
 else
    call read_dump(trim(dumpfile),time)      ! or read from dumpfile
 endif
!
!--change coordinate systems if necessary
!
 if (ifile.eq.0) call modify_dump
 print*,'geometry = ',geomsetup,geom
 if (geomsetup.ne.geom) stop 'unknown geometry'
 
 call check_setup       ! check for errors in the particle setup
!
!--setup additional quantities that are not done in setup
!
 alpha(1,:) = alphamin
 alpha(2,:) = alphaumin
 alpha(3,:) = alphaBmin
 divB = 0.
 curlB = 0.
 fmag = 0.
 if (iprterm.ne.11) psi = 0.
 sqrtg = 1.
 if (imhd.eq.0) then
    Bfield = 0.  ! zero mag field if turned off
    Bevol = 0.
 endif
!
!--set velocities to zero if damping is set
!
 if (damp.gt.tiny(damp)) then
    write(iprint,*) 'SETTING VELS TO ZERO FOR RELAXATION RUN'
    vel = 0.
 endif
!
!--set minimum density if using variable particle masses
!  and density iterations
!
 if (any(pmass(1:npart).ne.pmass(1)) .and. icty.eq.0 .and. ikernav.eq.3) then
    rhomin = 0.
    !!rhomin = minval(dens(1:npart))
    write(iprint,*) 'rhomin = ',rhomin
 else
    rhomin = 0.
    write(iprint,*) 'particle mass = ',pmass(1)
 endif
 h_min = 0.0
!
!--if using fixed particle boundaries, set them up
!
 if (any(ibound.eq.1)) call set_fixedbound
!
!--Set derivatives to zero until calculated
!      
 drhodt = 0.
 dudt = 0.
 dendt = 0.
 force = 0.
 dhdt = 0.
 daldt = 0.    
 dBevoldt = 0.
 dpsidt = 0.
 gradpsi = 0.
!
!--calculate the conservative quantities (rho, en, B/rho)
!  this also sets the smoothing length
!
 call primitive2conservative

! call check_neighbourlist

 do i=1,npart
    xin(:,i) = x(:,i)
    velin(:,i) = vel(:,i)
    Bevolin(:,i) = Bevol(:,i)
    rhoin(i) = rho(i)
    hhin(i) = hh(i)
    enin(i) = en(i)
    alphain(:,i) = alpha(:,i)
    psiin(i) = psi(i)
 enddo         
!
!--make sure ghost particle quantities are the same
!  (copy both conservative and primitive here)
!
 if (any(ibound.gt.1) .and. .not.any(ibound.eq.6)) then
    do i=npart+1,ntotal   
       j = ireal(i)
       call copy_particle(i,j)
    enddo
 endif
!
!--write second header to logfile/screen
!
 call write_header(2,infile,evfile,logfile)    
      
 return
!
!--error control
!
668   write(iprint,*) 'Initialise: Error opening ev file, exiting...'
      call quit 
669   write(iprint,*) 'Initialise: Error opening log file, exiting...'
      call quit       
      
end subroutine initialise
