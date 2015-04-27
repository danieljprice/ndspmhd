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

module infiles
 implicit none
 character(len=*), parameter :: disstype(3) = (/'viscosity   ','conductivity','resistivity '/)
 character(len=*), parameter :: disschar(3) = (/' ','u','B'/)

contains

!!-----------------------------------------------------------------
!! writes an input file
!!-----------------------------------------------------------------

subroutine write_infile(infile)
 use dimen_mhd
 use debug
 use loguns
 
 use artvi
 use eos
 use options
 use setup_params
 use timestep
 use xsph
 use infile_utils, only:write_inopt
!
!--define local variables
!      
 implicit none
 character(len=*), intent(in) :: infile
 integer :: i
!
!--allow for tracing flow
!      
 if (trace) write(iprint,*) ' entering subroutine write_infile'
              
 open(unit=iread,err=999,file=infile,status='replace',form='formatted')
  write(iread,"(a)") '# input file for ndspmhd, options not present assume default values'
  write(iread,"(/,a)") '# options affecting setup'
  call write_inopt(psep,'psep','particle separation',iread)

  write(iread,"(/,a)") '# options affecting timestepping/output'
  call write_inopt(tmax,'tmax','maximum time',iread)
  call write_inopt(tout,'tout','time between outputs',iread)
  call write_inopt(nmax,'nmax','maximum number of timesteps',iread)
  call write_inopt(nout,'nout','number of timesteps between outputs (-1=ignore)',iread)
  call write_inopt(C_cour,'C_cour','Courant factor on timestep',iread)
  call write_inopt(C_force,'C_force','Factor in force timestep condition',iread)
  call write_inopt(idumpghost,'idumpghost','dump ghost particles? (0: no 1: yes)',iread)

  write(iread,"(/,a)") '# options affecting density calculation'
  call write_inopt(icty,'icty','type of cty equation (0:direct sum 1:time deriv)',iread)
  call write_inopt(ndirect,'ndirect','perform direct sum every n timesteps (if icty=1)',iread)
  call write_inopt(maxdensits,'maxdensits','maximum number of density iterations',iread)
  call write_inopt(ihvar,'ihvar','type of variable smoothing length prediction',iread)
  call write_inopt(hfact,'hfact','h in units of mean particle spacing',iread)
  call write_inopt(tolh,'tolh','tolerance on h-rho iterations',iread)
  call write_inopt(ikernel,'ikernel','kernel type (0: cubic spline, 3:quintic)',iread)
  call write_inopt(ikernav,'ikernav','type of kernel averaging (1:average h, 2:average grad wab 3:springel/hernquist)',iread)
  call write_inopt(usenumdens,'usenumdens','Use number density formulation of gradh',iread)

  write(iread,"(/,a)") '# options affecting equation of state, energy equation'
  call write_inopt(iener,'iener','type of energy equation',iread)
  call write_inopt(gamma,'gamma','adiabatic index',iread)
  call write_inopt(polyk,'polyk','polytropic constant if iener=0',iread)

  write(iread,"(/,a)") '# options affecting momentum equation and shock-capturing terms'
  call write_inopt(iexternal_force,'iexternal_force','external force (1: toy star, 2:1/r^2 )',iread)
  call write_inopt(iav,'iav','type of artificial viscosity',iread)
  call write_inopt(alphamin,'alphamin','minimum alpha (viscosity)',iread)
  call write_inopt(alphaumin,'alphaumin','minimum alphau (conductivity)',iread)
  call write_inopt(alphaBmin,'alphaBmin','minimum alphaB (resistivity)',iread)
  call write_inopt(beta,'beta','beta in artificial viscosity',iread)
  do i=1,3
     call write_inopt(iavlim(i),'iavlim'//trim(disschar(i)),'use '//trim(disstype(i))//' switch',iread)
  enddo
  call write_inopt(avdecayconst,'avdecayconst','decay constant in av switch (0.1-0.2)',iread)
  call write_inopt(damp,'damp','artificial damping (0.0 or few percent)',iread)
  call write_inopt(dampz,'dampz','artificial damping in z (0.0 or few percent)',iread)
  call write_inopt(dampr,'dampr','artificial damping in r (0.0 or few percent)',iread)
  call write_inopt(ixsph,'ixsph','use XSPH (0:off 1:on)',iread)
  call write_inopt(xsphfac,'xsphfac','factor for XSPH correction',iread)
  call write_inopt(iprterm,'iprterm','type of pressure term (0:normal 1:pa+pb/rhoa*rhob 2:hernquist/katz',iread)

  write(iread,"(/,a)") '# options affecting magnetic fields'
  call write_inopt(imhd,'imhd','MHD (0:no 1-10:B/rho >10:B <0:A)',iread)
  call write_inopt(imagforce,'imagforce','MHD force type(1:vector 2:tensor)',iread)
  call write_inopt(idivbzero,'idivbzero','divergence correction method (0:none 1:projection 2: hyperbolic/parabolic)',iread)
  call write_inopt(psidecayfact,'psidecayfact','decay factor in hyperbolic/parabolic cleaning',iread)
  call write_inopt(iresist,'iresist','resistivity (0:off 1:explicit 2:implicit)',iread)
  call write_inopt(etamhd,'etamhd','eta for resistivity',iread)
  call write_inopt(iambipolar,'iambipolar','ambipolar diffusion',iread)

  write(iread,"(/,a)") '# options affecting dust'
  call write_inopt(idust,'idust','dust (0:off 1:one-f 2:two-f 3:diff-onef-1st 4:diff-onef-2ndderivs)',iread)
  call write_inopt(idrag_nature,'idrag_nature','drag type (0=none 1=const K 2=const ts 3=Epstein)',iread)
  !call write_inopt(idrag_structure,'idrag_structure','drag structure (1=default)',iread)
  call write_inopt(Kdrag,'Kdrag','drag coeff (idrag=1) or ts (idrag=2) or grain size in cm (idrag=3)',iread)
  call write_inopt(use_sqrtdustfrac,'use_sqrtdustfrac','evolve s=sqrt(rho*eps) instead of eps?',iread)

  write(iread,"(/,a)") '# options affecting self-gravity'
  call write_inopt(igravity,'igravity','self-gravity',iread)
  call write_inopt(hsoft,'hsoft','fixed softening length',iread)  

  write(iread,"(/,a)") '# options affecting physical viscosity'
  call write_inopt(ivisc,'ivisc','use physical viscosity?',iread)
  call write_inopt(shearvisc,'shearvisc','shear viscosity (nu)',iread)
  call write_inopt(bulkvisc,'bulkvisc','bulk viscosity (zeta)',iread)

  close(unit=iread)

 write(iprint,"(/,a)") ' input file '//trim(infile)//' created successfully'
 return
      
999   stop 'error creating input file, exiting...'
      
end subroutine write_infile

!!-----------------------------------------------------------------
!! reads parameters for the run from the input file
!!-----------------------------------------------------------------

subroutine read_infile(infile)
 use dimen_mhd
 use debug
 use loguns
 
 use artvi
 use eos
 use options
 use setup_params
 use timestep
 use xsph
 use infile_utils
!
!--define local variables
!      
 implicit none
 character(len=*), intent(in) :: infile 
 character(len=len(infile)+3) :: infilebak
 character(len=2*len(infile)+12) :: command
 character(len=1) :: ians
 logical :: iexist
 integer :: i,ierr,nerr,ierrold
 type(inopts), allocatable :: db(:)
!
!--allow for tracing flow
!      
 if (trace) write(iprint,*) ' entering subroutine read_infile'
              
 inquire(file=infile,exist=iexist)
 if (.not.iexist) then
    print "(' input file ',a20,' not found')",infile
    if (adjustl(infile(1:1)).ne.' ') then 
       write(*,*) ' would you like to create one with default options?'
       read*,ians
       if (ians.eq.'y'.or.ians.eq.'y') call write_infile(infile)
    endif
    stop 'exiting...'
 endif
 
  nerr = 0
  !print "(a)",'reading setup options from '//trim(infile)
  call open_db_from_file(db,infile,iread,ierr)
  if (ierr /= 0) then
     write (*,"(a)",advance='no') ' trying old format...'
     call read_infile_old(infile,ierrold)
     if (ierrold /= 0) then
        write(*,*) 'ERROR reading old format'
     else
        write(*,*) 'OK'     
     endif
  else
     ! particle setup
     call read_inopt(psep,'psep',db,errcount=nerr)

     ! timestepping
     call read_inopt(tmax,'tmax',db,errcount=nerr)
     call read_inopt(tout,'tout',db,errcount=nerr)
     call read_inopt(nmax,'nmax',db,errcount=nerr)
     call read_inopt(nout,'nout',db,errcount=nerr)
     call read_inopt(C_cour,'C_cour',db,errcount=nerr)
     call read_inopt(C_force,'C_force',db,errcount=nerr)
     call read_inopt(idumpghost,'idumpghost',db,errcount=nerr)

     ! equation of state
     call read_inopt(iener,'iener',db,errcount=nerr)
     call read_inopt(gamma,'gamma',db,errcount=nerr)
     call read_inopt(polyk,'polyk',db,errcount=nerr)

     ! density calculation
     call read_inopt(ikernel,'ikernel',db,errcount=nerr)
     call read_inopt(icty,'icty',db,errcount=nerr)
     call read_inopt(ndirect,'ndirect',db,errcount=nerr)
     call read_inopt(maxdensits,'maxdensits',db,errcount=nerr)
     call read_inopt(ikernav,'ikernav',db,errcount=nerr)
     call read_inopt(ihvar,'ihvar',db,errcount=nerr)
     call read_inopt(hfact,'hfact',db,errcount=nerr)
     call read_inopt(tolh,'tolh',db,errcount=nerr)

     ! momentum equation, dissipation
     call read_inopt(iprterm,'iprterm',db,errcount=nerr)
     call read_inopt(iexternal_force,'iexternal_force',db,errcount=nerr)
     call read_inopt(iav,'iav',db,errcount=nerr)
     call read_inopt(alphamin,'alphamin',db,errcount=nerr)
     call read_inopt(alphaumin,'alphaumin',db,errcount=nerr)
     call read_inopt(alphabmin,'alphaBmin',db,errcount=nerr)
     call read_inopt(beta,'beta',db,errcount=nerr)
     do i=1,3
        call read_inopt(iavlim(i),'iavlim'//trim(disschar(i)),db,iread)
     enddo
     call read_inopt(avdecayconst,'avdecayconst',db,errcount=nerr)
     call read_inopt(ixsph,'ixsph',db,errcount=nerr)
     call read_inopt(xsphfac,'xsphfac',db,errcount=nerr)
     call read_inopt(damp,'damp',db,errcount=nerr)
     call read_inopt(dampz,'dampz',db,errcount=nerr)
     call read_inopt(dampr,'dampr',db,errcount=nerr)

     ! mhd options
     call read_inopt(imhd,'imhd',db,errcount=nerr)
     if (imhd /= 0) then
        call read_inopt(imagforce,'imagforce',db,errcount=nerr)
        call read_inopt(idivbzero,'idivbzero',db,errcount=nerr)
        call read_inopt(psidecayfact,'psidecayfact',db,errcount=nerr)
        call read_inopt(iresist,'iresist',db,errcount=nerr)
        call read_inopt(etamhd,'etamhd',db,errcount=nerr)
        call read_inopt(iambipolar,'iambipolar',db,errcount=nerr)
     endif

     ! self-gravity
     call read_inopt(igravity,'igravity',db,errcount=nerr)
     if (igravity /= 0) then
        call read_inopt(hsoft,'hsoft',db,errcount=nerr)
     endif
     call read_inopt(usenumdens,'usenumdens',db,errcount=nerr)

     ! physical viscosity
     call read_inopt(ivisc,'ivisc',db,errcount=nerr)
     call read_inopt(shearvisc,'shearvisc',db,errcount=nerr)
     call read_inopt(bulkvisc,'bulkvisc',db,errcount=nerr)

     ! dust options
     call read_inopt(idust,'idust',db,errcount=nerr)
     if (idust /= 0) then
        call read_inopt(idrag_nature,'idrag_nature',db,errcount=nerr)
        !call read_inopt(idrag_structure,'idrag_structure',db,errcount=nerr)
        call read_inopt(Kdrag,'Kdrag',db,errcount=nerr)
        call read_inopt(use_sqrtdustfrac,'use_sqrtdustfrac',db,errcount=nerr)
     endif
  endif
  call close_db(db)

  if (nerr > 0) then  ! just overwrite existing file if new format
     print "(a)",' missing options in input file, re-writing...'
     call write_infile(infile)
     stop
  elseif (ierr /= 0) then ! save old format
     iexist = .true.
     i = 0
     do while(iexist)
        if (i > 1) then
           write(infilebak,"(a,i1)") trim(infile)//'.old',i
        else
           infilebak = trim(infile)//'.old'
        endif
        inquire(file=infilebak,exist=iexist)
        i = i + 1
     enddo
     command = 'mv '//trim(infile)//' '//trim(infilebak)
     print "(a)",' MOVING '//trim(infile)//'->'//trim(infilebak)
     call system(command)
     if (ierr /= 0) print "(a)",' RE-WRITING '//trim(infile)//' in new format'
     call write_infile(infile)
     stop
  endif
!
!--check options for possible errors
!      
 if (psep.lt.1.e-5) write(iprint,100) 'psep < 1.e-5'
 if (tout.gt.tmax) write(iprint,100) 'no output tout > tmax'
 if (nout.gt.nmax) write(iprint,100) 'no output nout > nmax'
 if (nout.eq.0) stop 'error in input: nout = 0'
 if (gamma.lt.1.) write(iprint,100) 'gamma < 1.0 '
 if (abs(gamma-1.).lt.1.e-3 .and. iener.ne.0) stop 'must use iener = 0 for isothermal eos'
 if ((iener.gt.0).and.(alphaumin.lt.0.).or.(alphabmin.lt.0.)) then
    write(iprint,100) 'alphaumin or alphabmin < 0.'
 elseif ((iener.eq.0).and.(polyk.lt.0.)) then
    write(iprint,100) 'polyk < 0.'      
 endif
 if ((iav.ne.0).and.(alphamin.lt.0. .or. beta.lt.0.) ) then
    write(iprint,100) 'av alpha or beta < 0.'      
 endif
 if ((iavlim(1).gt.0).and.(alphamin.ge.1.)) then
    write(iprint,100) 'using av limiter, but alphamin set > 1.0'
 endif
 if (any(iavlim.gt.0).and.((avdecayconst.le.0.01).or.(avdecayconst.gt.0.5))) then
    write(iprint,100) 'av decay constant not in range 0.01-0.5'
 endif     
 if ((ikernav.le.0).or.(ikernav.gt.3)) then
    write(iprint,100) 'kernel averaging not set (ikernav)'
 endif
 if ((hfact.le.1.0).or.(hfact.gt.2.0)) then
    write(iprint,100) 'hfact too low/high (1.0 < hfact < 2.0)'
 endif
 if (psidecayfact.lt.0.0) then
    write(iprint,100) 'psidecayfact < 0.0'
 endif
 if (tolh.lt.1.e-12) then
    write(iprint,100) 'tolh really, really tiny (probably zero)!!'
    stop
 endif
 if (iresist.lt.0 .or. iresist.gt.3) then
    write(iprint,100) 'invalid choice of resistivity formulation'
    stop
 endif
 if (etamhd.lt.0.) then
    write(iprint,100) 'eta < 0 in resistivity'
    stop
 endif
 if (shearvisc < 0.) then
    write(iprint,100) 'invalid choice of shear viscosity parameter'
    stop
 endif
 if (bulkvisc < 0.) then
    write(iprint,100) 'invalid choice of bulk viscosity parameter'
    stop
 endif
 if (idust < 0 .or. idust > 4) then
    write(iprint,100) 'invalid choice of dust formulation'
    stop
 endif
 if (iambipolar < 0 .or. iambipolar > 1) then
    write(iprint,100) 'invalid choice of ambipolar diffusion formulation'
    stop
 endif

100   format(/' read_infile: warning: ',a)
 return
      
end subroutine read_infile


!!-----------------------------------------------------------------
!! reads parameters for the run from the input file
!! to allow backwards compatibility with old input files
!!-----------------------------------------------------------------

subroutine read_infile_old(infile,ierr)
 use dimen_mhd
 use debug
 use loguns
 
 use artvi
 use eos
 use options
 use setup_params
 use timestep
 use xsph
!
!--define local variables
!
 implicit none
 character(len=*), intent(in)  :: infile
 integer,          intent(out) :: ierr
 character(len=1) :: ians
!
!--allow for tracing flow
!
 if (trace) write(iprint,*) ' entering subroutine read_infile_old'

 ierr = 0
 open(unit=iread,err=999,file=infile,status='old',form='formatted')
  read(iread,*,err=50,end=50) psep
  read(iread,*,err=50,end=50) tmax,tout,nmax,nout
  read(iread,*,err=50,end=50) gamma
  read(iread,*,err=50,end=50) iener,polyk
  read(iread,*,err=50,end=50) icty,ndirect,maxdensits
  read(iread,*,err=50,end=50) iprterm
  read(iread,*,err=50,end=50) iav,alphamin,alphaumin,alphabmin,beta
  read(iread,*,err=50,end=50) iavlim(:),avdecayconst
  read(iread,*,err=50,end=50) ikernav
  read(iread,*,err=50,end=50) ihvar,hfact,tolh
  read(iread,*,err=50,end=50) idumpghost
  read(iread,*,err=50,end=50) imhd,imagforce
  read(iread,*,err=50,end=50) idivbzero,psidecayfact
  read(iread,*,err=50,end=50) iresist,etamhd
  read(iread,*,err=50,end=50) ixsph,xsphfac
  read(iread,*,err=50,end=50) igravity,hsoft
  if (ndim.eq.3) then
     read(iread,*,err=50,end=50) damp,dampr,dampz
  elseif (ndim.eq.2) then
     read(iread,*,err=50,end=50) damp,dampr
  else
     read(iread,*,err=50,end=50) damp
  endif
  read(iread,*,err=50,end=50) ikernel
  read(iread,*,err=50,end=50) iexternal_force
  read(iread,*,err=50,end=50) C_Cour, C_force
  read(iread,*,err=50,end=50) usenumdens
  read(iread,*,err=50,end=50) idust,idrag_nature,idrag_structure,Kdrag,ismooth
  read(iread,*,err=50,end=50) ivisc,shearvisc,bulkvisc
 close(unit=iread)

 goto 55
50 continue
   close(unit=iread)
   ierr = 1
   return

55 continue
 return

!
! this should never happen (we check file exists before calling read_infile_old)
!
999   write(iprint,1000) infile
1000  format (' input file ',a20,' not found')
 ierr = 2
 if (adjustl(infile(1:1)).ne.' ') then 
    write(*,*) ' would you like to create one with default options?'
    read*,ians
    if (ians.eq.'y'.or.ians.eq.'y') call write_infile(infile)
 endif
      
end subroutine read_infile_old

end module infiles
