module infiles
 implicit none

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
 use kernels, only:ianticlump,eps,neps
!
!--define local variables
!      
 implicit none
 character(len=*), intent(in) :: infile
!
!--allow for tracing flow
!      
 if (trace) write(iprint,*) ' entering subroutine write_infile'
              
 open(unit=iread,err=999,file=infile,status='replace',form='formatted')
  write(iread,10) psep
  write(iread,20) tmax,tout,nmax,nout
  write(iread,30) gamma
  write(iread,40) iener,polyk
  write(iread,50) icty,ndirect,maxdensits
  write(iread,60) iprterm
  write(iread,70) iav,alphamin,alphaumin,alphabmin,beta
  write(iread,80) iavlim(:),avdecayconst
  write(iread,90) ikernav
  write(iread,100) ihvar,hfact,tolh
  write(iread,110) idumpghost
  write(iread,120) imhd,imagforce
  write(iread,130) idivbzero,psidecayfact
  write(iread,140) ianticlump,eps,neps
  write(iread,150) ixsph,xsphfac
  write(iread,160) igravity,hsoft
  write(iread,170) damp, dampz, dampr
  write(iread,180) ikernel
  write(iread,190) iexternal_force
  write(iread,200) C_cour, C_force
 close(unit=iread)

10 format(f14.10,22x,'! particle separation')
20 format(f7.3,2x,f7.3,2x,i9,1x,i5,3x,'! tmax, tout, nmax, nout')
30 format(f14.12,22x,'! gamma ')
40 format(i1,2x,f5.3,28x,'! type of energy equation, polyk(for iener=0)')
50 format(i1,2x,i9,2x,i4,18x,'! type of cty equation (0:direct sum 1:time deriv), ndirect, maxdensits')
60 format(i1,35x,'! type of pressure term (0:normal 1:pa+pb/rhoa*rhob 2:hernquist/katz )')
70 format(i1,2x,f5.3,2x,f5.3,2x,f5.3,2x,f5.3,7x,'! viscosity type, alpha(min), alphau(min), alphab(min), beta')
80 format(7x,i1,6x,i1,6x,i1,2x,f5.3,7x,'! use av, au, ab limiter, constant for this(0.1-0.2)')
90 format(i1,35x,'! type of kernel averaging (1:average h, 2:average grad wab 3:springel/hernquist)')
100 format(i1,2x,f5.3,2x,1pe8.3,18x,'! variable h, initial h factor, h tolerance')
110 format(i1,35x,'! dump ghost particles? (0: no 1: yes)')
120 format(i2,4x,i1,29x,'! magnetic field (0:off 1:on) and force algorithm(1:vector 2:tensor)')
130 format(i2,2x,f5.3,27x,'! divergence correction method (0:none 1:projection 2: hyperbolic/parabolic)')
140 format(i1,2x,f5.3,2x,i2,24x,'! anticlumping term (0:off 1:on), eps, power')
150 format(i1,2x,f5.3,28x,'! use xsph, parameter')
160 format(i1,2x,1pe9.3,24x,'! self-gravity, fixed softening length')
170 format(f7.4,1x,f7.4,1x,f7.4,13x,'! artificial damping (0.0 or few percent)')
180 format(i1,35x,'! kernel type (0: cubic spline, 3:quintic)')
190 format(i1,35x,'! external force (1: toy star, 2:1/r^2 )')
200 format(f7.3,2x,f7.3,2x,18x,'! C_cour, C_force')

 write(iprint,300) infile
300 format (' input file ',a20,' created successfully')

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
 use kernels, only:ianticlump,eps,neps
!
!--define local variables
!      
 implicit none
 character(len=*), intent(in) :: infile 
 character(len=len(infile)+3) :: infilenew   
 character(len=1) :: ians  
!
!--allow for tracing flow
!      
 if (trace) write(iprint,*) ' entering subroutine read_infile'
              
 open(unit=iread,err=999,file=infile,status='old',form='formatted')
  read(iread,*,err=50) psep
  read(iread,*,err=50) tmax,tout,nmax,nout
  read(iread,*,err=50) gamma
  read(iread,*,err=50) iener,polyk
  read(iread,*,err=50) icty,ndirect,maxdensits
  read(iread,*,err=50) iprterm
  read(iread,*,err=50) iav,alphamin,alphaumin,alphabmin,beta
  read(iread,*,err=50) iavlim(:),avdecayconst
  read(iread,*,err=50) ikernav
  read(iread,*,err=50) ihvar,hfact,tolh
  read(iread,*,err=50) idumpghost
  read(iread,*,err=50) imhd,imagforce
  read(iread,*,err=50) idivbzero,psidecayfact
  read(iread,*,err=50) ianticlump,eps,neps
  read(iread,*,err=50) ixsph,xsphfac
  read(iread,*,err=50) igravity,hsoft
  if (ndim.eq.3) then
     read(iread,*,err=50) damp,dampr,dampz
  elseif (ndim.eq.2) then
     read(iread,*,err=50) damp,dampr
  else
     read(iread,*,err=50) damp
  endif
  read(iread,*,err=50) ikernel
  read(iread,*,err=50) iexternal_force
  read(iread,*,err=50) C_Cour, C_force
 close(unit=iread)

 goto 55
50 continue
   close(unit=iread)
   infilenew = trim(infile)//'_new'
   write(iprint,*) 'error reading '//trim(infile)//': writing '//trim(infilenew)//' with current options'
   call write_infile(infilenew)
   stop
55 continue
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
 
100   format(/' read_infile: warning: ',a)
 return
            
999   write(iprint,1000) infile      
1000  format (' input file ',a20,' not found')
 if (adjustl(infile(1:1)).ne.' ') then 
    write(*,*) ' would you like to create one with default options?'
    read*,ians
    if (ians.eq.'y'.or.ians.eq.'y') call write_infile(infile)
 endif

 stop 'exiting...'
      
end subroutine read_infile

end module infiles
