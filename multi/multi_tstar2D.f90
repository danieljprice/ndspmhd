!!-----------------------------------------------------------------------
!! computes multiple infiles where one or several parameters are varied
!!
!! this one does toy star oscillations in 2D
!!-----------------------------------------------------------------------
program multirun
 use dimen_mhd
 use loguns
 use artvi
 use eos
 use options
 use setup_params
 use timestep
 use xsph
 use anticlumping
 
 use infiles
 implicit none
 integer :: i,j,nruns,jmode,mmode
 integer :: ierr
 character :: filename*115,infile*120,tstarfile*120,filenum*2,charnruns*3
 real :: h,c,a,sigma2,sigma,gamm1,period,omegasq
 real :: alpha,betatstar
 
 nruns = 0
 iread = 11
 call getarg(1,filename)
 call getarg(2,charnruns)
 
 read(charnruns,*,iostat=ierr) nruns
 
 if (filename.eq.'' .or. ierr.ne.0) then
  print*,'usage: multirun filename nruns'
  stop
 endif
 
 print*,' filename = ',trim(filename),' nruns = ',nruns
!
!--set default options
!
 call set_default_options

!
!--read the generic input file
! 
 print*,' reading multirun.in... '
 call read_infile('multirun.in')
 print*,' initial psep = ',psep
  
 jmode = 0
 mmode = 0 
 alpha = 1.0
 betatstar = 2.*pi
 
 do i=1,nruns
    if (i.ge.10) then
       filenum = achar(48+i/10)//achar(48+mod(i,10))
       infile = trim(filename)//filenum//'.in'
       tstarfile = trim(filename)//filenum//'.tstar2D'
    else
       filenum(1:1) = achar(48+mod(i,10))
       infile = trim(filename)//filenum(1:1)//'.in'
       tstarfile = trim(filename)//filenum(1:1)//'.tstar2D'
    endif
    print*,' writing input file ',trim(infile), ' psep = ',psep
    
    h = 1.0
    c = 1.0
    a = 0.05/sqrt(2.)
    !!jmode = jmode + 2
    mmode = mmode + 2

    gamm1 = gamma - 1.
    if (gamm1.lt.1.e-5) then
       stop 'error: gamma - 1 <= 0'
    endif

    omegasq = 1.0
    sigma2 = 0.5*omegasq*(gamm1)*((jmode+mmode)*(jmode+mmode + 2./gamm1) - mmode**2)
    if (sigma2.lt.1.e-5) then
       print*,'ERROR sigma2 <= 0 in perturbation'
       sigma = 1.e-1
    else
       sigma = sqrt(sigma2)
    endif
!
!--run for 10 periods with output every period
!
    period = 2.*pi/sigma
    tout = 0.25*period
    tmax = 10.*period
    print*,'jmode = ',jmode,' smode = ',mmode,' period = ',period
    
    print*,'writing ',trim(tstarfile),' with initial left/right states'
    open(unit=11,file=tstarfile,status='replace',form='formatted')
       write(11,*) h,c,a
       write(11,*) alpha,betatstar
       write(11,*) jmode,mmode
    close(unit=11)
    
    call write_infile(infile)

 enddo
 
end program multirun
