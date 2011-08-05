!!-----------------------------------------------------------------------
!! computes multiple infiles where one or several parameters are varied
!!
!! this varies h and the div b cleaning parameters
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
 use infiles
 implicit none
 integer :: i,j,nruns,iener_default
 integer :: ierr
 character :: filename*15,infile*20,filenum*2,charnruns*3,shkfile*20,hfile*25
 real :: hmin,hminlog,hmax,dhlog,hneighmin,hneighmax,hneigh,dhneigh
 
 nruns = 0
 iread = 11
 call getarg(1,filename)
 call getarg(2,charnruns)
 
 read(charnruns,*,iostat=ierr) nruns
 
 if (filename.eq.'' .or. ierr.ne.0) then
  print*,'usage: multirun filename nruns'
  stop
 endif
 
 print*,' filename = ',filename,' nruns = ',nruns 
  
!
!--set default options
!
 call set_default_options
 print*,'tolh = ',tolh
 print*,'hsoft = ',hsoft
!
!-- or read the generic input file
! 
 print*,' reading multirun.in... '
 call read_infile('multirun.in')
 print*,' initial psep = ',psep
 
 if (igravity.le.2) then
    hmin = 1.e-3
    hminlog = log10(hmin)
    hmax = 10.
    dhlog = (log10(hmax) - hminlog)/(nruns-1)
    hfile = 'hsoft_'//trim(filename)
 else
    hneighmin = 16.
    hminlog = log10(hneighmin)
    hneighmax = 1000.
    dhneigh = (log10(hneighmax) - log10(hneighmin))/(nruns-1)
    hfile = 'hfact_'//trim(filename)
 endif
 
 open(99,file=hfile,status='replace',form='formatted')
 do i=1,nruns
    if (i.ge.10) then
       filenum = achar(48+i/10)//achar(48+mod(i,10))
       infile = trim(filename)//filenum//'.in'
    else
       filenum(1:1) = achar(48+mod(i,10))
       infile = trim(filename)//filenum(1:1)//'.in'
    endif
    
    if (igravity.le.2) then
       hsoft = 10.**((i-1)*dhlog + hminlog)
       print*,' writing input file ',trim(infile), ' hsoft = ',hsoft
       write(99,*) hsoft,infile(1:len_trim(infile)-3)
    else
       hneigh = 10.**(hminlog + (i-1)*dhneigh)
       hfact = ((3.*hneigh)/(32.*pi))**(1./3.)
       print*,' writing input file ',trim(infile),' Nneigh = ',hneigh,' hfact = ',hfact
       write(99,*) hfact,infile(1:len_trim(infile)-3)
    endif
    
    call write_infile(infile)

 enddo
 close(99)
 
end program multirun
