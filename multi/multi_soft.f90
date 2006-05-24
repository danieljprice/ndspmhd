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
 use anticlumping
 
 use infiles
 implicit none
 integer :: i,j,nruns,iener_default
 integer :: ierr
 character :: filename*15,infile*20,filenum*2,charnruns*3,shkfile*20
 
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
!
!-- or read the generic input file
! 
 print*,' reading multirun.in... '
 call read_infile('multirun.in')
 print*,' initial psep = ',psep
 hmin = 1.e-3
 hmax = 10
 dhlog = (log10(hmax) - log10(hmin))/nruns
 
 do i=1,nruns
    if (i.ge.10) then
       filenum = achar(48+i/10)//achar(48+mod(i,10))
       infile = trim(filename)//filenum//'.in'
    else
       filenum(1:1) = achar(48+mod(i,10))
       infile = trim(filename)//filenum(1:1)//'.in'
    endif
    
    hsoft = 10
    print*,' writing input file ',trim(infile), ' psep = ',psep
    
    call write_infile(infile)

 enddo
 
end program multirun
