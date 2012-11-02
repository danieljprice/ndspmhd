!!-----------------------------------------------------------------------
!! writes multiple infiles where one or several parameters are varied
!! this version changes equation types etc
!!-----------------------------------------------------------------------
program multirun
 use dimen_mhd
 use loguns
 use artvi
 use options
 use setup_params
 use timestep
 use xsph
 use kernels, only:setkern,radkern
 
 use infiles
 implicit none
 integer :: i,j,nruns,nres
 integer :: ierr
 character :: filename*15,infile*30
 character(len=10) :: string
 character(len=6) :: charnruns,charhfac,chardhfac
 real :: C_cour_init,courmin,courmax,dcour,dhfac
 logical :: logspace
 
 nruns = 0
 iread = 11
 call getarg(1,filename)
 call getarg(2,charnruns)
 call getarg(3,charhfac)
 call getarg(4,chardhfac)
 
 read(charnruns,*,iostat=ierr) nruns
 
 if (filename.eq.'' .or. ierr.ne.0) then
  print*,'usage: multirun filename nruns start interval'
  stop
 endif
 print*,' filename = ',filename,' nruns = ',nruns 
!
!--set default options
!
 call set_default_options
!
!--read the generic input file
! 
 print*,' reading multirun.in... '
 call read_infile('multirun.in')

 if (command_argument_count().ge.4) then
    read(charhfac,*,iostat=ierr) hfact 
    if (ierr.ne.0) hfact = 0.75
    read(chardhfac,*,iostat=ierr) dhfac
    if (ierr.ne.0) dhfac = 0.025
 else
    hfact = 0.75
    dhfac = 0.025
 endif

 call setkern(ikernel,ndim,ierr)
 print*,' ikernel = ',ikernel
 print*,' radkernel = ',radkern
 hfact = hfact/(0.5*radkern)
 print*,' initial hfact = ',hfact

 do i=1,nruns
    hfact = hfact + dhfac
    write(string,"(f5.3)") hfact
    print*,' hfact = ',trim(string)
     
    infile = trim(filename)//trim(adjustl(string))//'.in'

    print*,' writing input file ',infile
    call write_infile(infile)
 enddo
 
end program multirun
