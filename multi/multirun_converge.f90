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
 
 use infiles
 implicit none
 integer :: i,j,nruns,nres
 integer :: ierr
 character :: filename*15,infile*30
 character(len=10) :: string
 character(len=2) :: charnruns
 real :: C_cour_init,courmin,courmax,dcour
 logical :: logspace
 
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
!--read the generic input file
! 
 print*,' reading multirun.in... '
 call read_infile('multirun.in')
 print*,' initial psep = ',psep

 do i=1,nruns
    psep = psep/2.
    print*,' psep = ',psep, 1./psep
    !write(string,"(i4.4)") int(1./psep)
    write(string,"(es8.2)") psep
     
    infile = trim(filename)//trim(adjustl(string))//'.in'

    print*,' writing input file ',infile
    call write_infile(infile)
 enddo
 
end program multirun
