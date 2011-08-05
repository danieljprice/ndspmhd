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
 real :: C_cour_init
 
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
 
 C_cour_init = C_cour
 C_cour = C_cour*sqrt(2.)
 do i=1,nruns
    C_cour = C_cour*1./sqrt(2.)
    dt = 1.078e-3*(psep/0.01)*(C_cour/C_cour_init)
    write(string,"(es10.3)") dt
     
    infile = trim(filename)//trim(adjustl(string))//'.in'

    print*,' writing input file ',infile
    call write_infile(infile)
 enddo
 
end program multirun
