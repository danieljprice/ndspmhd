!
!--utility to modify the contents of an output dump
!  (e.g. for restarting)
!
program modifydump
 use loguns
 use part
 use rates
 use options
 use derivB
 
 use dumpfiles
 use prompting
 implicit none
 integer :: i
 character(len=120) :: dumpfile,outfile
 real :: tfile
 logical :: ians,iexist
!
!--get filenames
!
 dumpfile = ' ' 
 call getarg(1,dumpfile)
 if (dumpfile == ' ') then
    call prompt('Enter name of dump to modify:',dumpfile)
 endif
 outfile = dumpfile
 call getarg(2,outfile)
 if (outfile == ' ') then
    call prompt('Enter name of output dump :',outfile)
 endif
!
!--check if output file exists
!
 inquire(file=outfile,exist=iexist)
 if (iexist) then
    ians = .true.
    call prompt(trim(outfile)//' already exists. Overwrite?',ians)
    if (.not.ians) stop
 endif
!
!--set unit numbers for i/o
!
 ireadf = 1
 idatfile = 2
 iprint = 6
!
!--set options required for output
!
 ians = .false.
 if (dumpfile(1:1) == 'm') ians = .true.
 call prompt('MHD dump?',ians) 
 imhd = 0 
 if (ians) imhd = 1
 
 call read_dump(dumpfile,tfile)
!
!--------- do modifications here ------------------
! 
 do i=1,npart 
    !vel(:,i) = x(:,i)
    !vel(1,i) = vel(1,i) - x(2,i)
    !vel(2,i) = vel(2,i) + x(1,i)
    !vel(:,i) = 0.
    Bfield(1:2,i) = 0.01
    Bfield(3,i) = 0.
 enddo
 tfile = 0.
 imhd = 1
 
!--------------------------------------------------
!
!--set unread variables to zero
!
 pr = 0.
 drhodt = 0.
 divB = 0.
 curlB = 0.
 psi = 0.
 idumpghost = 0
 
 print "(a)",'>> Writing modified dump to file '
 print *,' t = ',tfile,' file = ',trim(outfile)
 call write_dump(tfile,outfile)
 
end program modifydump
