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
 character(len=120) :: dumpfile,outfile
 real :: tfile,Bx,By,Bz
 logical :: ians,iexist
!
!--get filenames
!
 dumpfile = ' ' 
 call getarg(1,dumpfile)
 if (dumpfile == ' ') then
    call prompt('Enter name of dump to modify:',dumpfile)
 else
    print "(a)",'reading from '//trim(dumpfile)
 endif
 outfile = dumpfile
 call getarg(2,outfile)
 if (outfile == ' ') then
    call prompt('Enter name of output dump :',outfile)
 else
    print "(a)",'writing to '//trim(outfile)
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
 call prompt('Is the dump to read an MHD dump?',ians) 
 imhd = 0 
 if (ians) imhd = 1
 
 call read_dump(dumpfile,tfile,copysetup=.true.)
!
!--------- do modifications here ------------------
! 
 if (imhd.ne.0) then
    call prompt('Change magnetic field?',ians)
 else
    call prompt('Add magnetic field?',ians)
 endif
 if (ians) then
    call prompt('Enter Bx component:',Bx)
    if (ndimV.ge.2) call prompt('Enter By component:',By)
    if (ndimV.ge.3) call prompt('Enter Bz component:',Bz)
    print "(a)",'setting new field components...'
    Bfield(1,:) = Bx
    if (ndimV.ge.2) Bfield(2,:) = By
    if (ndimV.ge.3) Bfield(3,:) = Bz
    imhd = 1
 endif
 tfile = 0.
 
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
