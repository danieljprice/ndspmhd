program moddump
 use loguns
 use dumpfiles
 use prompting
 use options
 use part
 use derivB
 use rates, only:drhodt
 use eos, only:polyk
 use setup_params, only:geomsetup
 implicit none
 character(len=120) :: dumpfilein,dumpfileout
 logical :: iexist,iwrite,ians
 real :: tfile,t
 integer :: idash
!
!--set unit numbers for i/o
!
 ireadf = 1
 idatfile = 2
 iprint = 6
!
!--read command line arguments
! 
 call getarg(1,dumpfilein)
 if (dumpfilein == ' ') then
    call prompt('Enter name of input dump :',dumpfilein)
 endif 
 call getarg(2,dumpfileout)
 if (dumpfileout == ' ') then
    call prompt('Enter name of output dump :',dumpfileout)
 endif
!
!--set rootname from output filename
!
 idash = index(dumpfileout,'_')
 rootname = dumpfileout(1:idash-1)
 print*,'rootname = ',trim(rootname)
!
!--read input dump
!
 call read_dump(dumpfilein,tfile,copysetup=.true.)
 geom = geomsetup
!
!--perform modifications
!
 call prompt('Enter polytropic k if relevant',polyk,0.)

 call modify_dump
!
!--set new time
!
 print*
 t = tfile
 call prompt('Enter new time for dump',t)
!
!--turn MHD on/off
!
 ians = .false.
 if (dumpfilein(1:1) == 'm') ians = .true.
 call prompt('MHD dump?',ians) 
 imhd = 0 
 if (ians) imhd = 1
 
 inquire(file=dumpfileout,exist=iexist)
 iwrite = .true.
 if (iexist) then
    call prompt(trim(dumpfileout)//': file exists. Overwrite?',iwrite)
 endif

!
!--set unread variables to zero
!
 pr = 0.
 drhodt = 0.
 divB = 0.
 curlB = 0.
 psi = 0.
 idumpghost = 0
  
 if (iwrite) then 
    print "(a)",'>> Writing modified dump to file '//trim(dumpfileout)
    call write_dump(t,dumpfileout)
 else
    print*,' no output written '
 endif
 
 stop
end program moddump
