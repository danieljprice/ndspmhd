!----------------------------------------
!+
!  Utility program to sort a data file
!  (needed by genetic algorithm utility)
!+
!----------------------------------------
program sort
 use sortutils, only:indexx
 implicit none
 character(len=120) :: filename,line
 integer :: ierr,n,i,narg,ierr1
 integer, parameter :: maxk = 100
 real :: dat(4,maxk)
 integer :: index(maxk)
 
 narg = command_argument_count()
 if (narg /= 1) stop 'usage: sort filename'
 
 call get_command_argument(1,filename)
 
 open(unit=1,file=filename,status='old',iostat=ierr)
 if (ierr /= 0) then
    print*,'error opening '//trim(filename)
    stop
 endif
 n = 0
 ierr = 0
 do while (ierr==0)
   n = n + 1
   read(1,"(a)",iostat=ierr) line
   read(line,*,iostat=ierr1) dat(1:4,n)
   if (ierr /= 0) n = n - 1
 enddo
 close(unit=1)
 
 ! now sort
 call indexx(n,dat(2,:),index)
 
 ! now write
 do i=1,n
    write(*,"(i2.2,1(1x,es10.3),1x,i2)") nint(dat(1,index(i))),dat(2,index(i)),nint(dat(3,index(i)))
 enddo

end program sort
