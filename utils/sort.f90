!----------------------------------------
!+
!  Utility program to sort a data file
!  (needed by genetic algorithm utility)
!+
!----------------------------------------
program sort
 use sortutils, only:indexx
 implicit none
 character(len=120) :: filename
 integer :: ierr,n,i,narg
 integer, parameter :: maxk = 100
 real :: dat(3,maxk)
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
   read(1,*,iostat=ierr) dat(1:3,n)
   if (ierr /= 0) n = n - 1
 enddo
 close(unit=1)
 
 ! now sort
 call indexx(n,dat(2,:),index)
 
 ! now write
 do i=1,n
    write(*,"(i2.2,2(1x,es10.3))") nint(dat(1,index(i))),dat(2:3,index(i))
 enddo

end program sort
