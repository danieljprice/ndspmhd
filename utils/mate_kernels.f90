!------------------------------------------
! Dan's kernel mating program
! Written June 2016
! daniel.price@monash.edu
!------------------------------------------
! Combines two kernel tables in linear
! combination
!------------------------------------------
module mate_kernels
 use kernel_utils
 implicit none
 !integer, parameter :: dp = 8
 real(dp), parameter :: radkern = 2.
 real(dp), parameter :: radkern2 = radkern*radkern
 !real(dp), parameter :: pi = 3.1415926536
 
contains

 !------------------------------------------------------
 ! join two kernels to make a new "child"
 ! we take a linear combination with random coefficient
 !------------------------------------------------------
 subroutine mate2(ikern,wkern1,wkern2,wkern,f)
  integer,     intent(in)  :: ikern
  real(dp),    intent(in)  :: wkern1(0:ikern),wkern2(0:ikern)
  real(dp),    intent(out) :: wkern(0:ikern)
  real(dp),    intent(in)  :: f
  integer :: i
  
  do i=0,ikern
     wkern(i) = f*wkern1(i) + (1.-f)*wkern2(i)
  enddo

 end subroutine mate2

 !-----------------------------
 ! read kernel table from file
 !-----------------------------
 subroutine read_kernel(filename,ikern,wkern,ierr)
  character(len=*), intent(in) :: filename
  integer, intent(in)  :: ikern
  real(dp),    intent(out) :: wkern(0:ikern)
  integer, intent(out) :: ierr
  integer :: lu,i
  real(dp) :: c(3),dum
  
  open(newunit=lu,file=filename,status='old',iostat=ierr)
  if (ierr /= 0) return  
  read(lu,*) c(1:3)
  do i=0,ikern
     read(lu,*) dum,wkern(i)
  enddo
  close(lu)

 end subroutine read_kernel
 
 !----------------------------
 ! write kernel table to file
 !----------------------------
 subroutine write_kernel(filename,ikern,c,wkern,grkern,grgrkern,ierr)
  character(len=*), intent(in) :: filename
  integer, intent(in)  :: ikern
  real(dp),    intent(in)  :: wkern(0:ikern),grkern(0:ikern),grgrkern(0:ikern),c(3)
  integer, intent(out) :: ierr
  integer :: lu,i
  real(dp) :: q,dq2table,d2W

  dq2table = radkern2/real(ikern,kind=dp)  
  open(newunit=lu,file=filename,status='replace',iostat=ierr)
  if (ierr /= 0) return
  write(lu,*) c(1:3)
  do i=0,ikern
     q = sqrt(i*dq2table)
     if (q > 0.) then
        d2W = grgrkern(i) + 2.*grgrkern(i)/q
     else
        d2W = 0.
     endif
     write(lu,*) q,wkern(i),grkern(i),grgrkern(i),d2W
     !write(lu,*) q,wkern(i),grkern(i),grgrkern(i),d2W
  enddo
  close(lu)

 end subroutine write_kernel
 
 !-------------------------------------------------------
 ! read two kernels and breed them to create a "child"
 ! call the mutation routine to give a possible mutation
 !-------------------------------------------------------
 subroutine mate_pair(file1,file2,fileout,fac)
  integer, parameter :: ikern = 4000
  character(len=*), intent(in) :: file1,file2,fileout
  real,             intent(in) :: fac
  real(dp), dimension(0:ikern) :: wkern1,wkern2,wkern,grkern,grgrkern
  real(dp) :: c(3)
  integer :: ierr1,ierr2,ierrw

  call read_kernel(file1,ikern,wkern1,ierr1)
  call read_kernel(file2,ikern,wkern2,ierr2)
  if (ierr1 /= 0 .or. ierr2 /= 0) stop 'error reading kernel files'

  call mate2(ikern,wkern1,wkern2,wkern,fac)
  call diff(ikern,wkern,grkern,grgrkern,radkern2)
  call normalise(ikern,wkern,c,radkern2)
  wkern = wkern*c(3)
  c(:) = c(:)/c(3)
  call diff(ikern,wkern,grkern,grgrkern,radkern2)

  call write_kernel(fileout,ikern,c,wkern,grkern,grgrkern,ierrw)
  if (ierrw /= 0) print*,' ERROR during write'

 end subroutine mate_pair

end module mate_kernels

!----------------------------------------------
! driver program for mating procedure
!----------------------------------------------
program mate
 use mate_kernels
 implicit none
 integer :: nargs
 character(len=120) :: file1,file2,fileout,string
 real :: fac
 
 nargs = command_argument_count()
 if (nargs == 3) then
    call get_command_argument(1,file1)
    call get_command_argument(2,file2)
    call get_command_argument(3,string)
    read(string,*) fac
    fileout = 'kernel-new.dat'
    call mate_pair(file1,file2,fileout,fac)
    print "(5a)",trim(file1),' + ',trim(file2),' -> ',trim(fileout)
 elseif (nargs /= 1) then
    print*, 'usage: mate kernel.dat kernel2.dat epsilon'
    stop
 endif

end program mate
