!
! this module contains certain useful utilities which
! depend on system commands (in the system_commands module)
!
module system_utils
 use system_commands, only:get_environment
 implicit none
 public :: ienvironment,lenvironment,renvironment

contains
 !
 !--this routine returns an integer variable
 !  from an environment variable setting
 !
 !  if the errval argument is present then this is the
 !  value assigned when an error has occurred
 !  (default is zero)
 !
 integer function ienvironment(variable,errval)
    character(len=*), intent(in) :: variable
    character(len=30) :: string
    character(len=5) :: fmtstring
    integer, intent(in), optional :: errval
    integer :: ierr
    
    call get_environment(variable,string)
    if (len_trim(string).gt.0) then  
    !--use a formatted read - this is to avoid a compiler bug
    !  should in general be more robust anyway
       write(fmtstring,"(a,i2,a)",iostat=ierr) '(i',len_trim(string),')'
       read(string,fmtstring,iostat=ierr) ienvironment
    else
       ierr = 1
    endif    

    if (ierr /= 0) then
       if (present(errval)) then
          ienvironment = errval 
       else
          ienvironment = 0
       endif
    endif
 
 end function ienvironment

 !
 !--this routine returns a floating point (real) variable
 !  from an environment variable setting
 !
 !  if the errval argument is present then this is the
 !  value assigned when an error has occurred
 !  (default is zero)
 !
 real function renvironment(variable,errval)
    character(len=*), intent(in) :: variable
    character(len=30) :: string
    real, intent(in), optional :: errval
    integer :: ierr
    
    call get_environment(variable,string)
    read(string,*,iostat=ierr) renvironment
    if (ierr /= 0) then
       if (present(errval)) then
          renvironment = errval 
       else
          renvironment = 0.
       endif
    endif
 
 end function renvironment

 !
 !--this routine returns a logical variable
 !  from an environment variable setting
 !
 logical function lenvironment(variable)
    character(len=*), intent(in) :: variable
    character(len=30) :: string
    
    call get_environment(variable,string)
    if (string(1:1).eq.'y'.or.string(1:1).eq.'Y' &
    .or.string(1:1).eq.'t'.or.string(1:1).eq.'T' &
    .or.trim(string).eq.'on'.or.trim(string).eq.'ON') then
       lenvironment = .true.
    else
       lenvironment = .false.
    endif
 
 end function lenvironment

 !
 !--this routine returns an arbitrary number of 
 !  comma separated strings
 !
 subroutine envlist(variable,nlist,list)
    implicit none
    character(len=*), intent(in) :: variable
    integer, intent(out) :: nlist
    character(len=*), dimension(:), intent(out), optional :: list
    character(len=120) :: string
    character(len=10) :: dummy
    integer :: i1,i2,ierr
    logical :: notlistfull
    
    !--set list to blank strings if argument is present
    if (present(list)) then
       list = ' '
    endif
    
    !--get envlist from the environment
    call get_environment(variable,string)
    
    !--split the string on commas
    i1 = 1
    i2 = index(string,',')-1
    if (i2.eq.-1) i2 = len_trim(string)
    nlist = 0
    ierr = 0
    notlistfull = .true.

    !--for each comma separated string, add a list member
    do while(i2.ge.i1 .and. notlistfull .and. ierr.eq.0)
       nlist = nlist + 1
       !print*,'i1,i2,stringfrag= ',i1,i2,trim(string(i1:))
       if (present(list)) then
          read(string(i1:i2),"(a)",iostat=ierr) list(nlist)
          notlistfull = (nlist.lt.size(list))
       else
          read(string(i1:i2),"(a)",iostat=ierr) dummy
          notlistfull = .true.
       endif
       i1 = i2 + 2
       i2 = index(string(i1:),',')
       if (i2.eq.0) then
          i2 = len_trim(string)
       else
          i2 = i2 + i1 - 2
       endif
    enddo
    
    return
 end subroutine envlist
 
end module system_utils
