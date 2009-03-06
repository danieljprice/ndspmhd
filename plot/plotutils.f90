!---------------------------------------------------------------------------
! module containing application programming interfaces for basic
! plotting functions. The idea is to add more to this module to
! eventually use it to be able to change backends more easily.
!---------------------------------------------------------------------------
module plotutils
 implicit none
 public :: plotline,plotbins,formatreal
 
 private

contains

!
!  line plotting, with blanking
!
subroutine plotline(npts,xline,yline,blank)
 implicit none
 integer, intent(in) :: npts
 real, intent(in), dimension(:) :: xline,yline
 real, intent(in), optional :: blank
 integer :: i,nseg,istart
 
 if (present(blank)) then
    nseg = 0
    istart = 1
    !--plot line in segments, leaving blank segments where y=blank
    do i=1,npts
       if (abs(yline(i)-blank).lt.tiny(yline) .or. i.eq.npts) then
          if (nseg.gt.0) call pgline(nseg,xline(istart:istart+nseg),yline(istart:istart+nseg))
          istart = i+1
          nseg = 0
       else
          nseg = min(nseg + 1,npts-1)
       endif
    enddo
 else
    call pgline(npts,xline,yline) 
 endif

 return
end subroutine plotline

!
!  binned histogram plotting, with blanking
!
subroutine plotbins(nbins,xbins,ybins,blank)
 implicit none
 integer, intent(in) :: nbins
 real, intent(in), dimension(:) :: xbins,ybins
 real, intent(in), optional :: blank
 integer :: i,nseg,istart
 
 if (present(blank)) then
    nseg = 0
    istart = 1
    !--plot line in segments, leaving blank segments where y=blank
    do i=1,nbins
       if (abs(ybins(i)-blank).lt.tiny(ybins) .or. i.eq.nbins) then
          if (nseg.gt.0) call pgbin(nseg,xbins(istart:istart+nseg),ybins(istart:istart+nseg),.true.)
          istart = i+1
          nseg = 0
       else
          nseg = min(nseg + 1,nbins-1)
       endif
    enddo
 else
    call pgbin(nbins,xbins,ybins,.true.)
 endif
 
 return
end subroutine plotbins

!
!  formatting of real variables into strings (like PGNUMB)
!
subroutine formatreal(val,string,ierror)
 implicit none
 real, intent(in) :: val
 character(len=*), intent(out) :: string
 integer, intent(out), optional :: ierror
 integer :: ierr,i,idot
 logical :: nonzero

 if (abs(val).ge.1.d99) then
    write(string,"(1pe10.3)",iostat=ierr) val
 elseif (abs(val).lt.1.e-3 .or. abs(val).ge.1.e4) then
    write(string,"(1pe9.2)",iostat=ierr) val
 elseif (abs(val).lt.0.1) then
    write(string,"(f8.3)",iostat=ierr) val
 elseif (abs(val).ge.100.) then
    write(string,"(f8.0)",iostat=ierr) val
 else 
    write(string,"(f8.2)",iostat=ierr) val
 endif
 string = adjustl(trim(string))

 if (present(ierror)) ierror = ierr

 !
 !--strip trailing zeros after the decimal place
 !  (and the decimal place if it is the last character)
 !
 idot = index(string,'.')
 if (idot.gt.0) then
    nonzero = .false.
    do i = len_trim(string),idot,-1
       if (.not.nonzero .and. string(i:i).eq.'0') then
          string(i:i) = ' '
       elseif (.not.nonzero .and. string(i:i).eq.'.') then
          string(i:i) = ' '
          nonzero = .true.
       else
          nonzero = .true.
       endif
    enddo
 endif
 string = trim(string)

 return
end subroutine formatreal

end module plotutils
