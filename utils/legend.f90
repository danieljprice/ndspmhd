module legends

contains
!
!--draw a legend for different line styles
!  uses current line style and colour
!
subroutine legend(icall,text,hpos,vposin)
  implicit none
  integer, intent(inout) :: icall
  character(len=*), intent(in) :: text
  real, intent(in) :: vposin
  real, dimension(2) :: xline,yline
  real :: xch, ych, xmin, xmax, ymin, ymax
  real :: vspace, hpos, vpos

  call pgstbg(0)           ! opaque text to overwrite previous
!
!--set horizontal and vertical position and spacing
!  in units of the character height
!
  vspace = 1.5  ! (in units of character heights)
!  hpos = 0.4   ! distance from edge, in fraction of viewport
  vpos = vposin + (icall-1)*vspace  ! distance from top, in units of char height

  call pgqwin(xmin,xmax,ymin,ymax) ! query xmax, ymax
  call pgqcs(4,xch,ych) ! query character height in x and y units 

  yline = ymax - ((vpos - 0.5)*ych)
  xline(1) = xmin + hpos*(xmax-xmin)
  xline(2) = xline(1) + 3.*xch

!!--make up line style if > 5 calls (must match actual line drawn)
!   if (icall.eq.3) then
!     call pgpt(2,xline,yline,17) !mod(icall,5)+1)
!     call pgpt(1,0.5*(xline(1)+xline(2)),yline(1),17) !mod(icall,5)+1)
!     !!call pgline(2,xline,yline)            ! draw line segment
   if (icall.gt.0) then
     call pgline(2,xline,yline)            ! draw line segment
   endif
   icall = icall + 1
!
!--write text next to line segment
!  
  call PGTEXT(xline(2) + 0.5*xch,yline(1)-0.25*ych,trim(text))

!!  call pgmtxt('T',vpos,0.1,0.0,trim(text))  ! write text

  call pgstbg(-1) ! reset text background to transparent

end subroutine legend

end module legends
