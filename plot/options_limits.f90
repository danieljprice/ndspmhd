!----------------------------------------------------------------------
! submenu with options relating to plot limits
!----------------------------------------------------------------------
subroutine options_limits
 use settings_data
 use settings_limits
 use settings_page ! for iadapt
 use multiplot ! for itrans
 use prompting
 use limits
 use labels
 use transforms
 implicit none
 integer :: iaction,ipick,i,ierr
 real :: diff, mid
 
 iaction = 0
 if (iadapt) then
    print 10,iadapt,itrackpart,scalemax
 else
    print 10,iadapt,itrackpart,zoom
 endif
10 format(' 0) exit ',/,                 &
        ' 1) toggle adaptive/fixed limits  ( ',L1,' )   ',/,  &
        ' 2) set manual limits ',/,     &
        ' 3) xy limits track particle      ( ',i8,' )   ',/,   &
        ' 4) zoom in/out                   ( ',f4.2,' ) ',/,   &
        ' 5) apply transformations (log10,1/x) ',/, &
        ' 6) save current limits to file ',/, &
        ' 7) re-read limits file         ',/, &
        ' 8) reset limits for all plots  ')
 call prompt('enter option ',iaction,0,8)
!
!--limits
!
 select case(iaction)
!------------------------------------------------------------------------
 case(1)
    iadapt = .not.iadapt
    print*,'adaptive plot limits = ',iadapt
!------------------------------------------------------------------------
 case(2)
    ipick = 1
    iadapt = .false.
    do while (ipick.gt.0)
       ipick = 0
       call prompt('Enter plot number to set limits (0=finish)',ipick,0,numplot)
       if (ipick.gt.0) then
          call prompt('min ',lim(ipick,1))
          call prompt('max ',lim(ipick,2))
          print*,'>> limits set (min,max) = ',lim(ipick,1),lim(ipick,2)
       endif
    enddo
    return
!------------------------------------------------------------------------
 case(3)
    call prompt('Enter particle to track: ',itrackpart,0)
    print*,'tracking particle ',itrackpart
    if (itrackpart.gt.0) then
       do i=1,ndim
          call prompt('Enter offset for '//trim(label(ix(i)))//'min:', &
                      xminoffset_track(i))
          call prompt('Enter offset for '//trim(label(ix(i)))//'max :', &
                      xmaxoffset_track(i))
       enddo
    endif
!------------------------------------------------------------------------
 case(4)
    if (.not.iadapt) then
       call prompt('Enter zoom factor for fixed limits',zoom,0.0)
       do i=1,numplot
          diff = lim(i,2)- lim(i,1)
          mid = 0.5*(lim(i,1) + lim(i,2))
          lim(i,1) = mid - 0.5*zoom*diff
          lim(i,2) = mid + 0.5*zoom*diff
       enddo
    else
       call prompt('Enter scale factor (adaptive limits)',scalemax,0.0)
    endif
!------------------------------------------------------------------------
  case(5)
     ipick = 1
     do while (ipick.gt.0 .and. ipick.le.numplot)
        ipick = 0
        print*,'Enter plot number to apply transformation '
        call prompt('(0 = finish, -1 = set all) ',ipick)
        if (ipick.le.numplot .and. ipick.ne.0) then
           write(*,300) (i,trim(transform_label('x',i)), i=1,5)
300        format(1x,i1,') ',a)
           print*,'Enter transformation to apply (or a combination e.g. 321)'
           if (ipick.lt.0) then
              ipick = 0
              call prompt(' ',ipick,0)
              itrans(:) = ipick
              ipick = -99
           else
              call prompt(' ',itrans(ipick),0)
           endif
        endif
     enddo
     return
  case(6)
     call save_limits 
  case(7)
     call read_limits(ierr)
  case(8)
     if (ivegotdata) then
        call set_limits(nstart,n_end,1,numplot)
     else
        print*,'no data with which to set limits!!'
     endif
  end select
 
 return
end subroutine options_limits
