module interactive_routines
 implicit none
 public :: interactive_part,interactive_step
 private :: mvlegend,mvtitle,save_limits,save_rotation
 real, private :: xpt
 real, private :: ypt
 
contains
!
!--interactive tools on particle plots
!  allows user to change settings interactively
!
!  Arguments:
!
!  INPUT:
!   npart   : number of particles plotted
!   iplotx  : quantity plotted as x axis
!   iploty  : quantity plotted as y axis 
!   iplotz  : quantity to use in selecting particles
!   irender : quantity rendered
!   xcoords(npart) : x coordinates of particles
!   ycoords(npart) : y coordinates of particles
!   zcoords(npart) : z coordinates (or third quantity) of particles
!   hi(npart)      : smoothing lengths of particles
!   zmin, zmax     : range of z within which to plot particles
!   istep          : current step position
!   ilaststep      : position of last timestep
!
! CHANGEABLE:
!   icolourpart(npart) : flag indicating colour of particles
!   xmin, xmax, ymin, ymax : current plot limits
!   rendermin, rendermax   : current rendering limits
!   vecmax : maximum vector limits
!
! OUTPUT:
!   iadvance : integer telling the loop how to advance the timestep
!
subroutine interactive_part(npart,iplotx,iploty,iplotz,irender,ivecx,ivecy, &
  xcoords,ycoords,zcoords,hi,icolourpart,xmin,xmax,ymin,ymax, &
  rendermin,rendermax,vecmax,anglex,angley,anglez,ndim,x_sec,zslicepos,dzslice, &
  zobserver,dscreen,itrackpart,icolourscheme,iadvance,istep,ilaststep)
  implicit none
  integer, intent(in) :: npart,irender,ndim,iplotz,ivecx,ivecy,istep,ilaststep
  integer, intent(inout) :: iplotx,iploty,itrackpart,icolourscheme
  integer, intent(out) :: iadvance
  integer, dimension(npart), intent(inout) :: icolourpart
  real, dimension(npart), intent(in) :: xcoords,ycoords,zcoords,hi
  real, intent(inout) :: xmin,xmax,ymin,ymax,rendermin,rendermax,vecmax
  real, intent(inout) :: anglex,angley,anglez,zslicepos,dzslice,zobserver,dscreen
  logical, intent(in) :: x_sec
  real, parameter :: pi=3.141592653589
  integer :: i,iclosest,nc,ierr,ixsec
  integer :: nmarked,ncircpart,itrackparttemp
  integer, dimension(npart) :: icircpart
 !! real :: xpt,ypt
  real :: xpt2,ypt2,charheight
  real :: xptmin,xptmax,yptmin,yptmax,zptmin,zptmax
  real :: rmin,rr,gradient,yint,dx,dy,dr,anglerad
  real :: xlength,ylength,renderlength,renderpt,drender
  real, dimension(4) :: xline,yline
  character(len=1) :: char,char2
  character(len=20) :: string
  logical :: iexit, rotation

  call pgqinf('CURSOR',string,nc)
  if (string(1:nc).eq.'YES') then
     print*,'entering interactive mode...press h in plot window for help'
  else
     print*,'cannot enter interactive mode: device has no cursor'
     return
  endif
  char = 'A'
  xline = 0.
  yline = 0.
!  xpt = 0.
!  ypt = 0.
  xpt2 = 0.
  ypt2 = 0.
  nc = 0
  ncircpart = 0
  itrackparttemp = itrackpart
  iexit = .false.
  rotation = .false.
  if (iplotx.le.ndim .and. iploty.le.ndim .and. ndim.ge.2) rotation = .true.
  
  if (iplotz.gt.0 .and. x_sec) then
     zptmin = zslicepos - 0.5*dzslice
     zptmax = zslicepos + 0.5*dzslice
  else
  !--if not using z range, make it encompass all the particles
     zptmin = -huge(zptmin)
     zptmax = huge(zptmax)
  endif
  
  interactiveloop: do while (.not.iexit)
     call pgcurs(xpt,ypt,char)
     !
     !--exit if the device is not interactive
     !
     if (char.eq.achar(0)) return
  
     !
     !--find closest particle
     !  
     rmin = 1.e6
     do i=1,npart
        rr = (xcoords(i)-xpt)**2 + (ycoords(i)-ypt)**2
        if (rr.lt.rmin) then
           iclosest = i
           rmin = rr
        endif
     enddo
     
     select case(char)
     !
     !--particle plot stuff
     !
     case('p')
        print*,' closest particle = ',iclosest,'x = ',xcoords(iclosest),' y =',ycoords(iclosest)
        call pgnumb(iclosest,0,1,string,nc)
        call pgqch(charheight)
        call pgsch(2.0)
        call pgtext(xcoords(iclosest),ycoords(iclosest),string(1:nc))
        call pgsch(charheight)
     case('c','C')
        print*,'plotting circle of interaction on particle ',iclosest, &
               ' h = ',hi(iclosest)
        !--save settings for these
        ncircpart = ncircpart + 1
        if (ncircpart.gt.npart) then
           print*,'WARNING: ncircles > array limits, cannot save'
           ncircpart = npart
        else
           icircpart(ncircpart) = iclosest
        endif
        call pgsfs(2)
        call pgcirc(xcoords(iclosest),ycoords(iclosest),2.*hi(iclosest))
     case('t')
     !--track closest particle (must save to activate)
        if (iplotx.le.ndim .and. iploty.le.ndim) then
           if (itrackpart.ne.0 .and. itrackparttemp.eq.itrackpart) then
              itrackpart = 0
              itrackparttemp = 0
              print*,' particle tracking limits OFF'
           else
              itrackparttemp = iclosest
              print*,' limits set to track particle ',itrackparttemp
              print*,' save settings to activate '
           endif
        endif
     case('g')   ! draw a line between two points
        xline(2) = xpt
        yline(2) = ypt
        !--mark first point
        call pgpt1(xpt,ypt,4)
        !--select second point
        print*,' select another point (using left click or g) to plot line '
        call pgband(1,1,xline(2),yline(2),xline(3),yline(3),char2)
        !--draw line if left click or g
        select case(char2)
        case('A','g')
           !--mark second point
           call pgpt1(xline(3),yline(3),4)
           xlength = xline(3)-xline(2)
           if (abs(xlength).lt.tiny(xlength)) then
              xline(1) = xline(2)
              xline(4) = xline(2)
              yline(1) = ymin
              yline(4) = ymax
              print*,' error: gradient = infinite'
           elseif (xline(2).lt.xline(3)) then 
              xline(1) = xmin
              xline(4) = xmax
           else
              xline(1) = xmax
              xline(4) = xmin
           endif           
           if (abs(xlength).gt.tiny(xlength)) then
              ylength = yline(3)-yline(2)
              gradient = ylength/xlength
              yint = yline(3) - gradient*xline(3)
              dr = sqrt(xlength**2 + ylength**2)
              print*,' dr = ',dr,' dx = ',xlength,' dy = ',ylength
              print*,' gradient = ',gradient,' y intercept = ',yint
              yline(1) = gradient*xline(1) + yint
              yline(4) = gradient*xline(4) + yint
           endif
           !--plot line joining the two points
           call pgline(4,xline,yline)
        case default
           print*,' action cancelled'       
        end select
     !
     !--help
     !
     case('h')
        print*,'-------------- interactive mode commands --------------'
        print*,' left click (or A) : select region'
        print*,'  left click again to zoom on selection'
        if (irender.ne.0) then
           print*,'  select colour bar to change rendering limits'
        else
           print*,'  or press 1-9 to mark selected particles with colour 1-9'
        endif
        print*,' +   : zoom in by 10%'
        print*,' -(_): zoom out by 10(20)%'
        print*,' a   : (a)djust/reset plot limits to fit '
        print*,'      NB for these options cursor inside the plot changes both x and y,'
        print*,'       whereas cursor over a specific axis zooms on that axis only'
        if (irender.ne.0) then
           print*,'     (applies to colour bar if mouse is over colour bar)'
        endif
        print*,' l: (l)og / unlog axis (with cursor over the axis to change)'  
        if (irender.ne.0) then
           print*,'  (applies to colour bar if mouse is over colour bar)'
        endif
        print*,' o: re-centre plot on (o)rigin'
        print*,' r: (r)eplot current plot'
        print*,' p: label closest (p)article'
        if (iplotx.le.ndim .and. iploty.le.ndim) then
           print*,' t: t)rack closest particle/turn tracking off'
        endif
        print*,' c: plot (c)ircle of interaction for closest particle'        
        print*,' g: plot a line and find its g)radient'
        print*,' G: move le(G)end to current position'
        print*,' T: move (T)itle to current position'
        print*,' H: move vector plot legend H(ere)'
        if (irender.ge.0) then
           print*,' m: change colour m)ap to next'
           print*,' M: change colour M)ap to previous'
           print*,' i: i)nvert colour map'
        endif
        print*,' v: decrease arrow size on vector plots'
        print*,' V: increase arrow size on vector plots'
        if (rotation) then
           print*,' , .: rotate about z axis by +(-) 15 degrees'
           print*,' < >: rotate about z axis by +(-) 30 degrees'
           if (ndim.ge.3) then
              print*,'[ ]: rotate about x axis by +/- 15 degrees '
              print*,'{ }: rotate about x axis by +/- 30 degrees '
              print*,'/ \: rotate about y axis by +/- 15 degrees '
              print*,'? |: rotate about y axis by +/- 30 degrees '
              print*,' x: take cross section '
              if (iplotz.gt.0) then
                 print*,' u: move cross section/perspective position up (towards observer)'
                 print*,' U: move cross section/perspective position up more'
                 print*,' d: move cross section/perspective position down (away from observer)'
                 print*,' D: move cross section/perspective position down more'
              endif
           endif
        endif
        print*,' next timestep/plot   : space, n'
        print*,' previous timestep    : right click (or X), b'
        print*,' jump forward (back) by n timesteps  : 0,1,2,3..9 then left (right) click'
        print*,' h: (h)elp'
        print*,' s: (s)ave current settings for all steps'
        print*,' q,Q: (q)uit plotting'             
        print*,'-------------------------------------------------------'
     case('s','S')
        if (iplotx.le.ndim .and. iploty.le.ndim) itrackpart = itrackparttemp
        if (itrackpart.eq.0) then
           call save_limits(iplotx,xmin,xmax)
           call save_limits(iploty,ymin,ymax)
        else
           call save_limits_track(iplotx,xmin,xmax,xcoords(itrackpart))
           call save_limits_track(iploty,ymin,ymax,ycoords(itrackpart))
        endif
        if (irender.gt.0) call save_limits(irender,rendermin,rendermax)
        if (ivecx.gt.0 .and. ivecy.gt.0) then
           call save_limits(ivecx,-vecmax,vecmax)
           call save_limits(ivecy,-vecmax,vecmax)
        endif
        if (ncircpart.gt.0) call save_circles(ncircpart,icircpart)
        if (rotation) call save_rotation(ndim,anglex,angley,anglez)
        if (iplotz.gt.0) then
           if (x_sec) then
              call save_xsecpos(zslicepos)
           else
              call save_perspective(zobserver,dscreen)
           endif
        endif
        print*,'> plot settings saved <'
     !
     !--actions on left click
     !
     case('A') ! left click
        print*,'select area: '
        print*,'left click : zoom'
        !
        !--change colour bar limits
        !
        if (xpt.gt.xmax .and. irender.gt.0) then
           print*,'click to set rendering limits'
           call pgband(3,1,xpt,ypt,xpt2,ypt2,char2)
           if (char2 == 'A') then
              drender = (rendermax-rendermin)/(ymax-ymin)
              rendermax = rendermin + (max(ypt,ypt2)-ymin)*drender
              rendermin = rendermin + (min(ypt,ypt2)-ymin)*drender
              print*,'setting render min = ',rendermin
              print*,'setting render max = ',rendermax              
              iadvance = 0
              iexit = .true.
           endif
        !
        !--zoom or mark particles
        !
        else
           if (irender.le.0) then
              print*,'1-9 = mark selected particles with colours 1-9'
              print*,'0 = hide selected particles'
              print*,'p = plot selected particles only'
              print*,'c = plot circles of interaction on selected parts'
           endif
           call pgband(2,1,xpt,ypt,xpt2,ypt2,char2)
           xptmin = min(xpt,xpt2)
           xptmax = max(xpt,xpt2)
           yptmin = min(ypt,ypt2)
           yptmax = max(ypt,ypt2)
           print*,'xrange = ',xptmin,'->',xptmax
           print*,'yrange = ',yptmin,'->',yptmax
           if (iplotz.ne.0 .and. x_sec) then
              print*,'(zrange = ',zptmin,'->',zptmax,')'
           endif
           select case (char2)
           case('A')   ! zoom if another left click
              call pgsfs(2)
              call pgrect(xpt,xpt2,ypt,ypt2)
              xmin = xptmin
              xmax = xptmax
              ymin = yptmin
              ymax = yptmax
              iadvance = 0
              iexit = .true.
           case('0','1','2','3','4','5','6','7','8','9') ! mark particles
              if (irender.le.0) then
                 nmarked = 0
                 do i=1,npart
                    if ((xcoords(i).ge.xptmin .and. xcoords(i).le.xptmax) &
                    .and.(ycoords(i).ge.yptmin .and. ycoords(i).le.yptmax) &
                    .and.(zcoords(i).ge.zptmin .and. zcoords(i).le.zptmax)) then
                        read(char2,*,iostat=ierr) icolourpart(i)
                        if (ierr /=0) then
                           print*,'*** error marking particle' 
                           icolourpart(i) = 1
                        endif
                        nmarked = nmarked + 1
                    endif
                 enddo
                 print*,'marked ',nmarked,' particles in selected region'
              endif
              iadvance = 0
              iexit = .true.
           case('p') ! plot selected particles only
              if (irender.le.0) then
                 nmarked = 0
                 do i=1,npart
                    if ((xcoords(i).ge.xptmin .and. xcoords(i).le.xptmax) &
                    .and.(ycoords(i).ge.yptmin .and. ycoords(i).le.yptmax) &
                    .and.(zcoords(i).ge.zptmin .and. zcoords(i).le.zptmax)) then
                       nmarked = nmarked + 1
                       if (icolourpart(i).le.0) icolourpart(i) = 1
                    else
                       icolourpart(i) = 0
                    endif
                 enddo
                 print*,'plotting selected ',nmarked,' particles only'
              endif
              iadvance = 0
              iexit = .true.           
           case('c') ! set circles of interaction in marked region
              if (irender.le.0) then
                 ncircpart = 0
                 do i=1,npart
                    if ((xcoords(i).ge.xptmin .and. xcoords(i).le.xptmax) &
                    .and.(ycoords(i).ge.yptmin .and. ycoords(i).le.yptmax) &
                    .and.(zcoords(i).ge.zptmin .and. zcoords(i).le.zptmax)) then
                        ncircpart = ncircpart + 1
                        icircpart(ncircpart) = i
                        call pgsfs(2)
                        call pgcirc(xcoords(i),ycoords(i),2.*hi(i))
                    endif
                 enddo
                 print*,'set ',ncircpart,' circles of interaction in selected region'
              endif
           case default
              print*,' action cancelled'
           end select 
        endif   
     !
     !--zooming
     !
     case('-','_','+','o') ! zoom out by 10 or 20%
        xlength = xmax - xmin
        ylength = ymax - ymin
        renderlength = rendermax - rendermin
        select case(char)
        case('-')
           xlength = 1.1*xlength
           ylength = 1.1*ylength
           renderlength = 1.1*renderlength
        case('_')
           xlength = 1.2*xlength
           ylength = 1.2*ylength
           renderlength = 1.2*renderlength
        case('+')
           xlength = 0.9*xlength
           ylength = 0.9*ylength
           renderlength = 0.9*renderlength
        case('o') !--reset cursor to origin
           xpt = 0.
           ypt = 0.
        end select
        if (xpt.ge.xmin .and. xpt.le.xmax .and. ypt.le.ymax) then
           print*,'zooming on x axis'
           xmin = xpt - 0.5*xlength
           xmax = xpt + 0.5*xlength
           iadvance = 0
           iexit = .true.
        endif
        if (ypt.ge.ymin .and. ypt.le.ymax .and. xpt.le.xmax) then
           print*,'zooming on y axis'
           ymin = ypt - 0.5*ylength
           ymax = ypt + 0.5*ylength
           iadvance = 0
           iexit = .true.
        endif
        if (xpt.gt.xmax .and. irender.gt.0) then
           !--rendering zoom does not allow pan - renderpt is always centre of axis
           renderpt = 0.5*(rendermin + rendermax)
           rendermin = renderpt - 0.5*renderlength
           rendermax = renderpt + 0.5*renderlength
           iadvance = 0
           iexit = .true.
        endif
     case('a') ! reset plot limits
        if (xpt.gt.xmax .and. irender.gt.0) then
           call useadaptive
           iadvance = 0              ! that it should change the render limits
           iexit = .true.
        else
           if (xpt.ge.xmin .and. xpt.le.xmax .and. ypt.le.ymax) then
              print*,'resetting x limits'
              xmin = minval(xcoords)
              xmax = maxval(xcoords)
              iadvance = 0
              iexit = .true.
           endif
           if (ypt.ge.ymin .and. ypt.le.ymax .and. xpt.le.xmax) then
              print*,'resetting y limits'
              ymin = minval(ycoords)
              ymax = maxval(ycoords)
              iadvance = 0
              iexit = .true.
           endif
        endif
     !
     !--zoom in/out on vector plots (arrow size)
     !
     case('v')
        if (ivecx.gt.0 .and. ivecy.gt.0) then
           print*,'decreasing vector arrow size'
           vecmax = 1.2*vecmax
           iadvance = 0
           iexit = .true.
        endif
     case('V')
        if (ivecx.gt.0 .and. ivecy.gt.0) then
           print*,'increasing vector arrow size'
           vecmax = 0.8*vecmax
           iadvance = 0
           iexit = .true.
        endif
     !
     !--set/unset log axes
     !
     case('l')
        !
        !--change colour bar, y and x itrans between log / not logged
        !
        if (xpt.gt.xmax .and. irender.gt.0) then
           call change_itrans(irender,rendermin,rendermax)
           iadvance = 0
           iexit = .true.
        elseif (xpt.lt.xmin) then
           if (iploty.le.ndim .and. irender.gt.0) then
              print "(a)",'error: cannot log coordinate axes with rendering'
           else
              call change_itrans(iploty,ymin,ymax)
              iadvance = 0
              iexit = .true.
           endif
        elseif (ypt.lt.ymin) then
           if (iplotx.le.ndim .and. irender.gt.0) then
              print "(a)",'error: cannot log coordinate axes with rendering'
           else
              call change_itrans(iplotx,xmin,xmax)
              iadvance = 0
              iexit = .true.
           endif
        endif
     !
     !--rotation
     !
     case(',')
        if (rotation) then
           print*,'changing z rotation angle by -15 degrees...'
           anglez = anglez - 15.
           iadvance = 0
           iexit = .true.
        endif
     case('<')
        if (rotation) then
           print*,'changing z rotation angle by -30 degrees...'
           anglez = anglez - 30.
           iadvance = 0
           iexit = .true.
        endif
     case('.')
        if (rotation) then
           print*,'changing z rotation angle by 15 degrees...'
           anglez = anglez + 15.
           iadvance = 0
           iexit = .true.
        endif
     case('>')
        if (rotation) then
           print*,'changing z rotation angle by 30 degrees...'
           anglez = anglez + 30.
           iadvance = 0
           iexit = .true.
        endif
     case('/')
        if (rotation .and. ndim.ge.2) then
           print*,'changing y rotation angle by -15 degrees...'
           angley = angley - 15.
           iadvance = 0
           iexit = .true.
        endif
     case('?')
        if (rotation .and. ndim.ge.2) then
           print*,'changing y rotation angle by -30 degrees...'
           angley = angley - 30.
           iadvance = 0
           iexit = .true.
        endif
     case('\')
        if (rotation .and. ndim.ge.2) then
           print*,'changing y rotation angle by 15 degrees...'
           angley = angley + 15.
           iadvance = 0
           iexit = .true.
        endif
     case('|')
        if (rotation .and. ndim.ge.2) then
           print*,'changing y rotation angle by 30 degrees...'
           angley = angley + 30.
           iadvance = 0
           iexit = .true.
        endif
     case('[')
        if (rotation .and. ndim.ge.3) then
           print*,'changing x rotation angle by -15 degrees...'
           anglex = anglex - 15.
           iadvance = 0
           iexit = .true.
        endif
     case('{')
        if (rotation .and. ndim.ge.3) then
           print*,'changing x rotation angle by -30 degrees...'
           anglex = anglex - 30.
           iadvance = 0
           iexit = .true.
        endif
     case(']')
        if (rotation .and. ndim.ge.3) then
           print*,'changing x rotation angle by 15 degrees...'
           anglex = anglex + 15.
           iadvance = 0
           iexit = .true.
        endif
     case('}')
        if (rotation .and. ndim.ge.3) then
           print*,'changing x rotation angle by 30 degrees...'
           anglex = anglex + 30.
           iadvance = 0
           iexit = .true.
        endif
     !
     !--set cross section position
     !
     case('x')
        if (rotation .and. ndim.ge.3) then
           xline(1) = xpt
           yline(1) = ypt
           !--work out which is the third dimension
           do i=1,3
              if (i.ne.iplotx .and. i.ne.iploty) ixsec = i
           enddo
           print*,' select cross section position (using left click or x)'
           call pgband(1,1,xline(1),yline(1),xline(2),yline(2),char2)
           !--work out cross section if left click or x again
           select case(char2)
           case('A','x')
              !--plot the cross section line
              call pgline(2,xline(1:2),yline(1:2))
              !--work out angle with the x axis
              !  and offset of line from origin
              dx = xline(2) - xline(1)
              dy = yline(2) - yline(1)
              anglerad = ATAN2(dy,dx)
              select case(ixsec)
              case(1)
                 anglex = 180.*anglerad/pi + anglex
                 print*,'setting angle x = ',anglex
              case(2)
                 angley = 180.*anglerad/pi + angley
                 print*,'setting angle y = ',angley
              case(3)
                 anglez = 180.*anglerad/pi + anglez
                 print*,'setting angle z = ',anglez
              end select
              iploty = ixsec
              !--work out offset of cross section line
              ! y intercept
              yint = yline(2) - (dy/dx)*xline(2)
              zslicepos = yint/COS(anglerad)
              print*,'iploty = ',ixsec, ' xsecpos = ',zslicepos
              iadvance = 0
              iexit = .true.
           case default
              print*,' action cancelled'
           end select
        endif
     !
     !--cross sections
     !
     case('u') ! move cross section up by dxsec
        if (iplotz.gt.0 .and. ndim.eq.3) then
           if (x_sec) then
              print*,'shifting cross section position up by ',dzslice
              zslicepos = zslicepos + dzslice
              iadvance = 0
              iexit = .true.
           else
              print*,'shifting perspective position up ',dscreen
              zobserver = zobserver + dscreen
              iadvance = 0
              iexit = .true.           
           endif
        endif
     case('U') ! move cross section up by 2*dxsec
        if (iplotz.gt.0 .and. ndim.eq.3) then
           if (x_sec) then
              print*,'shifting cross section position up by ',2.*dzslice
              zslicepos = zslicepos + 2.*dzslice
              iadvance = 0
              iexit = .true.
           else
              print*,'shifting perspective position up by ',2.*dscreen
              zobserver = zobserver + 2.*dscreen
              iadvance = 0
              iexit = .true.           
           endif
        endif
     case('d') ! move cross section down by dxsec
        if (iplotz.gt.0 .and. ndim.eq.3) then
           if (x_sec) then
              print*,'shifting cross section position down by ',dzslice
              zslicepos = zslicepos - dzslice
              iadvance = 0
              iexit = .true.
           else
              print*,'shifting perspective position down by ',dscreen
              zobserver = zobserver - dscreen
              iadvance = 0
              iexit = .true.           
           endif
        endif     
     case('D') ! move cross section down by 2*dxsec
        if (iplotz.gt.0 .and. ndim.eq.3) then
           if (x_sec) then
              print*,'shifting cross section position down by ',2.*dzslice
              zslicepos = zslicepos - 2.*dzslice
              iadvance = 0
              iexit = .true.
           else
              print*,'shifting perspective position down by ',2.*dscreen
              zobserver = zobserver - 2.*dscreen
              iadvance = 0
              iexit = .true.           
           endif
           iadvance = 0
           iexit = .true.
        endif
     !
     !--general plot stuff
     !
     case('G') ! move legend here
        print*,'setting legend position to current location...'
        call mvlegend(xpt,ypt,xmin,xmax,ymax)
     case('T') ! move title here
        print*,'setting title position to current location...'
        call mvtitle(xpt,ypt,xmin,xmax,ymax)
     case('H') ! move vector legend here
        if (ivecx.gt.0 .and. ivecy.gt.0) then
           print*,'setting vector plot legend to current location...'
           call mvlegendvec(xpt,ypt,xmin,xmax,ymax)
        endif
     case('m') ! change colour map (next scheme)
        call change_colourmap(icolourscheme,1)
     case('M') ! change colour map (previous scheme)
        call change_colourmap(icolourscheme,-1)
     case('i') ! invert colour map
        icolourscheme = -icolourscheme
        call change_colourmap(icolourscheme,0)
     !
     !--timestepping
     !
     case('q','Q')
        iadvance = -666
        print*,'quitting...'
        iexit = .true.
     case('X','b','B') ! right click -> go back
        iadvance = -abs(iadvance)
        iexit = .true.
     case('r','R') ! replot
        iadvance = 0
        iexit = .true.
     case(' ','n','N') ! space
        iadvance = abs(iadvance)
        iexit = .true.
     case('0','1','2','3','4','5','6','7','8','9')
        read(char,*,iostat=ierr) iadvance
        if (ierr /=0) then
           print*,'*** internal error setting timestep jump' 
           iadvance = 1
        endif
        print*,' setting timestep jump = ',iadvance
     case(')')
        iadvance = 10
        print*,' setting timestep jump = ',iadvance
     !
     !--unknown
     !
     case default
        print*,' x, y = ',xpt,ypt,'; unknown option "',trim(char), '"'
     end select

     if (rotation) then
        if (anglez.ge.360.) anglez = anglez - 360.
        if (anglez.lt.0.) anglez = anglez + 360.
        if (ndim.gt.2) then
           if (angley.ge.360.) angley = angley - 360.
           if (angley.lt.0.) angley = angley + 360.
           if (anglex.ge.360.) anglex = anglex - 360.
           if (anglex.lt.0.) anglex = anglex + 360.        
        endif
     endif
     !
     !--do not let timestep go outside of bounds
     !  if we are at the first/last step, just print message and do nothing
     !  if iadvance trips over the bounds, jump to last/first step
     !
     if (iadvance.ne.-666 .and. iexit) then
        if (istep + iadvance .gt. ilaststep) then
           print "(1x,a)",'reached last timestep'
           if (ilaststep-istep .gt.0) then
              iadvance= ilaststep - istep
           else
              iexit = .false.
           endif
        elseif (istep + iadvance .lt. 1) then
           print "(1x,a)",'reached first timestep: can''t go back'
           if (1-istep .lt.0) then
              iadvance= 1 - istep
           else
              iexit = .false.
           endif
        endif
     endif

  enddo interactiveloop
  return
end subroutine interactive_part

!
! cut down version of interactive mode -> controls timestepping only
! used in powerspectrum / extra plots
!
subroutine interactive_step(iadvance,istep,ilaststep,xmin,xmax,ymin,ymax)
 implicit none
 integer, intent(inout) :: iadvance
 integer, intent(in) :: istep,ilaststep
 real, intent(inout) :: xmin,xmax,ymin,ymax
 integer :: nc,ierr
 real :: xpt,ypt,xpt2,ypt2
 real :: xlength, ylength
 character(len=1) :: char,char2
 character(len=5) :: string
 logical :: iexit
 
  call pgqinf('CURSOR',string,nc)
  if (string(1:nc).eq.'YES') then
     print*,'entering interactive mode...press h in plot window for help'
  else
     print*,'cannot enter interactive mode: device has no cursor'
     return
  endif
  char = 'A'
  xpt = 0.
  ypt = 0.
  iexit = .false.
  
  do while (.not.iexit)
     call pgcurs(xpt,ypt,char)
     !
     !--exit if the device is not interactive
     !
     if (char.eq.achar(0)) return
  
     print*,'x, y = ',xpt,ypt,' function = ',char
     
     select case(char)
     case('h')
        print*,'-------------- interactive mode commands --------------'
        print*,' select area and zoom : left click (or A)'
        print*,' zoom in by 10%       : +'
        print*,' zoom out by 10(20)%      : - (_)'
        print*,' (r)eplot current plot        : r'
        print*,' next timestep/plot   : space, n'
        print*,' previous timestep    : right click (or X), b'
        print*,' jump forward (back) by n timesteps  : 0,1,2,3..9 then left (right) click'
        print*,' G : move legend to current position'
        print*,' T : move title to current position'
        print*,' (h)elp                       : h'
        print*,' (q)uit plotting              : q, Q'             
        print*,'-------------------------------------------------------'
     !
     !--zoom
     !
     case('A') ! left click
        !
        !--draw rectangle from the point and reset the limits
        !
        print*,'select area: '
        print*,'left click : zoom'
        call pgband(2,1,xpt,ypt,xpt2,ypt2,char2)
        print*,xpt,ypt,xpt2,ypt2,char2
        select case (char2)
        case('A')   ! zoom if another left click
           call pgrect(xpt,xpt2,ypt,ypt2)
           xmin = min(xpt,xpt2)
           xmax = max(xpt,xpt2)
           ymin = min(ypt,ypt2)
           ymax = max(ypt,ypt2)
           iadvance = 0
           iexit = .true.
        case default
           print*,' action cancelled'
        end select    
     case('-','_') ! zoom out by 10 or 20%
        print*,'zooming out'
        xlength = xmax - xmin
        ylength = ymax - ymin
        select case(char)
        case('-')
           xlength = 1.1*xlength
           ylength = 1.1*ylength
        case('_')
           xlength = 1.2*xlength
           ylength = 1.2*ylength
        end select
        xmin = xpt - 0.5*xlength
        xmax = xpt + 0.5*xlength
        ymin = ypt - 0.5*ylength
        ymax = ypt + 0.5*ylength
        iadvance = 0
        iexit = .true.
     case('+') ! zoom in by 10%
        print*,'zooming in'
        xlength = xmax - xmin
        ylength = ymax - ymin
        xlength = 0.9*xlength
        ylength = 0.9*ylength
        xmin = xpt - 0.5*xlength
        xmax = xpt + 0.5*xlength
        ymin = ypt - 0.5*ylength
        ymax = ypt + 0.5*ylength
        iadvance = 0
        iexit = .true.
     !
     !--general plot stuff
     !
     case('G') ! move legend here
        print*,'setting legend position to current location...'
        call mvlegend(xpt,ypt,xmin,xmax,ymax)
     case('T') ! move title here
        print*,'setting title position to current location...'
        call mvtitle(xpt,ypt,xmin,xmax,ymax)
     !
     !--timestepping
     !
     case('q','Q')
        iadvance = -666
        print*,'quitting...'
        iexit = .true.
     case('X','b','B') ! right click -> go back
        iadvance = -abs(iadvance)
        iexit = .true.
     case('r','R') ! replot
        iadvance = 0
        iexit = .true.
     case(' ','n','N') ! space
        iadvance = abs(iadvance)
        iexit = .true.
     case('0','1','2','3','4','5','6','7','8','9')
        read(char,*,iostat=ierr) iadvance
        if (ierr /=0) then
           print*,'*** error setting timestep jump' 
           iadvance = 1
        endif
        print*,' setting timestep jump = ',iadvance
     case(')')
        iadvance = 10
        print*,' setting timestep jump = ',iadvance
     end select
     !
     !--do not let timestep go outside of bounds
     !  if we are at the first/last step, just print message and do nothing
     !  if iadvance trips over the bounds, jump to last/first step
     !
     if (iadvance.ne.-666 .and. iexit) then
        if (istep + iadvance .gt. ilaststep) then
           print "(1x,a)",'reached last timestep'
           if (ilaststep-istep .gt.0) then
              iadvance= ilaststep - istep
           else
              iexit = .false.
           endif
        elseif (istep + iadvance .lt. 1) then
           print "(1x,a)",'reached first timestep: can''t go back'
           if (1-istep .lt.0) then
              iadvance= 1 - istep
           else
              iexit = .false.
           endif
        endif
     endif
  enddo
  return
end subroutine interactive_step

!-----------------------------------------------------------
! These subroutines interface to the actual plot settings
!-----------------------------------------------------------

!
!--move the legend to the current position
!
subroutine mvlegend(xi,yi,xmin,xmax,ymax)
 use settings_page, only:hposlegend,vposlegend
 implicit none
 real, intent(in) :: xi,yi,xmin,xmax,ymax
 real :: xch,ych
 
 hposlegend = (xi - xmin)/(xmax-xmin)
 !--query character height in world coordinates
 call pgqcs(4,xch,ych)
 vposlegend = (ymax - yi)/ych
 print*,'hpos = ',hposlegend,' vpos = ',vposlegend
 
 return
end subroutine mvlegend
!
!--move the vector legend to the current position
!
subroutine mvlegendvec(xi,yi,xmin,xmax,ymax)
 use settings_vecplot, only:hposlegendvec,vposlegendvec
 implicit none
 real, intent(in) :: xi,yi,xmin,xmax,ymax
 real :: xch,ych
 
 hposlegendvec = (xi - xmin)/(xmax-xmin)
 !--query character height in world coordinates
 call pgqcs(4,xch,ych)
 vposlegendvec = (ymax - yi)/ych
 print*,'hpos = ',hposlegendvec,' vpos = ',vposlegendvec
 
 return
end subroutine mvlegendvec
!
!--move the title to the current position
!
subroutine mvtitle(xi,yi,xmin,xmax,ymax)
 use settings_page, only:hpostitle,vpostitle
 implicit none
 real, intent(in) :: xi,yi,xmin,xmax,ymax
 real :: xch,ych
 
 hpostitle = (xi - xmin)/(xmax-xmin)
 !--query character height in world coordinates
 call pgqcs(4,xch,ych)
 vpostitle = (yi - ymax)/ych
 print*,'hpos = ',hpostitle,' vpos = ',vpostitle
 
 return
end subroutine mvtitle

!
!--saves current plot limits
!
subroutine save_limits(iplot,xmin,xmax)
 use limits, only:lim
 use multiplot, only:itrans
 use settings_data, only:ndim
 use settings_limits, only:iadapt,iadaptcoords
 use transforms, only:transform_limits_inverse
 implicit none
 integer, intent(in) :: iplot
 real, intent(in) :: xmin,xmax
 real :: xmintemp,xmaxtemp
 
 if (itrans(iplot).ne.0) then
    xmintemp = xmin
    xmaxtemp = xmax
    call transform_limits_inverse(xmintemp,xmaxtemp,itrans(iplot))
    lim(iplot,1) = xmintemp
    lim(iplot,2) = xmaxtemp
 else
    lim(iplot,1) = xmin
    lim(iplot,2) = xmax
 endif
 !
 !--change appropriate plot limits to fixed (not adaptive)
 !
 if (iplot.le.ndim) then
    iadaptcoords = .false.
 else
    iadapt = .false.
 endif
 
 return
end subroutine save_limits
!
!--saves current plot limits for particle tracking
!
subroutine save_limits_track(iplot,xmin,xmax,xi)
 use multiplot, only:itrans
 use settings_data, only:ndim
 use settings_limits, only:xminoffset_track,xmaxoffset_track
 use transforms, only:transform_limits_inverse
 implicit none
 integer, intent(in) :: iplot
 real, intent(in) :: xmin,xmax,xi
 real :: xmintemp,xmaxtemp
 
 if (iplot.gt.ndim) then
    print*,'ERROR in save_limits_track: iplot>ndim'
    return
 elseif (itrans(iplot).ne.0) then
    xmintemp = xmin
    xmaxtemp = xmax
    call transform_limits_inverse(xmintemp,xmaxtemp,itrans(iplot))
    xminoffset_track(iplot) = abs(xi - xmintemp)
    xminoffset_track(iplot) = abs(xmaxtemp - xi)
 else
    xminoffset_track(iplot) = abs(xi - xmin)
    xmaxoffset_track(iplot) = abs(xmax - xi)
 endif
 
 return
end subroutine save_limits_track
!
!--toggles log/unlog
!  note this only changes a pure log transform: will not change combinations
!
subroutine change_itrans(iplot,xmin,xmax)
 use multiplot, only:itrans
 use settings_data, only:numplot
 use transforms, only:transform_limits,transform_limits_inverse
 implicit none
 integer, intent(in) :: iplot
 real, intent(inout) :: xmin, xmax
 
 if (iplot.le.numplot) then
    if (itrans(iplot).eq.1) then
       itrans(iplot) = 0
       !!--untransform the plot limits
       call transform_limits_inverse(xmin,xmax,1)
    else
       itrans(iplot) = 1
       !!--transform the plot limits
       call transform_limits(xmin,xmax,itrans(iplot))
    endif
 endif
 
end subroutine change_itrans

subroutine useadaptive
 use settings_limits, only:iadapt
 implicit none
 
 iadapt = .true.
 
end subroutine useadaptive
!
!--saves rotation options
!
subroutine save_rotation(ndim,anglexi,angleyi,anglezi)
 use settings_xsecrot, only:anglex,angley,anglez
 implicit none
 integer, intent(in) :: ndim
 real, intent(in) :: anglexi,angleyi,anglezi
 
 anglez = anglezi
 if (ndim.ge.3) then
    anglex = anglexi
    angley = angleyi
 endif
 
 return
end subroutine save_rotation

!
!--saves cross section position
!
subroutine save_xsecpos(xsecpos)
 use settings_xsecrot, only:xsecpos_nomulti
 implicit none
 real, intent(in) :: xsecpos
 
 xsecpos_nomulti = xsecpos
 
 return
end subroutine save_xsecpos

!
!--saves 3D perspective
!
subroutine save_perspective(zpos,dz)
 use settings_xsecrot, only:zobserver,dzscreenfromobserver
 implicit none
 real, intent(in) :: zpos,dz
 
 zobserver = zpos
 dzscreenfromobserver = dz
 
 return
end subroutine save_perspective

!
!--saves circles of interaction
!
subroutine save_circles(ncircpartset,icircpartset)
 use settings_part, only:ncircpart,icircpart
 implicit none
 integer, intent(in) :: ncircpartset
 integer, intent(in), dimension(:) :: icircpartset
 integer :: imax
 
 imax = min(size(icircpartset),size(icircpart),ncircpartset)
 ncircpart = imax
 icircpart(1:imax) = icircpartset(1:imax)
 print*,'saving ',imax,' circles of interaction'
 
end subroutine save_circles
!
!--change colour map
!
subroutine change_colourmap(imap,istep)
 use colours, only:colour_set,ncolourschemes
 implicit none
 integer, intent(inout) :: imap
 integer, intent(in) :: istep
 
 imap = imap + istep
 if (abs(imap).gt.ncolourschemes) imap = 1
 if (abs(imap).lt.1) imap = ncolourschemes
 call colour_set(imap)
 
end subroutine change_colourmap

end module interactive_routines
