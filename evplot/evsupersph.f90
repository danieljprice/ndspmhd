!
!--visualisation tool to plot the data contained in the .ev file
!  (ie. the time evolution of SPH conserved quantities and errors)
!  uses PGPLOT subroutines
!
program plotmeagraph
  implicit none
  integer :: ncol
  integer, parameter :: maxfile=21
  integer, parameter :: maxstep=100000
  integer, parameter :: maxcol=21	! (6)21 (non)MHD	maximum number of columns
  integer i,iprev,nfiles,ifile,ifilesteps
  integer mysteps,ipick,ipickx,nacross,ndown
  integer ichange, iongraph, ipt
  integer iplotx(maxcol),iploty(maxcol), nplots
  integer multiplotx(maxcol),multiploty(maxcol),nplotsmulti
  integer nstepsfile(maxfile)
  real evdata(maxstep,maxcol,maxfile),evplot(maxstep)
  real lim(maxcol,2)
  real hpos,vpos, freqmin, freqmax
  character, dimension(maxfile) :: rootname*20, legendtext*120
  character*23 :: filename
  character*24 :: title,label(maxcol)
  character*40 :: labely
  character*1 ans
  logical :: icycle, igetfreq
  print*,' Welcome to Dan''s supersphplotev 2004... '

  mysteps=maxstep
  ncol = 21
  icycle = .false.
  hpos = 0.4
  vpos = 8.0
  iongraph = 1
  igetfreq = .false.
  !
  !--get filename(s)
  !
  iprev = 1
  i = 1
  do while (rootname(iprev)(1:1).ne.' ' .and. i.le.maxfile)
     call getarg(i,rootname(i))
     !print*,i,rootname(i)
     iprev = i
     i = i + 1
     if (i.gt.maxfile .and. rootname(iprev)(1:1).ne.' ') then
        print*,'WARNING: number of files >= array size: setting nfiles = ',maxfile
        !print*,'press return to continue'
        !read*
     endif
  enddo
  if (i.gt.maxfile .and. rootname(maxfile)(1:1).ne.' ') then
     nfiles = maxfile
  else
     nfiles = iprev - 1
  endif
  print*,'number of files = ',nfiles

  if (rootname(1)(1:1).eq.' ') then
     nfiles = 1
     print*,' Enter filename : '
     read*,rootname(1)
  endif
!
!--read data from all files
!  
  do ifile=1,nfiles
     if (index(rootname(ifile),'.ev').eq.0) then
        filename = trim(rootname(ifile))//'.ev'
     else
        filename = trim(rootname(ifile))
     endif
     print*,' opening ',filename
     ifilesteps = maxstep
     mysteps = 0
     call readev(ifilesteps,evdata(1:ifilesteps,1:ncol,ifile),ncol,filename)
     mysteps = max(mysteps,ifilesteps)
     nstepsfile(ifile) = ifilesteps
     !
     !--pad array
     !
     !if (ifilesteps.lt.mysteps) then
     !   evdata(ifilesteps+1:mysteps,1:ncol,ifile) = 0.
     !endif
  enddo
!
!--read legend text from the legend file, if it exists
!
  if (nfiles.gt.1) then
     open(unit=50,file='legend',status='old',ERR=22)
     print*,'reading labels from legend file...'
     do ifile=1,nfiles  
        read(50,"(a)",ERR=21,END=20) legendtext(ifile)
     enddo
     close(unit=50)
     goto 23
20   continue
     print*,'error: end of file in legend file'
     legendtext(ifile:nfiles) = rootname(ifile:nfiles)
     close(unit=50)
     goto 23
21   continue
     print*,'error reading from legend file'
     close(unit=50)
22   continue
     legendtext(1:nfiles) = rootname(1:nfiles)
23   continue

     open(unit=51,file='legendpos',status='old',ERR=33)
     read(51,*,ERR=31,END=31) hpos,vpos,iongraph
31   continue
     close(51)     
33   continue
   endif
!
!--set plot limits
!
  print*,'setting plot limits'
  do i=1,ncol
     lim(i,1) = 1.e12
     lim(i,2) = -1.e12
     do ifile = 1,nfiles
        lim(i,1) = min(lim(i,1),minval(evdata(1:nstepsfile(ifile),i,ifile)))
        lim(i,2) = max(lim(i,2),maxval(evdata(1:nstepsfile(ifile),i,ifile)))
     enddo
     if (lim(i,1).eq.lim(i,2)) then
        print*,' equal plot limits column ',i
        lim(i,2) = lim(i,2)*1.05 
        lim(i,1) = lim(i,1) - 0.05*lim(i,2)
        if (lim(i,2).eq.0.0) lim(i,2) = 1.0
     endif
     !!if (i.gt.1) lim(i,2) = lim(i,2)*1.1
  enddo
  lim(1,1) = 0.0

!
!--open file for frequency output
!
  if (igetfreq) then
     print*,'opening ',trim(rootname(1))//'.freq for output'
     open(unit=8,file=trim(rootname(1))//'.freq',status='replace')
  endif

!------------------------------------------      
!  menu
200 continue

  label = 'crap    '  ! default name

  label(1) = 'time           '  ! name for each column of data
  label(2) = 'E_kinetic      '
  label(3) = 'E_internal     '
  label(4) = 'E_magnetic     '
  label(5) = 'E_total        '
  label(6) = 'Linear momentum'

  if (ncol.gt.6) then
     label(7) = 'Total flux     '
     label(8) = 'Cross helicity '
     label(9) = 'Plasma beta (min)'
     label(10) = 'Plasma beta (ave)'
     label(11) = 'Plasma beta (max)'
     label(12) = 'div B (average)  '
     label(13) = 'div B (maximum)  '
     label(14) = 'int (div B) dV   '
     label(15) = 'Fmag dot B/|Fmag| (av)'
     label(16) = 'Fmag dot B/|Fmag| (max)'
     label(17) = 'Fmag dot B/|force| (av)'
     label(18) = 'Fmag dot B/|force| (max)'            
     label(19) = 'omega_mhd (average)'
     label(20) = 'omega_mhd (max) '
     label(21) = '% particles with omega < 0.01 '
  endif

  title = ' '
  
  print*,' You may choose from a delectable sample of plots '
  print 12
  do i=1,ncol
     print 11,i,label(i)
  enddo
  print 12
  print 13,ncol+1,'Multiplot    '
  print 13,ncol+2,'Set multiplot    '
  print 12
  print 13,ncol+3,'Read new data'
  print 13,ncol+4,'Change number of timesteps read'
  print 13,ncol+5,'Get frequencies'
  print 13,ncol+6,'Adjust legend position'
  print 13,ncol+7,'Adjust plot limits'
  print 13,ncol+8,'Do auto update'          
  print 13,ncol+9,'Exit supersphplotev'
  print 12
  print*,'Enter y axis or selection: '
11 format(1x,i2,')',1x,a)
12 format(1x,45('-'))
13 format(1x,i2,')',1x,a)
14 format(1x,i2,')',1x,a1,a1,' projection')
  read *,ipick
  if (ipick.le.ncol .and. ipick.ge.1) then
     print*,'Enter x axis: '
     read*,ipickx
  endif
  if ((ipick.ge.ncol+8).or.(ipick.lt.1)) then
     close(8)
     stop
  endif

  if (ipick.eq.ncol+1) then
     if (nplotsmulti.eq.0) ipick = ncol+2
     nplots = nplotsmulti      
     iploty(1:nplots) = multiploty(1:nplots)
     iplotx(1:nplots) = multiplotx(1:nplots)
  else
     iploty(1) = ipick
     iplotx(1) = ipickx
     nplots = 1      
  endif

  if (ipick.eq.ncol+2) then
     print*,' Enter number of plots '
     read*,nplotsmulti
     nplots = nplotsmulti
     do i=1,nplotsmulti
        print*,' Enter y plot ',i
        read*,multiploty(i)
     enddo
     print*,' Enter x plot '
     read*,ipickx
     multiplotx(1:nplotsmulti) = ipickx
     goto 200
  elseif (ipick.eq.ncol+3) then
     nfiles = 1
     print*,'Enter new rootname:'
     read*,rootname(1)
     filename = trim(rootname(1))//'.ev'
     print*,'reading evolution file ',filename    	  
     call readev(mysteps,evdata(1:mysteps,1:ncol,1),ncol,filename)
     print*,'setting plot limits'
     do i=1,ncol
        lim(i,1) = minval(evdata(1:mysteps,i,1))
        lim(i,2) = maxval(evdata(1:mysteps,i,1))*1.05
     enddo
     lim(1,1) = 0.0
     goto 200
  elseif (ipick.eq.ncol+4) then
     print*,'Enter number of timesteps to read:'
     read*,mysteps
     print *,' Steps = ',mysteps
     goto 200
  elseif (ipick.eq.ncol+5) then
     igetfreq = .not.igetfreq
     print*,'frequencies = ',igetfreq
     goto 200
  elseif (ipick.eq.ncol+6) then
     print*,'Enter horizontal position as fraction of x axis:'
     read*,hpos
     print*,'Enter vertical offset from top of graph in character heights:'
     read*,vpos
     print*,'Enter number of plot on page to place legend'
     read*,iongraph
     open(unit=51,file='legendpos',status='replace',ERR=43)
         write(51,*) hpos,vpos,iongraph
     close(51)	 
43   continue     
     goto 200
  elseif (ipick.eq.ncol+7) then
     print*,'Enter plot number to change limits'
     read*,ichange
     if (ichange.gt.0 .and. ichange.le.ncol) then
        print*,' Enter ',trim(label(ichange)),' min:'
	read*,lim(ichange,1)
	print*,' Enter ',trim(label(ichange)),' max:'
        read*,lim(ichange,2)
     endif
     goto 200
  elseif (ipick.eq.ncol+8) then
     icycle = .not.icycle
     print *,' cycle = ',icycle
     goto 200	  
  endif
! ------------------------------------------------------------------------
! compute delta x/x (single files only)

  if (nfiles.eq.1) then
     evplot(1:2) = 0.
     do i=2,mysteps
        if (evdata(i,iploty(1),1).ne.0.0) then
           evplot(i) = abs(evdata(i,iploty(1),1)-evdata(i-1,iploty(1),1))  &
                /evdata(i,iploty(1),1)
        else
           evplot(i) = 0.0
        endif
     enddo
  endif

! ------------------------------------------------------------------------
! initialise PGPLOT
!
  if (nplots.eq.1) then
     if (nfiles.gt.1) then
        call PGBEGIN(0,'?',1,1)
     else
        call PGBEGIN(0,'?',1,2)
     endif
  else
     nacross = 1
     ndown = nplots/nacross
     !ndown = nplots/2
     !nacross = nplots/ndown
     if (ndown.eq.1 .and. nacross.gt.1) then
        call PGBEGIN(0,'?',nacross,ndown)
     else
        call PGBEGIN(0,'?',1,1) !uses pgtile
     endif
     if (nacross.eq.2 .and. ndown.eq.1) then
        call pgpap(11.7,0.5/sqrt(2.))
     endif
  endif
  !!if (nplots.gt.2) call PGSCH(2.0)
!
!--open frequency file for output
!
  if (igetfreq) then
     print*,'opening ',trim(rootname(1))//'.freq for output'
     open(unit=8,file=trim(rootname(1))//'.freq',status='replace')
  endif

!
!--plot graphs (icycle only sets up at the moment)
!
  if (icycle) then
     do i=1,nplots
        print*,i,' x = ',iplotx(i),' y = ',iploty(i)
        call PGENV(lim(iplotx(i),1),lim(iplotx(i),2)*2,  &
             lim(iploty(i),1)-2.*lim(iploty(i),1),   &
             1.5*lim(iploty(i),2),0,1)	    
        call PGLABEL(label(iplotx(i)),label(iploty(i)),title)
     enddo
  else
     do i=1,nplots
	!
        !--setup plotting page
        !
        if (all(iplotx(1:nplots).eq.iplotx(1)).and.ndown.gt.1 ) then
           
	   call pgsls(1)
	   call pgslw(3)
	   call danpgtile(i,nacross,ndown,  &
	       lim(iplotx(i),1),lim(iplotx(i),2),lim(iploty(i),1),lim(iploty(i),2), &
               label(iplotx(i)),label(iploty(i)),title,0,0)
        else
           !
           !--for one plot just use standard PGENV routine
           !
           call pgsls(1)
           if (nplots.gt.1) call pgsch(1.5) ! increase character height
           call pgpage
           call pgsvp(0.2,0.99,0.2,0.99)
           call pgswin(lim(iplotx(i),1),lim(iplotx(i),2), &
                lim(iploty(i),1),lim(iploty(i),2))
           call pgbox('bcnst',0.0,0,'1bvcnst',0.0,0)
           !call PGENV(lim(iplotx(i),1),lim(iplotx(i),2), &
           !     lim(iploty(i),1),lim(iploty(i),2),0,1)
            call pgmtxt('l',3.5,0.5,0.5,label(iploty(i)))
            call pglabel(label(iplotx(i)),' ',' ')
            !!call pglabel(label(iplotx(i)),label(iploty(i)),title)
        endif
        !
        !--draw the line
        !
        call pgslw(3)
        do ifile=1,nfiles
           call PGSLS(MOD(ifile-1,5)+1)  ! change line style between plots
           if (i.eq.iongraph .and. nfiles.gt.1) call legend(ifile,legendtext(ifile),hpos,vpos)
           !call PGSCI(ifile) ! or change line colour between plots
           
	   call PGLINE(nstepsfile(ifile),evdata(1:nstepsfile(ifile),iplotx(i),ifile), &
                evdata(1:nstepsfile(ifile),iploty(i),ifile))
           !
           !--if more than 5 files, plot points on the line to make a new
           !  line style
	   !
	   if (ifile.gt.5) then
	      call pgsls(1)
	      do ipt=1,nstepsfile(ifile)
	         if (mod(ipt,10).eq.0) then
	            call PGPT(1,evdata(ipt,iplotx(i),ifile), &
                         evdata(ipt,iploty(i),ifile),mod(ifile,5) + 1)
	         endif
	      enddo
	   endif
	   
	   !
           !--work out period of oscillation from spacing of minima/maxima
	   !
           if (igetfreq .and. iplotx(i).eq.1) then
	      print*,rootname(ifile)
	      call getmax(evdata(1:nstepsfile(ifile),iploty(i),ifile), &
	                  evdata(1:nstepsfile(ifile),iplotx(i),ifile), &
		          nstepsfile(ifile),freqmax,freqmin)
	   
 	   !
           !--output this to a file
           !
	      print*,'writing to frequency file'
	      write(8,*) freqmin,freqmax,rootname(ifile)
	   endif
	   
	   
	   
	enddo
        call pgslw(1)

     enddo
  endif
!
!--close frequency file
!
  if (igetfreq) close(unit=8)
!
!--plot graph(s)
!
  ans = ' '
  if (icycle) then
     do while(ans.ne.'q')
        mysteps = maxstep  
        do ifile=1,nfiles 
           if (index(rootname(ifile),'.ev').eq.0) then
              filename = trim(rootname(ifile))//'.ev'
           else
              filename = trim(rootname(ifile))
           endif
           print*,'reading file ',filename
           call readev(nstepsfile(ifile),evdata(1:nstepsfile(ifile),1:ncol,ifile),ncol,filename)
        enddo

        do i=1,nplots
           call PGPANL(1,i)
           do ifile=1,nfiles
              call PGSLS(MOD(ifile-1,5)+1)
              call PGLINE(nstepsfile(ifile),evdata(1:nstepsfile(ifile),iplotx(i),ifile), &
                   evdata(1:nstepsfile(ifile),iploty(i),ifile))
           enddo
        enddo
        
        print*,'type q to quit, any to update'    	
        read*,ans
     enddo
  endif
!
!--plot error between timesteps
!	    
  if (nplots.eq.1 .and. nfiles.eq.1) then
     call PGENV(lim(iplotx(1),1),lim(iplotx(1),2), &
          minval(evplot),maxval(evplot),0,1)
     labely = '| delta '//TRIM(label(iploty(1)))//'|/'//label(iploty(1))
     call PGLABEL(label(iplotx(1)),labely, title)
     
     call PGLINE(mysteps-1,evdata(2:mysteps,iplotx(1),1),evplot(2:mysteps))    	    
     
  endif
  call PGEND
  goto 200
  
! ------------------------------------------------------------------------      
  
999 continue  
    close(8)               
end program plotmeagraph

subroutine readev(nsteps,evdata,ncols,rootname)
  implicit none
  integer :: i
  integer, intent(in) :: ncols
  integer, intent(inout) :: nsteps
  real, dimension(nsteps,ncols), intent(out) :: evdata
  character*20, intent(in) :: rootname
  character evname*20
  logical iexist
  
  evname = rootname
  inquire (file=evname, exist=iexist)
  if (.not.iexist) then
     print*,evname,': file does not exist, exiting'
     stop
  endif
  open(unit=11,file=evname,status='old',form='formatted')
  do i=1,nsteps
     read(11,*,end=98) evdata(i,1:ncols)
     !	 print*,'t = ',evdata(i,1)
  enddo
  close(unit=11)
  
  print*,' WARNING: number of steps > array limits'
  goto 99 
98 continue
  close(unit=11)
  write(*,"(a)",advance="no") '>> end of file:'
  if (i.gt.1) then
     print*,' t=',evdata(i-1,1),' nsteps = ',i-1 
  else
     stop ' file empty!! no timesteps'
  endif
  nsteps = i-1     
99 continue
  
  return
end subroutine readev

!
! subroutine to work out the positions of maxima and minima in the
!  quantity plotted. Works out the period from the distance between two
!  minima or two maxima
!
!
subroutine getmax(datain,time,max,freq,freq2)
  implicit none
  integer i,max,nmax,nmin
  real datain(max),time(max),timeprev,timeprev2
  real freq,freq2,avper,avper2,period,omega,omega2,pi
  
  pi = 3.1415926536
  
  timeprev = 0.
  timeprev2 = 0.
  avper = 0.
  avper2 = 0.
  nmax = 0
  nmin = 0
  
  do i=2,max-1
     if ((datain(i).gt.datain(i-1)).and.(datain(i).gt.datain(i+1))) then
        
        nmax = nmax + 1
        if (mod(nmax,2).eq.0) then
           period = time(i)-timeprev            
           print 10,'max',time(i),period,1./period
           timeprev = time(i)
           avper = avper + period
        endif
            
     elseif ((datain(i).lt.datain(i-1)).and.   &
          (datain(i).lt.datain(i+1))) then
        
        nmin = nmin + 1
        if (mod(nmin,2).eq.1) then
           if (nmin.gt.1) then
              period = time(i)-timeprev2            
              print 10,'min',time(i),period,1./period
              avper2 = avper2 + period
           endif
           timeprev2 = time(i)	       
        endif
        
     endif
  enddo
10 format(a,' at t = ',f8.4,' period = ',f8.4,' freq = ',f8.4)
  
  nmin = nmin - 1
  avper = avper/(nmax/2)
  avper2 = avper2/(nmin/2)
  freq = 1./avper
  freq2 = 1./avper2
  omega = 2.*pi*freq
  omega2 = 2.*pi*freq2
  print 20,'average max ',avper,freq,omega,nmax
  print 20,'average min ',avper2,freq2,omega2,nmin
20 format(a,' period = ',f8.4,' freq = ',f8.4,' omega = ',f8.4,' using ',i3,' values')  
  
  return
end subroutine getmax

!
!--draw a legend for different line styles
!  uses current line style and colour
!
subroutine legend(icall,text,hpos,vposin)
  implicit none
  integer, intent(in) :: icall
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

  call pgline(2,xline,yline)            ! draw line segment
!!--make up line style if > 5 calls (must match actual line drawn)
  if (icall.gt.5) then
     call pgpt(2,xline,yline,mod(icall,5)+1)
     call pgpt1(1,0.5*(xline(1)+xline(2)),yline(1),mod(icall,5)+1)
  endif  
  call PGTEXT(xline(2) + 0.5*xch,yline(1)-0.25*ych,trim(text))

!!  call pgmtxt('T',vpos,0.1,0.0,trim(text))  ! write text

  call pgstbg(-1) ! reset text background to transparent

end subroutine legend
