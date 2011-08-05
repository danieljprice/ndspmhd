!
!--visualisation tool to plot the data contained in the .ev file
!  (ie. the time evolution of SPH conserved quantities and errors)
!  uses PGPLOT subroutines
!
program plotmeagraph
  use transforms
  use prompting
  use system_commands
  implicit none
  integer :: ncol
  integer, parameter :: maxfile=25
  integer, parameter :: maxstep=6000
  integer, parameter :: maxcol=22        ! (6)21 (non)MHD        maximum number of columns
  integer :: i,nfiles,ifile,ifilesteps
  integer :: mysteps,ipick,ipickx,nacross,ndown
  integer :: ihalf,iadjust,ierr
  integer :: ichange, iongraph, ipt
  integer :: nplots
  integer, dimension(maxcol) :: iplotx,iploty,itrans
  integer, dimension(maxcol) :: multiplotx,multiploty
  integer :: nplotsmulti
  integer, dimension(maxfile) :: nstepsfile
  real, dimension(maxstep,maxcol,maxfile) :: evdata
  real, dimension(maxstep) :: evplot,xplot,yplot
  real :: lim(maxcol,2)
  real :: xmin, xmax, ymin, ymax
  real :: hpos,vpos, freqmin, freqmax
  character, dimension(maxfile) :: rootname*120, legendtext*120
  character(len=120) :: title,label(maxcol)
  character(len=40) :: labely,text
  character(len=1) :: ans
  character(len=2) :: ioption
  logical :: icycle, igetfreq, isameXaxis, isameYaxis, ishowopts

  print*,' Welcome to Dan''s supersphplotev 2005... '

  mysteps=maxstep
  ncol = maxcol
  icycle = .false.
  hpos = 0.4
  vpos = 8.0
  iongraph = 1
  igetfreq = .false.
  ishowopts = .false.
  ipickx = 1
  do i=1,maxcol
     multiploty(i) = i+1
  enddo
  itrans = 0
  nplotsmulti = 1
  multiplotx(:) = 2
  
  !
  !--get filename(s)
  !
  i = 1
  nfiles = 0
  call get_number_arguments(nfiles)
  if (nfiles.gt.maxfile) then
     print "(a,i4)",'WARNING: number of files >= array size: setting nfiles = ',maxfile
     nfiles = maxfile
  endif
  do i=1,nfiles
     call get_argument(i,rootname(i))
  enddo
  print*,' number of files = ',nfiles

  if (nfiles.lt.1 .or. rootname(1)(1:1).eq.' ') then
     nfiles = 1
     call prompt(' Enter number of files to read : ',nfiles,1,maxfile)
     do i=1,nfiles
        call prompt(' Enter filename : ',rootname(i))
     enddo
  endif
!
!--read data from all files
!  
  do ifile=1,nfiles
     ifilesteps = maxstep
     mysteps = 0
     call readev(ifilesteps,evdata(1:ifilesteps,1:maxcol,ifile), &
                 maxcol,ncol,rootname(ifile))
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
  call read_legend_file
!
!--read defaults file
!
   open(unit=51,file='.evsupersph_defaults',status='old',iostat=ierr)
   if (ierr==0) then
      read(51,*,iostat=ierr) hpos,vpos,iongraph,igetfreq,nplotsmulti
      read(51,*,iostat=ierr) itrans(1:maxcol)
      read(51,*,iostat=ierr) multiplotx(1:maxcol)
      read(51,*,iostat=ierr) multiploty(1:maxcol)      
      if (ierr /= 0) then
         print*,'error reading default options from file'
      else 
         print*,'read default options from file'
      endif
      close(51)
   else
      print*,'defaults file not found'
   endif
!
!--set plot limits
!
  print*,'setting plot limits, ncol = ',ncol
  do i=1,ncol
     lim(i,1) = 1.e12
     lim(i,2) = -1.e12
     do ifile = 1,nfiles
        !!print*,ifile,nstepsfile(ifile),size(nstepsfile)
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

  title = ' '

!------------------------------------------      
!  menu
  menuloop: do

  !
  !--set column labels
  !
  label = 'crap    '  ! default name

  label(1) = 'time           '  ! name for each column of data
  label(2) = 'E_kinetic      '
  label(3) = 'E_internal     '
  label(4) = 'E_magnetic     '
  label(5) = 'E_pot          '
  label(6) = 'E_total        '
  label(7) = 'Linear momentum'

  if (ncol.gt.7) then
     label(8) = 'Total flux     '
     label(9) = 'Cross helicity '
     label(10) = 'Plasma beta (min)'
     label(11) = 'Plasma beta (ave)'
     label(12) = 'Plasma beta (max)'
     label(13) = 'div B (average)  '
     label(14) = 'div B (maximum)  '
     label(15) = 'int (div B) dV   '
     label(16) = 'Fmag dot B/|Fmag| (av)'
     label(17) = 'Fmag dot B/|Fmag| (max)'
     label(18) = 'Fmag dot B/|force| (av)'
     label(19) = 'Fmag dot B/|force| (max)'            
     label(20) = 'omega_mhd (average)'
     label(21) = 'omega_mhd (max) '
     label(22) = '% particles with omega < 0.01 '
  endif
  label(ncol) = 'nstep          '

  !
  !--overwrite column labels if columns file exists
  !
  call read_columns_file(ncol,label(1:ncol))

  if (ANY(itrans.ne.0)) then
     do i=1,ncol
        label(i) = transform_label(label(i),itrans(i))
     enddo
  endif
!
!--print data menu (in two columns)
!  
  print*,' You may choose from a delectable sample of plots '
  print 12
  ihalf = ncol/2                ! print in two columns
  iadjust = mod(ncol,2)
  print 11, (i,label(i)(1:20),ihalf+i+iadjust,label(ihalf+i+iadjust)(1:20),i=1,ihalf)
  if (iadjust.ne.0) then
     print 13, ihalf + iadjust,trim(label(ihalf+iadjust))
  endif
!
!--print menu options
!
  print 12
  print 18,ncol+1,'multiplot           ','m','set multiplot'
  print 12
  if (ishowopts) then
     print 14,'d','Read new data'
     print 14,'t','Change number of timesteps read'
     print 14,'f','Get frequencies'
     print 14,'g','Adjust legend / title options'
     print 14,'l','Adjust plot limits'
     print 14,'p','Set plot title'
     print 14,'a','Do auto update'          
     print 14,'h','toggle help'
     print 14,'s','Save settings'
     print 14,'q','Exit supersphplotev'
  else
     print*,' d(ata) t(imesteps) f(requencies) l(imits)'
     print*,' le(g)end a(uto update) h(elp) s(ave) q(uit)'
  endif
  print 12
11 format(1x,i2,')',1x,a20,1x,i2,')',1x,a)
12 format(1x,65('-'))
13 format(1x,i2,')',1x,a)
14 format(1x,a2,')',1x,a)
18 format(1x,i2,')',1x,a20,1x,a1,')',1x,a)
!
!--read option
!
  promptloop: do
  
  write(*,"(a)",ADVANCE="NO") 'Enter y axis or option: '
  read*,ioption
  read(ioption,*,iostat=ierr) ipick
  if (ierr /=0) ipick = 0

  if (ipick.le.ncol .and. ipick.ge.1) then
     call prompt('Enter x axis: ',ipickx,1,ncol)
  endif
  if ((ipick.gt.ncol+1).or.(ipick.lt.0)) then
     exit menuloop
  elseif (ipick.eq.ncol+1) then
     if (nplotsmulti.eq.0) then
        ipick = 0
        ioption = 'm'
        nplotsmulti = 1
     endif
     nplots = nplotsmulti      
     iploty(1:nplots) = multiploty(1:nplots)
     iplotx(1:nplots) = multiplotx(1:nplots)
  else
     iploty(1) = ipick
     iplotx(1) = ipickx
     nplots = 1      
  endif
!
!--menu options
!
  if (ipick.eq.0) then
  
  select case(ioption(1:1))
  case('m','M')
     call prompt('Enter number of plots ',nplotsmulti,1,maxcol)
     nplots = nplotsmulti
     do i=1,nplotsmulti
        write(text,"(a,i2)") 'Enter y plot ',i
        call prompt(text,multiploty(i),1,maxcol)
     enddo
     call prompt('Enter x plot ',ipickx,1,maxcol)
     multiplotx(1:nplotsmulti) = ipickx
     cycle menuloop
  case('d','D')
     nfiles = 1
     call prompt('Enter new rootname:',rootname(1))
     mysteps = maxstep           
     call readev(mysteps,evdata(1:mysteps,1:maxcol,1),maxcol,ncol,rootname(1))
     print*,'setting plot limits'
     do i=1,ncol
        lim(i,1) = minval(evdata(1:mysteps,i,1))
        lim(i,2) = maxval(evdata(1:mysteps,i,1))*1.05
     enddo
     lim(1,1) = 0.0
     cycle menuloop
  case('t','T')
     call prompt('Enter number of timesteps to read:',mysteps,1)
     print *,' Steps = ',mysteps
     cycle menuloop
  case('f','F')
     igetfreq = .not.igetfreq
     print*,'frequencies = ',igetfreq
     cycle menuloop
  case('g','G')
     call prompt('Enter horizontal position as fraction of x axis:',hpos,0.0,1.0)
     call prompt('Enter vertical offset from top of graph in character heights:',vpos)
     call prompt('Enter number of plot on page to place legend ',iongraph,1)
     call prompt('Enter plot title ',title)
     cycle menuloop
  case('l','L')
     ichange = 0
     call prompt('Enter plot number to change limits',ichange,0,ncol)
     if (ichange.gt.0) then
        call prompt(' Enter '//trim(label(ichange))//' min:',lim(ichange,1))
        call prompt(' Enter '//trim(label(ichange))//' max:',lim(ichange,2))
        call prompt(' Enter transformation (1=log,2=abs,3=sqrt,4=1/x):',itrans(ichange),0,4)
     endif
     cycle menuloop
!  case('r','R')
!     ichange = 0
!     call prompt('Enter plot number to apply transformation',ichange,0,ncol)
     
  case('s','S')
     open(unit=51,file='.evsupersph_defaults',status='replace',iostat=ierr)
     if (ierr /= 0) then
        print*,'*** error opening defaults file *** '
     else
        write(51,*,iostat=ierr) hpos,vpos,iongraph,igetfreq,nplotsmulti
        write(51,*,iostat=ierr) itrans(1:maxcol)
        write(51,*,iostat=ierr) multiplotx(1:maxcol)
        write(51,*,iostat=ierr) multiploty(1:maxcol)
        if (ierr /= 0) print*,'*** error writing to defaults file ***'
        close(51)
     endif
     cycle menuloop
  case('a','A')
     icycle = .not.icycle
     print *,' cycle = ',icycle
     cycle menuloop        
  case('h','H')
     ishowopts = .not.ishowopts
     cycle menuloop  
  case('q','Q')
     exit menuloop
  case default
     cycle promptloop
  end select
  
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
  isameXaxis = all(iplotx(1:nplots).eq.iplotx(1))
  isameYaxis = all(iploty(1:nplots).eq.iploty(1))

  if (nplots.eq.1) then
     if (nfiles.gt.1) then
        call PGBEGIN(0,'?',1,1)
        call PGSCH(1.0)
        nacross = 1
        ndown = 1
     else
        call PGBEGIN(0,'?',1,2)
        call PGSCH(2.0)
        nacross = 1
        ndown = 2
        isameYaxis = .false.
     endif
  else
     nacross = 1
     ndown = nplots/nacross
     !ndown = nplots/2
     !nacross = nplots/ndown
     if (isameXaxis .and. isameYaxis) then
        call PGBEGIN(0,'?',1,1) ! uses pgtile
        call PGSCH(1.0)
     else
        call PGBEGIN(0,'?',nacross,ndown)
        call PGSCH(2.0)
     endif
     if (nacross.eq.2 .and. ndown.eq.1) then
        call pgpap(11.7,0.5/sqrt(2.))
     endif
     
  endif
  call pgslw(3)
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
        !--apply transformations to plot limits
        !
        xmin = lim(iplotx(i),1)
        xmax = lim(iplotx(i),2)
        ymin = lim(iploty(i),1)
        ymax = lim(iploty(i),2)
        call transform_limits(xmin,xmax,itrans(iplotx(i)))
        call transform_limits(ymin,ymax,itrans(iploty(i)))
        !
        !--setup plotting page
        !
        if (nplots.ge.1 .and. isameXaxis .and. isameYaxis) then
           !
           !--tiled plots
           !
           call pgsls(1)
           call pgsci(1)
           call danpgtile(i,nacross,ndown,  &
               xmin,xmax,ymin,ymax, &
               label(iplotx(i)),label(iploty(i)),' ',0,0)
           !
           ! plot the title inside the plot boundaries
           !
           call pgmtxt('T',-1.5,0.2,1.0,title)

        else
           !
           !--non-tiled plots
           !
           call pgsls(1)
           call pgsci(1)
           if (nplots.gt.1) call pgsch(1.5) ! increase character height
           call pgpage
           call pgsvp(0.2,0.99,0.2,0.99)
           call pgswin(xmin,xmax,ymin,ymax)
           call pgbox('bcnst',0.0,0,'1bvcnst',0.0,0)
           !call PGENV(lim(iplotx(i),1),lim(iplotx(i),2), &
           !     lim(iploty(i),1),lim(iploty(i),2),0,1)
            call pgmtxt('l',4.5,0.5,0.5,label(iploty(i)))
            call pglabel(label(iplotx(i)),' ',' ')
            !!call pglabel(label(iplotx(i)),label(iploty(i)),title)
        endif
        !
        !--draw the line
        !
        do ifile=1,nfiles
           call pgsls(MOD(ifile-1,5)+1)  ! change line style between plots
           call pgsci(ifile+1)
           if (i.eq.iongraph .and. nfiles.gt.1) call legend(ifile,legendtext(ifile),hpos,vpos)
           !call PGSCI(ifile) ! or change line colour between plots
           
           !
           !--apply transformations to plot data
           !
           xplot(1:nstepsfile(ifile)) = evdata(1:nstepsfile(ifile),iplotx(i),ifile)
           yplot(1:nstepsfile(ifile)) = evdata(1:nstepsfile(ifile),iploty(i),ifile)

           call transform(xplot(1:nstepsfile(ifile)),itrans(iplotx(i)))
           call transform(yplot(1:nstepsfile(ifile)),itrans(iploty(i)))
           
           call PGLINE(nstepsfile(ifile),xplot(1:nstepsfile(ifile)), &
                yplot(1:nstepsfile(ifile)))
           !
           !--if more than 5 files, plot points on the line to make a new
           !  line style
           !
           if (ifile.gt.5) then
              call pgsls(1)
              do ipt=1,nstepsfile(ifile)
                 if (mod(ipt,10).eq.0) then
                    call PGPT(1,xplot(ipt),yplot(ipt),mod(ifile,5) + 1)
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
              write(8,*) freqmin,freqmax,trim(rootname(ifile))
           endif
           
           
           
        enddo

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
        do ifile=1,nfiles 
           mysteps = 0
           ifilesteps = maxstep
           call readev(ifilesteps,evdata(1:ifilesteps,1:maxcol,ifile), &
                maxcol,ncol,rootname(ifile))
           nstepsfile(ifile) = ifilesteps
        enddo

        do i=1,nplots
           do ifile=1,nfiles
              call pgsls(MOD(ifile-1,5)+1)
              call pgsci(ifile+1)
           !
           !--apply transformations to plot data
           !
           xplot(1:nstepsfile(ifile)) = evdata(1:nstepsfile(ifile),iplotx(i),ifile)
           yplot(1:nstepsfile(ifile)) = evdata(1:nstepsfile(ifile),iploty(i),ifile)

           call transform(xplot(1:nstepsfile(ifile)),itrans(iplotx(i)))
           call transform(yplot(1:nstepsfile(ifile)),itrans(iploty(i)))
           
           call PGLINE(nstepsfile(ifile),xplot(1:nstepsfile(ifile)), &
                yplot(1:nstepsfile(ifile)))

!              call PGLINE(nstepsfile(ifile),evdata(1:nstepsfile(ifile),iplotx(i),ifile), &
!                   evdata(1:nstepsfile(ifile),iploty(i),ifile))
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
     call pgsls(1)
     call pgsci(1)
     call PGENV(lim(iplotx(1),1),lim(iplotx(1),2), &
          minval(evplot),maxval(evplot),0,1)
     labely = '| delta '//TRIM(label(iploty(1)))//'|/'//label(iploty(1))
     call PGLABEL(label(iplotx(1)),labely, title)
     
     call PGLINE(mysteps-1,evdata(2:mysteps,iplotx(1),1),evplot(2:mysteps))                
     
  endif
  call PGEND
  
  cycle menuloop
  
  enddo promptloop
  enddo menuloop
  
! ------------------------------------------------------------------------      
  
  close(8)  
  
contains

!-------------------------------------------------------
! reads legend text from the legend file, if it exists
!-------------------------------------------------------

 subroutine read_legend_file
  implicit none

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
     print*,'legend file not found'
     legendtext(1:nfiles) = rootname(1:nfiles)
23   continue
   endif
 end subroutine read_legend_file



end program plotmeagraph


!-------------------------------------------------------
! reads column labels from columns file, if it exists
!-------------------------------------------------------

 subroutine read_columns_file(ncolumns,label)
  implicit none
  integer, intent(in) :: ncolumns
  character(len=*), dimension(ncolumns), intent(inout) :: label 
  integer :: icol,i

  if (ncolumns.gt.1) then
     open(unit=51,file='columns',status='old',ERR=22)
     print*,'reading column labels from columns file...'
     do icol=1,ncolumns
        read(51,"(a)",ERR=21,END=20) label(icol)
     enddo
     close(unit=51)
     goto 23
20   continue
     print*,'error: end of file in columns file'
     do i=icol,ncolumns
        write(label(i),*) 'column',i
     enddo
     close(unit=51)
     goto 23
21   continue
     print*,'error reading from columns file'
     close(unit=51)
22   continue
     print*,'columns file not found'
23   continue
   endif
   
 end subroutine read_columns_file

!-------------------------------------------------------------------------      
!  this subroutine reads the contents of the .ev file
!-------------------------------------------------------------------------      

subroutine readev(nsteps,evdata,maxcols,ncolumns,rootname)
  implicit none
  integer :: i
  integer, intent(in) :: maxcols
  integer, intent(inout) :: nsteps
  integer, intent(out) :: ncolumns
  real, dimension(nsteps,maxcols), intent(out) :: evdata
  character(len=*), intent(in) :: rootname
  character(len=len_trim(rootname)+3) :: evname,dummy
  logical :: iexist,redo
  integer :: ierr,j,nskip

!--initially try a filename with .ev appended
  if (index(rootname,'.ev').eq.0) then
     evname = trim(rootname)//'.ev'
  else
     evname = trim(rootname)
  endif
  
  inquire (file=evname, exist=iexist)
  if (.not.iexist) then
     evname = trim(rootname)
!--if no .ev file, try without the .ev
     inquire (file=evname, exist=iexist)
     if (.not.iexist) then
        print*,'***ERROR ',trim(evname),': file does not exist'
        nsteps = 0
        return
     endif
  endif
  print "(a,a)",' opening ',trim(evname)
  
  open(unit=11,file=evname,status='old',form='formatted')
  
  call get_ncolumns(11,ncolumns)
  ncolumns = ncolumns + 1
  if (ncolumns.gt.maxcols) then
     ncolumns = maxcols
     print*,'***WARNING: array too small: reading first ',ncolumns,' only'
  else
     print*,'ncolumns in file = ',ncolumns-1
  endif
  
  redo = .true.
  nskip = 0
  
  do while (redo)
     ierr = 0
     i = 1
     do while (i.le.nsteps .and. ierr.eq.0)
        evdata(i,ncolumns) = real(i)
        read(11,*,iostat=ierr) evdata(i,1:ncolumns-1)
        if (nskip.gt.0) then
           do j=1,nskip
              read(11,*,iostat=ierr) dummy
           enddo
        endif
        i = i + 1
     enddo

     if (i.ge.nsteps) then
        print*,'** WARNING: number of steps > array limits, resampling...'
        rewind(11)
        nskip = nskip + 1
        print*,'re-reading only every ',nskip+1,' steps'
     else
        i = i - 1
        if (ierr.gt.0) then
           write(*,"(a)",advance="no") '>> ERROR IN FILE '
        else
           write(*,"(a)",advance="no") '>> end of file:'
        endif
        redo = .false.
     endif
  
  enddo
  close(unit=11)
  
  if (i.gt.1) then
     print*,' t=',evdata(i-1,1),' nsteps = ',i-1 
  else
     print*,'** WARNING: file empty!! no timesteps'
  endif
  nsteps = i-1     
  
  return
end subroutine readev

!
! utility to work out number of columns of real numbers
! in an ascii output file
!
! file must already be open and at the start
! slightly ad-hoc but its the best way I could think of!
!
subroutine get_ncolumns(lunit,ncolumns)
 implicit none
 integer, intent(in) :: lunit
 integer, intent(out) :: ncolumns
 integer :: ierr,i,nblanklines
 character(len=2000) :: line
 real :: dummyreal(100)

 nblanklines = 0
 line = ' '
 ierr = 0
 do while (len_trim(line).eq.0 .and. ierr.eq.0)
    read(lunit,"(a)",iostat=ierr) line
    nblanklines = nblanklines + 1
 enddo
 if (ierr .ne.0 ) then
    ncolumns = 0
    return
 else
    if (nblanklines.gt.1) print*,'skipped ',nblanklines-1,' blank lines'
    rewind(lunit)
 endif
 dummyreal = -666.0
 
 ierr = 0
 read(line,*,iostat=ierr,end=10) (dummyreal(i),i=1,size(dummyreal))
 if (ierr /= 0) then
    print*,'*** ERROR: file does not contain real numbers'
    ncolumns = 0
    return
 endif
10 continue 

 i = 1
 ncolumns = 0
 do while(abs(dummyreal(i)+666.).gt.1.e-10)
    ncolumns = ncolumns + 1
    i = i + 1
    if (i.gt.size(dummyreal)) then
       print*,'*** ERROR: too many columns in file'
       return
    endif
 enddo

end subroutine get_ncolumns

!-------------------------------------------------------------------------
! subroutine to work out the positions of maxima and minima in the
!  quantity plotted. Works out the period from the distance between two
!  minima or two maxima
!-------------------------------------------------------------------------
subroutine getmax(datain,time,max,freq,freq2)
  implicit none
  integer, parameter :: noffset = 3 ! >1 gives extra checks
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
  
  do i=1+noffset,max-(noffset+1)
     if ((datain(i).gt.datain(i-1)).and.(datain(i).gt.datain(i+1))) then
        if ((datain(i).gt.datain(i-noffset)).and. &
            (datain(i).gt.datain(i+noffset))) then
           nmax = nmax + 1
           if (mod(nmax,2).eq.0) then
              period = time(i)-timeprev            
              print 10,'max',time(i),period,1./period
              timeprev = time(i)
              avper = avper + period
           endif
        else
           print*,'max period of ',time(i)-timeprev,' rejected'
        endif
            
     elseif ((datain(i).lt.datain(i-1)).and.   &
          (datain(i).lt.datain(i+1))) then
        if ((datain(i).lt.datain(i-noffset)).and. &
            (datain(i).lt.datain(i+noffset))) then        
           nmin = nmin + 1
           if (mod(nmin,2).eq.1) then
              if (nmin.gt.1) then
                 period = time(i)-timeprev2            
                 print 10,'min',time(i),period,1./period
                 avper2 = avper2 + period
              endif
              timeprev2 = time(i)       
           endif
        else
           print*,'min period of ',time(i)-timeprev2,' rejected'
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
  print 20,'average max ',avper,freq,omega,nmax/2
  print 20,'average min ',avper2,freq2,omega2,nmin/2
20 format(a,' period = ',f8.4,' freq = ',f8.4,' omega = ',f8.4,' using ',i3,' values')  
  
  return
end subroutine getmax

!-------------------------------------------------------------------------
!  draw a legend for different line styles
!  uses current line style and colour
!-------------------------------------------------------------------------
subroutine legend(icall,text,hpos,vposin)
  implicit none
  integer, intent(in) :: icall
  character(len=*), intent(in) :: text
  real, intent(in) :: vposin, hpos
  real, dimension(2) :: xline,yline
  real :: xch, ych, xmin, xmax, ymin, ymax
  real :: vspace, vpos

  call pgstbg(0)           ! opaque text to overwrite previous
!
!--set horizontal and vertical position and spacing
!  in units of the character height
!
  vspace = 1.5  ! (in units of character heights)
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
