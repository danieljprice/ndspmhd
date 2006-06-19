program plotit
 use legends
 implicit none
 integer, parameter :: maxrow = 1000, maxfile = 12
 character(len=120), dimension(maxfile) :: filename,text
 character(len=20) :: string
 real, dimension(maxrow) :: col1,col2
 integer :: imax,i,ifile,ilinestyle,nfirstplot,nfiles,ipower
 real :: xmin,xmax,ymin,ymax,dummy
  
 filename(1) = 'results_psoft100.txt'
 text(1) = 'fixed plummer softening (100)'
 filename(2) = 'results_ksoft100.txt'
 text(2) = 'fixed kernel softening (100)'
 filename(3) = 'results_psoftbig.txt'
 text(3) = 'fixed plummer softening'
 filename(4) = 'results_ksoftbig.txt'
 text(4) = 'fixed kernel softening'
 filename(5) = 'results_psoft10k.txt'
 text(5) = 'fixed plummer softening (10k)'
 filename(6) = 'results_ksoft10k.txt'
 text(6) = 'fixed kernel softening (10k)'

 filename(7) = 'results_ksoft100var.txt'
 text(7) = 'adaptive (100)'
 filename(8) = 'results_ksoft100cons.txt'
 text(8) = 'adaptive +term (100)'
 filename(9) = 'results_ksoftvar.txt'
 text(9) = 'adaptive kernel softening'
 filename(10) = 'results_ksoftcons.txt'
 text(10) = 'adaptive + extra term'
 filename(11) = 'results_ksoft10kvar.txt'
 text(11) = 'adaptive (10k)'
 filename(12) = 'results_ksoft10kcons.txt'
 text(12) = 'adaptive +term (10k)'

 nfirstplot = 6
 nfiles = 12
 
 call pgbegin(0,'?',1,1)
 call pgpap(5.5,1./sqrt(2.))
 xmin = 2.e-3
 xmax = 9.9
 ymin = 1.e-2
 ymax = 3.0
 
 do ifile=1,nfiles
    open(unit=1,file=filename(ifile),status='old')
    do i=1,maxrow
!     if (ifile.le.nfirstplot) then
        read(1,*,end=99) col1(i),col2(i)
!     else
!        read(1,*,end=99) dummy,col1(i),col2(i)     
!     endif
    enddo
99  continue
    imax = i-1
    print*,trim(filename(ifile)),': points = ',imax
    close(unit=1)
    
    call pgslw(2)
    call pgsch(1.5)
    if (mod(ifile,2).eq.0) then
       call pgsls(2)
    else
       call pgsls(1)
    endif
    if (ifile.eq.1) then
       ipower = 1
       call pgsvp(0.1,1.0,0.13,0.99)
       !call pgmtxt('T',0.5,0.5,0.5,'Plummer sphere')
       call pgsvp(0.1,0.55,0.13,0.99)
       call pgswin(log10(xmin),log10(xmax),log10(ymin),log10(ymax))
       call pgbox('BCSTLN',0.0,0,'BCSTLN',0.0,0)
       call pgmtxt('B',2.5,0.5,0.5,'softening length')
       call pgmtxt('L',2.5,0.5,0.5,'MASE')
       call legend_markers(1,1,1,1, &
           .false.,.true.,'fixed plummer',0.2,0.25)
       call legend_markers(2,1,1,2, &
           .false.,.true.,'fixed cubic spline',0.2,0.25)
    elseif (ifile.eq.nfirstplot+1) then
       ipower = 1
       xmin = 25
       xmax = 999.0
       call pgsvp(0.55,1.0,0.13,0.99)
       call pgswin(log10(xmin),log10(xmax),log10(ymin),log10(ymax))
       call pgbox('BCSTLN',0.0,0,'BCSTL',0.0,0)
       call pgmtxt('B',2.5,0.5,0.5,'Number of Neighbours')
       call legend_markers(1,1,1,1, &
           .false.,.true.,'adaptive cubic spline',0.1,0.25)
       call legend_markers(2,1,1,2, &
           .false.,.true.,'adaptive + extra term',0.1,0.25)
    endif
    i = 1
    do while (col1(i).lt.xmin .or. col2(i).gt.col2(i+1))
       i = i + 1  
    enddo
    if (mod(ifile,2).eq.0 .and. ifile.le.nfirstplot) then
       ipower = ipower + 1
       write(string,"('10\u'i1)") ipower
       call pgsch(1.3)
       call pgtext(log10(col1(i))+0.25,log10(col2(i)),trim(string))
       call pgsch(1.5)
    endif
   ! call pgsci(ifile)

    call pgline(imax,log10(col1(1:imax)),log10(col2(1:imax)))
    call pgqls(ilinestyle)
!    if (ifile.le.nfirstplot) then
!       call legend_markers(ifile,ifile,1,ilinestyle, &
!           .false.,.true.,text(ifile),0.2,1.0)
!    else
!       call legend_markers(ifile-nfirstplot,ifile,1,ilinestyle, &
!           .false.,.true.,text(ifile),0.2,1.0)
!    endif
    
    if (.false.) then
       call pgsls(3)
       col1(1) = xmin
       col1(2) = xmax
       col2(1) = 0.06755937113374015
       col2(2) = col2(1)
       call pgline(2,log10(col1(1:2)),log10(col2(1:2)))
       call pgsls(4)
       col2(1) = 0.08886530004574236
       col2(2) = col2(1)
       call pgline(2,log10(col1(1:2)),log10(col2(1:2)))
    endif
    call pgsci(1)
 
 enddo

 call pgend


end program plotit
