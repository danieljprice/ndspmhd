!
! program to plot to content of the .freq file
!
program plotfreq
 implicit none
 integer, parameter :: max = 100
 integer :: i,npts,jmode,mmode,ierr,nacross,ndown
 real, dimension(max) :: rmode,freq,freq2,freq3,freqexact
 character(len=30) :: filename,rootname
 real, parameter :: pi=3.1415926536
 real :: omegasq,sigma2,gamm1,xpos,ypos,xch,ych

 call getarg(1,rootname)
 if (rootname(1:1).EQ.' ') stop 'usage: plotfreq filename'
 print*,' rootname = ',rootname
 if (index(rootname,'.freq').eq.0) then
    filename = trim(rootname)//'.freq'
 else
    filename = trim(rootname)
 endif
 freq = 1.
 freq2 = 1.
 freq3 = 1.
 
 open(1,file=filename,status='old',form='formatted')
 jmode = 0
 mmode = 2
 do i=1,max
    !!read(1,*,end=10,iostat=ierr) freq(i),freq2(i),freq3(i)
    read(1,*,end=10,iostat=ierr) freq2(i)
    !!freq(i) = 1./(freq(i)*period(i))
    
    jmode = jmode + 2
    if (jmode.eq.20) then
       jmode = 2
       mmode = mmode + 2
    endif
    rmode(i) = real(jmode)

    omegasq = 1.0
    gamm1 = 1.0
    sigma2 = 0.5*omegasq*(gamm1)*((jmode+mmode)*(jmode+mmode + 2./gamm1) - mmode**2)
    freqexact(i) = sqrt(sigma2)/(2.*pi)
    print*,i,freq(i),freq2(i),0.5*freq3(i),'period = ',1./freq(i),1./freq2(i),2./freq3(i),' exact = ',1./freqexact(i)
 enddo
10 continue
 close(unit=1)
 npts = i-1

 call pgbegin(0,'?',1,1)
!! call pgsch(1.5)
 call pgslw(2)
 
 nacross = 1
 ndown = 1
! call danpgtile(1,nacross,ndown,  &
!                0.,maxval(rmode(1:npts))+2.,0.,maxval(freqexact(1:npts)), &
!                'radial mode (j)','frequency',' ',0,1)
 !!call pglab(' radial mode ',' frequency ',' ')
 !!call pgpt(npts,rmode(1:npts),freq(1:npts),17)
! call pgpt(npts,rmode(1:npts),freqexact(1:npts),2)

 call pgsls(1)
! call pgline(npts,rmode(1:npts),freqexact(1:npts))
 !!call pgsci(2)
! call pgpt(npts,rmode(1:npts),0.5*freq2(1:npts),4)
! call pgpt(npts,rmode(1:npts),0.5*freq3(1:npts),17)

 call pgsls(4)
! call pgline(npts,rmode(1:npts),0.5*freq2(1:npts))
! call pgsls(2)
! call pgline(npts,rmode(1:npts),0.5*freq3(1:npts))

!
!--plot legend
! 
! call pgqcs(4,xch,ych)
! xpos = rmode(1)
! ypos = freqexact(npts-1) - 0.5*ych

!--exact
! call pgpt1(xpos,ypos+0.25*ych,2)
! call pgtext(xpos+xch,ypos,'exact')

! ypos = ypos - 1.25*ych
! call pgpt1(xpos,ypos+0.25*ych,4)
! call pgtext(xpos+xch,ypos,'equal mass particles ')

! ypos = ypos - 1.25*ych
! call pgpt1(xpos,ypos+0.25*ych,17)
! call pgtext(xpos+xch,ypos,'variable particle masses')

!
!--residuals
!
!! call pgpanl(2,1)
 call pgsls(1)
 call danpgtile(1,nacross,ndown,  &
                0.,maxval(rmode(1:npts))+2.,-0.15,0.15, &
                'radial mode (j)','(frequency-exact)/exact',' ',0,1)

 freq(1:npts) = (0.5*freq2(1:npts) - freqexact(1:npts))/freqexact(1:npts)
 call pgpt(npts,rmode(1:npts),freq(1:npts),4)
 call pgsls(4)
! call pgline(npts,rmode(1:npts),freq(1:npts))

! freq(1:npts) = (0.5*freq3(1:npts) - freqexact(1:npts))/freqexact(1:npts)
! call pgpt(npts,rmode(1:npts),freq(1:npts),17)
! call pgsls(2)
! call pgline(npts,rmode(1:npts),freq(1:npts))
 call pgsls(1)

 call pgend

 print*,'writing to output.freq'
 open(unit=1,file='output.freq',status='replace',form='formatted')
  do i=1,npts
     write(1,*) int(rmode(i)),freqexact(i),0.5*freq2(i),0.5*freq3(i), &
                (0.5*freq2(i)-freqexact(i))/freqexact(i), &
                (0.5*freq3(i)-freqexact(i))/freqexact(i)
  enddo
 close(unit=1)

end program plotfreq
