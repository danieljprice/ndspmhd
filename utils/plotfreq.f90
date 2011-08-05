!
! program to plot to content of the .freq file
!
program plotfreq
 implicit none
 integer, parameter :: max = 100
 integer :: i,npts,jmode,mmode,ierr
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
 do i=1,max
    read(1,*,end=10,iostat=ierr) freq(i),freq2(i),freq3(i)
    !!freq(i) = 1./(freq(i)*period(i))
    rmode(i) = real(i)*2

    omegasq = 1.0
    jmode = int(rmode(i))
    mmode = 0
    gamm1 = 1.0
    sigma2 = 0.5*omegasq*(gamm1)*((jmode+mmode)*(jmode+mmode + 2./gamm1) - mmode**2)
    freqexact(i) = sqrt(sigma2)/(2.*pi)
    print*,i,freq(i),freq2(i),0.5*freq3(i),'period = ',1./freq(i),1./freq2(i),2./freq3(i),' exact = ',1./freqexact(i)
 enddo
10 continue
 close(unit=1)
 npts = i-1

 call pgbegin(0,'?',1,1)
 call pgslw(2)
 call pgenv(0.,maxval(rmode(1:npts)),0.,maxval(freqexact(1:npts)),0,1)
 call pglab(' radial mode ',' frequency ',' ')
 !!call pgpt(npts,rmode(1:npts),freq(1:npts),17)
 call pgpt(npts,rmode(1:npts),freqexact(1:npts),2)
 call pgsls(4)
 call pgslw(1)
 call pgline(npts,rmode(1:npts),freqexact(1:npts))
 !!call pgsci(2)
 !!call pgpt(npts,rmode(1:npts),freq2(1:npts),4)
 call pgpt(npts,rmode(1:npts),0.5*freq3(1:npts),17)
 
!
!--plot legend
! 
 xpos = rmode(1)
 ypos = freqexact(npts-1)
 call pgqcs(4,xch,ych)
!--exact
 call pgpt1(xpos,ypos+0.25*ych,2)
 call pgtext(xpos+xch,ypos,'exact')

 ypos = ypos - 1.25*ych
 call pgpt1(xpos,ypos+0.25*ych,17)
 call pgtext(xpos+xch,ypos,'1000 particles (variable mass)')


 call pgend

end program plotfreq
