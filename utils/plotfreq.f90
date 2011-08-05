!
! program to plot to content of the .freq file
!
program plotfreq
 implicit none
 integer, parameter :: max = 100
 integer :: i,npts
 real, dimension(max) :: hfact,freq
 character(len=30) :: filename,rootname

 call getarg(1,rootname)
 if (rootname(1:1).EQ.' ') stop 'usage: plotfreq filename'
 print*,' rootname = ',rootname
 if (index(rootname,'.freq').eq.0) then
    filename = trim(rootname)//'.freq'
 else
    filename = trim(rootname)
 endif
 
 open(1,file=filename)
 do i=1,max
    read(1,*,end=10) hfact(i),freq(i)
    print*,i,hfact(i),freq(i)
 enddo
10 continue
 close(1)
 npts = i-1

 call pgbegin(0,'?',1,1)
 call pgenv(1.0,2.0,minval(freq(1:npts)),maxval(freq(1:npts)),0,1)
 call pglab(' hfact ',' numerical sound speed ',' ')
 call pgline(npts,hfact(1:npts),freq(1:npts))
 call pgend

end program plotfreq
