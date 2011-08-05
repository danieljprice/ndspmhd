subroutine kernelstability1D(iplot,nacross,ndown)
  use kernel
  implicit none
  integer, parameter :: ny = 100, nkx = 100, npart = 20, ncont = 40
  integer, intent(in) :: iplot, nacross, ndown
  integer :: i,j,ipart,mm,pp,nc
  real, parameter :: pi = 3.1415926536
  real, dimension(nkx,ny) :: dat
  real, dimension(ny) :: yaxis
  real, dimension(nkx) :: kxarray
  real, dimension(npart) :: x
  real, dimension(ncont) :: levels
  real, dimension(6) :: trans
  real :: xmin, psep, pmass, rhozero, cs2, R
  real :: ymin, ymax, dy, dkx, kxmin, kxmax, h, kx
  real :: datmin, datmax, dcont, omegasq, omegasq1D
  real :: charheight
  character(len=5) :: string,labely
  character(len=50) :: filename, text
  logical :: negstress, contours
  
  cs2 = 1.0
  rhozero = 1.0
  psep = 1.0
  pmass = rhozero*psep
!
!--set various options
!  
  contours = .false.     ! plot whole dispersion relation or just kx=0
  negstress = .false.    ! plot vs h or R
  R = 1.0       ! R=1 gives usual hydrodynamics, R < 0 gives negative stress
  h = 1.2*psep   ! value of smoothing length
!
!--set up an array of particle positions to sum over
!
  xmin = -npart/2*psep
  do i=1,npart
     x(i) = xmin + i*psep
  enddo
  ipart = npart/2
!
!--set up an array of smoothing length values for the y axis
!
  if (negstress) then   ! y axis is negative stress parameter R
     ymin = -1.
     ymax = 1.
     labely = 'R'
  else                  ! y axis is h
     ymin = 0.5*psep
     ymax = 2.0*psep
     labely = 'h'
  endif
!
!--set up y axis
!
  dy = (ymax - ymin)/REAL(ny)
  do i=1,ny
     yaxis(i) = (ymin + (i-1)*dy)
  enddo 
!
!--set up an array of wavenumber values
!
  kxmin = 0.0
  kxmax = 2.*pi
  dkx = (kxmax - kxmin)/REAL(nkx)
  do i=1,nkx
     kxarray(i) = kxmin + i*dkx
  enddo
!
!--calculates the 1D dispersion relation for SPH (no av)
!  
  do j=1,ny
     if (negstress) then  ! y axis is R
	R = yaxis(j)
     else                 ! y axis is h
        h = yaxis(j)
     endif
     do i=1,nkx
        kx = kxarray(i)/psep
	omegasq = omegasq1D(h,kxarray(i),cs2,pmass,rhozero,R,x,npart,ipart)
        dat(i,j) = omegasq/(kx**2*cs2)
!!       print*,i,j,' kx = ',kx,' h = ',h,' R = ',R,' dat = ',dat(i,j)
     enddo
  enddo
!
!--now plot this using PGPLOT
!
  if (contours .or. negstress) then
!!  call pgbegin(0,'?',1,1)
!!  call pgsch(1.2)
     call danpgtile(iplot,nacross,ndown,kxmin,kxmax,ymin,ymax-0.001, &  ! tiled plots
                 'kx',TRIM(labely),TRIM(kernelname),0,0)
!  call pgenv(kxmin,kxmax,ymin,ymax,0,0)            ! use this for movie
!  call pglabel('kx',TRIM(labely),TRIM(kernelname)) ! use this for movie

     trans(1) = kxmin - 0.5*dkx               ! this is for the pgimag call
     trans(2) = dkx                  ! see help for pgimag/pggray/pgcont
     trans(3) = 0.0
     trans(4) = ymin - 0.5*dy
     trans(5) = 0.0
     trans(6) = dy
!
!--set contour levels
! 
     datmin = minval(dat)
     datmax = maxval(dat)
     print*,'datmax = ',datmax,datmin  !!maxval(dat)*100,nint(maxval(dat)*100)
     dcont = (datmax-datmin)/real(ncont)   ! even contour levels
     do i=1,ncont
        levels(i) = datmin + real(i)*dcont
     enddo
     call pgcont(dat,nkx,ny,1,nkx,1,ny,levels,ncont,trans)
!
!--label levels
!
     do i=1,ncont
        write(string,"(f5.2)") levels(i)
        !!print*,'level = ',levels(i),string
        call pgqch(charheight) ! query character height
        call pgsch(0.5*charheight)   ! shrink character height
        call pgconl(dat,nkx,ny,1,nkx,1,ny,levels(i),trans,TRIM(string),35,14)
        call pgsch(charheight) ! restore character height
     enddo
  
  else
!
!--plot just the line along kx = 0
!
     call pgsch(1.0)
     call danpgtile(iplot,nacross,ndown, &
          1.0,1.99,0.8,1.11, &  ! tiled plots
          'h','cs',TRIM(kernelname),0,0)
     call pgline(ny,yaxis(1:ny),sqrt(dat(1,1:ny)))
     !
     !--plot line from file
     !
     call pgsls(2)
     if (iplot.eq.1) then
        open(unit=11,file='legend',status='old')
20      continue
        do i=1,4
           call getarg(i,filename)
	   call pgsls(i)
           if (filename(1:1).ne.' ') then
	      call plotfreq(filename)
	      read(11,"(a)",end=21) text
21            continue
	      print*,trim(text),filename
              call legend(i,trim(text),0.6,25.0)
	   endif
	enddo
        close(unit=11)
     endif

     call pgsch(charheight)
     call pgsls(1)
  endif

end subroutine kernelstability1D

real function omegasq1D(hh,kx,cs2,pmass,rhozero,RR,x,npart,ipart)
  use kernel
  implicit none
  integer :: i,index,index1, ipart, npart
  real, dimension(npart) :: x
  real :: hh,kx,cs2,pmass,rhozero,RR
  real :: sum1, sum2
  real :: dxx, dgrwdx, dgrgrwdx, gradW, gradgradW
  real :: dx, q2

  sum1 = 0.
  sum2 = 0.

  do i=1,npart
     dx = abs(x(i) - x(ipart))
!
!--find nearest index in kernel table
! 
     q2 = (dx/hh)**2
     index = INT(q2*ddq2table)
     index1 = index + 1
     IF (index.GT.ikern) index = ikern
     IF (index1.GT.ikern) index1 = ikern
     dxx = q2 - index*dq2table
     dgrwdx = (grwij(index1)-grwij(index))*ddq2table
     dgrgrwdx = (grgrwij(index1)-grgrwij(index))*ddq2table
     gradW = (grwij(index) + dgrwdx*dxx)/hh**2
     gradgradW = (grgrwij(index) + dgrgrwdx*dxx)/hh**3 
     
     sum1 = sum1 + (1. - COS(kx*dx))*gradgradW  ! times unit vector in x dir (=1)
     sum2 = sum2 + SIN(kx*dx)*gradW
  enddo
  omegasq1D = (2.*cs2*pmass/rhozero)*RR*sum1 + (1.-2.*RR)*cs2*(pmass/rhozero*sum2)**2

end function omegasq1D

!
! subroutine to plot to content of the .freq file
!
subroutine plotfreq(filename)
 implicit none
 integer, parameter :: max = 100
 integer :: i,npts
 real, dimension(max) :: hfact,freq
 character(len=50) :: filename
 
 if (index(filename,'.freq').eq.0) then
    filename = trim(filename)//'.freq'
 endif
 
 print*,'opening ',trim(filename)
 
 open(1,file=filename,status='old',err=20)
 do i=1,max
    read(1,*,end=10) hfact(i),freq(i)
    !!print*,i,hfact(i),freq(i)
 enddo
10 continue
 close(1)
 npts = i-1

 if (filename(1:3).eq.'hfi') then
    call pgpt(npts,hfact(1:npts),freq(1:npts),17)
 else
    call pgline(npts,hfact(1:npts),freq(1:npts))
 endif
! call pgpt(npts,hfact(1:npts),freq(1:npts),17)
20 continue

 return
end subroutine plotfreq
