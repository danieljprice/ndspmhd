subroutine kernelstability1D
  use kernel
  implicit none
  integer, parameter :: nh = 100, nkx = 100, npart = 20, ncont = 40
  integer :: i,j,ipart,mm,pp,nc,icall
  real, parameter :: pi = 3.1415926536
  real, dimension(nh,nkx) :: dat
  real, dimension(nh) :: harray
  real, dimension(nkx) :: kxarray
  real, dimension(npart) :: x
  real, dimension(ncont) :: levels
  real, dimension(6) :: trans
  real :: xmin, psep, pmass, rhozero, cs2
  real :: hmin, hmax, dh, dkx, kxmin, kxmax, h, kx
  real :: datmin, datmax, dcont, omegasq, omegasq1D
  character(len=5) :: string
  save icall
  
  icall = icall + 1
  !print*,'icall = ',icall
  cs2 = 1.0
  rhozero = 1.0
  psep = 1.0
  pmass = rhozero*psep
!
!--set up an array of particle positions to sum over
!
  xmin = -npart/2*psep
  do i=1,npart
     x(i) = xmin + i*psep
  enddo
  ipart = npart/2
!
!--set up an array of smoothing length values
!
  hmin = 0.5
  hmax = 2.0
  dh = (hmax - hmin)/REAL(nh)
  do i=1,nh
     harray(i) = hmin + (i-1)*dh
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
  do j=1,nh
     h = harray(j)*psep
     do i=1,nkx
        kx = kxarray(i)/psep
	omegasq = omegasq1D(h,kxarray(i),cs2,pmass,rhozero,x,npart,ipart)
        dat(i,j) = 0.5*omegasq/(kx**2*cs2)
!       print*,i,j,' kx = ',kx,' h = ',harray(j),h,' dat = ',dat(i,j)
     enddo
  enddo
!
!--now plot this using PGPLOT
!
!!  call pgbegin(0,'?',1,1)
!  call pgsch(1.2)
  call danpgtile(icall,3,2,kxmin,kxmax,hmin,hmax-0.001,'kx','h',TRIM(kernelname),0)
!  call pgenv(kxmin,kxmax,hmin,hmax,0,1)
!  call pglabel('kx','h',TRIM(kernelname))
  call pgsci(1)
  trans(1) = kxmin - 0.5*dkx               ! this is for the pgimag call
  trans(2) = dkx                  ! see help for pgimag/pggray/pgcont
  trans(3) = 0.0
  trans(4) = hmin - 0.5*dh
  trans(5) = 0.0
  trans(6) = dh
!
!--set contour levels
! 
  datmin = 0.0  !!minval(dat)
  datmax = int(maxval(dat)) + 1
  !!print*,'datmax = ',datmax,maxval(dat)*100,nint(maxval(dat)*100)
  dcont = (datmax-datmin)/real(ncont)   ! even contour levels
  do i=1,ncont
     levels(i) = datmin + real(i)*dcont
  enddo
  call pgcons(dat,nkx,nh,1,nkx,1,nh,levels,ncont,trans)
!
!--label levels
!
  do i=1,ncont

     MM=nint(levels(i)*100)
     PP=nint(log10(levels(i))-log10(levels(i)*100))
     call pgnumb(MM,PP,1,string,nc)
     !!print*,'level = ',levels(i),string(1:nc)
     call pgsch(0.5) ! character height
     call pgconl(dat,nkx,nh,1,nkx,1,nh,levels(i),trans,string(1:nc),35,20)
     call pgsch(0.6)
  enddo
!!  call pgend

end subroutine kernelstability1D

real function omegasq1D(hh,kx,cs2,pmass,rhozero,x,npart,ipart)
  use kernel
  implicit none
  integer :: i,index,index1, ipart, npart
  real, dimension(npart) :: x
  real :: hh,kx,cs2,pmass,rhozero
  real :: sum1, sum2
  real :: dxx, dgrwdx, dgrgrwdx, gradW, gradgradW
  real :: dx, q2

  sum1 = 0.
  sum2 = 0.

  do i=1,npart
     dx = x(i) - x(ipart)
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
  omegasq1D = (2.*cs2*pmass/rhozero)*sum1 - cs2*(pmass/rhozero*sum2)**2

end function omegasq1D
