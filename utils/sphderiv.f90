program sphderiv
  !
  !--this program is to see how well SPH estimates derivatives for
  !  various situations
  !
  use kernel
  use loguns
  use options
  implicit none
  integer, parameter :: npart = 16
  real, parameter :: pi = 3.1415926536
  integer :: i,j
  real, dimension(npart) :: x,rho,h,pmass,Apart
  real, dimension(npart) :: gradApart1,gradApart2,gradApart3
  real :: psep,xmin,xmax,dx,dr,hfact,rhozero,rholeft,rhoright,totmass
  real :: q2,wab,grkern
  real :: psepright
  real, external :: A, gradA
  real :: Amin,Amax,xminpart,xmaxpart
  !
  !--set parameters
  !
  hfact = 1.2
  !
  !--setup kernel tables
  !
  ikernel = 0
  iprint = 6
  call setkern
  print*,'ndim = ',ndim
  !
  !--setup particles
  !
  xmin = -1.0
  xmax = 1.0  !2.*pi
  psep = (xmax - xmin)/REAL(npart-6)
  xminpart = xmin - 3.*psep
  xmaxpart = xmax + 3.*psep
  rhozero = 1.0
  rholeft = rhozero
  rhoright = rhozero/2.
  psepright = psep*rholeft/rhoright
  totmass = rhozero*(xmax-xmin)

  do i=1,npart
     if (i.eq.1) then
        x(i) = xminpart
        rho(i) = rholeft  ! just to estimate h
     elseif (x(i-1).lt.0.) then
        x(i) = x(i-1) + psep
        rho(i) = rholeft ! just to estimate h
     else
        x(i) = x(i-1) + psepright
        rho(i) = rhoright  ! just to estimate h
     endif
     pmass(i) = totmass/REAL(npart-6)
     h(i) = hfact*(pmass(i)/rho(i))
     Apart(i) = A(x(i))
  enddo

  print*,'pmass = ',pmass(5)
  !
  !--calculate density on central particles
  !
  do i=3,npart-3
     rho(i) = 0.
     do j=1,npart
        dx = x(i)-x(j)
        q2 = dx*dx/h(i)**2
        call interpolate_kernel(q2,wab,grkern)
        rho(i) = rho(i) + pmass(j)*wab/h(i)
     enddo
  enddo
  !
  !--copy density and h to edge particles
  !
  do i=1,3
     rho(i) = rho(4)
  enddo
  do i=npart-2,npart
     rho(i) = rho(npart-3)
  enddo
  print*,'rho = ',rho
  !
  !--now estimate SPH derivative estimates of function
  !
  gradApart1 = 0.
  gradApart2 = 0.
  gradApart3 = 0.

  do i=3,npart-3
     do j=1,npart
        if (j.ne.i) then
           dx = x(i)-x(j)
           dr = dx / abs(dx)
           q2 = dx*dx/h(i)**2
           call interpolate_kernel(q2,wab,grkern)
           gradApart1(i) = gradApart1(i) &
                + pmass(j)*Apart(j)/rho(j)*dr*grkern/h(i)**2
           gradApart2(i) = gradApart2(i) &
                + 1./rho(i)*pmass(j)*(Apart(j) - Apart(i))*dr*grkern/h(i)**2
           !!print*,i,j,Apart(i),Apart(j),q2,grkern
           gradApart3(i) = gradApart3(i) &
                + rho(i)*pmass(j)*(Apart(j)/rho(j)**2 + Apart(i)/rho(i)**2)*dr*grkern/h(i)**2

        endif
     enddo
  enddo  
  !
  !--plot using PGPLOT
  !
  Amin = min(minval(Apart),minval(gradApart1),minval(gradApart2),minval(gradApart3)) - 0.1
  Amax = max(maxval(Apart),maxval(gradApart1),maxval(gradApart2),maxval(gradApart3)) + 0.1
  call pgbegin(0,'?',1,1)
  call pgenv(xmin,xmax,Amin,Amax,0,0)
  call pglab('x','A','SPH derivatives')
  !
  !--plot density
  !
  print*,'plotting density'
  call pgpt(npart,x,rho,3)
  !
  !--plot function
  !
  call pgfunx(A,100,xmin,xmax,1)
  call pgpt(npart,x,Apart,17)
 
  !
  !--plot derivative
  !
  print*,'gradApart1 = ',gradApart1
  call pgsls(2)
  call pgfunx(gradA,100,xmin,xmax,1)
  call pgsls(1)
  print*,'crap version'
  call pgpt(npart,x,gradApart1,6)
  call pgline(npart,x,gradApart1)
  read*
  print*,'const exact'
  call pgpt(npart,x,gradApart2,5)
  call pgline(npart,x,gradApart2)
  read*
  print*,'mom cons'
  call pgpt(npart,x,gradApart3,4)
  call pgline(npart,x,gradApart3)
  read*
  call pgend

end program sphderiv

!
!--function to plot
!
real function A(x)
  implicit none
  real :: x

!
!--sin
!
!  A = SIN(x)
!
!--step function
!
!  if (x < 0) then
!     A = 0.1
!  else
     A = 1.0
!  endif

end function A
!
!--gradient of function to plot
!
real function gradA(x)
  implicit none
  real :: x

!
!--sin
!
!  gradA = COS(x)
!
!--step function
!
  gradA = 0.

end function gradA
