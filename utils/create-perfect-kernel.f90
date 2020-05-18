module kfuncs
 implicit none
 real, parameter :: pi = 4.d0*atan(1.d0)
 real :: r1 = 0.

contains

 ! kernel function defined by Fourier transform
 real function wk(k,r)
  real, intent(in) :: k,r
  real :: x1,x2

  !wk = sinc(k) !sinc(0.5*k)**4*k*besj0(k*r)
  x1 = k
  x2 = sqrt(2.d0)*k
  wk = sinc(x1)*pi**2/(pi**2 - x1**2)*sinc(x2)*pi**2/(pi**2 - x2**2)*k*besj0(k*r)
 
 end function wk

 ! kernel function defined by Fourier transform
 real function wk1(k)
  real, intent(in) :: k
  real :: x1,x2,x3

  !wk = sinc(k) !sinc(0.5*k)**4*k*besj0(k*r)
  x1 = k*0.5
  x2 = sqrt(2.)*k*0.5
  x3 = sqrt(5.)*k*0.5
  wk1 = sinck(x1)*sinck(x2)*sinck(x3)*k*besj0(k*r1)
  !wk1 = sinc(x1)*sinc(x2)*sinc(x3)*k*besj0(k*r1)
 
 end function wk1

 ! kernel function defined by Fourier transform
 real function m4(k)
  real, intent(in) :: k

  m4 = sinc(0.5*k)**2*cos(k*r1)
 
 end function m4

 ! Sinc function
 real function sinck(k)
  real, intent(in) :: k

  sinck = sinc(k)*pi**2/(pi**2 - k**2)

 end function sinck

 ! Sinc function
 real function sinc(x)
  real, intent(in) :: x

  if (abs(x) > epsilon(x)) then
     sinc = sin(x)/x
  else
     sinc = 1.
  endif

 end function sinc

 ! Sinc function
 real function wksinc(x,r)
  real, intent(in) :: x,r

  if (abs(x) > epsilon(x)) then
     wksinc = sin(x)/x
  else
     wksinc = 1.
  endif

 end function wksinc
 
 ! exact 3D kernel
 real function w3D(x)
  real, intent(in) :: x

  if (x < 0.5*(1. + sqrt(2.) - sqrt(3.))) then
     w3d = 0.5*(-3. + sqrt(2.) + sqrt(3.) + sqrt(6.) - 2.*x**2)
  elseif (x < 0.5*(1. - sqrt(2.) + sqrt(3.))) then
     w3d = 0.25*(-3. + 3.*sqrt(2.) + sqrt(3.) + sqrt(6.) - 2.*x &
           - 2.*sqrt(2.)*x + 2.*sqrt(3.)*x - 2.*x**2)  
  elseif (x < 0.5*(-1. + sqrt(2.) + sqrt(3.))) then
     w3d = 0.5*(sqrt(2.) + sqrt(3.) - 2.*x)
  elseif (x < 0.5*(1. + sqrt(2.) + sqrt(3.))) then
     w3d = 0.25*(3. + sqrt(2.) + sqrt(3.) + sqrt(6.) &
          - 2.*x - 2.*sqrt(2.)*x - 2.*sqrt(3.)*x  + 2.*x**2)
  else
     w3d = 0.
  endif
  w3d = w3d / ( sqrt(6.) * pi )
  
 end function w3d

 real function integrate_kernel(w,r,kmin,kmax,npts) result(wint)
  real, external :: w
  real, intent(in) :: r,kmin,kmax
  integer, intent(in) :: npts
  integer :: i
  real :: k,dk,wfunc,wprev
 
  dk = (kmax-kmin)/(npts-1)
  wint = 0.
  wprev = 0.
  do i=1,npts
     k = kmin + (i-1)*dk
     wfunc = w(k,r)
     ! trapezoidal rule
     wint  = wint + 0.5*(wfunc + wprev)*dk
     wprev = wfunc
  enddo
 
 end function integrate_kernel

end module kfuncs

program create_perfect_kernel
 use kernels,      only:ikern,write_kernel_to_file,dq2table
 use kernel_utils, only:diff,dp,normalise
 use kfuncs
 use gaussm3
 implicit none
 integer :: ndim,i,nrom
 integer, parameter :: lu = 66
 real(dp) :: q2,q,romberg_trap
 real(dp), dimension(0:ikern) :: wkern,grkern,grgrkern
 real(dp), dimension(3) :: cnormkd
 real(dp) :: radkern,radkern2
 
 radkern = 0.5*(1. + sqrt(2.) + sqrt(3.))
 !radkern = 0.5 + 0.5*sqrt(2.) + 0.5*sqrt(5.) !0.5 + 0.5*sqrt(2.) + 0.5*sqrt(5.) 
 !radkern = 1._dp + sqrt(2._dp)
 print*,' Rkern = ',radkern
 radkern2 = radkern*radkern
 dq2table = radkern2/real(ikern)
 ndim = 2
 print*,pi,2.*integrate_kernel(wksinc,0.,0.,1000.,100000)
 print*,'with Gauss = ',2.*qgauss(sinc,0.,20000.,10000)
 print*,'with romberg = ',2.*romberg_trap(0.,20000.,sinc,1.e-12,nrom)
 print "(a,i1,a,i5)",'generating exact kernel in ',ndim,'D: table size=',ikern
 
 r1 = 0.
! print*,'wk(0),k=10 = ',qgauss(wk1,0.,10.,1000),romberg_trap(0.,10.,wk1,1.e-8,nrom)
! print*,'wk(0),k=100 = ',qgauss(wk1,0.,100.,10000),romberg_trap(0.,100.,wk1,1.e-12,nrom)
! print*,'wk(0),k=1000 = ',qgauss(wk1,0.,1000.,1000),romberg_trap(0.,1000.,wk1,1.e-12,nrom)
! print*,nrom
 cnormkd(2) = 1./(2.*pi)
 call test_kernel()
 !stop
! print*,'wk(0),k=10000 = ',qgauss(wk1,0.,2000.,1000),romberg_trap(0.,2000.,wk1,1.e-12,nrom)
 
! !$omp parallel do default(none) shared(dq2table,wkern) private(i,q2,q)
 do i=0,ikern
    q2 = i*dq2table
    q  = sqrt(q2)
    r1 = q
    wkern(i) = w3d(q) !qgauss(wk1,0.,1000.,1000)
    !wkern(i) = romberg_trap(0.,100.,wk1,1.e-12,nrom)
    !if (mod(i,100).eq.0) print*,i,'nrom =',nrom
    !wkern(i) = integrate_kernel(wk,q,0.,100.,100000)
 enddo
! !$omp end parallel do
 !call create_kernel(qtable,wkern)
 call diff(ikern,wkern,grkern,grgrkern,radkern2)
 call normalise(ikern,wkern,cnormkd,radkern2)
 print*,' got norm = ',1./cnormkd(2),' should be ',2.*pi
 !cnormkd(2) = 1./(2.*pi)
 call write_kernel_to_file(200,lu,cnormkd,wkern,grkern,grgrkern)

contains

 subroutine test_kernel()
  integer, parameter :: nx = 10
  real, dimension(nx**3) :: x,y,z
  real :: dx,rmin,r2,q2,xi,yi,zi,hi,mi,rhoi
  real :: wij,wij1,wij2,wij3,wij4,wij5,kmax
  real :: tol,rhogi
  integer :: i,j,k,n,ip,nrom,nneigh
  
  dx = 1./nx  ! unit square
  rmin = huge(rmin)
  n = 0
  ip = 0
  do k=1,nx
     zi = (k-1)*dx
     do j=1,nx
        yi = (j-1)*dx
        do i=1,nx
           xi = (i-1)*dx
           n = n + 1
           x(n) = xi
           y(n) = yi
           z(n) = zi
           r2 = (xi-0.5)**2 + (yi-0.5)**2 + (zi - 0.5)**2
           if (r2 < rmin) then
              rmin = r2
              ip = n
           endif
        enddo
     enddo
  enddo
  hi = dx

  print*,'testing 3D density estimate ' !,ip,' x,y = ',x(ip),y(ip)
  print*,' kmax = ',2.*pi/dx
  print*,0.5*(1. + sqrt(2.) - sqrt(3.))
  print*,0.5*(1. - sqrt(2.) + sqrt(3.))
  print*,0.5*(-1. + sqrt(2.) + sqrt(3.))
  print*,0.5*(1. + sqrt(2.) + sqrt(3.))

  xi = x(ip)
  yi = y(ip)
  zi = z(ip)
  mi = 1.*dx*dx*dx
  rhoi = 0.
  rhogi = 0.
  nneigh = 0
  do j=1,n
     r2 = (xi-x(j))**2 + (yi - y(j))**2 + (zi - z(j))**2
     q2 = r2/dx**2
     if (q2 <= radkern2) then
        nneigh = nneigh + 1
        r1 = sqrt(q2)
        wij = w3D(r1)
        !wij = wij + qgauss(wk1,10.,10000.,4000)
        !wij = romberg_trap(0.,1000.,wk1,1.e-13,nrom)
        !wij = wij + romberg_trap(100.,100000.,wk1,1.e-12,nrom)
        !print*,j,r1,wij,cnormkd(2)
        print*,rhoi,mi*wij/hi**3,r1
        rhoi = rhoi + mi*wij/hi**3
        rhogi = rhogi + 1./(pi*sqrt(pi))*mi*exp(-q2)/hi**3
     endif
  enddo
  print*,' GOT rhoi = ',rhoi,' nneigh = ',nneigh
  print*,' truncated gaussian rhogi = ',rhogi

  print*,'testing 2D density estimate ' !,ip,' x,y = ',x(ip),y(ip)
  print*,' kmax = ',2.*pi/dx
  xi = x(ip)
  yi = y(ip)
  hi = dx
  mi = 1.*dx*dx
  rhoi = 0.
  rhogi = 0.
  nneigh = 0
  do j=1,nx*nx
     r2 = (xi-x(j))**2 + (yi - y(j))**2
     q2 = r2/dx**2
     if (q2 <= radkern2) then
        nneigh = nneigh + 1
        r1 = sqrt(q2)
        wij = qgauss(wk1,0.,10000.,10000)
        !wij = wij + qgauss(wk1,10.,10000.,4000)
        !wij = romberg_trap(0.,1000.,wk1,1.e-13,nrom)
        !wij = wij + romberg_trap(100.,100000.,wk1,1.e-12,nrom)
        !print*,j,r1,wij,cnormkd(2)
        print*,rhoi,cnormkd(2)*mi*wij/hi**2,r1
        rhoi = rhoi + cnormkd(2)*mi*wij/hi**2
        rhogi = rhogi + 1./pi*mi*exp(-q2)/hi**2
     endif
  enddo
  print*,' GOT rhoi = ',rhoi,' nneigh = ',nneigh
  print*,' truncated gaussian rhogi = ',rhogi
  
  kmax = 100000. !20.*pi/dx !10000.
  tol = 1.e-10
  r1 = 0.
  wij1 = romberg_trap(0.,kmax,wk1,tol,nrom)
  
  r1 = 1.
  wij2 = romberg_trap(0.,kmax,wk1,tol,nrom)

  r1 = sqrt(2.)
  wij3 = romberg_trap(0.,kmax,wk1,tol,nrom)

  r1 = 2.
  wij4 = romberg_trap(0.,kmax,wk1,tol,nrom)

  r1 = sqrt(5.)
  wij5 = romberg_trap(0.,kmax,wk1,tol,nrom)

  print*,'sum = ',(1.*wij1 + 4.*wij2 + 4.*wij3 + 4.*wij4 + 8.*wij5)*cnormkd(2)*mi/hi**2

  ! 1D particles in line
  print*,'1D test: m4 '
  ip = nx/2 + 1
  xi = x(ip)
  rhoi = 0.
  mi = 1.*dx
  do j=1,nx
     r2 = (xi - x(j))**2
     q2 = r2/dx**2
     if (q2 <= radkern2) then
        r1 = sqrt(q2)
        wij = 2.*romberg_trap(0.,100000.,m4,1.e-13,nrom) !/(2.*pi)
        rhoi = rhoi + cnormkd(2)*mi*wij/hi
     endif
  enddo
  print*,' GOT rhoi = ',rhoi
 
 end subroutine test_kernel
end program create_perfect_kernel
