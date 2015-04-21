!-------------------------------------------------------
!--module to compute csplines for n>=4
!--uses the Hankel transform of order 0 of a nice kernel
!--rlim, eps and tol finely tuned to have beautiful results !
!-------------------------------------------------------
module csplinekernels
 implicit none
 public :: getcsplinekernel,getcsplinekernelder

 private

contains
!------------------------------------------------------
!  Subroutine to compute the kernel and its derivatives
!------------------------------------------------------
subroutine getcsplinekernel(nc,radkern,q,w,gw,ggw)
 implicit none
 integer, intent(in) :: nc
 real, intent(in)    :: q,radkern
 real, intent(out)   :: w
 real, intent(out)   :: gw
 real, intent(out)   :: ggw
 real                :: absq
 real, parameter     :: rlim = 2.d2
 real, external :: besselstuff,dbesselstuff,ddbesselstuff,ddzbesselstuff

 absq = abs(q)
 if (abs(q).lt.radkern) then
    call qromb_twoparam(besselstuff,nc,q,0.,rlim,w)
    if (w.lt.0.) w = 0.
    call qromb_twoparam(dbesselstuff,nc,q,0.,rlim,gw)
    gw = -gw
    if (gw.gt.0.) gw = 0.
    if (q.gt.tiny(0.)) then
       call qromb_twoparam(ddbesselstuff,nc,q,0.,rlim,ggw)
       ggw = - ggw - gw/q
    else
       call qromb_twoparam(ddzbesselstuff,nc,q,0.,rlim,ggw)
       ggw = -0.5*ggw
    endif
 else
    w   = 0.
    gw  = 0.
    ggw = 0.
 endif

end subroutine getcsplinekernel

subroutine getcsplinekernelder(nc1,nc2,secondz,radkern,q,w,gw,ggw)
 implicit none
 integer, intent(in) :: nc1,nc2
 real, intent(in)    :: q,radkern,secondz
 real, intent(out)   :: w
 real, intent(out)   :: gw
 real, intent(out)   :: ggw
 real                :: absq
 real, parameter     :: rlim = 2.d2
 real, external :: besselstuffder,dbesselstuffder
 real, external :: ddbesselstuffder,ddzbesselstuffder

 absq = abs(q)
 if (abs(q).lt.radkern) then
    call qromb_twoparamder(besselstuffder,nc1,nc2,secondz,q,0.,rlim,w)
    if (w.lt.0.) w = 0.
    call qromb_twoparamder(dbesselstuffder,nc1,nc2,secondz,q,0.,rlim,gw)
    gw = -gw
    if (gw.gt.0.) gw = 0.
    if (q.gt.tiny(0.)) then
       call qromb_twoparamder(ddbesselstuffder,nc1,nc2,secondz,q,0.,rlim,ggw)
       ggw = - ggw - gw/q
    else
       call qromb_twoparamder(ddzbesselstuffder,nc1,nc2,secondz,q,0.,rlim,ggw)
       ggw = -0.5*ggw
    endif
 else
    w   = 0.
    gw  = 0.
    ggw = 0.
 endif

end subroutine getcsplinekernelder

end module

!------------------------------------------------------------------
!--utilities for the standard csplines
!------------------------------------------------------------------

!------------------------------------------------------------------
     SUBROUTINE qromb_twoparam(func,nc,q,a,b,ss)
     implicit none
     INTEGER JMAX,JMAXP,K,KM,nc
     REAL a,b,func,q,ss,EPS
     EXTERNAL func
     PARAMETER (EPS=1.d-8, JMAX=30, JMAXP=JMAX+1, K=5, KM=K-1)
     INTEGER j
     REAL dss,h(JMAXP),s(JMAXP)
     real, parameter :: tol = 1.d-9
     h(1)=1.
     do 11 j=1,JMAX
       call trapzd_twoparam(func,nc,q,a,b,s(j),j)
       if (j.ge.K) then
         call polint_two(h(j-KM),s(j-KM),K,0.,ss,dss)
  !--lines added to hack some spurious breaks------------
         if (j.ne.K .and. j.ne.K+1) then !added
            if (abs(dss).le.EPS*abs(ss)) return
            if (abs(dss).lt.tol) return  !added
  !            print*,abs(dss),abs(ss)
         endif
  !-------------------------------------------------------         
       endif
       s(j+1)=s(j)
       h(j+1)=0.25*h(j)
11    continue
     print*,'too many steps in qromb_twoparam'
     end subroutine qromb_twoparam

!------------------------------------------------------------------
     SUBROUTINE polint_two(xa,ya,n,x,y,dy)
     implicit none
     INTEGER n,NMAX
     REAL dy,x,y,xa(n),ya(n)
     PARAMETER (NMAX=10)
     INTEGER i,m,ns
     REAL den,dif,dift,ho,hp,w,c(NMAX),d(NMAX)
     ns=1
     dif=abs(x-xa(1))
     do i=1,n
       dift=abs(x-xa(i))
       if (dift.lt.dif) then
         ns=i
         dif=dift
       endif
       c(i)=ya(i)
       d(i)=ya(i)
     enddo
     y=ya(ns)
     ns=ns-1
     do m=1,n-1
       do i=1,n-m
         ho=xa(i)-x
         hp=xa(i+m)-x
         w=c(i+1)-d(i)
         den=ho-hp
         if(den.eq.0.) print*,'failure in polint_two'
         den=w/den
         d(i)=hp*den
         c(i)=ho*den
      enddo
       if (2*ns.lt.n-m)then
         dy=c(ns+1)
       else
         dy=d(ns)
         ns=ns-1
       endif
       y=y+dy
     enddo
     return
     end subroutine polint_two

!------------------------------------------------------------------
     SUBROUTINE trapzd_twoparam(func,nc,q,a,b,s,n)
     implicit none
     INTEGER n,nc
     REAL a,b,s,func,q
     INTEGER it,j
     REAL del,sum,tnm,x
     if (n.eq.1) then
       s=0.5*(b-a)*(func(a,q,nc)+func(b,q,nc))
     else
       it=2**(n-2)
       tnm=it
       del=(b-a)/tnm
       x=a+0.5*del
       sum=0.
       do j=1,it
         sum=sum+func(x,q,nc)
         x=x+del
       enddo
       s=0.5*(s+(b-a)*sum/tnm)
     endif
     return
     END subroutine trapzd_twoparam

!------------------------------------------------------------------
     real function besselstuff(x,q,nc)
     implicit none
     integer, intent(in) :: nc
     real, intent(in)    :: x,q
     real :: ker,rinter1
     real, parameter    :: pi = 3.141592653589

     if (x .le. tiny(0.)) then
        ker = 1.
     else
        ker = 2.*bessel_j1(0.5*x)/(0.5*x)
     endif
     rinter1 = 1./(2.*pi)*ker**nc
     besselstuff = rinter1*x*bessel_j0(q*x)

     end function besselstuff
!------------------------------------------------------------------
     real function dbesselstuff(x,q,nc)    
     implicit none
     integer, intent(in) :: nc
     real, intent(in)    :: x,q
     real :: ker,rinter1
     real, parameter    :: pi = 3.141592653589

     if (x .le. tiny(0.)) then
        ker = 1.
     else
        ker = 2.*bessel_j1(0.5*x)/(0.5*x)
     endif
     rinter1 = 1./(2.*pi)*ker**nc
     dbesselstuff = rinter1*x*x*bessel_j1(q*x)

     end function dbesselstuff
 !------------------------------------------------------------------
     real function ddbesselstuff(x,q,nc)  
     implicit none
     integer, intent(in) :: nc
     real, intent(in)    :: x,q
     real :: ker,rinter1
     real, parameter    :: pi = 3.141592653589

     if (x .le. tiny(0.)) then
        ker = 1.
     else
        ker = 2.*bessel_j1(0.5*x)/(0.5*x)
     endif
     rinter1 = 1./(2.*pi)*ker**nc
     ddbesselstuff = rinter1*x*x*x*bessel_j0(q*x)

     end function ddbesselstuff    
 !------------------------------------------------------------------
     real function ddzbesselstuff(x,q,nc)   
     implicit none
     integer, intent(in) :: nc
     real, intent(in)    :: x,q
     real :: ker,rinter1
     real, parameter    :: pi = 3.141592653589

     if (x .le. tiny(0.)) then
        ker = 1.
     else
        ker = 2.*bessel_j1(0.5*x)/(0.5*x)
     endif
     rinter1 = 1./(2.*pi)*ker**nc
     ddzbesselstuff = rinter1*x*x*x

     end function ddzbesselstuff

!------------------------------------------------------------------
!--utilities for the csplines derivatives
!------------------------------------------------------------------

!------------------------------------------------------------------
     SUBROUTINE qromb_twoparamder(func,nc1,nc2,secondz,q,a,b,ss)
     implicit none
     INTEGER JMAX,JMAXP,K,KM,nc1,nc2
     REAL a,b,func,q,ss,EPS,secondz
     EXTERNAL func
     PARAMETER (EPS=1.d-8, JMAX=30, JMAXP=JMAX+1, K=5, KM=K-1)
     INTEGER j
     REAL dss,h(JMAXP),s(JMAXP)
     real, parameter :: tol = 1.d-9
     h(1)=1.
     do 11 j=1,JMAX
       call trapzd_twoparamder(func,nc1,nc2,secondz,q,a,b,s(j),j)
       if (j.ge.K) then
         call polint_two(h(j-KM),s(j-KM),K,0.,ss,dss)
  !--lines added to hack some spurious breaks------------
         if (j.ne.K .and. j.ne.K+1) then !added
            if (abs(dss).le.EPS*abs(ss)) return
            if (abs(dss).lt.tol) return  !added
  !            print*,abs(dss),abs(ss)
         endif
  !-------------------------------------------------------         
       endif
       s(j+1)=s(j)
       h(j+1)=0.25*h(j)
11    continue
     print*,'too many steps in qromb_twoparam'
     end subroutine qromb_twoparamder

!------------------------------------------------------------------
     SUBROUTINE trapzd_twoparamder(func,nc1,nc2,secondz,q,a,b,s,n)
     implicit none
     INTEGER n,nc1,nc2
     REAL a,b,s,func,q,secondz
     INTEGER it,j
     REAL del,sum,tnm,x
     if (n.eq.1) then
       s=0.5*(b-a)*(func(a,q,nc1,nc2,secondz)+func(b,q,nc1,nc2,secondz))
     else
       it=2**(n-2)
       tnm=it
       del=(b-a)/tnm
       x=a+0.5*del
       sum=0.
       do j=1,it
         sum=sum+func(x,q,nc1,nc2,secondz)
         x=x+del
       enddo
       s=0.5*(s+(b-a)*sum/tnm)
     endif
     return
     END subroutine trapzd_twoparamder

!------------------------------------------------------------------
     real function besselstuffder(x,q,nc1,nc2,secondz)
     implicit none
     integer, intent(in) :: nc1,nc2
     real, intent(in)    :: x,q,secondz
     real :: ker1,ker2,rinter1
     real, parameter    :: pi    = 3.141592653589

     if (x .le. tiny(0.)) then
        ker1 = 1.
        ker2 = 1.        
     else
        ker1 = 2.*bessel_j1(0.5*x)/(0.5*x)
        ker2 = 2.*bessel_j1(0.5*secondz*x)/(0.5*secondz*x)
     endif
     rinter1 = 1./(2.*pi)*ker1**nc1*ker2**nc2
     besselstuffder = rinter1*x*bessel_j0(q*x)

     end function besselstuffder
!------------------------------------------------------------------
     real function dbesselstuffder(x,q,nc1,nc2,secondz)    
     implicit none
     integer, intent(in) :: nc1,nc2
     real, intent(in)    :: x,q,secondz
     real :: ker1,ker2,rinter1
     real, parameter    :: pi    = 3.141592653589
     
     if (x .le. tiny(0.)) then
        ker1 = 1.
        ker2 = 1.
     else
        ker1 = 2.*bessel_j1(0.5*x)/(0.5*x)
        ker2 = 2.*bessel_j1(0.5*secondz*x)/(0.5*secondz*x)
     endif
     rinter1 = 1./(2.*pi)*ker1**nc1*ker2**nc2
     dbesselstuffder = rinter1*x*x*bessel_j1(q*x)

     end function dbesselstuffder
 !------------------------------------------------------------------
     real function ddbesselstuffder(x,q,nc1,nc2,secondz)  
     implicit none
     integer, intent(in) :: nc1,nc2
     real, intent(in)    :: x,q,secondz
     real :: ker1,ker2,rinter1
     real, parameter    :: pi    = 3.141592653589    

     if (x .le. tiny(0.)) then
        ker1 = 1.
        ker2 = 1.        
     else
        ker1 = 2.*bessel_j1(0.5*x)/(0.5*x)
        ker2 = 2.*bessel_j1(0.5*secondz*x)/(0.5*secondz*x)        
     endif
     rinter1 = 1./(2.*pi)*ker1**nc1*ker2**nc2
     ddbesselstuffder = rinter1*x*x*x*bessel_j0(q*x)

     end function ddbesselstuffder    
 !------------------------------------------------------------------
     real function ddzbesselstuffder(x,q,nc1,nc2,secondz)   
     implicit none
     integer, intent(in) :: nc1,nc2
     real, intent(in)    :: x,q,secondz
     real :: ker1,ker2,rinter1
     real, parameter    :: pi    = 3.141592653589

     if (x .le. tiny(0.)) then
        ker1 = 1.
        ker2 = 1.        
     else
        ker1 = 2.*bessel_j1(0.5*x)/(0.5*x)
        ker2 = 2.*bessel_j1(0.5*secondz*x)/(0.5*secondz*x)         
     endif
     rinter1 = 1./(2.*pi)*ker1**nc1*ker2**nc2
     ddzbesselstuffder = rinter1*x*x*x

     end function ddzbesselstuffder
