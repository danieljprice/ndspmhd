module erbskernels
 implicit none
 public :: geterbskernel1,geterbskernel2

 private

contains
!------------------------------------------------------
!  Subroutine to compute the kernel and its derivatives
!------------------------------------------------------
subroutine geterbskernel1(q,w,gw,ggw)
 implicit none
 real, intent(in)    :: q
 real, intent(out)   :: w
 real, intent(out)   :: gw
 real, intent(out)   :: ggw
 real :: absq,rinter1,rinter2,rinter3,rinter4
 real, external :: psierbs
 real, parameter :: normerbs = 0.6034501612189381

 absq = abs(q)
 if (absq .ge. 2) then
    w   = 0.
    gw  = 0.
    ggw = 0.       
 elseif(absq .le. tiny(0.)) then
    w   = 0.5
    gw  = 0.
    ggw = 0.
 else
    rinter1 = absq-1.
    rinter2 = absq*(2.-absq)
    rinter3 = -rinter1*rinter1/rinter2
    gw      = -0.25/normerbs*exp(rinter3)     
    rinter4 = -2*rinter1/(rinter2*rinter2)
    ggw     = rinter4*gw  
    call qromb(psierbs,absq,2.,w)
    w = 0.25*w/normerbs
    if (q.lt.0.) then
       gw = -gw
    endif  
 endif

end subroutine geterbskernel1

subroutine geterbskernel2(q,w,gw,ggw)
 implicit none
real, intent(in)    :: q
real, intent(out)   :: w
real, intent(out)   :: gw
real, intent(out)   :: ggw
real :: absq
real, external :: psisupererbs,intsupererbs
real, parameter :: normerbs = 0.6034501612189381

 absq = abs(q)
 ggw = psisupererbs(absq)
 gw  = intsupererbs(absq)
 if (q.lt.0.) then
    gw = -gw
 endif
 if (absq.ge.2) then
    w   = 0.
 elseif(absq.le.tiny(0.)) then
    w = 0.5
 else
    call qromb(intsupererbs,absq,2.,w)
    w = -w
 endif 

end subroutine geterbskernel2

end module

!------------------------------------------------------------------
     SUBROUTINE qromb(func,a,b,ss)
     INTEGER JMAX,JMAXP,K,KM
     REAL a,b,func,ss,EPS
     EXTERNAL func
     PARAMETER (EPS=5.e-6, JMAX=20, JMAXP=JMAX+1, K=5, KM=K-1)
     INTEGER j
     REAL dss,h(JMAXP),s(JMAXP)
     h(1)=1.
     do 11 j=1,JMAX
       call trapzd(func,a,b,s(j),j)
       if (j.ge.K) then
         call polint(h(j-KM),s(j-KM),K,0.,ss,dss)
         if (abs(dss).le.EPS*abs(ss)) return
       endif
       s(j+1)=s(j)
       h(j+1)=0.25*h(j)
11    continue
     print*,'too many steps in qromb'
     end subroutine qromb

!------------------------------------------------------------------
     SUBROUTINE polint(xa,ya,n,x,y,dy)
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
         if(den.eq.0.) print*,'failure in polint'
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
     end subroutine polint

!------------------------------------------------------------------
     SUBROUTINE trapzd(func,a,b,s,n)
     INTEGER n
     REAL a,b,s,func
     INTEGER it,j
     REAL del,sum,tnm,x
     if (n.eq.1) then
       s=0.5*(b-a)*(func(a)+func(b))
     else
       it=2**(n-2)
       tnm=it
       del=(b-a)/tnm
       x=a+0.5*del
       sum=0.
       do j=1,it
         sum=sum+func(x)
         x=x+del
       enddo
       s=0.5*(s+(b-a)*sum/tnm)
     endif
     return
     END subroutine trapzd

!------------------------------------------------------------------
     real function psierbs(x)
     implicit none
     real :: x
     real :: rinter1,rinter2,rinter3

      if (x .le. tiny(0.)) then
         psierbs = 0.
      elseif (x .ge. 2.) then
         psierbs = 0.       
      else
         rinter1 = x - 1.
         rinter2 = x*(2.-x)
         rinter3 = - rinter1*rinter1/rinter2
         psierbs = exp(rinter3)
      endif

     end function psierbs
!------------------------------------------------------------------
     real function psisupererbs(x)
     implicit none
     real    :: x,xint
     integer :: fsym
     real    :: rinter1,rinter2,rinter3
     real,parameter:: normerbs = 0.6034501612189381

     if (x .lt. 1.) then
        xint = 2.-x
        fsym = -1
     else   
        xint = x
        fsym = 1
     endif

     if(abs(xint - 1.) .le. tiny(0.)) then  
        psisupererbs = 0.
     elseif (xint .ge. 2.) then
        psisupererbs = 0.       
     else   
        rinter1 = xint-1.5
        rinter2 = (xint - 1.)*(2.-xint)
        rinter3 = - rinter1*rinter1/rinter2
        psisupererbs = fsym*exp(rinter3)     
     endif
        psisupererbs = 0.5*psisupererbs/normerbs

     end function psisupererbs

!------------------------------------------------------------------
     real function intsupererbs(x)
     implicit none
     real :: x,xint,xtest,psisupererbs

     if (x.lt.1.) then
        xint = 2.- x
     else
        xint = x
     endif

     if (abs(xint - 1.) .le. tiny(0.)) then
        intsupererbs = -0.5
     elseif (xint .ge. 2.) then
        intsupererbs = 0.
     else
        xtest =psisupererbs(xint)
        call qromb(psisupererbs,xint,2.,intsupererbs)
        intsupererbs = -intsupererbs      
      endif
     end function intsupererbs

