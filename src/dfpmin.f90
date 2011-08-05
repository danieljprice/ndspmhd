!-------------------------------------------------------------------------
!
! Given a starting point p(1:ndim) that is a vector of length n
!  a minimisation is performed on a function func using its
!  gradient as calculated by the routine dfunc
!
! Convergence requirement on zeroing the gradient is input as gtol
! Returned quantities are:
!  p(1:ndim) (location of the minimum)
!  iter      (number of iterations that were performed)
!  fret      (minimum value of the function)
!
! The routine lnsrch is called to perform approximate line minimisations
!
! Parameters:
! NMAX is the maximum anticipated value of ndim
! ITMAX is the maximum allowed number of iterations
! STPMX is the scaled maximum step length allowed in line searches
! TOLX is the convergence criterion on x values
!
! From Numerical Recipes in Fortran, Press et al (1986)
!
! This version coded/modified slightly by 
! Daniel Price, University of Exeter, Jan 2008
!
!-------------------------------------------------------------------------
subroutine dfpmin(p,ndim,gtol,iter,fret,func,dfunc)
 implicit none
 integer, intent(in) :: ndim
 integer, intent(out) :: iter
 real, dimension(ndim), intent(inout) :: p
 real, intent(in) :: gtol
 real, intent(out) :: fret
 real, external :: func
 external :: dfunc
 integer, parameter :: NMAX = 3, ITMAX = 200
 real, parameter :: STPMX=100., eps=3.e-8, TOLX=4.*eps

 integer :: i,its,j
 logical :: check
 real :: den,fac,fad,fae,fp,stpmax,sum,sumdg,sumxi,temp,test
 real, dimension(NMAX) :: dg,g,hdg,pnew,xi
 real, dimension(NMAX,NMAX) :: hessin
 
 !--calculate starting function value and gradient
 fp = func(p)
 call dfunc(p,g)
 
 sum = 0.
 do i=1,ndim
    do j=1,ndim
       hessin(i,j) = 0.
    enddo
    hessin(i,i) = 1.
    xi(i) = -g(i)
    sum = sum + p(i)**2
 enddo
 stpmax=STPMX*max(sqrt(sum),float(ndim))
 
 !--main loop over the iterations
 do its=1,ITMAX
    iter = its
    call lnsrch(ndim,p,fp,g,xi,pnew,fret,stpmax,check,func)
    !--save function evaluation for next line search
    fp = fret
    do i=1,ndim
       xi(i) = pnew(i)-p(i) ! update line direction
       p(i) = pnew(i)       ! and the current point
    enddo
    test = 0.
    do i=1,ndim
       temp = abs(xi(i))/max(abs(p(i)),1.)
       if (temp.gt.test) test = temp
    enddo
    if (test.lt.TOLX) return
    do i=1,ndim
       dg(i) = g(i)
    enddo
    
    !--get new gradient
    call dfunc(p,g)
    !--test for convergence on new gradient
    test = 0.
    den=max(fret,1.)
    do i=1,ndim
       temp= abs(g(i))*max(abs(p(i)),1.)/den
       if (temp.gt.test) test = temp
    enddo
    if (test.lt.gtol) return
    
    !--compute difference of gradients
    do i=1,ndim
       dg(i) = g(i) - dg(i)
    enddo
    !--and difference times current matrix
    do i=1,ndim
       hdg(i) = 0.
       do j=1,ndim
          hdg(i) = hdg(i) + hessin(i,j)*dg(j)
       enddo
    enddo
    
    !--calculate dot products for the denominators
    fac = 0.
    fae = 0.
    sumdg=0.
    sumxi=0.
    do i=1,ndim
       fac = fac+dg(i)*xi(i)
       fae = fae+dg(i)*hdg(i)
       sumdg = sumdg + dg(i)**2
       sumxi = sumxi + xi(i)**2
    enddo
    
    !--skip update if fac not sufficiently positive
    if (fac*fac.gt.EPS*sumdg*sumxi) then
       fac=1./fac
       fad=1./fae
       do i=1,ndim ! the vector that makes BFGS different from DFP
          dg(i) = fac*xi(i) - fad*hdg(i)
       enddo
       do i=1,ndim
          do j=1,ndim
             hessin(i,j) = hessin(i,j) + fac*xi(i)*xi(j) &
                         - fad*hdg(i)*hdg(j) + fae*dg(i)*dg(j)
          enddo
       enddo
    endif
    
    !--now calculate the next direction to go
    do i=1,ndim
       xi(i) = 0.
       do j=1,ndim
          xi(i) = xi(i) - hessin(i,j)*g(j)
       enddo
    enddo
 enddo ! go back for another iteration

 print*,'WARNING! too many iterations in dfpmin'
 return
end subroutine dfpmin

!-------------------------------------------------------------------------
!
! Given an n dimensional point xold(1:ndim), the value of
! the function and gradient there, fold and g(1:ndim),
! finds a new point x(1:ndim) along the direction p from
! xold where the function func has decreased "sufficiently".
!
! The new function value is returned in f
!
! stpmax is an input quantity that limits the length of the
! steps so that you do not try to evaluate the function in
! regions where it is undefined or subject to overflow.
!
! p is usually the Newton direction. The output quantity
! check is false on a normal exit. It is true when x is too
! close to xold. In a minimisation algorithm, this usually
! signals convergence and can be ignored. However, in a 
! zero-finding algorithm, the calling program should check
! whether the convergence is spurious
! 
! Parameters:
! ALF ensures sufficient decrease in function value
! TOLX is the convergence criterion on delta x
!
! From Numerical Recipes in Fortran, Press et al (1986)
!
! This version coded/modified slightly by 
! Daniel Price, University of Exeter, Jan 2008
!
!-------------------------------------------------------------------------

subroutine lnsrch(ndim,xold,fold,g,p,x,f,stpmax,check,func)
 implicit none
 integer, intent(in) :: ndim
 logical, intent(out) :: check
 real, intent(in) :: fold,stpmax
 real, intent(out) :: f
 real, dimension(ndim), intent(inout) :: p
 real, dimension(ndim), intent(in) :: g,xold
 real, dimension(ndim), intent(out) :: x
 real, parameter :: ALF=1.e-4, TOLX=1.e-7
 real, external :: func

 integer :: i
 real :: a,alam,alam2,alamin,b,disc,f2,fold2,rhs1,rhs2,slope
 real :: sum,temp,test,tmplam
 
 check = .false.
 sum = 0.
 do i=1,ndim
    sum = sum + p(i)*p(i)
 enddo
 sum = sqrt(sum)
 
 !--scale if attempted step is too big
 if (sum.gt.stpmax) then
    do i=1,ndim
       p(i) = p(i)*stpmax/sum
    enddo
 endif
 
 slope = 0.
 do i=1,ndim
    slope = slope + g(i)*p(i)
 enddo
 
 test = 0.
 do i=1,ndim
    temp = abs(p(i))/max(abs(xold(i)),1.)
    if (temp.gt.test) test = temp
 enddo
 alamin = TOLX/test
 alam = 1.     !--always try full Newton step first
 
 alam2 = 0. ! unnecessary, but for caution
 f2 = 0.
 fold2 = 0.
!--iteration loop
1 continue
     do i=1,ndim
        x(i) = xold(i) + alam*p(i)
     enddo
     f = func(x)
     if (alam.lt.alamin) then
        do i=1,ndim
           x(i) = xold(i)
        enddo
        check = .true.
        return
     elseif (f.le.fold+ALF*alam*slope) then ! sufficient function decrease
        return
     else ! backtrack
        if (abs(alam - 1.).lt.tiny(alam)) then ! first time
           tmplam = -slope/(2.*(f-fold-slope))
        else
           rhs1 = f - fold - alam*slope
           rhs2 = f2 - fold2 - alam2*slope
           a = (rhs1/alam**2 - rhs2/alam2**2)/(alam-alam2)
           b = (-alam2*rhs1/alam**2 + alam*rhs2/alam2**2)/(alam-alam2)
           if (abs(a).lt.tiny(a)) then
              tmplam = -slope/(2.*b)
           else
              disc = b*b - 3.*a*slope
              if (disc.lt.0.) print*,' roundoff problem in lnsrch'
              tmplam = (-b + sqrt(disc))/(3.*a)
           endif
           if (tmplam.gt.0.5*alam) tmplam = 0.5*alam  ! lambda < 0.5 lambda1
        endif
     endif
     alam2 = alam
     f2 = f
     fold2 = fold
     alam = max(tmplam,0.1*alam)    ! lambda > 0.1 lambda 1
  goto 1
end subroutine lnsrch
