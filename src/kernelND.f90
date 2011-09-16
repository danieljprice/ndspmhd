!-----------------------------------------------------------------------------
! This module contains everything needed for the SPH kernel
!
! contains:
!  setkern : sets up kernel tables
!  interpolate_kernel : interpolation function from kernel tables
!  interpolate_kernels : interpolation function from kernel tables
!  interpolate_softening : interpolation function for softening kernel tables
!
!-----------------------------------------------------------------------------

module kernels
 implicit none
 integer, parameter :: ikern=4000    ! dimensions of kernel table
 integer :: ianticlump,neps
 real, parameter, private :: pi = 3.141592653589
 real, dimension(0:ikern) :: wij,grwij,wijalt,grwijalt
 real, dimension(0:ikern) :: wijdrag,grwijdrag,grgrwijdrag
 real :: dq2table,ddq2table,radkern2,radkern,eps
!--these variables for force softening only
 real, dimension(0:ikern) :: potensoft,fsoft,dphidh
!--these variables needed for plotting and analysis only (not in rates etc)
 real, dimension(0:ikern) :: grgrwij,grgrwijalt
 character(len=100) :: kernelname,kernelnamealt,kernelnamedrag
 
 public  :: wij,grwij,grgrwij
 public  :: setkern,interpolate_kernel,interpolate_kernels,interpolate_softening
 public  :: setkerndrag, interpolate_kerneldrag
 private :: setkerntable

contains

!-----------------------------------------------------------------
! This is the interface routine (public) -- calls setkern once only
!
!-----------------------------------------------------------------
subroutine setkern(ikernel,ndim,ierr)
 implicit none
 integer, intent(in)  :: ikernel, ndim
 integer, intent(out) :: ierr
!
!--setup kernel tables for primary kernel
!
 call setkerntable(ikernel,ndim,wij,grwij,grgrwij,kernelname,ierr)
 
end subroutine setkern

!----------------------------------------------------------------------
! This is the interface routine (public) -- calls setkerndrag once only
!
!----------------------------------------------------------------------
subroutine setkerndrag(ikernel,ndim,ierr)
 implicit none
 integer, intent(in)  :: ikernel, ndim
 integer, intent(out) :: ierr
!
!--setup kernel tables for primary kernel
!
 call setkerntable(ikernel,ndim,wijdrag,grwijdrag,grgrwijdrag,kernelnamedrag,ierr)
 
end subroutine setkerndrag

!-----------------------------------------------------------------
! This is another interface routine (public) -- calls setkern twice
! for both usual kernel and alternative kernel
!
!-----------------------------------------------------------------
subroutine setkernels(ikernel,ikernelalt,ndim,ierr1,ierr2)
 implicit none
 integer, intent(in)  :: ikernel,ikernelalt, ndim
 integer, intent(out) :: ierr1,ierr2
!
!--setup kernel tables for primary kernel
!
 call setkerntable(ikernel,ndim,wij,grwij,grgrwij,kernelname,ierr1)
!
!--setup kernel tables for alternative kernel
! 
 call setkerntable(ikernelalt,ndim,wijalt,grwijalt,grgrwijalt,kernelnamealt,ierr2)
 
end subroutine setkernels

!-----------------------------------------------------------------
! Sets up the tables for the kernel
! Returns kernel, and derivative.
!
! Default kernel is the cubic spline, but I have experimented
! with lots more.
!
!-----------------------------------------------------------------
subroutine setkerntable(ikernel,ndim,wkern,grwkern,grgrwkern,kernellabel,ierr)
 implicit none         !  define local variables
 integer, intent(in) :: ikernel, ndim
 real, intent(out), dimension(0:ikern) :: wkern,grwkern,grgrwkern
 character(len=*), intent(out) :: kernellabel
 integer, intent(out) :: ierr
 integer :: i,j,npower,n
 real :: q,q2,q4,q3,q5,q6,q7,cnormk,cnormkaniso
 real :: term1,term2,term3,term4,term
 real :: dterm1,dterm2,dterm3,dterm4
 real :: ddterm1,ddterm2,ddterm3,ddterm4,w0
 real :: alpha,beta,gamma,a,b,c,d,e,f,u,u2,qs,wdenom,wint

 cnormk = 0.0
 wkern = 0.
 grwkern = 0.
 grgrwkern = 0.
 fsoft = 1.
 potensoft = 0.
 ierr = 0

 select case(ikernel)

  case(2)
!
!--quartic spline
!  
    kernellabel = 'M_5 quartic'    

    radkern = 2.5
    radkern2 = radkern*radkern
    dq2table = radkern2/real(ikern)
    select case(ndim)
      case(1)
         cnormk = 1./24.
      case(2)
         cnormk = 96./(1199*pi)
      case(3)
         cnormk = 1./(20*pi)
    end select  
    do i=0,ikern
       q2 = i*dq2table
       q = sqrt(q2)
       if (q.lt.0.5) then
          wkern(i) = (2.5-q)**4 - 5.*(1.5-q)**4 + 10.*(0.5-q)**4
          grwkern(i) = -4.*((2.5-q)**3 - 5.*(1.5-q)**3 + 10*(0.5-q)**4)
          grgrwkern(i) = 12.*((2.5-q)**2 - 5.*(1.5-q)**2 + 10*(0.5-q)**2)
       elseif (q.lt.1.5) then
          wkern(i) = (2.5-q)**4 - 5.*(1.5-q)**4
          grwkern(i) = -4.*((2.5-q)**3 - 5.*(1.5-q)**3)
          grgrwkern(i) = 12.*((2.5-q)**2 - 5.*(1.5-q)**2)   
       elseif (q.lt.2.5) then
          wkern(i) = (2.5-q)**4
          grwkern(i) = -4.*((2.5-q)**3)
          grgrwkern(i) = 12.*((2.5-q)**2)
       else
          wkern(i) = 0.
          grwkern(i) = 0.
          grgrwkern(i) = 0.
       endif
    enddo 

  case(3)
!
!--this is the m_6 quintic spline (see e.g. morris 1996, phd thesis)
!
    kernellabel = 'M_6 quintic'  
    radkern = 3.0
    radkern2 = radkern*radkern
    dq2table = radkern2/real(ikern)
    select case(ndim)
      case(1) 
       cnormk = 1./120.
      case(2)
       cnormk = 7./(478*pi)
      case(3)
       cnormk = 1./(120.*pi)
    end select
    do i=0,ikern         
       q2 = i*dq2table
       q4 = q2*q2
       q = sqrt(q2)
       term1 = -5.*(3.-q)**4.
       if (q.lt.1.0) then
          wkern(i) = 66.-60.*q2 + 30.*q4 - 10.*q4*q
          grwkern(i) = term1 + 30*(2.-q)**4. - 75.*(1.-q)**4.
          grgrwkern(i) = 20.*(3.-q)**3. - 120.*(2.-q)**3. + 300.*(1.-q)**3.
       elseif ((q.ge.1.0).and.(q.lt.2.0)) then
          wkern(i) = (3.-q)**5. - 6.*(2.-q)**5.
          grwkern(i) = term1 + 30*(2.-q)**4.
          grgrwkern(i) = 20.*(3.-q)**3. - 120.*(2.-q)**3.
       elseif ((q.ge.2.0).and.(q.le.3.0)) then
          wkern(i) = (3.-q)**5.
          grwkern(i) = term1
          grgrwkern(i) = 20.*(3.-q)**3.
       else
          wkern(i) = 0.0
          grwkern(i) = 0.0
          grgrwkern(i) = 0.
       endif
    enddo

  case(5,6,7)
!
!--this is the do-it-yourself quintic kernel (general class of quintic splines)
!
    radkern = 2.0
    radkern2 = radkern*radkern
    dq2table = radkern2/real(ikern)

    gamma = 0.
    if (ikernel.eq.5) then
       beta = 0.5
       alpha = 1.7   !!1.4
       kernellabel = 'New quintic (1)'    
    elseif (ikernel.eq.6) then
       beta = 0.7
       alpha = 1.5
       kernellabel = 'New quintic (2)'    
    else
    !--match to cubic spline, ie w''(0) = -2
      kernellabel = 'Cubic-like quintic' 
      beta = 0.85
      !print*,' enter beta'
      !read*,beta
      print*,'beta = ',beta, ' calculating alpha'
      term1 = (beta+2.)*(beta**3 - 2.*beta**2 - 4.*beta + 128)
      alpha = -0.5*(4.*beta + beta**2 + 4. - sqrt(term1))/(beta + 2.)
   
      term2 = (beta-2.)*(beta+2.)*(beta-alpha)*(beta+alpha)
      dterm2 = (alpha+2)*(-alpha**2 + beta**2 - 4.)
      q = 2.*alpha*(-alpha**2 - 2.*alpha + beta**2 - 4. + sqrt(term2))/dterm2
      print*,' 3rd derivative zero at q = ',q    
    endif
    c =0.
    a =  (-radkern**4 + (radkern**2 + c*gamma**2)*beta**2)      &
        /(alpha**2*(alpha**2-beta**2))      
    b = -(radkern**4 + a*alpha**4 + c*gamma**4)/(beta**4)
    print*,'matching points = ',beta,alpha
    select case(ndim)
      case(1)
        cnormk = 3./(a*alpha**6 + b*beta**6 + c*gamma**6 + radkern**6)   ! for radkern = 2 and 1d
        print*,'1d cnormk = ',cnormk,' a,b = ',a,b
      case(2)
        cnormk = 42./(2.*pi*(a*alpha**7 + b*beta**7 + c*gamma**7 + radkern**7))
        print*,'2d cnormk = ',cnormk,' a,b = ',a,b,beta,alpha
      case default
       write(*,666)
       ierr = 1
       return
       !stop  
    end select
  
    do i=0,ikern         
      q2 = i*dq2table
      q = sqrt(q2)
      term1 = (radkern-q)**5
      term2 = (alpha-q)**5
      term3 = (beta-q)**5
      term4 = (gamma-q)**5
      dterm1 = -5*(radkern-q)**4
      dterm2 = -5*(alpha-q)**4
      dterm3 = -5*(beta-q)**4
      dterm4 = -5*(gamma-q)**4
      ddterm1 = 20*(radkern-q)**3
      ddterm2 = 20*(alpha-q)**3
      ddterm3 = 20*(beta-q)**3
      ddterm4 = 20*(gamma-q)**3
      if (q.lt.gamma) then
         wkern(i) = term1 + a*term2  + b*term3 + c*term4
         grwkern(i) = dterm1 + a*dterm2 + b*dterm3 + c*dterm4
         grgrwkern(i) = ddterm1 + a*ddterm2 + b*ddterm3 + c*ddterm4   
      elseif ((q.ge.gamma).and.(q.lt.beta)) then
         wkern(i) = term1 + a*term2  + b*term3
         grwkern(i) = dterm1 + a*dterm2 + b*dterm3
         grgrwkern(i) = ddterm1 + a*ddterm2 + b*ddterm3       
      elseif ((q.ge.beta).and.(q.lt.alpha)) then
         wkern(i) = term1 + a*term2
         grwkern(i) = dterm1 + a*dterm2
         grgrwkern(i) = ddterm1 + a*ddterm2
      elseif ((q.ge.alpha).and.(q.lt.radkern)) then
         wkern(i) = term1
         grwkern(i) = dterm1
         grgrwkern(i) = ddterm1
      else
         wkern(i) = 0.0
         grwkern(i) = 0.0
         grgrwkern(i) = 0.
      endif
    enddo

  case (8)
!
!--(1-r^2)^3
!   
    npower = 6
    write(kernellabel,"(a,i1)") '(1-r^2)^',npower     

    radkern = 2.0
    radkern2 = radkern*radkern
    dq2table = radkern2/real(ikern)
    select case(ndim)
      case(1)
         select case(npower)
         case(1)
            cnormk = 0.5*3./(2.*radkern**3)
         case(2)
            cnormk = 0.5*15./(8.*radkern**5)
         case(3)
            cnormk = 0.5*35./(16.*radkern**7)
         case(4)
            cnormk = 0.5*315./(128.*radkern**9)
         case(5)
            cnormk = 0.5*693./(256.*radkern**11)
         case(6)
            cnormk = 0.5*3003./(1024.*radkern**13)
         case(7)
            cnormk = 0.5*6435./(2048.*radkern**15)    
         case(8)
            cnormk = 0.5*109395./(32768.*radkern**17)    
         case default
            cnormk = 0.
         end select    
      case(2)
         cnormk = 3./(64.*pi)
      case(3)
         cnormk = 105./(4096.*pi)
    end select  
    do i=0,ikern
       q2 = i*dq2table
       q = sqrt(q2)
       if (q.lt.radkern) then
          wkern(i) = (radkern**2-q2)**npower
          grwkern(i) = -2.*npower*q*(radkern**2-q2)**(npower-1)
          grgrwkern(i) = (4.*npower*(npower-1)*q2*(radkern**2-q2)**(npower-2) &
                      - 2.*npower*(radkern**2-q2)**(npower-1))
       else
          wkern(i) = 0.
          grwkern(i) = 0.
          grgrwkern(i) = 0.
       endif
    enddo
    
   case (9)
!
!--(1-r)^n - a peaked kernel (deriv non-zero at origin) - truly awful
!   
    npower = 4
    write(kernellabel,"(a,i1)") '(2-r)^',npower     

    radkern = 2.0
    radkern2 = radkern*radkern
    dq2table = radkern2/real(ikern)
    select case(ndim)
      case(1)
         cnormk = 0.5*(npower+1)/radkern**(npower+1)
      case(2,3)
         write(*,666)
         ierr = 1
         return
         !stop 'normalisation const not defined in kernel'
    end select  
    do i=0,ikern
       q2 = i*dq2table
       q = sqrt(q2)
       if (q.lt.radkern) then
          wkern(i) = (radkern-q)**npower
          grwkern(i) = -npower*(radkern-q)**(npower-1)
          grgrwkern(i) = npower*(npower-1)*(radkern-q)**(npower-2)
       else
          wkern(i) = 0.
          grwkern(i) = 0.
          grgrwkern(i) = 0.
       endif
    enddo   

  case(10)
!
!--gaussian
!  
    kernellabel = 'Gaussian'    

    radkern = 10.0
    radkern2 = radkern*radkern
    dq2table = radkern2/real(ikern)
    select case(ndim)
      case(1)
         cnormk = 1./sqrt(pi)
      case(2)
         cnormk = 1./pi
      case(3)
         cnormk = 1./(pi*sqrt(pi))
    end select  
    do i=0,ikern
       q2 = i*dq2table
       q = sqrt(q2)
       if (q.lt.radkern) then
          wkern(i) = exp(-q2)
          grwkern(i) = -2.*q*wkern(i)
          grgrwkern(i) = -2.*q*grwkern(i) - 2.*wkern(i)
       else
          wkern(i) = 0.
          grwkern(i) = 0.
          grgrwkern(i) = 0.
       endif
    enddo

  case(11)
!
!--this is the usual spline based kernel modified for r/h < 2/3 to
!   prevent particles from clumping (see thomas & couchman '92)
!      
    kernellabel = 'Thomas & Couchman anti-clumping'
    radkern = 2.0      ! interaction radius of kernel
    radkern2 = radkern*radkern
    dq2table = radkern*radkern/real(ikern)
    select case(ndim)
      case(1)
        cnormk = 2./3. 
!        cnormk = 54./85. ! normalisation constant
      case(2)
        cnormk = 10./(7.*pi)
      case(3)
        cnormk = 1./pi
    end select
   
    do i=0,ikern
       q2 = i*dq2table
       q = sqrt(q2)
       if (q.lt.2./3.) then
          wkern(i) = 1. - 1.5*q2 + 0.75*q*q2
!          wkern(i) = 11./9. - q
          grwkern(i) = -1.
          grgrwkern(i) = 0.
       elseif (q.le.1.0) then
          wkern(i) = 1. - 1.5*q2 + 0.75*q*q2
          grwkern(i) = -3.*q+ 2.25*q2
          grgrwkern(i) = -3. + 4.5*q
       elseif (q.le.2.0) then
          wkern(i) = 0.25*(2.-q)**3.
          grwkern(i) = -0.75*(2.-q)**2.
          grgrwkern(i) = 1.5*(2.-q)
       else
          wkern(i) = 0.0
          grwkern(i) = 0.0
          grgrwkern(i) = 0.
       endif
    enddo
  
  case(12)
!
!--this is the squashed quintic spline from bonet & kulesegaram
!
    kernellabel = 'BK squashed quintic spline'    
   
    radkern = 2.0
    radkern2 = radkern*radkern
    dq2table = radkern2/real(ikern)
    select case(ndim)
      case(1)
       cnormk = 1./16.
      case default
       write(*,666)
       ierr = 1
       return
    end select
    do i=0,ikern
      q2 = i*dq2table
      q = sqrt(q2)
      if (q.lt.1.0) then
         wkern(i) = (2.-q)**5 - 16.*(1.-q)**5
         grwkern(i) = -5.*(2.-q)**4 + 80.*(1.-q)**4.
         grgrwkern(i) = 20.*(2.-q)**3 - 320.*(1.-q)**3.
      elseif ((q.ge.1.0).and.(q.le.2.0)) then
         wkern(i) = (2.-q)**5
         grwkern(i) = -5.*(2.-q)**4
         grgrwkern(i) = 20.*(2.-q)**3
      else
         wkern(i) = 0.0
         grwkern(i) = 0.0
         grgrwkern(i) = 0.
      endif
    enddo  
    
  case(13)
!
!--this is the lucy kernel (to 2h not to 1h)
!
    kernellabel = 'Lucy kernel'    
   
    radkern = 2.0
    radkern2 = radkern*radkern
    dq2table = radkern2/real(ikern)
    select case(ndim)
      case(1)
       cnormk = 5./8. ! normalisation is probably wrong
      case(2)
       cnormk = 5./(4.*pi)
      case(3)
       cnormk = 105/(128.*pi)
      case default
       write(*,666)
       ierr = 1
       return
      end select
    do i=0,ikern
      q2 = i*dq2table
      q = sqrt(q2)
      if (q.lt.2.0) then
         wkern(i) = (1.+1.5*q)*(1.-0.5*q)**3
         grwkern(i) = 1.5*(1.-0.5*q)**3 - 1.5*(1.+1.5*q)*(1.-0.5*q)**2
         grgrwkern(i) = -4.5*(1.-0.5*q)**2 + 1.5*(1.+1.5*q)*(1.-0.5*q)
      else
         wkern(i) = 0.0
         grwkern(i) = 0.0
         grgrwkern(i) = 0.
      endif
    enddo

  case(14)
!
!--this is a modification of the cubic spline
!
    kernellabel = 'Peaked cubic spline'    
   
    radkern = 2.0
    radkern2 = radkern*radkern
    dq2table = radkern2/real(ikern)
    select case(ndim)
      case(1)
       cnormk = 1./8.
      case(2)
       cnormk = 5./(16.*pi)
      case(3)
       cnormk = 15./(64.*pi)
      case default
       write(*,666)
       ierr = 1
       return
    end select
    do i=0,ikern
      q2 = i*dq2table
      q = sqrt(q2)
      if (q.lt.2.0) then
         wkern(i) = (2.-q)**3
         grwkern(i) = -3.*(2.-q)**2
         grgrwkern(i) = 6.*(2.-q)
      else
         wkern(i) = 0.0
         grwkern(i) = 0.0
         grgrwkern(i) = 0.
      endif
    enddo
    
  case(15)
!
!--this is a modification of the cubic spline
!
    kernellabel = 'Peaked cubic spline 2'    
   
    radkern = 2.0
    radkern2 = radkern*radkern
    dq2table = radkern2/real(ikern)
    select case(ndim)
      case(1)
       cnormk = 4./7.
      case(2)
       cnormk = 4./(3.*pi)
      case(3)
       cnormk = 30./(31.*pi)
      case default
       write(*,666)
       ierr = 1
       return
    end select
    do i=0,ikern
      q2 = i*dq2table
      q = sqrt(q2)
      if (q.lt.1.0) then
         wkern(i) = 0.25*(2.-q)**3 - 0.5*(1.-q)**3
         grwkern(i) = -0.75*(2.-q)**2 + 1.5*(1.-q)**2
         grgrwkern(i) = 1.5*q
      elseif(q.lt.2.0) then
         wkern(i) = 0.25*(2.-q)**3
         grwkern(i) = -0.75*(2.-q)**2
         grgrwkern(i) = 3. - 1.5*q
      else
         wkern(i) = 0.0
         grwkern(i) = 0.0
         grgrwkern(i) = 0.
      endif
    enddo

  case(16)
!
!--another version of the peaked cubic spline
!
    alpha = 1.0
    print*,'enter alpha,'
    read*,alpha
    write(kernellabel,"(a,f6.2)") 'peaked cubic spline alpha = ',alpha
!!    kernellabel = 'peaked cubic spline 3'    
   
    radkern = 2.0
    radkern2 = radkern*radkern
    dq2table = radkern2/real(ikern)
    select case(ndim)
      case(1)
       cnormk = 1./(1. + 0.25/alpha + 1.5*alpha - alpha**2)
      case(2)
       cnormk = -20.*alpha/(pi*(10.*alpha**3 - 20.*alpha**2 - alpha - 4.))
      case(3)
       cnormk = -30.*alpha/(pi*(20.*alpha**3 - 45.*alpha**2 + 4.*alpha - 10.))
      case default
       write(*,666)
       ierr = 1
       return
    end select
    do i=0,ikern
      q2 = i*dq2table
      q = sqrt(q2)
      if (q.lt.alpha) then
         wkern(i) = 0.25*(2.-q)**3 - 0.5/alpha*(alpha-q)**3
         grwkern(i) = -0.75*(2.-q)**2 + 1.5/alpha*(alpha-q)**2
         grgrwkern(i) = 1.5*(2.-q) - 3./alpha*(alpha-q)
      elseif (q.lt.2.0) then
         wkern(i) = 0.25*(2.-q)**3
         grwkern(i) = -0.75*(2.-q)**2
         grgrwkern(i) = 1.5*(2.-q)
      else
         wkern(i) = 0.0
         grwkern(i) = 0.0
         grgrwkern(i) = 0.
      endif
    enddo
    
  case(17)
!
!--cosine**n kernel
!
    kernellabel = 'cosine**n'
   
    radkern = 2.0
    radkern2 = radkern*radkern
    dq2table = radkern2/real(ikern)
    n = 5
    !print*,'enter n'
    !read*,n
    select case(ndim)
      case(1)
       select case(n)
       case(3)
          cnormk = 1./1.694930
       case(4)
          cnormk = 1./1.5
       case(5)
          cnormk = 1./1.359516
       case(6)
          cnormk = 1./1.249811
       case(7)
          cnormk = 1./1.163532
       case(8)
          cnormk = 1./1.095270
       case default
          cnormk = 1.
       end select
      case default
       write(*,666)
       ierr = 1
       return
    end select

    do i=0,ikern
      q2 = i*dq2table
      q = sqrt(q2)
      if (q.lt.radkern) then
         wkern(i) = cos(0.25*pi*q)**n !!sin(0.5*pi*q)/q
         grwkern(i) = -0.25*pi*n*(cos(0.25*pi*q))**(n-1)*sin(0.25*pi*q)
         grgrwkern(i) = 1./16.*pi*pi*n* &
           ((n-1)*(cos(0.25*pi*q))**(n-2)*sin(0.25*pi*q)**2 - &
            cos(0.25*pi*q)**n)
      else
         wkern(i) = 0.0
         grwkern(i) = 0.0
         grgrwkern(i) = 0.
      endif
    enddo    

  case(18)
!
!--Dirichlet kernel (sin(x)/x)
!
    kernellabel = 'Dirichlet'
   
    radkern = 2.0
    radkern2 = radkern*radkern
    dq2table = radkern2/real(ikern)
    select case(ndim)
      case(1)
        cnormk = 1./5.194046
      case default
       write(*,666)
       ierr = 1
       return
    end select

    do i=0,ikern
      q2 = i*dq2table
      q = sqrt(q2)
      if (q.lt.radkern) then
         if (q.gt.0.) then
            wkern(i) = sin(0.75*pi*q)/q + 0.5
            grwkern(i) = 0.75*pi*cos(0.75*pi*q)/q - sin(0.75*pi*q)/q2
         else
            wkern(i) = 0.75*pi + 0.5
            grwkern(i) = 0.
         endif
         grgrwkern(i) = 0.
      else
         wkern(i) = 0.0
         grwkern(i) = 0.0
         grgrwkern(i) = 0.
      endif
    enddo    

  case(19)
!
!--sin**2/x**2
!
    kernellabel = 'sin^2/x^2'
   
    radkern = 2.0
    radkern2 = radkern*radkern
    dq2table = radkern2/real(ikern)
    select case(ndim)
      case(1)
        cnormk = 1./4.35
      case default
       write(*,666)
       ierr = 1
       return
    end select

    do i=0,ikern
      q2 = i*dq2table
      q = sqrt(q2)
      if (q.lt.radkern) then
         term = 0.5*pi*q
         if (q.gt.0.) then
            wkern(i) = (sin(term)**2)/q2
            grwkern(i) = (pi*sin(term)*cos(term))/q2 - 2.*wkern(i)/q
            grgrwkern(i) = 0.5*pi*pi*(cos(term)**2 - sin(term)**2)/q2 &
                         - 2.*pi*sin(term)*cos(term)/(q2*q) &
                         + 6.*sin(term)**2/(q2*q2)
         else
            wkern(i) = 0.5*pi
            grwkern(i) = 0.
            grgrwkern(i) = 0.
         endif
      else
         wkern(i) = 0.0
         grwkern(i) = 0.0
         grgrwkern(i) = 0.
      endif
    enddo

  case(20)
!
!--Jackson-Feyer de la Vallee Poussin Kernel
!
    kernellabel = '[sin(q)/q]**4'
   
    radkern = 2.0
    radkern2 = radkern*radkern
    dq2table = radkern2/real(ikern)
    select case(ndim)
      case(1)
        cnormk = 1./8.097925
      case default
       write(*,666)
       ierr = 1
       return
    end select

    do i=0,ikern
      q2 = i*dq2table
      q = sqrt(q2)
      if (q.lt.radkern) then
         term = 0.5*pi*q
         if (q.gt.0.) then
            wkern(i) = (sin(term)**4)/(q2*q2)
            grwkern(i) = 2./(q2*q2)*(sin(term)**3*cos(term)*pi - 2.*sin(term)**4/q)
            grgrwkern(i) = 1./(q2*q2)*(3.*sin(term)**2*cos(term)**2*pi**2 &
                        -16.*sin(term)**3*cos(term)*pi/q - sin(term)**4*pi**2 &
                        + 20.*(sin(term)**4)/q2)
         else
            wkern(i) = (0.5*pi)**4
            grwkern(i) = 0.
            grgrwkern(i) = 0.
         endif
      else
         wkern(i) = 0.0
         grwkern(i) = 0.0
         grgrwkern(i) = 0.
      endif
    enddo
    grgrwkern(0) = grgrwkern(1)

  case(21)
!
!--Jackson kernel
!
    kernellabel = 'Jackson kernel'
   
    radkern = 2.0
    radkern2 = radkern*radkern
    dq2table = radkern2/real(ikern)
    select case(ndim)
      case(1)
        cnormk = 1./24.02428 !!/113.3185 !!/136.77!!/24.02428
      case default
       write(*,666)
       ierr = 1
       return
    end select

    npower = 2
    do i=0,ikern
      q2 = i*dq2table
      q = sqrt(q2)
      if (q.lt.radkern) then
         term = 0.5*pi*q
         if (q.gt.0.) then
            wkern(i) = (sin(term)/sin(term/npower))**4
            grwkern(i) = (2.*pi/sin(term/npower)**4)*(sin(term)**3*cos(term)  &
                     - (sin(term)**4)*cos(term/npower)/(npower*sin(term/npower)))
            grgrwkern(i) = pi*pi/sin(term/npower)**4* &
               (3.*sin(term)**2*cos(term)**2 &
              - 8.*sin(term)**3*cos(term)*cos(term/npower)/(npower*sin(term/npower)) &
              - sin(term)**4  &
              + 5.*sin(term)**4*cos(term/npower)**2/(npower*sin(term/npower))**2 &
              + sin(term)**4/(npower**2))
         else
            wkern(i) = (npower)**4
            grwkern(i) = 0.
            grgrwkern(i) = 0.
         endif
      else
         wkern(i) = 0.0
         grwkern(i) = 0.0
         grgrwkern(i) = 0.
      endif
    enddo
    grgrwkern(0) = grgrwkern(1)    

  case(22)
!
!--Q-gaussian from Fulk & Quinn 1996
!  
    kernellabel = 'Q-Gaussian'    

    radkern = 2.0
    radkern2 = radkern*radkern
    dq2table = radkern2/real(ikern)
    select case(ndim)
      case(1)
         cnormk = 0.643998 ! 0.7764 from Mathematica
      case default
       write(*,666)
       ierr = 1
       return
    end select  
    do i=0,ikern
       q2 = i*dq2table
       q = sqrt(q2)
       if (q.lt.radkern) then
          wkern(i) = (1.-0.25*q2)*exp(-q2)
          grwkern(i) = 0.5*q*(q2 - 5.)*exp(-q2)
          grgrwkern(i) = -0.5*exp(-q2)*(2.*q2*q2 - 13.*q2 + 5.)
       else
          wkern(i) = 0.
          grwkern(i) = 0.
          grgrwkern(i) = 0.
       endif
    enddo

  case(23)
!
!--this is the m_6 quintic spline (S. Rosswog)
!
    kernellabel = 'Linear Core M_6 quintic'  
    radkern = 3.0
    radkern2 = radkern*radkern
    dq2table = radkern2/real(ikern)
    w0 = 80.7721808180
    qs = 0.7592984881
    a  = 55.2040173925
    select case(ndim)
      case(1) 
!       cnormk = 1./120.
       cnormk = 0.0079541232
      case(2)
       !cnormk = 7./(478*pi)
       cnormk = 0.0046023237
      case(3)
       !cnormk = 1./(120.*pi)
       cnormk = 0.0026427689
    end select
    do i=0,ikern         
       q2 = i*dq2table
       q4 = q2*q2
       q = sqrt(q2)
       term1 = -5.*(3.-q)**4.
       if (q.lt.qs) then
          wkern(i) = w0 - a*q
          grwkern(i) = -a
          grgrwkern(i) = 0.       
       elseif (q.lt.1.0) then
          wkern(i) = 66.-60.*q2 + 30.*q4 - 10.*q4*q
          grwkern(i) = term1 + 30*(2.-q)**4. - 75.*(1.-q)**4.
          grgrwkern(i) = 20.*(3.-q)**3. - 120.*(2.-q)**3. + 300.*(1.-q)**3.
       elseif ((q.ge.1.0).and.(q.lt.2.0)) then
          wkern(i) = (3.-q)**5. - 6.*(2.-q)**5.
          grwkern(i) = term1 + 30*(2.-q)**4.
          grgrwkern(i) = 20.*(3.-q)**3. - 120.*(2.-q)**3.
       elseif ((q.ge.2.0).and.(q.le.3.0)) then
          wkern(i) = (3.-q)**5.
          grwkern(i) = term1
          grgrwkern(i) = 20.*(3.-q)**3.
       else
          wkern(i) = 0.0
          grwkern(i) = 0.0
          grgrwkern(i) = 0.
       endif
    enddo

  case(24)
!
!--this is the m_6 quintic spline (S. Rosswog)
!
    kernellabel = 'Quartic Core M_6 quintic'  
    radkern = 3.0
    radkern2 = radkern*radkern
    dq2table = radkern2/real(ikern)
    qs = 0.7592984881
    a  = 11.017537
    b  = -38.111922
    c  = -16.619585
    d  = 69.785768
    select case(ndim)
      case(1) 
       cnormk = 8.245880e-3
      case(2)
       cnormk = 4.649647e-3
      case(3)
       cnormk = 2.650839e-3
    end select
    do i=0,ikern
       q2 = i*dq2table
       q4 = q2*q2
       q = sqrt(q2)
       term1 = -5.*(3.-q)**4
       if (q.lt.qs) then
          wkern(i) = a*q4 + b*q2 + c*q + d
          grwkern(i) = 4.*a*q2*q + 2.*b*q + c
          grgrwkern(i) = 12.*a*q2 + 2.*b
       elseif (q.lt.1.0) then
          wkern(i) = 66.-60.*q2 + 30.*q4 - 10.*q4*q
          grwkern(i) = term1 + 30*(2.-q)**4. - 75.*(1.-q)**4.
          grgrwkern(i) = 20.*(3.-q)**3. - 120.*(2.-q)**3. + 300.*(1.-q)**3.
       elseif ((q.ge.1.0).and.(q.lt.2.0)) then
          wkern(i) = (3.-q)**5. - 6.*(2.-q)**5.
          grwkern(i) = term1 + 30*(2.-q)**4.
          grgrwkern(i) = 20.*(3.-q)**3. - 120.*(2.-q)**3.
       elseif ((q.ge.2.0).and.(q.le.3.0)) then
          wkern(i) = (3.-q)**5.
          grwkern(i) = term1
          grgrwkern(i) = 20.*(3.-q)**3.
       else
          wkern(i) = 0.0
          grwkern(i) = 0.0
          grgrwkern(i) = 0.
       endif
    enddo

  case(25)
!
!--Linear Quartic kernel from Valcke et al. (2010), as used by Rosswog (2010)
!  
    kernellabel = 'Linear-Quartic'

    radkern = 2.0
    radkern2 = radkern*radkern
    dq2table = radkern2/real(ikern)
    qs = 0.3
    alpha = 1./(qs**3 - 3.*qs**2 + 3.*qs - 1.)
    a = 0.5*alpha
    b = -alpha*(1. + qs)
    c = 3.*alpha*qs
    d = -alpha*(3.*qs - 1.)
    e = alpha*(2.*qs - 1.)/2.
    f = a*qs**4 + b*qs**3 + c*qs*qs + d*qs + e + qs
    select case(ndim)
      case(1)
         cnormk = 2.235
      case(2)
         cnormk = 2.962
      case(3)
!         cnormk = 4.*pi*((f/3.*qs**3 - 0.25*qs**4) &
!                + a/7.*qs**7 + b/6.*qs**6 + c/5.*qs**5 + 0.25*d*qs**4 + e/3.*qs**3)
!         print*,cnormk,3.947
         cnormk = 3.947
      case default
       write(*,666)
       ierr = 1
       return
    end select  
    cnormk = cnormk/radkern**ndim
    do i=0,ikern
       q2 = i*dq2table
       q = sqrt(q2)
       u = 0.5*q
       u2 = u*u
       if (u.lt.qs) then
          wkern(i) = f - u
          grwkern(i) = -0.5
          grgrwkern(i) = 0.
       elseif (u.le.1.) then
          wkern(i) = a*u2*u2 + b*u2*u + c*u2 + d*u + e
          grwkern(i) = 2.*a*u2*u + 1.5*b*u2 + 1.*c*u + 0.5*d
          grgrwkern(i) = 3.*a*u2 + 1.5*b*u + 0.5*c
       else
          wkern(i) = 0.
          grwkern(i) = 0.
          grgrwkern(i) = 0.
       endif
    enddo

  case(30)
!
!--these are the Ferrer's spheres from Dehnen (2001)
!
    kernellabel = 'Ferrers n=3 sphere'    
   
    radkern = 2.0
    radkern2 = radkern*radkern
    dq2table = radkern2/real(ikern)
    select case(ndim)
      case(1)
       cnormk = 35./4096.
      case(2)
       cnormk = 1./(64.*pi)
      case(3)
       cnormk = 315./(32768.*pi)
      case default
       write(*,666)
       ierr = 1
       return
    end select
    do i=0,ikern
      q2 = i*dq2table
      q = sqrt(q2)
      if (q.lt.2.0) then
         wkern(i) = (4.-q2)**3
         grwkern(i) = -6.*q*(4.-q2)**2
         grgrwkern(i) = 24.*q2*(4.-q2) - 6.*(4.-q2)**2
      else
         wkern(i) = 0.0
         grwkern(i) = 0.0
         grgrwkern(i) = 0.
      endif
    enddo
    
  case(31)
!
!--these are the Ferrer's spheres from dehnen (2001)
!
    kernellabel = 'Ferrers n=6 sphere'    
   
    radkern = 2.0
    radkern2 = radkern*radkern
    dq2table = radkern2/real(ikern)
    select case(ndim)
      case(1)
       cnormk = 3003./4096.
      case(2)
       cnormk = 7./(4.*pi)
      case(3)
       cnormk = 45045./(32768.*pi)
      case default
       write(*,666)
       ierr = 1
       return
    end select
    do i=0,ikern
      q2 = i*dq2table
      q = sqrt(q2)
      if (q.lt.2.0) then
         term = 1.-0.25*q2
         wkern(i) = term**6
         grwkern(i) = -3.*q*term**5
         grgrwkern(i) = 7.5*q2*term**4 - 3.*term**5
      else
         wkern(i) = 0.0
         grwkern(i) = 0.0
         grgrwkern(i) = 0.
      endif
    enddo

  case(32)
!
!--these are the Ferrer's spheres from dehnen (2001)
!
    kernellabel = 'Ferrers n=7 sphere'    
   
    radkern = 2.0
    radkern2 = radkern*radkern
    dq2table = radkern2/real(ikern)
    select case(ndim)
      case(1)
       cnormk = 6435./8192.
      case(2)
       cnormk = 2./pi
      case(3)
       cnormk = 109395./(65536.*pi)
      case default
       write(*,666)
       ierr = 1
       return
    end select
    do i=0,ikern
      q2 = i*dq2table
      q = sqrt(q2)
      if (q.lt.2.0) then
         term = 1.-0.25*q2
         wkern(i) = term**7
         grwkern(i) = -3.5*q*term**6
         grgrwkern(i) = 0.5*21.*q2*term**5 - 3.5*term**6
      else
         wkern(i) = 0.0
         grwkern(i) = 0.0
         grgrwkern(i) = 0.
      endif
    enddo

  case(33)
!
!--these are the Ferrer's spheres from dehnen (2001)
!
    kernellabel = 'Ferrers n=8 sphere'    
   
    radkern = 2.0
    radkern2 = radkern*radkern
    dq2table = radkern2/real(ikern)
    select case(ndim)
      case(1)
       cnormk = 109395./131072.
      case(2)
       cnormk = 9./(4.*pi)
      case(3)
       cnormk = 2078505./(1048576.*pi)
      case default
       write(*,666)
       ierr = 1
       return
    end select
    do i=0,ikern
      q2 = i*dq2table
      q = sqrt(q2)
      if (q.lt.2.0) then
         term = 1.-0.25*q2
         wkern(i) = term**8
         grwkern(i) = -4.*q*term**7
         grgrwkern(i) = 14.*q2*term**6 - 4.*term**7
      else
         wkern(i) = 0.0
         grwkern(i) = 0.0
         grgrwkern(i) = 0.
      endif
    enddo

  case(34)
!
!--Cauchy kernel
!
    kernellabel = 'Cauchy'    
   
    radkern = 20.0
    radkern2 = radkern*radkern
    dq2table = radkern2/real(ikern)
    select case(ndim)
      case(1)
       cnormk = 1./pi
      case default
       write(*,666)
       ierr = 1
       return
    end select
    do i=0,ikern
      q2 = i*dq2table
      q = sqrt(q2)
      if (q.lt.radkern) then
         wkern(i) = 1./(1. + q2)
         grwkern(i) = -1./(1. + q2)**2*2.*q
         grgrwkern(i) = -2./(1. + q2)**2 + 4.*q/(1. + q2)**3*2.*q
      else
         wkern(i) = 0.0
         grwkern(i) = 0.0
         grgrwkern(i) = 0.
      endif
    enddo

  case(41)
!
!--particle splitting kernel (cubic spline gradient  * q)
!   
    kernellabel = 'particle splitting'    
  
    radkern = 2.0      ! interaction radius of kernel
    radkern2 = radkern*radkern
    dq2table = radkern*radkern/real(ikern)    
    select case(ndim)
      case(1)
        cnormk = 0.66666666666
      case(2)
        cnormk = 10./(7.*pi)
      case(3)
        cnormk = 1./(3.*pi)
    end select
!
!--setup kernel table
!   
    do i=0,ikern
       q2 = i*dq2table
       q = sqrt(q2)
       q4 = q2*q2
       if (q.lt.1.0) then
          wkern(i) = q*(3.*q - 2.25*q2)
          grwkern(i) = 6.*q - 3.*2.25*q2
          grgrwkern(i) = 6. - 6.*2.25*q
       elseif ((q.ge.1.0).and.(q.le.2.0)) then
          wkern(i) = q*0.75*(2.-q)**2
          grwkern(i) = 0.75*(2.-q)**2 - 1.5*q*(2.-q)
          grgrwkern(i) = -1.5*(2.-q) - 1.5*(2.-q) + 1.5*q
       else
          wkern(i) = 0.0
          grwkern(i) = 0.0
          grgrwkern(i) = 0.
       endif
    enddo

  case(42)
!
!--Double hump kernel from Fulk & Quinn (1996)
!   
    kernellabel = 'double cubic spline'
  
    radkern = 2.0      ! interaction radius of kernel
    radkern2 = radkern*radkern
    dq2table = radkern*radkern/real(ikern)    
    select case(ndim)
      case(1)
        cnormk = 2.0
      case(2)
        cnormk = 70./(31.*pi)
      case(3)
        cnormk = 10./(9.*pi)
    end select
!
!--setup kernel table
!   
    do i=0,ikern
       q2 = i*dq2table
       q = sqrt(q2)
       q4 = q2*q2
       if (q.lt.1.0) then
          wkern(i) = q2 - 1.5*q4 + 0.75*q4*q
          grwkern(i) = 2.*q - 6.*q2*q + 3.75*q4
          grgrwkern(i) = 2. - 18.*q2 + 15.*q2*q
       elseif ((q.ge.1.0).and.(q.le.2.0)) then
          wkern(i) = 0.25*q2*(2.-q)**3
          grwkern(i) = 0.5*q*(2.-q)**3 - 0.75*q2*(2.-q)**2
          grgrwkern(i) = 0.5*(2.-q)**3 - 3.*q*(2.-q)**2 + 1.5*q2*(2.-q)
       else
          wkern(i) = 0.0
          grwkern(i) = 0.0
          grgrwkern(i) = 0.
       endif
    enddo

  case(43)
!
!--double hump version of the m_6 quintic spline
!
    kernellabel = 'double hump M_6 quintic'  
    radkern = 3.0
    radkern2 = radkern*radkern
    dq2table = radkern2/real(ikern)
    select case(ndim)
      case(1) 
       cnormk = 1./60.
      case(2)
       cnormk = 42./(2771.*pi)
      case(3)
       cnormk = 1./(168.*pi)
    end select
    do i=0,ikern         
       q2 = i*dq2table
       q4 = q2*q2
       q = sqrt(q2)
       if (q.lt.1.0) then
          wkern(i) = q2*(66.-60.*q2 + 30.*q4 - 10.*q4*q)
          grwkern(i) = 2.*q*(-35.*q4*q + 90.*q4 - 120.*q2 + 66.)
          grgrwkern(i) = -12.*(35.*q4*q - 75.*q4 + 60.*q2 - 11.)
       elseif ((q.ge.1.0).and.(q.lt.2.0)) then
          wkern(i) = q2*((3.-q)**5. - 6.*(2.-q)**5.)
          grwkern(i) = q*(35.*q4*q - 270.*q4 + 750.*q2*q - 840.*q2 + 225.*q + 102)
          grgrwkern(i) = 6.*(35.*q4*q - 225.*q4 + 500.*q2*q - 420.*q2 + 75.*q + 17)
       elseif ((q.ge.2.0).and.(q.le.3.0)) then
          wkern(i) = q2*(3.-q)**5.
          grwkern(i) = -q*(7.*q - 6.)*(3.-q)**4
          grgrwkern(i) = 6.*(3.-q)**3*(7.*q2 - 12.*q + 3.)
       else
          wkern(i) = 0.0
          grwkern(i) = 0.0
          grgrwkern(i) = 0.
       endif
    enddo

  case(44)
!
!--Another double hump kernel from Fulk & Quinn (1996)
!   
    kernellabel = 'double hump Lucy/Wendland 4'    
  
    radkern = 2.0      ! interaction radius of kernel
    radkern2 = radkern*radkern
    dq2table = radkern*radkern/real(ikern)    
    select case(ndim)
      case(1)
        cnormk = 0.102539
      case(2)
        cnormk = 10./(7.*pi)
      case(3)
        cnormk = 1./(3.*6.5*pi)
    end select
!
!--setup kernel table
!   
    do i=0,ikern
       q2 = i*dq2table
       q = sqrt(q2)
       q4 = q2*q2
       if (q.lt.2.0) then
          term = (2.+3.*q)*(2.-q)**3
          wkern(i) = q2*term
          grwkern(i) = 2.*q*term + q2*(3.*(2.-q)**3 - 3.*(2.+3.*q)*(2.-q)**2)
          grgrwkern(i) = 2.*term + 2.*q*(3.*(2.-q)**3 - 3.*(2.+3.*q)*(2.-q)**2) &
                       + q2*(-18.*(2.-q)**2 + 6.*(2.+3.*q)*(2.-q))&
                       + 2.*q*(3.*(2.-q)**3 - 3.*(2.+3.*q)*(2.-q)**2)
       else
          wkern(i) = 0.0
          grwkern(i) = 0.0
          grgrwkern(i) = 0.
       endif
    enddo

  case(45)
!
!--double hump Gaussian
!   
    kernellabel = 'double Gaussian'    
  
    radkern = 10.0      ! interaction radius of kernel
    radkern2 = radkern*radkern
    dq2table = radkern*radkern/real(ikern)    
    cnormk = 2./real(ndim)*pi**(-ndim/2.)
!
!--setup kernel table
!   
    do i=0,ikern
       q2 = i*dq2table
       q = sqrt(q2)
       q4 = q2*q2
       if (q.lt.radkern) then
          wkern(i) = q2*exp(-q2)
          grwkern(i) = -2.*exp(-q2)*q*(q2 - 1.)
          grgrwkern(i) = 2.*exp(-q2)*(2.*q2*q2 - 5.*q2 + 1.)
       else
          wkern(i) = 0.0
          grwkern(i) = 0.0
          grgrwkern(i) = 0.
       endif
    enddo

  case(46)
!
!--particle splitting kernel (cubic spline gradient  * q)
!   
    kernellabel = 'cubic spline gradient-as-kernel'    
  
    radkern = 2.0      ! interaction radius of kernel
    radkern2 = radkern*radkern
    dq2table = radkern*radkern/real(ikern)    
    select case(ndim)
      case(1)
        cnormk = 0.5
      case(2)
        cnormk = 10./(7.*pi)
      case(3)
        cnormk = 1./(3.*pi)/0.9245
    end select
!
!--setup kernel table
!   
    do i=0,ikern
       q2 = i*dq2table
       q = sqrt(q2)
       q4 = q2*q2
       if (q.lt.1.0) then
          wkern(i) = (3.*q - 2.25*q2)
          grwkern(i) = 3. - 4.5*q
          grgrwkern(i) = -4.5
       elseif ((q.ge.1.0).and.(q.le.2.0)) then
          wkern(i) = 0.75*(2.-q)**2
          grwkern(i) = -1.5*(2.-q)
          grgrwkern(i) = 1.5
       else
          wkern(i) = 0.0
          grwkern(i) = 0.0
          grgrwkern(i) = 0.
       endif
    enddo

  case(51)
!
!--cubic spline with vanishing second moment (using monaghan 1992 method)
!
    kernellabel = 'Higher order cubic spline'    
   
    radkern = 2.0
    radkern2 = radkern*radkern
    dq2table = radkern2/real(ikern)
    select case(ndim)
      case(1)
       a = 0.9
       cnormk = 2./(3.*a - 1.)
      case(3)
       a = 85./63.
       cnormk = 1./pi * (630./283.)
      case default
       write(*,666)
       ierr = 1
       return
    end select
    do i=0,ikern
      q2 = i*dq2table
      q = sqrt(q2)
      if (q.lt.1.0) then
         wkern(i) = (a - q2)*(1. - 1.5*q2 + 0.75*q2*q)
         grwkern(i) = -2.*q*(1. - 1.5*q2 + 0.75*q2*q) + (a -q2)*(-3.*q + 9./4.*q2)
         grgrwkern(i) = -2. + 18.*q2 - 15.*q2*q - 3.*a + 4.5*a*q
      elseif (q.lt.2.0) then
         wkern(i) = 0.25*(a - q2)*(2.-q)**3
         grwkern(i) = -0.5*q*(2.-q)**3 - 0.75*(a-q2)*(2.-q)**2
         grgrwkern(i) = -4. + 3.*a + 18.*q - 18.*q2 - 1.5*a*q + 5.*q2*q
      else
         wkern(i) = 0.0
         grwkern(i) = 0.0
         grgrwkern(i) = 0.
      endif
    enddo

  case(52)
!
!--quartic spline with vanishing second moment (using monaghan 1992 method)
!
    kernellabel = 'Higher order quartic spline'    
   
    radkern = 2.5
    radkern2 = radkern*radkern
    dq2table = radkern2/real(ikern)
    select case(ndim)
      case(1)
       cnormk = 1.
      case(3)
       cnormk = 1./pi
      case default
       write(*,666)
       ierr = 1
       return
    end select
    do i=0,ikern
      q2 = i*dq2table
      q = sqrt(q2)
      q3 = q*q2
      q4 = q2*q2
      q5 = q3*q2
      q6 = q3*q4
      q7 = q6*q
      if (q.lt.0.5) then
         wkern(i)     = 1./40.*(575./8. - 105.*q2 + 54.*q4)
         grwkern(i)   = 1./40.*(-210.*q + 216.*q2*q)
         fsoft(i)     = 1./10.*(575./24.*q3 - 21.*q5 + 54./7.*q7)
         potensoft(i) = 1./10.*(575./48.*q2 - 21./4.*q4 + 9./7.*q6 - 1199./64.)
         dphidh(i)    = 1./10.*(1199./64. - 575./16.*q2 + 105./4.*q4 - 9.*q6)
      elseif (q.lt.1.5) then
         wkern(i)   =  1./40.*(275./4. + 30.*q - 210.*q2 + 160.*q3 - 36.*q4)
         grwkern(i) = 1./40.*(30. - 420.*q + 480.*q2 - 144.*q3)
         fsoft(i)   = 1./10.*(275./12.*q3 + 15./2.*q4 - 42.*q5 + 80./3.*q6 - 36./7.*q7 + 1./672.)
         potensoft(i) = 1./10.*(275./24.*q2 + 5./2.*q3 - 21./2.*q4 + 16./3.*q5 - 6./7.*q6 - 599./32.)
         dphidh(i)    = 1./10.*(599./32. - 275./8.*q2 - 10.*q3 + 105./2.*q4 - 32.*q5 + 6.*q6)
      elseif (q.lt.2.5) then
         wkern(i)     = 1./40.*(3125./16. - 375.*q + 525./2.*q2 -  80.*q3 + 9.*q4)
         grwkern(i)   = 1./40.*(-375. + 525.*q - 240.*q2 + 36.*q3)
         fsoft(i)     = 1./10.*(3125./48.*q3 - 375./4.*q4 + 105./2.*q5 - 40./3.*q6 + 9./7.*q7 - 2185./1344.)
         potensoft(i) = 1./10.*(3125./96.*q2 - 125./4.*q3 + 105./8.*q4 - 8./3.*q5 + 3./14.*q6 - 3125./128.)
         dphidh(i)    = 1./10.*(3125./128. - 3125./32.*q2 + 125.*q3 - 525./8.*q4 + 16.*q5 - 3./2.*q6)      
      else
         wkern(i) = 0.0
         grwkern(i) = 0.0
         grgrwkern(i) = 0.
      endif
    enddo
    
   case(53)
!
!--this is a modification of the cubic spline
!
    kernellabel = 'High order peaked cubic spline 2'    
   
    radkern = 2.0
    radkern2 = radkern*radkern
    dq2table = radkern2/real(ikern)
    select case(ndim)
      case(1)
       a = 381./434.
       cnormk = 1./(7.*a/4. - 31./60.)
      case(2)
       cnormk = 4./(3.*pi)
      case(3)
       cnormk = 30./(31.*pi)/(-0.877)
      case default
       write(*,666)
       ierr = 1
       return
    end select
    do i=0,ikern
      q2 = i*dq2table
      q = sqrt(q2)
      if (q.lt.1.0) then
         wkern(i) = (a - q2)*(0.25*(2.-q)**3 - 0.5*(1.-q)**3)
         grwkern(i) = -2.*q*(0.25*(2.-q)**3 - 0.5*(1.-q)**3) &
                  + (a - q2)*(-0.75*(2.-q)**2 + 1.5*(1.-q)**2)
         grgrwkern(i) = -2.*(0.25*(2.-q)**3 - 0.5*(1.-q)**3) &
                      -4.*q*(-0.75*(2.-q)**2 + 1.5*(1.-q)**2) &
                      +(a - q2)*(1.5*(2.-q) - 3.*(1.-q))
      elseif(q.lt.2.0) then
         wkern(i) = (a - q2)*(0.25*(2.-q)**3)
         grwkern(i) = -2.*q*(0.25*(2.-q)**3) + (a - q2)*(-0.75*(2.-q)**2)
         grgrwkern(i) = -2.*(0.25*(2.-q)**3) -4.*q*(-0.75*(2.-q)**2) &
                      + (a - q2)*(1.5*(2.-q))
      else
         wkern(i) = 0.0
         grwkern(i) = 0.0
         grgrwkern(i) = 0.
      endif
    enddo

  case(58)
!
!--Super Gaussian (Monaghan & Gingold 1983)
!   
    kernellabel = '1D Super Gaussian'
  
    radkern = 2.0      ! interaction radius of kernel
    radkern2 = radkern*radkern
    dq2table = radkern*radkern/real(ikern)    
    select case(ndim)
      case(1)
        cnormk = 1./sqrt(pi)
      case default
        write(*,666)
        ierr = 1
        return
    end select
!
!--setup kernel table
!   
    do i=0,ikern
       q2 = i*dq2table
       q = sqrt(q2)
       q4 = q2*q2
       if (q2.lt.radkern2) then
          wkern(i) = exp(-q2)*(1.5 - q2)
          grwkern(i) = exp(-q2)*q*(2.*q2 - 5.)
          grgrwkern(i) = exp(-q2)*(-4.*q4 + 16.*q2 - 5.)
       else
          wkern(i) = 0.0
          grwkern(i) = 0.0
          grgrwkern(i) = 0.
       endif
    enddo

  case(60)
!
!--this is the m_6 quintic spline (see e.g. morris 1996, phd thesis)
!
    kernellabel = 'Y'''' (M_6 quintic)'  
    radkern = 3.0
    radkern2 = radkern*radkern
    dq2table = radkern2/real(ikern)
    select case(ndim)
      case(1) 
       cnormk = 1./120.
      case(2)
       cnormk = 7./(478*pi)
      case(3)
       cnormk = 1./(120.*pi)
    end select
    do i=0,ikern         
       q2 = i*dq2table
       q4 = q2*q2
       q = sqrt(q2)
       term1 = -5.*(3.-q)**4
       if (q.lt.1.0) then
          wkern(i) = 66.-60.*q2 + 30.*q4 - 10.*q4*q
          grwkern(i) = term1 + 30.*(2.-q)**4 - 75.*(1.-q)**4
          if (q.lt.epsilon(q)) then
             grgrwkern(i) = 240.
          else
             grgrwkern(i) = -2*grwkern(i)/q
          endif
       elseif ((q.ge.1.0).and.(q.lt.2.0)) then
          wkern(i) = (3.-q)**5. - 6.*(2.-q)**5.
          grwkern(i) = term1 + 30*(2.-q)**4.
          grgrwkern(i) = -2.*grwkern(i)/q
       elseif ((q.ge.2.0).and.(q.le.3.0)) then
          wkern(i) = (3.-q)**5.
          grwkern(i) = term1
          grgrwkern(i) = -2.*grwkern(i)/q
       else
          wkern(i) = 0.0
          grwkern(i) = 0.0
          grgrwkern(i) = 0.
       endif
    enddo


  case(61)
!
!--default is cubic spline (see monaghan 1992; monaghan & lattanzio 1985)
!   
    kernellabel = 'Y'''' (M_4 cubic)'    
  
    radkern = 2.0      ! interaction radius of kernel
    radkern2 = radkern*radkern
    dq2table = radkern*radkern/real(ikern)    
    select case(ndim)
      case(1)
        cnormk = 0.66666666666
      case(2)
        cnormk = 10./(7.*pi)
      case(3)
        cnormk = 1./pi
    end select
!
!--setup kernel table
!   
    do i=0,ikern
       q2 = i*dq2table
       q = sqrt(q2)
       q4 = q2*q2
       !
       ! potential must be divided by h
       ! force must be divided by h^2
       !
       if (q.lt.1.0) then
          wkern(i) = 3.*q2 - 0.75*q2*q - 6.*log(2.)*q + 2.
          grwkern(i) = 6.*q - 9./4.*q2 - 6.*log(2.)
          grgrwkern(i) = -2.*(-3. + 2.25*q)
       elseif ((q.ge.1.0).and.(q.le.2.0)) then
          wkern(i) = 0.25*q2*q - 3.*q2 + 6.*q*log(q) + 3.*q - 6.*log(2.)*q + 4.
          grwkern(i) = 0.75*q2 - 6.*q + 6.*log(q) +9. - 6.*log(2.)
          grgrwkern(i) = -2.*(-0.75*(2.-q)**2)/q
       else
          wkern(i) = 0.0
          grwkern(i) = 0.0
          grgrwkern(i) = 0.
       endif
    enddo

  case(62)
!
!--Wendland's radial functions of minimal degree
!   
    kernellabel = 'Wendland 3D kernel of degree 2'    
  
    radkern = 2.0      ! interaction radius of kernel
    radkern2 = radkern*radkern
    dq2table = radkern*radkern/real(ikern)    
    select case(ndim)
      case(1)
        cnormk = 0.75
      case(2)
        cnormk = 7./(4.*pi)
      case(3)
        cnormk = 21./(16.*pi)
    end select
!
!--setup kernel table
!   
    do i=0,ikern
       q2 = i*dq2table
       q = sqrt(q2)
       q4 = q2*q2
       !
       ! potential must be divided by h
       ! force must be divided by h^2
       !
       if (q.lt.2.0) then
          wkern(i) = (1.-0.5*q)**4*(2.*q + 1.)
          grwkern(i) = -5.*q + 7.5*q2 - 15./4.*q2*q + 5./8.*q**4
          grgrwkern(i) = -5. + 15.*q - 45./4.*q2 + 2.5*q2*q
       else
          wkern(i) = 0.0
          grwkern(i) = 0.0
          grgrwkern(i) = 0.
       endif
    enddo

  case(63)
!
!--Wendland's radial functions of minimal degree
!   
    kernellabel = 'Wendland 3D kernel of degree 4'    
  
    radkern = 2.0      ! interaction radius of kernel
    radkern2 = radkern*radkern
    dq2table = radkern*radkern/real(ikern)    
    select case(ndim)
      case(1)
        cnormk = 9./32.
      case(2)
        cnormk = 3./(4.*pi)
      case(3)
        cnormk = 165./(256.*pi)
    end select
!
!--setup kernel table
!   
    do i=0,ikern
       q2 = i*dq2table
       q = sqrt(q2)
       q4 = q2*q2
       !
       ! potential must be divided by h
       ! force must be divided by h^2
       !
       if (q.lt.2.0) then
          term = 1.-0.5*q
          wkern(i) = term**6*(35./4.*q2 + 9.*q + 3.)
          grwkern(i) = -3.*term**5*(35./4.*q2 + 9.*q + 3.) + term**6*(0.5*35.*q + 9.)
          grgrwkern(i) = 7.5*term**4*(35./4.*q2 + 9.*q + 3.) - 6.*term**5*(0.5*35.*q + 9.) + 0.5*35.*term**6
       else
          wkern(i) = 0.0
          grwkern(i) = 0.0
          grgrwkern(i) = 0.
       endif
    enddo

  case(64)
!
!--Wendland's radial functions of minimal degree
!   
    kernellabel = 'Wendland 1D kernel of degree 2'    
  
    radkern = 2.0      ! interaction radius of kernel
    radkern2 = radkern*radkern
    dq2table = radkern*radkern/real(ikern)    
    select case(ndim)
      case(1)
        cnormk = 5./8.
      case(2)
        cnormk = 5./(4.*pi)
      case(3)
        cnormk = 105./(128.*pi)
    end select
!
!--setup kernel table
!   
    do i=0,ikern
       q2 = i*dq2table
       q = sqrt(q2)
       q4 = q2*q2
       !
       ! potential must be divided by h
       ! force must be divided by h^2
       !
       if (q.lt.2.0) then
          term = 1.-0.5*q
          wkern(i) = term**3*(1.5*q + 1.)
          grwkern(i) = -3.*q + 3.*q2 - 0.75*q2*q
          grgrwkern(i) = -3. + 6.*q - 9./4.*q2
       else
          wkern(i) = 0.0
          grwkern(i) = 0.0
          grgrwkern(i) = 0.
       endif
    enddo

  case(65)
!
!--Wendland's radial functions of minimal degree
!   
    kernellabel = 'Wendland 1D kernel of degree 4'    
  
    radkern = 2.0      ! interaction radius of kernel
    radkern2 = radkern*radkern
    dq2table = radkern*radkern/real(ikern)    
    select case(ndim)
      case(1)
        cnormk = 0.75
      case(2)
        cnormk = 9./(5.*pi)
      case(3)
        cnormk = 45./(32.*pi)
    end select
!
!--setup kernel table
!   
    do i=0,ikern
       q2 = i*dq2table
       q = sqrt(q2)
       q4 = q2*q2
       !
       ! potential must be divided by h
       ! force must be divided by h^2
       !
       if (q.lt.2.0) then
          term = 1.-0.5*q
          wkern(i) = term**5*(2.*q2 + 2.5*q + 1.)
          grwkern(i) = -2.5*term**4*(2.*q2 + 2.5*q + 1.) + term**5*(4.*q + 2.5)
          grgrwkern(i) = 105./4.*q2 -3.5 -35.*q2*q + 525./32.*q4 - 21./8.*q4*q
       else
          wkern(i) = 0.0
          grwkern(i) = 0.0
          grgrwkern(i) = 0.
       endif
    enddo

  case(66)
!
!--Bessel functions
!   
    kernellabel = 'Bessel function kernel'    
  
    radkern = 2.0      ! interaction radius of kernel
    radkern2 = radkern*radkern
    dq2table = radkern*radkern/real(ikern)    
    select case(ndim)
      case(1)
        cnormk = 1.
      case(2)
        cnormk = 9./(5.*pi)
      case(3)
        cnormk = 45./(32.*pi)
    end select
!
!--setup kernel table
!   
    do i=0,ikern
       q2 = i*dq2table
       q = sqrt(q2)
       q4 = q2*q2
       !
       ! potential must be divided by h
       ! force must be divided by h^2
       !
       if (q.lt.2.0) then
          wkern(i) = 0. !sqrt(q)*(besj0(2.*sqrt(2.*q))  + besj1(2.*sqrt(2.*q)))
          grwkern(i) = 0.
          grgrwkern(i) = 0.
       else
          wkern(i) = 0.0
          grwkern(i) = 0.0
          grgrwkern(i) = 0.
       endif
    enddo

  case(67)
!
!--second derivative using ferrers n=6
!   
    kernellabel = 'second derivative kernel (ferrers n=6)'    
  
    radkern = 2.0      ! interaction radius of kernel
    radkern2 = radkern*radkern
    dq2table = radkern*radkern/real(ikern)    
    select case(ndim)
      case(1)
        cnormk = 45045./8192.
      case(2,3)
       write(*,666)
       ierr = 1
       return
      end select
!
!--setup kernel table
!   
    do i=0,ikern
       q2 = i*dq2table
       q = sqrt(q2)
       q4 = q2*q2
       !
       ! potential must be divided by h
       ! force must be divided by h^2
       !
       if (q.lt.2.0) then
          wkern(i) = 0.
          grwkern(i) = 0.
          grgrwkern(i) = (1.-0.25*q2)**6
       else
          wkern(i) = 0.0
          grwkern(i) = 0.0
          grgrwkern(i) = 0.
       endif
    enddo

  case(68)
!
!--second derivative using the cubic spline
!   
    kernellabel = 'second derivative kernel (cubic)'    
  
    radkern = 2.0      ! interaction radius of kernel
    radkern2 = radkern*radkern
    dq2table = radkern*radkern/real(ikern)    
    select case(ndim)
      case(1)
        cnormk = 4.
      case(2,3)
       write(*,666)
       ierr = 1
       return
    end select
!
!--setup kernel table
!   
    do i=0,ikern
       q2 = i*dq2table
       q = sqrt(q2)
       q4 = q2*q2
       !
       ! potential must be divided by h
       ! force must be divided by h^2
       !
       if (q.lt.1.0) then
          wkern(i) = 0.
          grwkern(i) = 0.
          grgrwkern(i) = 1. - 1.5*q2 + 0.75*q*q2
       elseif (q.le.2.0) then
          wkern(i) = 0.0
          grwkern(i) = 0.0
          grgrwkern(i) = 0.25*(2.-q)**3
       else
          wkern(i) = 0.0
          grwkern(i) = 0.0
          grgrwkern(i) = 0.
       endif
    enddo

  case (69)
!
!--hacked cubic spline
!   
    kernellabel = 'hacked M_4 cubic'    
  
    radkern = 2.0      ! interaction radius of kernel
    radkern2 = radkern*radkern
    dq2table = radkern*radkern/real(ikern)    
    select case(ndim)
      case(1)
        cnormk = 2./5.
      case(2)
        cnormk = 10./(9.*pi)
      case(3)
        cnormk = 15./(17.*pi)
    end select
!
!--setup kernel table
!   
    do i=0,ikern
       q2 = i*dq2table
       q = sqrt(q2)
       q4 = q2*q2
       !
       ! potential must be divided by h
       ! force must be divided by h^2
       !
       if (q.lt.1.0) then
          wkern(i) = 0.25*(2.-q)**3 + (1.-q)**3
          grwkern(i) = -0.75*(2.-q)**2 - 3.*(1.-q)**2
          grgrwkern(i) = 1.5*(2.-q) + 6.*(1.-q)
       elseif ((q.ge.1.0).and.(q.le.2.0)) then
          wkern(i) = 0.25*(2.-q)**3
          grwkern(i) = -0.75*(2.-q)**2
          grgrwkern(i) = 1.5*(2.-q)
       else
          potensoft(i) = -1./q
          fsoft(i) = 1.0
          dphidh(i) = 0.
          wkern(i) = 0.0
          grwkern(i) = 0.0
          grgrwkern(i) = 0.
       endif
    enddo

  case default  
!
!--default is cubic spline (see monaghan 1992; monaghan & lattanzio 1985)
!   
    kernellabel = 'M_4 cubic spline'    
  
    radkern = 2.0      ! interaction radius of kernel
    radkern2 = radkern*radkern
    dq2table = radkern*radkern/real(ikern)    
    select case(ndim)
      case(1)
        cnormk = 0.66666666666
      case(2)
        cnormk = 10./(7.*pi)
      case(3)
        cnormk = 1./pi
    end select
!
!--setup kernel table
!   
    do i=0,ikern
       q2 = i*dq2table
       q = sqrt(q2)
       q4 = q2*q2
       !
       ! potential must be divided by h
       ! force must be divided by h^2
       !
       if (q.lt.1.0) then
          potensoft(i) = 2./3.*q2 - 0.3*q4 + 0.1*q4*q - 1.4
          fsoft(i) = 4./3.*q - 1.2*q2*q + 0.5*q4
          dphidh(i) = -2.*q**2 + 1.5*q4 - 0.6*q4*q + 1.4
          wkern(i) = 1. - 1.5*q2 + 0.75*q*q2
          grwkern(i) = -3.*q+ 2.25*q2
          grgrwkern(i) = -3. + 4.5*q
       elseif ((q.ge.1.0).and.(q.le.2.0)) then
          potensoft(i) = 4./3.*q2 - q2*q + 0.3*q4 - q4*q/30. - 1.6 + 1./(15.*q)
          fsoft(i) = 8./3.*q - 3.*q2 + 1.2*q2*q - q4/6. - 1./(15.*q2)
          dphidh(i) = -4.*q2 + 4.*q2*q - 1.5*q4 + 0.2*q4*q + 1.6
          wkern(i) = 0.25*(2.-q)**3
          grwkern(i) = -0.75*(2.-q)**2
          grgrwkern(i) = 1.5*(2.-q)
       else
          potensoft(i) = -1./q
          fsoft(i) = 1.0
          dphidh(i) = 0.
          wkern(i) = 0.0
          grwkern(i) = 0.0
          grgrwkern(i) = 0.
       endif
    enddo

 end select

!
!--calculate modified kernel to use on anisotropic forces
!
 if (ianticlump.eq.1) then
!
!--this is joe's standard anticlumping term
!
    !!j = nint((1./hfact)**2/dq2table)
    j = nint((1./1.5)**2/dq2table)
    wdenom = wkern(j)

    select case(neps) ! integral of w^{n+1} dv
    case(3)
       wint = 122559./160160.
    case(4)
       wint = 200267./292864.
    case(5)
       wint = 825643615./1324331008.
    end select   
    cnormkaniso = cnormk
    !!cnormkaniso = 1./(1./cnormk + eps/(neps+1.)*wint/wdenom**neps)
    
    do i=0,ikern
       q2 = i*dq2table
       q = sqrt(q2)
       grwijalt(i) = cnormkaniso*grwkern(i)*(1. - 0.5*eps*(wkern(i)/wdenom)**neps)
       wijalt(i) = cnormkaniso*wkern(i)*(1. - 0.5*eps/(neps+1.)*(wkern(i)/wdenom)**neps) 
       grgrwijalt(i) = cnormkaniso*(grgrwkern(i)*  &
                  (1. - 0.5*eps*(wkern(i)/wdenom)**neps) &
                 -0.5*eps*neps*wkern(i)**(neps-1.)/wdenom**neps*grwkern(i)*grwkern(i))   
    enddo

    print*,'cnormk = ',cnormkaniso,eps,neps
    kernellabel=trim(kernellabel)//' & anti-clumping cubic spline' 
 
 elseif (ianticlump.eq.2) then
!
!--this is a modified version of the cubic spline
!
    beta = eps
    print*,'enter beta (0.6-0.99) , currently = ',beta
    !!read*,beta
    a = -0.25*(2.-beta)**2/(1.-beta)**2
    if (beta.eq.0.) then
       b = 0.
    else
       b = -(a+1.)/beta**2
    endif
    cnormkaniso = cnormk
    !cnormkaniso = 0.5/(1. + 0.25*(a + b*beta**4))
    do i=0,ikern
       q2 = i*dq2table
       q = sqrt(q2)
       if (q.lt.beta) then
          wijalt(i) = cnormkaniso*(0.25*(2.-q)**3 + a*(1.-q)**3 + b*(beta-q)**3)
          grwijalt(i) = -3.*cnormkaniso*(0.25*(2.-q)**2 + a*(1.-q)**2 + b*(beta-q)**2)
          grgrwijalt(i) = 6.*cnormkaniso*(0.25*(2.-q) + a*(1.-q) + b*(beta-q))       
       elseif(q.lt.1.) then
          wijalt(i) = cnormkaniso*(0.25*(2.-q)**3 + a*(1.-q)**3)
          grwijalt(i) = -3.*cnormkaniso*(0.25*(2.-q)**2 + a*(1.-q)**2)
          grgrwijalt(i) = 6.*cnormkaniso*(0.25*(2.-q) + a*(1.-q))                
       else
          wijalt(i) = cnormkaniso*wkern(i)
          grwijalt(i) = cnormkaniso*grwkern(i)
          grgrwijalt(i) = cnormkaniso*grgrwkern(i)
       endif

    enddo

    print*,'cnormk = ',cnormkaniso,eps,neps
    kernellabel=trim(kernellabel)//' & modified cubic spline'   
 
 endif
!
!--normalise kernel
!
 wkern = cnormk*wkern
 grwkern = cnormk*grwkern
 grgrwkern = cnormk*grgrwkern
 
666 format(/,'ERROR!!! normalisation constant not defined in kernel',/)
!
!--the variable ddq2table is used in interpolate_kernel
!
 ddq2table = 1./dq2table 

end subroutine setkerntable

!!----------------------------------------------------------------------
!! function to interpolate linearly from kernel tables
!! returns kernel and derivative given q^2 = (r_a-r_b)^2/h^2
!!
!! must then divide returned w, grad w by h^ndim, h^ndim+1 respectively
!!----------------------------------------------------------------------

subroutine interpolate_kernel(q2,w,gradw)
 implicit none
 integer :: index,index1
 real, intent(in) :: q2
 real, intent(out) :: w,gradw
 real :: dxx,dwdx,dgrwdx
!
!--find nearest index in kernel table
! 
 index = int(q2*ddq2table)
 index1 = index + 1
 if (index.gt.ikern .or. index.lt.0) index = ikern
 if (index1.gt.ikern .or. index1.lt.0) index1 = ikern
!
!--find increment from index point to actual value of q2
!
 dxx = q2 - index*dq2table
!
!--calculate slope for w, gradw and interpolate for each
! 
 dwdx =  (wij(index1)-wij(index))*ddq2table
 w = (wij(index)+ dwdx*dxx)

 dgrwdx =  (grwij(index1)-grwij(index))*ddq2table
 gradw = (grwij(index)+ dgrwdx*dxx)
 
end subroutine interpolate_kernel

subroutine interpolate_kerneldrag(q2,w)
 implicit none
 integer :: index,index1
 real, intent(in) :: q2
 real, intent(out) :: w
 real :: dxx,dwdx
!
!--find nearest index in kernel table
! 
 index = int(q2*ddq2table)
 index1 = index + 1
 if (index.gt.ikern .or. index.lt.0) index = ikern
 if (index1.gt.ikern .or. index1.lt.0) index1 = ikern
!
!--find increment from index point to actual value of q2
!
 dxx = q2 - index*dq2table
!
!--calculate slope for w, gradw and interpolate for each
! 
 dwdx =  (wijdrag(index1)-wijdrag(index))*ddq2table
 w = (wijdrag(index)+ dwdx*dxx)
 
end subroutine interpolate_kerneldrag

!!----------------------------------------------------------------------
!! same but for kernel *and* modified kernel in anticlumping term
!!----------------------------------------------------------------------
subroutine interpolate_kernels(q2,w,gradw,gradwalt,gradgradwalt)
 implicit none
 integer :: index,index1
 real, intent(in) :: q2
 real, intent(out) :: w,gradw,gradwalt,gradgradwalt
 real :: dxx,dwdx,dgrwdx,dgrwaltdx,dgrgrwaltdx
!
!--find nearest index in kernel table
! 
 index = int(q2*ddq2table)
 index1 = index + 1
 if (index.gt.ikern .or. index.lt.0) index = ikern
 if (index1.gt.ikern .or. index1.lt.0) index1 = ikern
!
!--find increment from index point to actual value of q2
!
 dxx = q2 - index*dq2table
!
!--calculate slope for w, gradw, waniso, gradwaniso
!  and interpolate for each
!
 w = wij(index)
 dwdx =  (wij(index1)-w)*ddq2table
 w = w + dwdx*dxx

 gradw = grwij(index)
 dgrwdx =  (grwij(index1)-gradw)*ddq2table
 gradw = gradw + dgrwdx*dxx
!
!--interpolate for alternative kernel and derivative
!
! walt = wijalt(index)
! dwaltdx =  (wijalt(index1)-walt)*ddq2table
! walt = walt + dwaltdx*dxx

 gradwalt = grwijalt(index)
 dgrwaltdx =  (grwijalt(index1)-gradwalt)*ddq2table
 gradwalt = gradwalt + dgrwaltdx*dxx

 gradgradwalt = grgrwijalt(index)
 dgrgrwaltdx =  (grgrwijalt(index1)-gradgradwalt)*ddq2table
 gradgradwalt = gradgradwalt + dgrgrwaltdx*dxx
 
end subroutine interpolate_kernels

!!----------------------------------------------------------------------
!! function to interpolate linearly from kernel tables
!! returns kernel and derivative given q^2 = (r_a-r_b)^2/h^2
!!
!! must then divide returned w, grad w by h^ndim, h^ndim+1 respectively
!!----------------------------------------------------------------------

subroutine interpolate_softening(q2,phi,force,gradw)
 implicit none
 integer :: index,index1
 real, intent(in) :: q2
 real, intent(out) :: phi,force,gradw
 real :: dxx,dphidx,dfdx,dgrwdx
!
!--find nearest index in kernel table
! 
 index = int(q2*ddq2table)
 index1 = index + 1
 if (index.gt.ikern .or. index.lt.0) index = ikern
 if (index1.gt.ikern .or. index1.lt.0) index1 = ikern
!
!--find increment from index point to actual value of q2
!
 dxx = q2 - index*dq2table
!
!--calculate slope for phi, force
! 
 dphidx =  (potensoft(index1)-potensoft(index))*ddq2table
 phi = (potensoft(index)+ dphidx*dxx)

 dfdx =  (fsoft(index1)-fsoft(index))*ddq2table
 force = (fsoft(index)+ dfdx*dxx)

 dgrwdx = (grwij(index1)-grwij(index))*ddq2table
 gradw = (grwij(index)+ dgrwdx*dxx)
 
end subroutine interpolate_softening

!!----------------------------------------------------------------------
!! function to interpolate linearly from kernel tables
!! returns kernel and derivative and dphidh given q^2 = (r_a-r_b)^2/h^2
!!
!! must then divide returned w, grad w by h^ndim, h^ndim+1 respectively
!!----------------------------------------------------------------------

subroutine interpolate_kernel_soft(q2,w,gradw,dphidhi)
 implicit none
 integer :: index,index1
 real, intent(in) :: q2
 real, intent(out) :: w,gradw,dphidhi
 real :: dxx,dwdx,dgrwdx,dpotdx
!
!--find nearest index in kernel table
! 
 index = int(q2*ddq2table)
 index1 = index + 1
 if (index.gt.ikern .or. index.lt.0) index = ikern
 if (index1.gt.ikern .or. index1.lt.0) index1 = ikern
!
!--find increment from index point to actual value of q2
!
 dxx = q2 - index*dq2table
!
!--linear interpolation
! 
 dwdx =  (wij(index1)-wij(index))*ddq2table
 w = (wij(index)+ dwdx*dxx)

 dgrwdx =  (grwij(index1)-grwij(index))*ddq2table
 gradw = (grwij(index)+ dgrwdx*dxx)

 dpotdx =  (dphidh(index1)-dphidh(index))*ddq2table
 dphidhi = (dphidh(index) + dpotdx*dxx)
 
end subroutine interpolate_kernel_soft

!!----------------------------------------------------------------------
!! function to interpolate linearly from kernel tables
!! returns kernel and second derivative given q^2 = (r_a-r_b)^2/h^2
!! (required in the new densityiterate routine)
!!
!! must then divide returned w, grad grad w by h^ndim, h^ndim+2 respectively
!!----------------------------------------------------------------------

subroutine interpolate_kernels_dens(q2,w,gradw,gradgradw,walt,gradwalt)
 implicit none
 integer :: index,index1
 real, intent(in) :: q2
 real, intent(out) :: w,gradw,gradgradw,walt,gradwalt
 real :: dxx,dwdx,dgrwdx,dgrgrwdx,dwaltdx,dgrwaltdx
!
!--find nearest index in kernel table
! 
 index = int(q2*ddq2table)
 index1 = index + 1
 if (index.gt.ikern .or. index.lt.0) index = ikern
 if (index1.gt.ikern .or. index1.lt.0) index1 = ikern
!
!--find increment from index point to actual value of q2
!
 dxx = q2 - index*dq2table
!
!--calculate slope for w, gradw, gradgradw and interpolate for each
! 
 dwdx =  (wij(index1)-wij(index))*ddq2table
 w = (wij(index)+ dwdx*dxx)

 dgrwdx =  (grwij(index1)-grwij(index))*ddq2table
 gradw = (grwij(index)+ dgrwdx*dxx)

 dgrgrwdx =  (grgrwij(index1)-grgrwij(index))*ddq2table
 gradgradw = (grgrwij(index)+ dgrgrwdx*dxx)
!
!--interpolate for alternative kernel and derivative
!
 walt = wijalt(index)
 dwaltdx =  (wijalt(index1)-walt)*ddq2table
 walt = walt + dwaltdx*dxx

 gradwalt = grwijalt(index)
 dgrwaltdx =  (grwijalt(index1)-gradwalt)*ddq2table
 gradwalt = gradwalt + dgrwaltdx*dxx
 
end subroutine interpolate_kernels_dens

!!----------------------------------------------------------------------
!! kernels used in calculating the curl in get_curl.f90
!!----------------------------------------------------------------------
subroutine interpolate_kernel_curl(q2,gradwalt,gradgradwalt)
 implicit none
 integer :: index,index1
 real, intent(in) :: q2
 real, intent(out) :: gradwalt,gradgradwalt
 real :: dxx,dgrwaltdx,dgrgrwaltdx
!
!--find nearest index in kernel table
! 
 index = int(q2*ddq2table)
 index1 = index + 1
 if (index.gt.ikern .or. index.lt.0) index = ikern
 if (index1.gt.ikern .or. index1.lt.0) index1 = ikern
!
!--find increment from index point to actual value of q2
!
 dxx = q2 - index*dq2table
!
!--calculate slope for w, gradw, waniso, gradwaniso
!  and interpolate for each
! 
 gradwalt = grwijalt(index)
 dgrwaltdx =  (grwijalt(index1)-gradwalt)*ddq2table
 gradwalt = gradwalt + dgrwaltdx*dxx

 gradgradwalt = grgrwijalt(index)
 dgrgrwaltdx =  (grgrwijalt(index1)-gradgradwalt)*ddq2table
 gradgradwalt = gradgradwalt + dgrgrwaltdx*dxx
 
end subroutine interpolate_kernel_curl

end module kernels
