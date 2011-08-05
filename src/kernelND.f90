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
 real, dimension(0:ikern) :: wij,grwij,wijaniso,grwijaniso
 real :: dq2table,ddq2table,radkern2,radkern,eps
!--these variables for force softening only
 real, dimension(0:ikern) :: potensoft,fsoft
!--these variables needed for plotting and analysis only (not in rates etc)
 real, dimension(0:ikern) :: grgrwij,grgrwijaniso
 character(len=100) :: kernelname

contains

!-----------------------------------------------------------------
! Sets up the tables for the kernel
! Returns kernel, and derivative.
!
! Default kernel is the cubic spline, but I have experimented
! with lots more.
!
!-----------------------------------------------------------------
subroutine setkern(ikernel,ndim)
 implicit none         !  define local variables
 integer, intent(in) :: ikernel, ndim
 integer :: i,j,npower
 real :: q,q2,q4,cnormk,cnormkaniso
 real :: term1,term2,term3,term4
 real :: dterm1,dterm2,dterm3,dterm4
 real :: ddterm1,ddterm2,ddterm3,ddterm4
 real :: alpha,beta,gamma,a,b,c,wdenom,wint

 cnormk = 0.0
 wij = 0.
 grwij = 0.
 grgrwij = 0.
 fsoft = 1.
 potensoft = 0.

 select case(ikernel)

  case(2)
!
!--quartic spline
!  
    kernelname = 'quartic spline'    

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
          wij(i) = (2.5-q)**4 - 5.*(1.5-q)**4 + 10.*(0.5-q)**4
          grwij(i) = -4.*((2.5-q)**3 - 5.*(1.5-q)**3 + 10*(0.5-q)**4)
          grgrwij(i) = 12.*((2.5-q)**2 - 5.*(1.5-q)**2 + 10*(0.5-q)**2)
       elseif (q.lt.1.5) then
          wij(i) = (2.5-q)**4 - 5.*(1.5-q)**4
          grwij(i) = -4.*((2.5-q)**3 - 5.*(1.5-q)**3)
          grgrwij(i) = 12.*((2.5-q)**2 - 5.*(1.5-q)**2)   
       elseif (q.lt.2.5) then
          wij(i) = (2.5-q)**4
          grwij(i) = -4.*((2.5-q)**3)
          grgrwij(i) = 12.*((2.5-q)**2)
       else
          wij(i) = 0.
          grwij(i) = 0.
          grgrwij(i) = 0.
       endif
    enddo 

  case(3)
!
!--this is the m_6 quintic spline (see e.g. morris 1996, phd thesis)
!
    kernelname = 'quintic spline'  
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
          wij(i) = 66.-60.*q2 + 30.*q4 - 10.*q4*q
          grwij(i) = term1 + 30*(2.-q)**4. - 75.*(1.-q)**4.
          grgrwij(i) = 20.*(3.-q)**3. - 120.*(2.-q)**3. + 300.*(1.-q)**3.
       elseif ((q.ge.1.0).and.(q.lt.2.0)) then
          wij(i) = (3.-q)**5. - 6.*(2.-q)**5.
          grwij(i) = term1 + 30*(2.-q)**4.
          grgrwij(i) = 20.*(3.-q)**3. - 120.*(2.-q)**3.
       elseif ((q.ge.2.0).and.(q.le.3.0)) then
          wij(i) = (3.-q)**5.
          grwij(i) = term1
          grgrwij(i) = 20.*(3.-q)**3.
       else
          wij(i) = 0.0
          grwij(i) = 0.0
          grgrwij(i) = 0.
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
       kernelname = 'New quintic (1)'    
    elseif (ikernel.eq.6) then
       beta = 0.7
       alpha = 1.5
       kernelname = 'New quintic (2)'    
    else
    !--match to cubic spline, ie w''(0) = -2
      kernelname = 'Cubic-like quintic' 
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
       stop  
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
         wij(i) = term1 + a*term2  + b*term3 + c*term4
         grwij(i) = dterm1 + a*dterm2 + b*dterm3 + c*dterm4
         grgrwij(i) = ddterm1 + a*ddterm2 + b*ddterm3 + c*ddterm4   
      elseif ((q.ge.gamma).and.(q.lt.beta)) then
         wij(i) = term1 + a*term2  + b*term3
         grwij(i) = dterm1 + a*dterm2 + b*dterm3
         grgrwij(i) = ddterm1 + a*ddterm2 + b*ddterm3       
      elseif ((q.ge.beta).and.(q.lt.alpha)) then
         wij(i) = term1 + a*term2
         grwij(i) = dterm1 + a*dterm2
         grgrwij(i) = ddterm1 + a*ddterm2
      elseif ((q.ge.alpha).and.(q.lt.radkern)) then
         wij(i) = term1
         grwij(i) = dterm1
         grgrwij(i) = ddterm1
      else
         wij(i) = 0.0
         grwij(i) = 0.0
         grgrwij(i) = 0.
      endif
    enddo

  case (8)
!
!--(1-r^2)^3
!   
    npower = 6
    write(kernelname,"(a,i1)") '(1-r^2)^',npower     

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
          wij(i) = (radkern**2-q2)**npower
          grwij(i) = -2.*npower*q*(radkern**2-q2)**(npower-1)
          grgrwij(i) = (4.*npower*(npower-1)*q2*(radkern**2-q2)**(npower-2) &
                      - 2.*npower*(radkern**2-q2)**(npower-1))
       else
          wij(i) = 0.
          grwij(i) = 0.
          grgrwij(i) = 0.
       endif
    enddo
    
   case (9)
!
!--(1-r)^n - a peaked kernel (deriv non-zero at origin) - truly awful
!   
    npower = 4
    write(kernelname,"(a,i1)") '(2-r)^',npower     

    radkern = 2.0
    radkern2 = radkern*radkern
    dq2table = radkern2/real(ikern)
    select case(ndim)
      case(1)
         cnormk = 0.5*(npower+1)/radkern**(npower+1)
      case(2,3)
         stop 'normalisation const not defined in kernel'
    end select  
    do i=0,ikern
       q2 = i*dq2table
       q = sqrt(q2)
       if (q.lt.radkern) then
          wij(i) = (radkern-q)**npower
          grwij(i) = -npower*(radkern-q)**(npower-1)
          grgrwij(i) = npower*(npower-1)*(radkern-q)**(npower-2)
       else
          wij(i) = 0.
          grwij(i) = 0.
          grgrwij(i) = 0.
       endif
    enddo   

  case(10)
!
!--gaussian
!  
    kernelname = 'Gaussian'    

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
          wij(i) = exp(-q2)
          grwij(i) = -2.*q*wij(i)
          grgrwij(i) = -2.*q*grwij(i) - 2.*wij(i)
       else
          wij(i) = 0.
          grwij(i) = 0.
          grgrwij(i) = 0.
       endif
    enddo

  case(11)
!
!--this is the usual spline based kernel modified for r/h < 2/3 to
!   prevent particles from clumping (see thomas & couchman '92)
!      
    kernelname = 'Thomas & Couchman anti-clumping'
    radkern = 2.0      ! interaction radius of kernel
    radkern2 = radkern*radkern
    dq2table = radkern*radkern/real(ikern)
    select case(ndim)
      case(1)
        cnormk = 0.66666666666   ! normalisation constant
      case(2)
        cnormk = 10./(7.*pi)
      case(3)
        cnormk = 1./pi
    end select
   
    do i=0,ikern
       q2 = i*dq2table
       q = sqrt(q2)
       if (q.lt.1.0) then
          wij(i) = 1. - 1.5*q2 + 0.75*q*q2
          if (q.lt.2./3.) then
             grwij(i) = -1.
             grgrwij(i) = 0.
          else
             grwij(i) = -3.*q+ 2.25*q2
             grgrwij(i) = -3. + 4.5*q
          endif
       elseif ((q.ge.1.0).and.(q.le.2.0)) then
          wij(i) = 0.25*(2.-q)**3.
          grwij(i) = -0.75*(2.-q)**2.
          grgrwij(i) = 1.5*(2.-q)
       else
          wij(i) = 0.0
          grwij(i) = 0.0
          grgrwij(i) = 0.
       endif
    enddo
  
  case(12)
!
!--this is the squashed quintic spline from bonet & kulesegaram
!
    kernelname = 'BK squashed quintic spline'    
   
    radkern = 2.0
    radkern2 = radkern*radkern
    dq2table = radkern2/real(ikern)
    select case(ndim)
      case(1)
       cnormk = 1./16.
      case default
       write(*,666)
       stop
    end select
    do i=0,ikern
      q2 = i*dq2table
      q = sqrt(q2)
      if (q.lt.1.0) then
         wij(i) = (2.-q)**5 - 16.*(1.-q)**5
         grwij(i) = -5.*(2.-q)**4 + 80.*(1.-q)**4.
         grgrwij(i) = 20.*(2.-q)**3 - 320.*(1.-q)**3.
      elseif ((q.ge.1.0).and.(q.le.2.0)) then
         wij(i) = (2.-q)**5
         grwij(i) = -5.*(2.-q)**4
         grgrwij(i) = 20.*(2.-q)**3
      else
         wij(i) = 0.0
         grwij(i) = 0.0
         grgrwij(i) = 0.
      endif
    enddo  
    
  case(13)
!
!--this is the lucy kernel (to 2h not to 1h)
!
    kernelname = 'Lucy kernel'    
   
    radkern = 2.0
    radkern2 = radkern*radkern
    dq2table = radkern2/real(ikern)
    select case(ndim)
      case(1)
       cnormk = 3./5. ! normalisation is probably wrong
      case(2)
       cnormk = 1/(pi)
      case(3)
       cnormk = 105/(16.*pi)
      case default
       write(*,666)
       stop
    end select
    do i=0,ikern
      q2 = i*dq2table
      q = sqrt(q2)
      if (q.lt.2.0) then
         wij(i) = (1.+1.5*q)*(1.-0.5*q)**3
         grwij(i) = (1.5*(1.-0.5*q)**3 - 1.5*(1.+1.5*q)*(1.-0.5*q)**2)
         !--note second deriv is wrong
         grgrwij(i) = 6666*(-1.5*(1.-0.5*q)**2 + 1.5*(1.+0.5*q)*(1.-0.5*q))
      else
         wij(i) = 0.0
         grwij(i) = 0.0
         grgrwij(i) = 0.
      endif
    enddo

  case(14)
!
!--this is a modification of the cubic spline
!
    kernelname = 'Peaked cubic spline'    
   
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
       stop
    end select
    do i=0,ikern
      q2 = i*dq2table
      q = sqrt(q2)
      if (q.lt.2.0) then
         wij(i) = (2.-q)**3
         grwij(i) = -3.*(2.-q)**2
         grgrwij(i) = 6.*(2.-q)
      else
         wij(i) = 0.0
         grwij(i) = 0.0
         grgrwij(i) = 0.
      endif
    enddo
    
  case(15)
!
!--this is a modification of the cubic spline
!
    kernelname = 'Peaked cubic spline 2'    
   
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
       stop
    end select
    do i=0,ikern
      q2 = i*dq2table
      q = sqrt(q2)
      if (q.lt.1.0) then
         wij(i) = 0.25*(2.-q)**3 - 0.5*(1.-q)**3
         grwij(i) = -0.75*(2.-q)**2 + 1.5*(1.-q)**2
         grgrwij(i) = 1.5*q
      elseif(q.lt.2.0) then
         wij(i) = 0.25*(2.-q)**3
         grwij(i) = -0.75*(2.-q)**2
         grgrwij(i) = 3. - 1.5*q
      else
         wij(i) = 0.0
         grwij(i) = 0.0
         grgrwij(i) = 0.
      endif
    enddo

  case(16)
!
!--another version of the peaked cubic spline
!
    alpha = 1.0
    print*,'enter alpha,'
    read*,alpha
    write(kernelname,"(a,f6.2)") 'peaked cubic spline alpha = ',alpha
!!    kernelname = 'peaked cubic spline 3'    
   
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
       stop
    end select
    do i=0,ikern
      q2 = i*dq2table
      q = sqrt(q2)
      if (q.lt.alpha) then
         wij(i) = 0.25*(2.-q)**3 - 0.5/alpha*(alpha-q)**3
         grwij(i) = -0.75*(2.-q)**2 + 1.5/alpha*(alpha-q)**2
         grgrwij(i) = 1.5*(2.-q) - 3./alpha*(alpha-q)
      elseif (q.lt.2.0) then
         wij(i) = 0.25*(2.-q)**3
         grwij(i) = -0.75*(2.-q)**2
         grgrwij(i) = 1.5*(2.-q)
      else
         wij(i) = 0.0
         grwij(i) = 0.0
         grgrwij(i) = 0.
      endif
    enddo
    
  case(17)
!
!--exponential kernel
!
    kernelname = 'Exponential'
   
    radkern = 10.0
    radkern2 = radkern*radkern
    dq2table = radkern2/real(ikern)
    select case(ndim)
      case(1)
       cnormk = 10./(6.*pi)
      case default
       write(*,666)
       stop
    end select
    do i=0,ikern
      q2 = i*dq2table
      q = sqrt(q2)
      if (q.lt.radkern) then
         wij(i) = exp(-q)
         grwij(i) = -wij(i)
         grgrwij(i) = -grwij(i)
      else
         wij(i) = 0.0
         grwij(i) = 0.0
         grgrwij(i) = 0.
      endif
    enddo    
  case(21)
!
!--cubic spline with vanishing second moment (using monaghan 1992 method)
!
    kernelname = 'higher order cubic spline'    
   
    radkern = 2.0
    radkern2 = radkern*radkern
    dq2table = radkern2/real(ikern)
    select case(ndim)
      case(1)
       a = 0.9
       cnormk = 2./(3.*a - 1.)
      case default
       write(*,666)
       stop
    end select
    do i=0,ikern
      q2 = i*dq2table
      q = sqrt(q2)
      if (q.lt.1.0) then
         wij(i) = (a - q2)*(1. - 1.5*q2 + 0.75*q2*q)
         grwij(i) = -2.*q*(1. - 1.5*q2 + 0.75*q2*q) + (a -q2)*(-3.*q + 9./4.*q2)
         grgrwij(i) = -2. + 18.*q2 - 15.*q2*q - 3.*a + 4.5*a*q
      elseif (q.lt.2.0) then
         wij(i) = 0.25*(a - q2)*(2.-q)**3
         grwij(i) = -0.5*q*(2.-q)**3 - 0.75*(a-q2)*(2.-q)**2
         grgrwij(i) = -4. + 3.*a + 18.*q - 18.*q2 - 1.5*a*q + 5.*q2*q
      else
         wij(i) = 0.0
         grwij(i) = 0.0
         grgrwij(i) = 0.
      endif
    enddo
    
   case(23)
!
!--this is a modification of the cubic spline
!
    kernelname = 'High order peaked cubic spline 2'    
   
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
       cnormk = 30./(31.*pi)
      case default
       write(*,666)
       stop
    end select
    do i=0,ikern
      q2 = i*dq2table
      q = sqrt(q2)
      if (q.lt.1.0) then
         wij(i) = (a - q2)*(0.25*(2.-q)**3 - 0.5*(1.-q)**3)
         grwij(i) = -2.*q*(0.25*(2.-q)**3 - 0.5*(1.-q)**3) &
                  + (a - q2)*(-0.75*(2.-q)**2 + 1.5*(1.-q)**2)
         grgrwij(i) = -2.*(0.25*(2.-q)**3 - 0.5*(1.-q)**3) &
                      -4.*q*(-0.75*(2.-q)**2 + 1.5*(1.-q)**2) &
                      +(a - q2)*(1.5*(2.-q) - 3.*(1.-q))
      elseif(q.lt.2.0) then
         wij(i) = (a - q2)*(0.25*(2.-q)**3)
         grwij(i) = -2.*q*(0.25*(2.-q)**3) + (a - q2)*(-0.75*(2.-q)**2)
         grgrwij(i) = -2.*(0.25*(2.-q)**3) -4.*q*(-0.75*(2.-q)**2) &
                      + (a - q2)*(1.5*(2.-q))
      else
         wij(i) = 0.0
         grwij(i) = 0.0
         grgrwij(i) = 0.
      endif
    enddo
    
  case(30)
!
!--these are the Ferrer's spheres from Dehnen (2001)
!
    kernelname = 'Ferrers n=3 sphere'    
   
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
       stop
    end select
    do i=0,ikern
      q2 = i*dq2table
      q = sqrt(q2)
      if (q.lt.2.0) then
         wij(i) = (4.-q2)**3
         grwij(i) = -6.*q*(4.-q2)**2
         grgrwij(i) = 24.*q2*(4.-q2) - 6.*(4.-q2)**2
      else
         wij(i) = 0.0
         grwij(i) = 0.0
         grgrwij(i) = 0.
      endif
    enddo
    
  case(31)
!
!--these are the Ferrer's spheres from dehnen (2001)
!
    kernelname = 'Ferrers n=6 sphere'    
   
    radkern = 2.0
    radkern2 = radkern*radkern
    dq2table = radkern2/real(ikern)
    select case(ndim)
      case(1)
       cnormk = 3003./16777216.
      case(2)
       cnormk = 7./(16384.*pi)
      case(3)
       cnormk = 45045./(134217728.*pi)
      case default
       write(*,666)
       stop
    end select
    do i=0,ikern
      q2 = i*dq2table
      q = sqrt(q2)
      if (q.lt.2.0) then
         wij(i) = (4.-q2)**6
         grwij(i) = -12.*q*(4.-q2)**5
         grgrwij(i) = 120.*q2*(4.-q2)**4 - 12.*(4.-q2)**5
      else
         wij(i) = 0.0
         grwij(i) = 0.0
         grgrwij(i) = 0.
      endif
    enddo
  case(101)
!
!--this is the potential and force corresponding to the cubic spline
!      
    kernelname = 'Potential'    
   
    radkern = 3.5
    radkern2 = radkern*radkern
    dq2table = radkern2/real(ikern)
    select case(ndim)
      case(3)
       cnormk = 1.
      case default
       write(*,666)
       stop
    end select
    do i=0,ikern
      q2 = i*dq2table
      q = sqrt(q2)
      q4 = q2*q2
      !
      ! note that the potential does not need to be divided by r
      ! (tabulated like this to prevent numerical divergences)
      ! however the force must be divided by 1/r^2
      !
      if (q.lt.1.0) then 
         potensoft(i) = 2./3.*q2 - 0.3*q4 + 0.1*q4*q - 1.4
         fsoft(i) = 4./3.*q2*q - 6./5.*q4*q + 0.5*q4*q2
      elseif (q.lt.2.0) then
         potensoft(i) = 4./3.*q2 - q2*q + 0.3*q4 - q4*q/30. - 1.6 + 1./(15.*q)
         fsoft(i) = 8./3.*q2*q - 3.*q4 + 6./5.*q4*q - q4*q2/6. - 1./15.
      else
         potensoft(i) = -1./q
         fsoft(i) = 1.0
      endif
    enddo
    
  case(102)
!
!--this is the force softening kernel corresponding to the cubic spline
!      
    kernelname = 'force'    
   
    radkern = 3.5
    radkern2 = radkern*radkern
    dq2table = radkern2/real(ikern)
    select case(ndim)
      case(3)
       cnormk = 1.
      case default
       write(*,666)
       stop
    end select
    do i=0,ikern
      q2 = i*dq2table
      q = sqrt(q2)
      q4 = q2*q2
      if (q.lt.1.0) then
         !!wij(i) = 2./3.*q2*q - 0.3*q4*q + 0.1*q4*q2 - 1.4*q
         wij(i) = 4./3.*q2*q - 6./5.*q4*q + 0.5*q4*q2
      elseif (q.lt.2.0) then
         !!wij(i) = 4./3.*q2*q - q4 + 0.3*q4*q - q4*q2/30. - 1.6*q + 1./(15.)
         wij(i) = 8./3.*q2*q - 3.*q4 + 6./5.*q4*q - q4*q2/6. - 1./15.
      else
         !!wij(i) = -1.
         wij(i) = 1.0
         grgrwij(i) = 0.
      endif
    enddo

  case default  
!
!--default is cubic spline (see monaghan 1992; monaghan & lattanzio 1985)
!   
    kernelname = 'cubic spline'    
  
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
       ! note that the potential does not need to be divided by r
       ! (tabulated like this to prevent numerical divergences)
       ! however the force must be divided by 1/r^2
       !
       if (q.lt.1.0) then
          potensoft(i) = 2./3.*q2 - 0.3*q4 + 0.1*q4*q - 1.4
          fsoft(i) = 4./3.*q2*q - 6./5.*q4*q + 0.5*q4*q2
          wij(i) = 1. - 1.5*q2 + 0.75*q*q2
          grwij(i) = -3.*q+ 2.25*q2
          grgrwij(i) = -3. + 4.5*q
       elseif ((q.ge.1.0).and.(q.le.2.0)) then
          potensoft(i) = 4./3.*q2 - q2*q + 0.3*q4 - q4*q/30. - 1.6 + 1./(15.*q)
          fsoft(i) = 8./3.*q2*q - 3.*q4 + 6./5.*q4*q - q4*q2/6. - 1./15.
          wij(i) = 0.25*(2.-q)**3.
          grwij(i) = -0.75*(2.-q)**2.
          grgrwij(i) = 1.5*(2.-q)
       else
          potensoft(i) = -1./q
          fsoft(i) = 1.0
          wij(i) = 0.0
          grwij(i) = 0.0
          grgrwij(i) = 0.
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
    wdenom = wij(j)

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
       grwijaniso(i) = cnormkaniso*grwij(i)*(1. - 0.5*eps*(wij(i)/wdenom)**neps)
       wijaniso(i) = cnormkaniso*wij(i)*(1. - 0.5*eps/(neps+1.)*(wij(i)/wdenom)**neps) 
       grgrwijaniso(i) = cnormkaniso*(grgrwij(i)*  &
                  (1. - 0.5*eps*(wij(i)/wdenom)**neps) &
                 - 0.5*eps*neps*wij(i)**(neps-1.)/wdenom**neps*grwij(i)*grwij(i))   
    enddo

    print*,'cnormk = ',cnormkaniso,eps,neps
    kernelname=trim(kernelname)//' & anti-clumping cubic spline' 
 
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
          wijaniso(i) = cnormkaniso*(0.25*(2.-q)**3 + a*(1.-q)**3 + b*(beta-q)**3)
          grwijaniso(i) = -3.*cnormkaniso*(0.25*(2.-q)**2 + a*(1.-q)**2 + b*(beta-q)**2)
          grgrwijaniso(i) = 6.*cnormkaniso*(0.25*(2.-q) + a*(1.-q) + b*(beta-q))       
       elseif(q.lt.1.) then
          wijaniso(i) = cnormkaniso*(0.25*(2.-q)**3 + a*(1.-q)**3)
          grwijaniso(i) = -3.*cnormkaniso*(0.25*(2.-q)**2 + a*(1.-q)**2)
          grgrwijaniso(i) = 6.*cnormkaniso*(0.25*(2.-q) + a*(1.-q))                
       else
          wijaniso(i) = cnormkaniso*wij(i)
          grwijaniso(i) = cnormkaniso*grwij(i)
          grgrwijaniso(i) = cnormkaniso*grgrwij(i)
       endif

    enddo

    print*,'cnormk = ',cnormkaniso,eps,neps
    kernelname=trim(kernelname)//' & modified cubic spline'   
 
 endif
!
!--normalise kernel
!
 wij = cnormk*wij
 grwij = cnormk*grwij
 grgrwij = cnormk*grgrwij

666 format(/,'ERROR!!! normalisation constant not defined in kernel',/)
!
!--the variable ddq2table is used in interpolate_kernel
!
 ddq2table = 1./dq2table 

end subroutine setkern

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
!--calculate slope for w, gradw
! 
 dwdx =  (wij(index1)-wij(index))*ddq2table
 dgrwdx =  (grwij(index1)-grwij(index))*ddq2table
!
!--interpolate for kernel and derivative
!
 w = (wij(index)+ dwdx*dxx)
 gradw = (grwij(index)+ dgrwdx*dxx)
 
end subroutine interpolate_kernel

!!----------------------------------------------------------------------
!! same but for kernal *and* modified kernel in anticlumping term
!!----------------------------------------------------------------------
subroutine interpolate_kernels(q2,w,gradw,waniso,gradwaniso)
 implicit none
 integer :: index,index1
 real, intent(in) :: q2
 real, intent(out) :: w,gradw,waniso,gradwaniso
 real :: dxx,dwdx,dgrwdx,dwanisodx,dgrwanisodx
!
!--find nearest index in kernel table
! 
 index = int(q2*ddq2table)
 index1 = index + 1
 if (index.gt.ikern .or. index.lt.0) index = ikern
 if (index1.gt.ikern .or. index.lt.0) index1 = ikern
!
!--find increment from index point to actual value of q2
!
 dxx = q2 - index*dq2table
!
!--calculate slope for w, gradw, waniso, gradwaniso
! 
 dwdx =  (wij(index1)-wij(index))*ddq2table
 dgrwdx =  (grwij(index1)-grwij(index))*ddq2table
 dwanisodx =  (wijaniso(index1)-wijaniso(index))*ddq2table
 dgrwanisodx =  (grwijaniso(index1)-grwijaniso(index))*ddq2table
!
!--interpolate for kernel and derivative
!
 w = (wij(index)+ dwdx*dxx)
 gradw = (grwij(index)+ dgrwdx*dxx)
!
!--interpolate for anticlumping kernel and derivative
!
 waniso = (wijaniso(index)+ dwanisodx*dxx)
 gradwaniso = (grwijaniso(index)+ dgrwanisodx*dxx)
 
end subroutine interpolate_kernels

!!----------------------------------------------------------------------
!! function to interpolate linearly from kernel tables
!! returns kernel and derivative given q^2 = (r_a-r_b)^2/h^2
!!
!! must then divide returned w, grad w by h^ndim, h^ndim+1 respectively
!!----------------------------------------------------------------------

subroutine interpolate_softening(q2,phi,force)
 implicit none
 integer :: index,index1
 real, intent(in) :: q2
 real, intent(out) :: phi,force
 real :: dxx,dphidx,dfdx
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
 dfdx =  (fsoft(index1)-fsoft(index))*ddq2table
!
!--interpolate for potential and force
!
 phi = (potensoft(index)+ dphidx*dxx)
 force = (fsoft(index)+ dfdx*dxx)
 
end subroutine interpolate_softening

end module kernels
