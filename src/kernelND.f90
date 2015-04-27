!------------------------------------------------------------------------------!
! NDSPMHD: A Smoothed Particle (Magneto)Hydrodynamics code for (astrophysical) !
! fluid dynamics simulations in 1, 2 and 3 spatial dimensions.                 !
!                                                                              !
! (c) 2002-2014 Daniel Price                                                   !
!                                                                              !
! http://users.monash.edu.au/~dprice/ndspmhd                                   !
! daniel.price@monash.edu -or- dprice@cantab.net (forwards to current address) !
!                                                                              !
!  NDSPMHD comes with ABSOLUTELY NO WARRANTY.                                  !
!  This is free software; and you are welcome to redistribute                  !
!  it under the terms of the GNU General Public License                        !
!  (see LICENSE file for details) and the provision that                       !
!  this notice remains intact. If you modify this file, please                 !
!  note section 2a) of the GPLv2 states that:                                  !
!                                                                              !
!  a) You must cause the modified files to carry prominent notices             !
!     stating that you changed the files and the date of any change.           !
!                                                                              !
!  ChangeLog:                                                                  !
!------------------------------------------------------------------------------!

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
 public  :: wijdrag,grwijdrag,grgrwijdrag
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
 
 radkern = 2. ! default value for kernel radius
!
!--setup kernel tables for primary kernel
!
 call setkerntable(ikernel,ndim,wij,grwij,grgrwij,kernelname,ierr)
 
end subroutine setkern

!----------------------------------------------------------------------
! This is the interface routine (public) -- calls setkerndrag once only
!
!----------------------------------------------------------------------
subroutine setkerndrag(ikerneldrag,ndim,ierr)
 implicit none
 integer, intent(in)  :: ikerneldrag, ndim
 integer, intent(out) :: ierr
!
!--setup kernel tables for drag kernel
!
 call setkerntable(ikerneldrag,ndim,wijdrag,grwijdrag,grgrwijdrag,kernelnamedrag,ierr)
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

 radkern = 2. ! default value for kernel radius
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
 real :: q,q2,q4,q3,q5,q6,q7,q8,cnormk,cnormkaniso
 real :: term1,term2,term3,term4,term
 real :: dterm1,dterm2,dterm3,dterm4
 real :: ddterm1,ddterm2,ddterm3,ddterm4,w0
 real :: alpha,beta,gamma,a,b,c,d,e,f,u,u2,qs,wdenom,wint
 integer, parameter :: lu = 55

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
!--M5 quartic (auto-generated by kernels.py)
!
    kernellabel = 'M_5 quartic' 

    radkern = max(radkern,  2.5)
    radkern2 = radkern*radkern
    dq2table = radkern2/real(ikern)
    print*,' setting up quartic with radkern = ',radkern
    select case(ndim)
      case(1)
         cnormk = 1./24.
      case(2)
         cnormk = 96./(1199.*pi)
      case(3)
         cnormk = 0.05/pi
    end select
    do i=0,ikern
       q2 = i*dq2table
       q = sqrt(q2)
       if (q < 0.5) then
          wkern(i)     = 6.0*q2*q2 - 15.0*q2 + 14.375 
          grwkern(i)   = q*(24.0*q2 - 30.0) 
          grgrwkern(i) = 72.0*q2 - 30.0 
          fsoft(i)     = q*(144.*q2*q2 - 504.*q2 + 805.)/840. 
          potensoft(i) = q2*q2*q2/35. - 3.*q2*q2/20. + 23.*q2/48. - 1199./960. 
          dphidh(i)    = -0.2*q2*q2*q2 + 0.75*q2*q2 - 1.4375*q2 + 1.24895833333333 
       elseif (q < 1.5) then
          wkern(i)     = -5.0*(-q + 1.5)**4 + (-q + 2.5)**4 
          grwkern(i)   = -16.0*q2*q + 60.0*q2 - 60.0*q + 5.0 
          grgrwkern(i) = -48.0*q2 + 120.0*q - 60.0 
          fsoft(i)     = (-768.*q**7 + 4480.*q2*q2*q2 - 8064.*q2*q2*q + 1680.*q2*q2 + 6160.*q2*q &
                         + 1.)/(6720.*q2) 
          potensoft(i) = (-128.*q**7 + 896.*q2*q2*q2 - 2016.*q2*q2*q + 560.*q2*q2 + 3080.*q2*q - &
                         8386.*q - 1.)/(6720.*q) 
          dphidh(i)    = 2.*q2*q2*q2/15. - 4.*q2*q2*q/5. + 3.*q2*q2/2. - q2*q/3. - 11.*q2/8. + &
                         599./480. 
       elseif (q < 2.5) then
          wkern(i)     = (-q + 2.5)**4 
          grwkern(i)   = -4.0*(-q + 2.5)**3 
          grgrwkern(i) = 12.0*q2 - 60.0*q + 75.0 
          fsoft(i)     = (384*q**7 - 4480.*q2*q2*q2 + 20160.*q2*q2*q - 42000.*q2*q2 + &
                         35000.*q2*q - 2185.)/(13440.*q2) 
          potensoft(i) = (64*q**7 - 896.*q2*q2*q2 + 5040.*q2*q2*q - 14000.*q2*q2 + 17500.*q2*q - &
                         21875.*q + 2185.)/(13440.*q) 
          dphidh(i)    = -q2*q2*q2/30. + 2.*q2*q2*q/5. - 15.*q2*q2/8. + 25.*q2*q/6. - &
                         125.*q2/32. + 625./384. 
       else
          wkern(i)     = 0.0 
          grwkern(i)   = 0.0 
          grgrwkern(i) = 0.0 
          fsoft(i)     = 1./q2 
          potensoft(i) = -1.0/q 
          dphidh(i)    = 0.0 
       endif
    enddo

  case(3)
!
!--this is the m_6 quintic spline (see e.g. morris 1996, phd thesis)
!
    kernellabel = 'M_6 quintic'  
    radkern = max(radkern,  3.0)
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
          grwkern(i) = term1 + 30*(2.-q)**4 - 75.*(1.-q)**4
          grgrwkern(i) = 20.*(3.-q)**3 - 120.*(2.-q)**3 + 300.*(1.-q)**3
       elseif ((q.ge.1.0).and.(q.lt.2.0)) then
          wkern(i) = (3.-q)**5 - 6.*(2.-q)**5
          grwkern(i) = term1 + 30*(2.-q)**4
          grgrwkern(i) = 20.*(3.-q)**3 - 120.*(2.-q)**3
       elseif ((q.ge.2.0).and.(q.le.3.0)) then
          wkern(i) = (3.-q)**5
          grwkern(i) = term1
          grgrwkern(i) = 20.*(3.-q)**3
       else
          wkern(i) = 0.0
          grwkern(i) = 0.0
          grgrwkern(i) = 0.
       endif
    enddo

  case(4)
!   
!--Hexic M7 kernel, used to test the streaming instability
!
    kernellabel = 'M_7 hexic'    

    radkern = max(radkern,  3.5)
    radkern2 = radkern*radkern
    dq2table = radkern2/real(ikern)
    select case(ndim)
      case(1)
         cnormk = 1./720.  !--=1/(7!):OK
      case(2)
         cnormk = 256./(113149.*pi)
      case(3)
         cnormk = 1./(840*pi)
    end select  
    do i=0,ikern 
       q2 = i*dq2table
       q  = sqrt(q2)
       q3 = q*q2
       q4 = q2*q2
       q5 = q3*q2
       q6 = q3*q3
       if (q.lt.0.5) then
          wkern(i) = -20.*q6 + 105.*q4 - 288.75*q2+367.9375
          grwkern(i) = -120.*q5+420.*q3-577.5*q
          grgrwkern(i) =  -600.*q4+1260.*q2-577.5
       elseif ((q.ge.0.5) .and. (q.lt.1.5)) then
          wkern(i) =15.*q6-105.*q5+236.25*q4-87.5*q3 &
                    -255.9375*q2-6.5625*q+368.484375
          grwkern(i) = 90.*q5-525.*q4+945.*q3 &
                      -262.5*q2-511.875*q-6.5625
          grgrwkern(i) = 450.*q4-2100.*q3+2835.*q2-525.*q-511.875   
       elseif ((q.ge.1.5).and. (q.lt.2.5)) then
          wkern(i)=-6.*q6+84.*q5-472.5*q4+1330.*q3 &
                   -1850.625*q2+950.25*q+129.28125
          grwkern(i) =  -36.*q5+420.*q4-1890.*q3 &
                       +3990.*q2-3701.25*q+950.25
          grgrwkern(i) =  -180.*q4+1680.*q3-5670.*q2+7980.*q-3701.25
       elseif ((q.ge.2.5).and. (q.lt.3.5)) then
          wkern(i) =q6-21.*q5+183.75*q4-857.5*q3 &
                    +2250.9375*q2-3151.3125*q+1838.265625
          grwkern(i) = 6.*q5-105.*q4+735.*q3-2572.5*q2 &
                      +4501.875*q-3151.3125
          grgrwkern(i) = 30.*q4-420.*q3+2205.*q2-5145.*q+4501.875 
       else
          wkern(i) = 0.
          grwkern(i) = 0.
          grgrwkern(i) = 0.
       endif
    enddo

  case(5)
!   
!--Heptic M8 kernel, used to test the streaming instability
!
    kernellabel = 'M_8 heptic'    

    radkern = max(radkern,  4.0)
    radkern2 = radkern*radkern
    dq2table = radkern2/real(ikern)
    select case(ndim)
      case(1)
         cnormk = 1./5040.  !--=1/(7!):OK
      case(2)
         cnormk = 9./(29740.*pi)
      case(3)
         cnormk = 1./(6720.*pi)
    end select  
    do i=0,ikern 
       q2 = i*dq2table
       q  = sqrt(q2)
       q3 = q*q2
       q4 = q2*q2
       q5 = q3*q2
       q6 = q3*q3
       q7 = q6*q
       if (q.lt.1.0) then
          wkern(i) = 35.*q7 - 140.*q6 + 560.*q4 &
                     - 1680.*q2+ 2416.
          grwkern(i) = 245.*q6 - 840.*q5 + 2240.*q3 - 3360.*q
          grgrwkern(i) =  1470.*q5  - 4200.*q4  &
                        + 6720.*q2 - 3360.
       elseif ((q.ge.1.0) .and. (q.lt.2.0)) then
          wkern(i) = -21.*q7 + 252.*q6 - 1176.*q5 + 2520.*q4 &
                     - 1960.*q3  - 504.*q2  - 392.*q + 2472.
          grwkern(i) = -147.*q6 + 1512.*q5 - 5880.*q4 &
                       + 10080.*q3 - 5880.*q2 - 1008.*q - 392.
          grgrwkern(i) = -882.*q5 + 7560.*q4 &
                        - 23520.*q3 +30240.*q2 -11760.*q - 1008.   
       elseif ((q.ge.2.0).and. (q.lt.3.0)) then
          wkern(i) = 7.*q7 - 140.*q6 + 1176.*q5 - 5320.*q4 &
                     + 13720.*q3 - 19320.*q2 + 12152.*q - 1112.
          grwkern(i) =  49.*q6 -  840.*q5 + 5880.*q4 &
                       - 21280.*q3 + 41160.*q2 - 38640.*q + 12152.
          grgrwkern(i) =  294.*q5 - 4200.*q4 &
                        + 23520.*q3 - 63840.*q2 + 82320.*q - 38640.
       elseif ((q.ge.3.0).and. (q.lt.4.0)) then
          wkern(i) = -q7 + 28.*q6 - 336.*q5 + 2240.*q4 &
                     - 8960.*q3 + 21504.*q2 - 28672.*q + 16384.
          grwkern(i) = -7.*q6 + 168.*q5 - 1680.*q4 &
                      + 8960.*q3 - 26880.*q2 + 43008.*q - 28672.
          grgrwkern(i) = -42.*q5 +  840.*q4 &
                          - 6720.*q3 + 26880.*q2 - 53760.*q + 43008.  
       else
          wkern(i) = 0.
          grwkern(i) = 0.
          grgrwkern(i) = 0.
       endif
    enddo

  case(6,7)
!
!--this is the do-it-yourself quintic kernel (general class of quintic splines)
!
    radkern = max(radkern,  2.0)
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

    radkern = max(radkern,  2.0)
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

 case(9)
!
!--Better cubic (auto-generated by kernels.py)
!
    kernellabel = 'Better cubic' 

    radkern = max(radkern,  2.0)
    radkern2 = radkern*radkern
    dq2table = radkern2/real(ikern)
    select case(ndim)
      case(1)
         cnormk = 1.
      case(2)
         cnormk = 240./(119.*pi)
      case(3)
         cnormk = 4./(3.*pi)
    end select
    do i=0,ikern
       q2 = i*dq2table
       q4 = q2*q2
       q6 = q4*q2
       q8 = q4*q4
       q = sqrt(q2)
       if (q < 1.0) then
          wkern(i)     = 0.375*q2*q - 0.8125*q2 + 0.625 
          grwkern(i)   = q*(1.125*q - 1.625) 
          grgrwkern(i) = 2.25*q - 1.625 
          fsoft(i)     = q*(15.*q2*q - 39.*q2 + 50.)/45. 
          potensoft(i) = q4*q/15. - 13.*q4/60. + 5.*q2/9. - 119./90. 
          dphidh(i)    = -2.*q4*q/5. + 13.*q4/12. - 5.*q2/3. + 119./90. 
       elseif (q < 2.0) then
          wkern(i)     = -0.125*q2*q + 0.8125*q2 - 1.75*q + 1.25 
          grwkern(i)   = -0.375*q2 + 1.625*q - 1.75 
          grgrwkern(i) = -0.75*q + 1.625 
          fsoft(i)     = (-5.*q6 + 39.*q4*q - 105.*q4 + 100.*q2*q - 3.)/(45.*q2) 
          potensoft(i) = (-4.*q6 + 39.*q4*q - 140.*q4 + 200.*q2*q - 272.*q + 12.)/(180.*q) 
          dphidh(i)    = 2.*q4*q/15. - 13.*q4/12. + 28.*q2*q/9. - 10.*q2/3. + 68./45. 
       else
          wkern(i)     = 0.0 
          grwkern(i)   = 0.0 
          grgrwkern(i) = 0.0 
          fsoft(i)     = 1./q2 
          potensoft(i) = -1.0/q 
          dphidh(i)    = 0.0 
       endif
    enddo

  case(10)
!
!--gaussian
!  
    kernellabel = 'Gaussian'    

    radkern = max(radkern,  10.0)
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
    radkern = max(radkern,  2.0)     ! interaction radius of kernel
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
   
    radkern = max(radkern,  2.0)
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

  case (13)
!
!--(1-r)^n - a peaked kernel (deriv non-zero at origin) - truly awful
!   
    npower = 4
    write(kernellabel,"(a,i1)") '(2-r)^',npower     

    radkern = max(radkern,  2.0)
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

  case(14)
!
!--this is a modification of the cubic spline
!
    kernellabel = 'Peaked cubic spline'    
   
    radkern = max(radkern,  2.0)
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
   
    radkern = max(radkern,  2.0)
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
    !print*,'enter alpha,'
    !read*,alpha
    write(kernellabel,"(a,f6.2)") 'peaked cubic spline alpha = ',alpha
!!    kernellabel = 'peaked cubic spline 3'    
   
    radkern = max(radkern,  2.0)
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
   
    radkern = max(radkern,  2.0)
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
!--Dirichlet kernel (sin(x)/x) (auto-generated by kernels.py)
!
    kernellabel = 'Dirichlet [sin(q)/q]'

    radkern = max(radkern,  2.0)
    radkern2 = radkern*radkern
    dq2table = radkern2/real(ikern)
    select case(ndim)
      case(1)
         cnormk = 1./3.74
         !cnormk = 0.5/(-Integral(sin(0.5*pi*q)/q, (q, 0.0)) + Integral(sin(0.5*pi*q)/q, (q, 2.0)))
      case(2)
         cnormk = 1./8.
      case(3)
         cnormk = 1./16.
    end select
    do i=0,ikern
       q2 = i*dq2table
       q4 = q2*q2
       q6 = q4*q2
       q8 = q4*q4
       q = sqrt(q2)
       if (q < epsilon(q)) then
          wkern(i)     = 0.5*pi
          grwkern(i)   = 0. 
          grgrwkern(i) = 0.
       elseif (q < 2.0) then
          wkern(i)     = sin(0.5*pi*q)/q 
          grwkern(i)   = (0.5*pi*q*cos(0.5*pi*q) - 1.0*sin(0.5*pi*q))/q2 
          grgrwkern(i) = (-0.25*pi**2*q2*sin(0.5*pi*q) - 1.0*pi*q*cos(0.5*pi*q) + &
                         2.0*sin(0.5*pi*q))/q2*q 
          fsoft(i)     = (0.5*pi*q*cos(0.5*pi*q) - 1.0*sin(0.5*pi*q))/q2 
          potensoft(i) = (0.5*pi*q*cos(0.5*pi*q) - 1.0*sin(0.5*pi*q))/q2 
          dphidh(i)    = (0.5*pi*q*cos(0.5*pi*q) - 1.0*sin(0.5*pi*q))/q2 
       else
          wkern(i)     = 0.0 
          grwkern(i)   = 0.0 
          grgrwkern(i) = 0.0 
          fsoft(i)     = 0.0 
          potensoft(i) = 0.0 
          dphidh(i)    = 0.0 
       endif
    enddo

  case(19)
!
!--sin**2/x**2
!
    kernellabel = '[sin(q)/q]**2' 
   
    radkern = max(radkern,  2.0)
    radkern2 = radkern*radkern
    dq2table = radkern2/real(ikern)
    select case(ndim)
      case(1)
         cnormk = 1./4.35 !0.5/(-Integral(sin(0.5*pi*q)**2/q**2, (q, 0.0)) + Integral(sin(0.5*pi*q)**2/q**2, (q, 2.0)))
      case(2)
         cnormk = 1.   !1.0/(-Integral(2.0*pi*sin(0.5*pi*q)**2/q, (q, 0.0)) + Integral(2.0*pi*sin(0.5*pi*q)**2/q, (q, 2.0)))
      case(3)
         cnormk = 0.25/pi
    end select

    do i=0,ikern
      q2 = i*dq2table
      q4 = q2*q2
      q = sqrt(q2)
      if (q.lt.radkern) then
         term = 0.5*pi*q
         if (q.gt.0.) then
            wkern(i)     = sin(0.5*pi*q)**2/q2 
            grwkern(i)   = (1.0*pi*q*cos(0.5*pi*q) - 2.0*sin(0.5*pi*q))*sin(0.5*pi*q)/q2*q 
            grgrwkern(i) = (-1.0*pi**2*q2*sin(0.5*pi*q)**2 + 0.5*pi**2*q2 - &
                           4.0*pi*q*sin(0.5*pi*q)*cos(0.5*pi*q) + 6.0*sin(0.5*pi*q)**2)/q4 
            fsoft(i)     = (1.0*pi*q*cos(0.5*pi*q) - 2.0*sin(0.5*pi*q))*sin(0.5*pi*q)/q2*q 
            potensoft(i) = (1.0*pi*q*cos(0.5*pi*q) - 2.0*sin(0.5*pi*q))*sin(0.5*pi*q)/q2*q 
            dphidh(i)    = (1.0*pi*q*cos(0.5*pi*q) - 2.0*sin(0.5*pi*q))*sin(0.5*pi*q)/q2*q 
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
   
    radkern = max(radkern,  2.0)
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
   
    radkern = max(radkern,  2.0)
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

    radkern = max(radkern,  2.0)
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
    radkern = max(radkern,  3.0)
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
    radkern = max(radkern,  3.0)
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

    radkern = max(radkern,  2.0)
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

 case(26)
!
!--A peaked cubic (auto-generated by kernels.py)
!
    kernellabel = 'A peaked cubic'

    radkern = max(radkern,  2.0)
    radkern2 = radkern*radkern
    dq2table = radkern2/real(ikern)
    select case(ndim)
      case(1)
         cnormk = 1./7.
      case(2)
         cnormk = 1./(3.*pi)
      case(3)
         cnormk = 15./(62.*pi)
    end select
    do i=0,ikern
       q2 = i*dq2table
       q4 = q2*q2
       q6 = q4*q2
       q8 = q4*q4
       q = sqrt(q2)
       if (q < 1.0) then
          wkern(i)     = 1.0*q2*q - 6.0*q + 6.0 
          grwkern(i)   = 3.0*q2 - 6.0 
          grgrwkern(i) = 6.0*q 
          fsoft(i)     = 5.*q*(q2*q - 9.*q + 12.)/31. 
          potensoft(i) = q4*q/31. - 15.*q2*q/31. + 30.*q2/31. - 45./31. 
          dphidh(i)    = -6.*q4*q/31. + 60.*q2*q/31. - 90.*q2/31. + 45./31. 
       elseif (q < 2.0) then
          wkern(i)     = (-q + 2.0)**3 
          grwkern(i)   = -3.0*q2 + 12.0*q - 12.0 
          grgrwkern(i) = -6.0*q + 12.0 
          fsoft(i)     = (-5.*q6 + 36.*q4*q - 90.*q4 + 80.*q2*q - 1.)/(31.*q2) 
          potensoft(i) = (-q6 + 9.*q4*q - 30.*q4 + 40.*q2*q - 48.*q + 1.)/(31.*q) 
          dphidh(i)    = 6.*q4*q/31. - 45.*q4/31. + 120.*q2*q/31. - 120.*q2/31. + 48./31. 
       else
          wkern(i)     = 0.0 
          grwkern(i)   = 0.0 
          grgrwkern(i) = 0.0 
          fsoft(i)     = 1./q2 
          potensoft(i) = -1.0/q 
          dphidh(i)    = 0.0 
       endif
    enddo

  case(30)
!
!--these are the Ferrer's spheres from Dehnen (2001)
!
    kernellabel = 'Ferrers n=3 sphere'    
   
    radkern = max(radkern,  2.0)
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
   
    radkern = max(radkern,  2.0)
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
   
    radkern = max(radkern,  2.0)
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
   
    radkern = max(radkern,  2.0)
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
   
    radkern = max(radkern,  20.0)
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
    
  case(35)
!  
!--Bump function
!
    kernellabel = 'Bump kernel'    

    radkern = max(radkern,  2.0)
    radkern2 = radkern*radkern
    dq2table = radkern2/real(ikern)
    select case(ndim)
      case(1)
         cnormk = 1./2.481250006
      case(2)
         cnormk = 1./6.505988642
      case(3)
         cnormk = 1./14.85711424
    end select  
    do i=0,ikern 
       q2    = i*dq2table
       alpha    = 1./(4.-q2)
       beta     = exp(-alpha)
       q        = sqrt(q2)
       if (q2.lt.4.) then    
          wkern    = beta
          grwkern  = -2.*q*beta*alpha**2
          grgrwkern = 2.*beta*alpha**4*(-6.*q2 + 3*q2*q2-16.)
       else
          wkern = 0.
          grwkern = 0.
          grgrwkern = 0.
       endif
    enddo

  case(41)
!
!--particle splitting kernel (cubic spline gradient  * q)
!   
    kernellabel = 'particle splitting'    
  
    radkern = max(radkern,  2.0)      ! interaction radius of kernel
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
  
    radkern = max(radkern,  2.0)      ! interaction radius of kernel
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
    radkern = max(radkern,  3.0)
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
  
    radkern = max(radkern,  2.0)      ! interaction radius of kernel
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
  
    radkern = max(radkern,  10.0)      ! interaction radius of kernel
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
  
    radkern = max(radkern,  2.0)      ! interaction radius of kernel
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

  case(47)
!
!--Double-hump-on-steroids M_4 cubic (auto-generated by kernels.py)
!
    kernellabel = 'Double-hump-on-steroids M_4 cubic' 

    radkern = 2.0
    radkern2 = radkern*radkern
    dq2table = radkern2/real(ikern)
    select case(ndim)
      case(1)
         cnormk = 70./31.
      case(2)
         cnormk = 20./(9.*pi)
      case(3)
         cnormk = 126./(127.*pi)
    end select
    do i=0,ikern
       q2 = i*dq2table
       q = sqrt(q2)
       q3 = q2*q
       q4 = q2*q2
       q5 = q4*q
       q6 = q4*q2
       q7 = q6*q
       q8 = q4*q4
       if (q < 1.0) then
          wkern(i)     = q3*(-0.25*(q - 2.0)**3 + (q - 1.0)**3) 
          grwkern(i)   = q2*(4.5*q3 - 7.5*q2 + 3.0) 
          grgrwkern(i) = q*(22.5*q3 - 30.0*q2 + 6.0) 
          fsoft(i)     = 21.*q4*(4.*q3 - 9.*q2 + 8.)/254. 
          potensoft(i) = 21.*q8/508. - 27.*q7/254. + 84.*q5/635. - 567./635. 
          dphidh(i)    = -189.*q8/508. + 108.*q7/127. - 504.*q5/635. + 567./635. 
       elseif (q < 2.0) then
          wkern(i)     = -0.25*q3*(q - 2.0)**3 
          grwkern(i)   = q2*(-1.5*q + 1.5)*(q - 2.0)**2 
          grgrwkern(i) = q*(-7.5*q3 + 30.0*q2 - 36.0*q + 12.0) 
          fsoft(i)     = (-28.*q8*q + 189.*q8 - 432.*q7 + 336.*q6 - 2.)/(254.*q2) 
          potensoft(i) = (q*(-35.*q8 + 270.*q7 - 720.*q6 + 672.*q5 - 2304.) + 20.)/(2540.*q) 
          dphidh(i)    = 63.*q8/508. - 108.*q7/127. + 252.*q6/127. - 1008.*q5/635. + 576./635. 
       else
          wkern(i)     = 0.0 
          grwkern(i)   = 0.0 
          grgrwkern(i) = 0.0 
          fsoft(i)     = q**(-2.0) 
          potensoft(i) = -1.0/q 
          dphidh(i)    = 0.0 
       endif
    enddo

  case(48)
!
!--Double-hump-on-overdrive M_4 cubic (auto-generated by kernels.py)
!
    kernellabel = 'Double-hump-on-overdrive M_4 cubic' 

    radkern = 2.0
    radkern2 = radkern*radkern
    dq2table = radkern2/real(ikern)
    select case(ndim)
      case(1)
         cnormk = 252./127.
      case(2)
         cnormk = 28./(17.*pi)
      case(3)
         cnormk = 330./(511.*pi)
    end select
    do i=0,ikern
       q2 = i*dq2table
       q4 = q2*q2
       q6 = q4*q2
       q8 = q4*q4
       q = sqrt(q2)
       if (q < 1.0) then
          wkern(i)     = q4*q*(-0.25*(q - 2.0)**3 + (q - 1.0)**3) 
          grwkern(i)   = q4*(6.0*q2*q - 10.5*q2 + 5.0) 
          grgrwkern(i) = q2*q*(42.0*q2*q - 63.0*q2 + 20.0) 
          fsoft(i)     = 3.*q6*(30.*q2*q - 66.*q2 + 55.)/511. 
          potensoft(i) = 9.*q6*q4/511. - 22.*q8*q/511. + 165.*q6*q/3577. - 2805./3577. 
          dphidh(i)    = -99.*q6*q4/511. + 220.*q8*q/511. - 1320.*q6*q/3577. + 2805./3577. 
       elseif (q < 2.0) then
          wkern(i)     = -0.25*q4*q*(q - 2.0)**3 
          grwkern(i)   = q4*(-2.0*q + 2.5)*(q - 2.0)**2 
          grgrwkern(i) = q2*q*(-14.0*q2*q + 63.0*q2 - 90.0*q + 40.0) 
          fsoft(i)     = (-30.*q6*q4*q + 198.*q6*q4 - 440.*q8*q + 330.*q8 - 1.)/(511.*q2) 
          potensoft(i) = (q*(-21.*q6*q4 + 154.*q8*q - 385.*q8 + 330.*q6*q - 2816.) + &
                         7.)/(3577.*q) 
          dphidh(i)    = 33.*q6*q4/511. - 220.*q8*q/511. + 495.*q8/511. - 2640.*q6*q/3577. + &
                         2816./3577. 
       else
          wkern(i)     = 0.0 
          grwkern(i)   = 0.0 
          grgrwkern(i) = 0.0 
          fsoft(i)     = q**(-2.0) 
          potensoft(i) = -1.0/q 
          dphidh(i)    = 0.0 
       endif
    enddo

  case(51)
!
!--cubic spline with vanishing second moment (using monaghan 1992 method)
!
    kernellabel = 'Higher order cubic spline'    
   
    radkern = max(radkern,  2.0)
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
   
    radkern = max(radkern,  2.5)
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
      q6 = q3*q3
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
   
    radkern = max(radkern,  2.0)
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
  
    radkern = max(radkern,  2.0)      ! interaction radius of kernel
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

  case(59)
!
!--default is cubic spline (see monaghan 1992; monaghan & lattanzio 1985)
!   
    kernellabel = 'Y'''' (M_4 cubic)'    
  
    radkern = max(radkern,  2.0)      ! interaction radius of kernel
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

  case(60)
!
!--this is the m_6 quintic spline (see e.g. morris 1996, phd thesis)
!
    kernellabel = 'Y'''' (M_6 quintic)'  
    radkern = max(radkern,  3.0)
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
!--Wendland 3D kernel of degree 2 (auto-generated by kernels.py)
!
    kernellabel = 'Wendland 2/3D C^2' 

    radkern = max(radkern,  2.0)
    radkern2 = radkern*radkern
    dq2table = radkern2/real(ikern)
    select case(ndim)
      case(1)
         cnormk = 3./4.
      case(2)
         cnormk = 1.75/pi
      case(3)
         cnormk = 1.3125/pi
    end select
    do i=0,ikern
       q2 = i*dq2table
       q4 = q2*q2
       q6 = q4*q2
       q8 = q4*q4
       q = sqrt(q2)
       if (q < 2.0) then
          wkern(i)     = (-0.5*q + 1.0)**4*(2.0*q + 1.0) 
          grwkern(i)   = q*(0.625*q2*q - 3.75*q2 + 7.5*q - 5.0) 
          grgrwkern(i) = 2.5*q2*q - 11.25*q2 + 15.0*q - 5.0 
          fsoft(i)     = q*(0.625*q2*q - 3.75*q2 + 7.5*q - 5.0) 
          potensoft(i) = q*(0.625*q2*q - 3.75*q2 + 7.5*q - 5.0) 
          dphidh(i)    = q*(0.625*q2*q - 3.75*q2 + 7.5*q - 5.0) 
       else
          wkern(i)     = 0.0 
          grwkern(i)   = 0.0 
          grgrwkern(i) = 0.0 
          fsoft(i)     = 0.0 
          potensoft(i) = 0.0 
          dphidh(i)    = 0.0 
       endif
    enddo

  case(62)
!
!--Wendland 3D kernel of degree 4 (auto-generated by kernels.py)
!
    kernellabel = 'Wendland 2/3D C^4' 

    radkern = max(radkern,  2.0)
    radkern2 = radkern*radkern
    dq2table = radkern2/real(ikern)
    select case(ndim)
      case(1)
         cnormk = 27./32.
      case(2)
         cnormk = 2.25/pi
      case(3)
         cnormk = 495./(256.*pi)
    end select
    do i=0,ikern
       q2 = i*dq2table
       q4 = q2*q2
       q6 = q4*q2
       q8 = q4*q4
       q = sqrt(q2)
       if (q < 2.0) then
          wkern(i)     = (-q/2. + 1.)**6*(35.*q2/12. + 3.*q + 1.) 
          grwkern(i)   = (-q/2. + 1.)**6*(35.*q/6. + 3.) - 3.*(-q/2. + 1.)**5*(35.*q2/12. + 3.*q &
                         + 1.) 
          grgrwkern(i) = 2.55208333333333*q6 - 21.0*q4*q + 65.625*q4 - 93.3333333333333*q2*q + &
                         52.5*q2 - 4.66666666666667 
          fsoft(i)     = (-q/2. + 1.)**6*(35.*q/6. + 3.) - 3.*(-q/2. + 1.)**5*(35.*q2/12. + 3.*q &
                         + 1.) 
          potensoft(i) = (-q/2. + 1.)**6*(35.*q/6. + 3.) - 3.*(-q/2. + 1.)**5*(35.*q2/12. + 3.*q &
                         + 1.) 
          dphidh(i)    = (-q/2. + 1.)**6*(35.*q/6. + 3.) - 3.*(-q/2. + 1.)**5*(35.*q2/12. + 3.*q &
                         + 1.) 
       else
          wkern(i)     = 0.0 
          grwkern(i)   = 0.0 
          grgrwkern(i) = 0.0 
          fsoft(i)     = 0.0 
          potensoft(i) = 0.0 
          dphidh(i)    = 0.0 
       endif
    enddo

  case(63)
!
!--Wendland 3D kernel of degree 6 (auto-generated by kernels.py)
!
    kernellabel = 'Wendland 2/3D C^6' 

    radkern = max(radkern,  2.0)
    radkern2 = radkern*radkern
    dq2table = radkern2/real(ikern)
    select case(ndim)
      case(1)
         cnormk = 15./16.
      case(2)
         cnormk = 39./(14.*pi)
      case(3)
         cnormk = 1365./(512.*pi)
    end select
    do i=0,ikern
       q2 = i*dq2table
       q4 = q2*q2
       q6 = q4*q2
       q8 = q4*q4
       q = sqrt(q2)
       if (q < 2.0) then
          wkern(i)     = (-0.5*q + 1.0)**8*(4.0*q2*q + 6.25*q2 + 4.0*q + 1.0) 
          grwkern(i)   = (-q/2. + 1.)**8*(12.*q2 + 25.*q/2. + 4.) - 4.*(-q/2. + 1.)**7*(4.*q2*q &
                         + 25.*q2/4. + 4.*q + 1.) 
          grgrwkern(i) = 1.71875*q8*q - 20.302734375*q8 + 99.0*q6*q - 252.65625*q6 + 346.5*q4*q &
                         - 216.5625*q4 + 49.5*q2 - 5.5 
          fsoft(i)     = (-q/2. + 1.)**8*(12.*q2 + 25.*q/2. + 4.) - 4.*(-q/2. + 1.)**7*(4.*q2*q &
                         + 25.*q2/4. + 4.*q + 1.) 
          potensoft(i) = (-q/2. + 1.)**8*(12.*q2 + 25.*q/2. + 4.) - 4.*(-q/2. + 1.)**7*(4.*q2*q &
                         + 25.*q2/4. + 4.*q + 1.) 
          dphidh(i)    = (-q/2. + 1.)**8*(12.*q2 + 25.*q/2. + 4.) - 4.*(-q/2. + 1.)**7*(4.*q2*q &
                         + 25.*q2/4. + 4.*q + 1.)
       else
          wkern(i)     = 0.0 
          grwkern(i)   = 0.0 
          grgrwkern(i) = 0.0 
          fsoft(i)     = 0.0 
          potensoft(i) = 0.0 
          dphidh(i)    = 0.0 
       endif
    enddo

  case(64)
!
!--Wendland 1D kernel of degree 2 (auto-generated by kernels.py)
!
    kernellabel = 'Wendland 1D C^2 [Lucy kernel]' 

    radkern = max(radkern,  2.0)
    radkern2 = radkern*radkern
    dq2table = radkern2/real(ikern)
    select case(ndim)
      case(1)
         cnormk = 5./8.
      case(2)
         cnormk = 1.25/pi
      case(3)
         cnormk = 105./(128.*pi)
    end select
    do i=0,ikern
       q2 = i*dq2table
       q4 = q2*q2
       q6 = q4*q2
       q = sqrt(q2)
       if (q < 2.0) then
          wkern(i)     = -0.1875*q4 + 1.0*q2*q - 1.5*q2 + 1.0 
          grwkern(i)   = q*(-0.75*q2 + 3.0*q - 3.0) 
          grgrwkern(i) = -2.25*q2 + 6.0*q - 3.0 
          fsoft(i)     = q*(-45.*q4 + 280.*q2*q - 504.*q2 + 560.)/512. 
          potensoft(i) = -15.*q6/1024. + 7.*q4*q/64. - 63.*q4/256. + 35.*q2/64. - 21./16. 
          dphidh(i)    = 105.*q6/1024. - 21.*q4*q/32. + 315.*q4/256. - 105.*q2/64. + 21./16. 
       else
          wkern(i)     = 0.0 
          grwkern(i)   = 0.0 
          grgrwkern(i) = 0.0 
          fsoft(i)     = 1./q2 
          potensoft(i) = -1.0/q 
          dphidh(i)    = 0.0 
       endif
    enddo

  case(65)
!
!--Wendland 1D kernel of degree 4 (auto-generated by kernels.py)
!
    kernellabel = 'Wendland 1D C^4' 

    radkern = max(radkern,  2.0)
    radkern2 = radkern*radkern
    dq2table = radkern2/real(ikern)
    select case(ndim)
      case(1)
         cnormk = 3./4.
      case(2)
         cnormk = 1.8/pi
      case(3)
         cnormk = 45./(32.*pi)
    end select
    do i=0,ikern
       q2 = i*dq2table
       q4 = q2*q2
       q6 = q4*q2
       q = sqrt(q2)
       if (q < 2.0) then
          wkern(i)     = (-q/2. + 1.)**5*(2.*q2 + 5.*q/2. + 1.) 
          grwkern(i)   = q*(-0.4375*q4*q + 3.28125*q4 - 8.75*q2*q + 8.75*q2 - 3.5) 
          grgrwkern(i) = -2.625*q4*q + 16.40625*q4 - 35.0*q2*q + 26.25*q2 - 3.5 
          fsoft(i)     = q*(-18.*q6*q + 175.*q6 - 630.*q4*q + 900.*q4 - 1008.*q2 + 960.)/512. 
          potensoft(i) = -q4*q4*q/256. + 175.*q4*q4/4096. - 45.*q6*q/256. + 75.*q6/256. - &
                         63.*q4/128. + 15.*q2/16. - 25./16. 
          dphidh(i)    = 5.*q4*q4*q/128. - 1575.*q4*q4/4096. + 45.*q6*q/32. - 525.*q6/256. + &
                         315.*q4/128. - 45.*q2/16. + 25./16. 
       else
          wkern(i)     = 0.0 
          grwkern(i)   = 0.0 
          grgrwkern(i) = 0.0 
          fsoft(i)     = 1./q2 
          potensoft(i) = -1.0/q 
          dphidh(i)    = 0.0 
       endif
    enddo

  case(66)
!
!--Wendland 1D kernel of degree 6 (auto-generated by kernels.py)
!
    kernellabel = 'Wendland 1D C^6' 

    radkern = max(radkern,  2.0)
    radkern2 = radkern*radkern
    dq2table = radkern2/real(ikern)
    select case(ndim)
      case(1)
         cnormk = 55./64.
      case(2)
         cnormk = 33./(14.*pi)
      case(3)
         cnormk = 2145./(1024.*pi)
    end select
    do i=0,ikern
       q2 = i*dq2table
       q4 = q2*q2
       q6 = q4*q2
       q8 = q4*q4
       q = sqrt(q2)
       if (q < 2.0) then
          wkern(i)     = (-0.5*q + 1.0)**7*(2.625*q2*q + 4.75*q2 + 3.5*q + 1.0) 
          grwkern(i)   = q*(-0.205078125*q8 + 2.25*q6*q - 9.84375*q6 + 21.0*q4*q - 19.6875*q4 + &
                         10.5*q2 - 4.5) 
          grgrwkern(i) = -1.845703125*q8 + 18.0*q6*q - 68.90625*q6 + 126.0*q4*q - 98.4375*q4 + &
                         31.5*q2 - 4.5 
          fsoft(i)     = q*(-3465.*q6*q4 + 45760.*q8*q - 245700.*q8 + 658944.*q6*q - 800800.*q6 &
                         + 823680.*q4 - 988416.*q2 + 732160.)/262144. 
          potensoft(i) = -1155.*q6*q6/1048576. + 65.*q6*q4*q/4096. - 12285.*q6*q4/131072. + &
                         143.*q8*q/512. - 25025.*q8/65536. + 2145.*q6/4096. - &
                         3861.*q4/4096. + 715.*q2/512. - 455./256. 
          dphidh(i)    = 15015.*q6*q6/1048576. - 195.*q6*q4*q/1024. + 135135.*q6*q4/131072. - &
                         715.*q8*q/256. + 225225.*q8/65536. - 15015.*q6/4096. + &
                         19305.*q4/4096. - 2145.*q2/512. + 455./256. 
       else
          wkern(i)     = 0.0 
          grwkern(i)   = 0.0 
          grgrwkern(i) = 0.0 
          fsoft(i)     = 1./q2 
          potensoft(i) = -1.0/q 
          dphidh(i)    = 0.0 
       endif
    enddo

  case(67)
!
!--Bessel functions
!   
    kernellabel = 'Bessel function kernel'    
  
    radkern = max(radkern,  2.0)      ! interaction radius of kernel
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

  case(68)
!
!--second derivative using the cubic spline
!   
    kernellabel = 'second derivative kernel (cubic)'    
  
    radkern = max(radkern,  2.0)      ! interaction radius of kernel
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
  
    radkern = max(radkern,  2.0)      ! interaction radius of kernel
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

  case(70)
!
!--Double-hump Wendland 2/3D C^2 (auto-generated by kernels.py)
!
    kernellabel = 'Double-hump Wendland 2/3D C^2' 

    radkern = 2.0
    radkern2 = radkern*radkern
    dq2table = radkern2/real(ikern)
    select case(ndim)
      case(1)
         cnormk = 21./8.
      case(2)
         cnormk = 3.15/pi
      case(3)
         cnormk = 105./(64.*pi)
    end select
    do i=0,ikern
       q2 = i*dq2table
       q4 = q2*q2
       q6 = q4*q2
       q8 = q4*q4
       q = sqrt(q2)
       if (q < 2.0) then
          wkern(i)     = q2*(0.5*q - 1.0)**4*(2.0*q + 1.0) 
          grwkern(i)   = q*(0.875*q4*q - 5.625*q4 + 12.5*q2*q - 10.0*q2 + 2.0) 
          grgrwkern(i) = 5.25*q4*q - 28.125*q4 + 50.0*q2*q - 30.0*q2 + 2.0 
          fsoft(i)     = q2*q*(21.*q4*q - 175.*q4 + 525.*q2*q - 600.*q2 + 336.)/256. 
          potensoft(i) = 7.*q8*q/768. - 175.*q8/2048. + 75.*q6*q/256. - 25.*q6/64. + 21.*q4/64. &
                         - 25./24. 
          dphidh(i)    = -35.*q8*q/384. + 1575.*q8/2048. - 75.*q6*q/32. + 175.*q6/64. - &
                         105.*q4/64. + 25./24. 
       else
          wkern(i)     = 0.0 
          grwkern(i)   = 0.0 
          grgrwkern(i) = 0.0 
          fsoft(i)     = q**(-2.0) 
          potensoft(i) = -1.0/q 
          dphidh(i)    = 0.0 
       endif
    enddo

  case(71)
!
!--Double-hump Wendland 2/3D C^4 (auto-generated by kernels.py)
!
    kernellabel = 'Double-hump Wendland 2/3D C^4' 

    radkern = 2.0
    radkern2 = radkern*radkern
    dq2table = radkern2/real(ikern)
    select case(ndim)
      case(1)
         cnormk = 495./128.
      case(2)
         cnormk = 297./(56.*pi)
      case(3)
         cnormk = 6435./(2048.*pi)
    end select
    do i=0,ikern
       q2 = i*dq2table
       q4 = q2*q2
       q6 = q4*q2
       q8 = q4*q4
       q = sqrt(q2)
       if (q < 2.0) then
          wkern(i)     = q2*(-q/2. + 1.)**6*(35.*q2/12. + 3.*q + 1.) 
          grwkern(i)   = q2*(-q/2. + 1.)**6*(35.*q/6. + 3.) - 3.*q2*(-q/2. + 1.)**5*(35.*q2/12. &
                         + 3.*q + 1.) + 2.*q*(-q/2. + 1.)**6*(35.*q2/12. + 3.*q + 1.) 
          grgrwkern(i) = 4.1015625*q8 - 36.0*q6*q + 122.5*q6 - 196.0*q4*q + 131.25*q4 + &
                         7.105427357601e-15*q2*q - 28.0*q2 + 2.0 
          fsoft(i)     = q2*q*(5775.*q8 - 68640.*q6*q + 327600.*q6 - 768768.*q4*q + 800800.*q4 - &
                         549120.*q2 + 329472.)/131072. 
          potensoft(i) = 1925.*q6*q6/524288. - 195.*q6*q4*q/4096. + 4095.*q6*q4/16384. - &
                         1001.*q8*q/1536. + 25025.*q8/32768. - 715.*q6/1024. + &
                         1287.*q4/2048. - 455./384. 
          dphidh(i)    = -25025.*q6*q6/524288. + 585.*q6*q4*q/1024. - 45045.*q6*q4/16384. + &
                         5005.*q8*q/768. - 225225.*q8/32768. + 5005.*q6/1024. - &
                         6435.*q4/2048. + 455./384. 
       else
          wkern(i)     = 0.0 
          grwkern(i)   = 0.0 
          grgrwkern(i) = 0.0 
          fsoft(i)     = q**(-2.0) 
          potensoft(i) = -1.0/q 
          dphidh(i)    = 0.0 
       endif
    enddo

  case(72)
!
!--Double-hump Wendland 2/3D C^6 (auto-generated by kernels.py)
!
    kernellabel = 'Double-hump Wendland 2/3D C^6' 

    radkern = 2.0
    radkern2 = radkern*radkern
    dq2table = radkern2/real(ikern)
    select case(ndim)
      case(1)
         cnormk = 1365./256.
      case(2)
         cnormk = 8.125/pi
      case(3)
         cnormk = 1365./(256.*pi)
    end select
    do i=0,ikern
       q2 = i*dq2table
       q4 = q2*q2
       q6 = q4*q2
       q8 = q4*q4
       q = sqrt(q2)
       if (q < 2.0) then
          wkern(i)     = q2*(0.5*q - 1.0)**8*(4.0*q2*q + 6.25*q2 + 4.0*q + 1.0) 
          grwkern(i)   = q*(0.203125*q6*q4*q - 2.70703125*q6*q4 + 15.125*q8*q - 45.1171875*q8 + &
                         74.25*q6*q - 57.75*q6 + 24.75*q4 - 11.0*q2 + 2.0) 
          grgrwkern(i) = 2.4375*q6*q4*q - 29.77734375*q6*q4 + 151.25*q8*q - 406.0546875*q8 + &
                         594.0*q6*q - 404.25*q6 + 123.75*q4 - 33.0*q2 + 2.0 
          fsoft(i)     = q2*q*(1365.*q6*q4*q - 21021.*q6*q4 + 137280.*q8*q - 485100.*q8 + &
                         960960.*q6*q - 917280.*q6 + 640640.*q4 - 549120.*q2 + &
                         279552.)/65536. 
          potensoft(i) = 91.*q**15/65536. - 3003.*q**14/131072. + 165.*q**13/1024. - &
                         40425.*q6*q6/65536. + 1365.*q6*q4*q/1024. - 5733.*q6*q4/4096. + &
                         5005.*q8/4096. - 715.*q6/512. + 273.*q4/256. - 21./16. 
          dphidh(i)    = -91.*q**15/4096. + 45045.*q**14/131072. - 1155.*q**13/512. + &
                         525525.*q6*q6/65536. - 4095.*q6*q4*q/256. + 63063.*q6*q4/4096. - &
                         45045.*q8/4096. + 5005.*q6/512. - 1365.*q4/256. + 21./16. 
       else
          wkern(i)     = 0.0 
          grwkern(i)   = 0.0 
          grgrwkern(i) = 0.0 
          fsoft(i)     = q**(-2.0) 
          potensoft(i) = -1.0/q 
          dphidh(i)    = 0.0 
       endif
    enddo

 case(91)
!
!--integrated M4 (auto-generated by kernels.py)
!
    kernellabel = 'Integrated M_4' 

    radkern = max(radkern,  2.0)
    radkern2 = radkern*radkern
    dq2table = radkern2/real(ikern)
    select case(ndim)
      case(1)
         cnormk = 2.
      case(2)
         cnormk = 140./(31.*pi)
      case(3)
         cnormk = 10./(3.*pi)
    end select
    do i=0,ikern
       q2 = i*dq2table
       q = sqrt(q2)
       if (q < 1.0) then
          wkern(i)     = -0.15*q2*q2*q + 0.375*q2*q2 - 0.5*q2 + 0.35 
          grwkern(i)   = q*(-0.75*q2*q + 1.5*q2 - 1.0) 
          grgrwkern(i) = -3.0*q2*q + 4.5*q2 - 1.0 
          fsoft(i)     = q*(-63.*q2*q2*q + 180.*q2*q2 - 336.*q2 + 392.)/252. 
          potensoft(i) = -q**7/28. + 5.*q2*q2*q2/42. - q2*q2/3. + 7.*q2/9. - 31./21. 
          dphidh(i)    = 2.*q**7/7. - 5.*q2*q2*q2/6. + 5.*q2*q2/3. - 7.*q2/3. + 31./21. 
       elseif (q < 2.0) then
          wkern(i)     = 0.05*q2*q2*q - 0.375*q2*q2 + 1.0*q2*q - 1.0*q2 + 0.4 
          grwkern(i)   = q*(0.25*q2*q - 1.5*q2 + 3.0*q - 2.0) 
          grgrwkern(i) = 1.0*q2*q - 4.5*q2 + 6.0*q - 2.0 
          fsoft(i)     = (21*q**8 - 180.*q**7 + 560.*q2*q2*q2 - 672.*q2*q2*q + 448.*q2*q - &
                         4.)/(252.*q2) 
          potensoft(i) = (3*q**8 - 30.*q**7 + 112.*q2*q2*q2 - 168.*q2*q2*q + 224.*q2*q - 384.*q &
                         + 4.)/(252.*q) 
          dphidh(i)    = -2.*q**7/21. + 5.*q2*q2*q2/6. - 8.*q2*q2*q/3. + 10.*q2*q2/3. - 8.*q2/3. &
                         + 32./21. 
       else
          wkern(i)     = 0.0 
          grwkern(i)   = 0.0 
          grgrwkern(i) = 0.0 
          fsoft(i)     = 1./q2 
          potensoft(i) = -1.0/q 
          dphidh(i)    = 0.0 
       endif
    enddo

 case(92)
!
!--integrated M5 (auto-generated by kernels.py)
!
    kernellabel = 'Integrated M_5'
    radkern = max(radkern,  2.5)
    radkern2 = radkern*radkern
    dq2table = radkern2/real(ikern)
    select case(ndim)
      case(1)
         cnormk = 1./10.
      case(2)
         cnormk = 7168./(35783.*pi)
      case(3)
         cnormk = 3./(23.*pi)
    end select
    do i=0,ikern
       q2 = i*dq2table
       q = sqrt(q2)
       if (q < 0.5) then
          wkern(i)     = -q2*q2*q2 + 15.*q2*q2/4. - 115.*q2/16. + 1199./192. 
          grwkern(i)   = q*(-6.0*q2*q2 + 15.0*q2 - 14.375) 
          grgrwkern(i) = -30.0*q2*q2 + 45.0*q2 - 14.375 
          fsoft(i)     = q*(-448.*q2*q2*q2 + 2160.*q2*q2 - 5796.*q2 + 8393.)/7728. 
          potensoft(i) = -q**8/138. + 15.*q2*q2*q2/322. - 3.*q2*q2/16. + 1199.*q2/2208. - &
                         107349./82432. 
          dphidh(i)    = 3.*q**8/46. - 15.*q2*q2*q2/46. + 15.*q2*q2/16. - 1199.*q2/736. + &
                         107349./82432. 
       elseif (q < 1.5) then
          wkern(i)     = 2.*q2*q2*q2/3. - 4.*q2*q2*q + 15.*q2*q2/2. - 5.*q2*q/3. - 55.*q2/8. + &
                         599./96. 
          grwkern(i)   = q*(4.0*q2*q2 - 20.0*q2*q + 30.0*q2 - 5.0*q - 13.75) 
          grgrwkern(i) = 20.0*q2*q2 - 80.0*q2*q + 90.0*q2 - 10.0*q - 13.75 
          fsoft(i)     = (7168*q**9 - 48384.*q**8 + 103680.*q**7 - 26880.*q2*q2*q2 - &
                         133056.*q2*q2*q + 201264.*q2*q + 1.)/(185472.*q2) 
          potensoft(i) = (1792*q**9 - 13824.*q**8 + 34560.*q**7 - 10752.*q2*q2*q2 - &
                         66528.*q2*q2*q + 201264.*q2*q - 483057.*q - 2.)/(370944.*q) 
          dphidh(i)    = -q**8/23. + 48.*q**7/161. - 15.*q2*q2*q2/23. + 4.*q2*q2*q/23. + &
                         165.*q2*q2/184. - 599.*q2/368. + 53673./41216. 
       elseif (q < 2.5) then
          wkern(i)     = -q2*q2*q2/6. + 2.*q2*q2*q - 75.*q2*q2/8. + 125.*q2*q/6. - 625.*q2/32. + &
                         3125./384. 
          grwkern(i)   = q*(-1.0*q2*q2 + 10.0*q2*q - 37.5*q2 + 62.5*q - 39.0625) 
          grgrwkern(i) = -5.0*q2*q2 + 40.0*q2*q - 112.5*q2 + 125.0*q - 39.0625 
          fsoft(i)     = (-3584.*q**9 + 48384.*q**8 - 259200.*q**7 + 672000.*q2*q2*q2 - &
                         756000.*q2*q2*q + 525000.*q2*q - 19681.)/(370944.*q2) 
          potensoft(i) = (-1792.*q**9 + 27648.*q**8 - 172800.*q**7 + 537600.*q2*q2*q2 - &
                         756000.*q2*q2*q + 1050000.*q2*q - 2109375.*q + &
                         78724.)/(1483776.*q) 
          dphidh(i)    = q**8/92. - 24.*q**7/161. + 75.*q2*q2*q2/92. - 50.*q2*q2*q/23. + &
                         1875.*q2*q2/736. - 3125.*q2/1472. + 234375./164864. 
       else
          wkern(i)     = 0.0 
          grwkern(i)   = 0.0 
          grgrwkern(i) = 0.0 
          fsoft(i)     = 1./q2 
          potensoft(i) = -1.0/q 
          dphidh(i)    = 0.0 
       endif
    enddo

 case(93)
!
!--integrated M6 (auto-generated by kernels.py)
!
    kernellabel = 'Integrated M_6' 

    radkern = max(radkern,  3.0)
    radkern2 = radkern*radkern
    dq2table = radkern2/real(ikern)
    select case(ndim)
      case(1)
         cnormk = 1./60.
      case(2)
         cnormk = 84./(2771.*pi)
      case(3)
         cnormk = 1./(56.*pi)
    end select
    do i=0,ikern
       q2 = i*dq2table
       q = sqrt(q2)
       if (q < 1.0) then
          wkern(i)     = 10.*q**7/7. - 5.*q2*q2*q2 + 15.*q2*q2 - 33.*q2 + 239./7. 
          grwkern(i)   = q*(10.0*q2*q2*q - 30.0*q2*q2 + 60.0*q2 - 66.0) 
          grgrwkern(i) = 60.0*q2*q2*q - 150.0*q2*q2 + 180.0*q2 - 66.0 
          fsoft(i)     = q*(45.*q**7 - 175.*q2*q2*q2 + 675.*q2*q2 - 2079.*q2 + 3585.)/4410. 
          potensoft(i) = q**9/882. - 5.*q**8/1008. + 5.*q2*q2*q2/196. - 33.*q2*q2/280. + &
                         239.*q2/588. - 2771./2352. 
          dphidh(i)    = -5.*q**9/441. + 5.*q**8/112. - 5.*q2*q2*q2/28. + 33.*q2*q2/56. - &
                         239.*q2/196. + 2771./2352. 
       elseif (q < 2.0) then
          wkern(i)     = -5.*q**7/7. + 15.*q2*q2*q2/2. - 30.*q2*q2*q + 105.*q2*q2/2. - 25.*q2*q &
                         - 51.*q2/2. + 473./14. 
          grwkern(i)   = q*(-5.0*q2*q2*q + 45.0*q2*q2 - 150.0*q2*q + 210.0*q2 - 75.0*q - 51.0) 
          grgrwkern(i) = -30.0*q2*q2*q + 225.0*q2*q2 - 600.0*q2*q + 630.0*q2 - 150.0*q - 51.0 
          fsoft(i)     = (-90.*q**10 + 1050.*q**9 - 4725.*q**8 + 9450.*q**7 - 5250.*q2*q2*q2 - &
                         6426.*q2*q2*q + 14190.*q2*q + 5.)/(17640.*q2) 
          potensoft(i) = (-40.*q**10 + 525.*q**9 - 2700.*q**8 + 6300.*q**7 - 4200.*q2*q2*q2 - &
                         6426.*q2*q2*q + 28380.*q2*q - 83055.*q - 20.)/(70560.*q) 
          dphidh(i)    = 5.*q**9/882. - 15.*q**8/224. + 15.*q**7/49. - 5.*q2*q2*q2/8. + &
                         5.*q2*q2*q/14. + 51.*q2*q2/112. - 473.*q2/392. + 113./96. 
       elseif (q < 3.0) then
          wkern(i)     = q**7/7. - 5.*q2*q2*q2/2. + 18.*q2*q2*q - 135.*q2*q2/2. + 135.*q2*q - &
                         243.*q2/2. + 729./14. 
          grwkern(i)   = q*(1.0*q2*q2*q - 15.0*q2*q2 + 90.0*q2*q - 270.0*q2 + 405.0*q - 243.0) 
          grgrwkern(i) = 6.0*q2*q2*q - 75.0*q2*q2 + 360.0*q2*q - 810.0*q2 + 810.0*q - 243.0 
          fsoft(i)     = (18*q**10 - 350.*q**9 + 2835.*q**8 - 12150.*q**7 + 28350.*q2*q2*q2 - &
                         30618.*q2*q2*q + 21870.*q2*q - 2043.)/(17640.*q2) 
          potensoft(i) = (8*q**10 - 175.*q**9 + 1620.*q**8 - 8100.*q**7 + 22680.*q2*q2*q2 - &
                         30618.*q2*q2*q + 43740.*q2*q - 98415.*q + 8172.)/(70560.*q) 
          dphidh(i)    = -q**9/882. + 5.*q**8/224. - 9.*q**7/49. + 45.*q2*q2*q2/56. - &
                         27.*q2*q2*q/14. + 243.*q2*q2/112. - 729.*q2/392. + 2187./1568. 
       else
          wkern(i)     = 0.0 
          grwkern(i)   = 0.0 
          grgrwkern(i) = 0.0 
          fsoft(i)     = 1./q2 
          potensoft(i) = -1.0/q 
          dphidh(i)    = 0.0 
       endif
    enddo

 case(94)
!
!--twice integrated M4 (auto-generated by kernels.py)
!
    kernellabel = 'Twice integrated M_4' 

    radkern = max(radkern,  2.0)
    radkern2 = radkern*radkern
    dq2table = radkern2/real(ikern)
    select case(ndim)
      case(1)
         cnormk = 20./3.
      case(2)
         cnormk = 2016./(127.*pi)
      case(3)
         cnormk = 210./(17.*pi)
    end select
    do i=0,ikern
       q2 = i*dq2table
       q = sqrt(q2)
       if (q < 1.0) then
          wkern(i)     = 3.*q**7/140. - q2*q2*q2/16. + q2*q2/8. - 7.*q2/40. + 31./280. 
          grwkern(i)   = q*(0.15*q2*q2*q - 0.375*q2*q2 + 0.5*q2 - 0.35) 
          grgrwkern(i) = 0.9*q2*q2*q - 1.875*q2*q2 + 1.5*q2 - 0.35 
          fsoft(i)     = q*(54.*q**7 - 175.*q2*q2*q2 + 450.*q2*q2 - 882.*q2 + 930.)/510. 
          potensoft(i) = q**9/85. - 35.*q**8/816. + 5.*q2*q2*q2/34. - 147.*q2*q2/340. + &
                         31.*q2/34. - 635./408. 
          dphidh(i)    = -2.*q**9/17. + 105.*q**8/272. - 35.*q2*q2*q2/34. + 147.*q2*q2/68. - &
                         93.*q2/34. + 635./408. 
       elseif (q < 2.0) then
          wkern(i)     = -q**7/140. + q2*q2*q2/16. - q2*q2*q/5. + q2*q2/4. - q2/5. + 4./35. 
          grwkern(i)   = q*(-0.05*q2*q2*q + 0.375*q2*q2 - 1.0*q2*q + 1.0*q2 - 0.4) 
          grgrwkern(i) = -0.3*q2*q2*q + 1.875*q2*q2 - 4.0*q2*q + 3.0*q2 - 0.4 
          fsoft(i)     = (-18.*q**10 + 175.*q**9 - 630.*q**8 + 900.*q**7 - 1008.*q2*q2*q + &
                         960.*q2*q - 2.)/(510.*q2) 
          potensoft(i) = (-16.*q**10 + 175.*q**9 - 720.*q**8 + 1200.*q**7 - 2016.*q2*q2*q + &
                         3840.*q2*q - 6400.*q + 16.)/(4080.*q) 
          dphidh(i)    = 2.*q**9/51. - 105.*q**8/272. + 24.*q**7/17. - 35.*q2*q2*q2/17. + &
                         42.*q2*q2/17. - 48.*q2/17. + 80./51. 
       else
          wkern(i)     = 0.0 
          grwkern(i)   = 0.0 
          grgrwkern(i) = 0.0 
          fsoft(i)     = 1./q2 
          potensoft(i) = -1.0/q 
          dphidh(i)    = 0.0 
       endif
    enddo

 case(95)
!
!--twice integrated M5 (auto-generated by kernels.py)
!
    kernellabel = 'Twice integrated M_5' 

    radkern = max(radkern,  2.5)
    radkern2 = radkern*radkern
    dq2table = radkern2/real(ikern)
    select case(ndim)
      case(1)
         cnormk = 6./23.
      case(2)
         cnormk = 516096./(947039.*pi)
      case(3)
         cnormk = 84./(227.*pi)
    end select
    do i=0,ikern
       q2 = i*dq2table
       q = sqrt(q2)
       if (q < 0.5) then
          wkern(i)     = q**8/8. - 5.*q2*q2*q2/8. + 115.*q2*q2/64. - 1199.*q2/384. + &
                         35783./14336. 
          grwkern(i)   = q**7 - 15.*q2*q2*q/4. + 115.*q2*q/16. - 1199.*q/192. 
          grgrwkern(i) = 7.0*q2*q2*q2 - 18.75*q2*q2 + 21.5625*q2 - 6.24479166666667 
          fsoft(i)     = q*(80640.*q**8 - 492800.*q2*q2*q2 + 1821600.*q2*q2 - 4431504.*q2 + &
                         5904195.)/4794240. 
          potensoft(i) = 21.*q**10/12485. - 35.*q**8/2724. + 115.*q2*q2*q2/1816. - &
                         8393.*q2*q2/36320. + 35783.*q2/58112. - 947039./697344. 
          dphidh(i)    = -21.*q**10/1135. + 105.*q**8/908. - 805.*q2*q2*q2/1816. + &
                         8393.*q2*q2/7264. - 107349.*q2/58112. + 947039./697344. 
       elseif (q < 1.5) then
          wkern(i)     = -q**8/12. + 4.*q**7/7. - 5.*q2*q2*q2/4. + q2*q2*q/3. + 55.*q2*q2/32. - &
                         599.*q2/192. + 17891./7168. 
          grwkern(i)   = -2.*q**7/3. + 4.*q2*q2*q2 - 15.*q2*q2*q/2. + 5.*q2*q2/3. + 55.*q2*q/8. &
                         - 599.*q/96. 
          grgrwkern(i) = -14.*q2*q2*q2/3. + 24.*q2*q2*q - 75.*q2*q2/2. + 20.*q2*q/3. + &
                         165.*q2/8. - 599./96. 
          fsoft(i)     = (-53760.*q**11 + 405504.*q**10 - 985600.*q**9 + 295680.*q**8 + &
                         1742400.*q**7 - 4427808.*q2*q2*q + 5904030.*q2*q + &
                         1.)/(4794240.*q2) 
          potensoft(i) = (-21504.*q**11 + 180224.*q**10 - 492800.*q**9 + 168960.*q**8 + &
                         1161600.*q**7 - 4427808.*q2*q2*q + 11808060.*q2*q - 26043545.*q - &
                         4.)/(19176960.*q) 
          dphidh(i)    = 14.*q**10/1135. - 64.*q**9/681. + 105.*q**8/454. - 16.*q**7/227. - &
                         385.*q2*q2*q2/908. + 4193.*q2*q2/3632. - 53673.*q2/29056. + &
                         473519./348672. 
       elseif (q < 2.5) then
          wkern(i)     = q**8/48. - 2.*q**7/7. + 25.*q2*q2*q2/16. - 25.*q2*q2*q/6. + &
                         625.*q2*q2/128. - 3125.*q2/768. + 78125./28672. 
          grwkern(i)   = q**7/6. - 2.*q2*q2*q2 + 75.*q2*q2*q/8. - 125.*q2*q2/6. + 625.*q2*q/32. &
                         - 3125.*q/384. 
          grgrwkern(i) = 7.*q2*q2*q2/6. - 12.*q2*q2*q + 375.*q2*q2/8. - 250.*q2*q/3. + &
                         1875.*q2/32. - 3125./384. 
          fsoft(i)     = (26880*q**11 - 405504.*q**10 + 2464000.*q**9 - 7392000.*q**8 + &
                         9900000.*q**7 - 11550000.*q2*q2*q + 12890625.*q2*q - &
                         177145.)/(9588480.*q2) 
          potensoft(i) = (21504*q**11 - 360448.*q**10 + 2464000.*q**9 - 8448000.*q**8 + &
                         13200000.*q**7 - 23100000.*q2*q2*q + 51562500.*q2*q - &
                         107421875.*q + 1417160.)/(76707840.*q) 
          dphidh(i)    = -7.*q**10/2270. + 32.*q**9/681. - 525.*q**8/1816. + 200.*q**7/227. - &
                         4375.*q2*q2*q2/3632. + 21875.*q2*q2/14528. - 234375.*q2/116224. + &
                         1953125./1394688. 
       else
          wkern(i)     = 0.0 
          grwkern(i)   = 0.0 
          grgrwkern(i) = 0.0 
          fsoft(i)     = 1./q2 
          potensoft(i) = -1.0/q 
          dphidh(i)    = 0.0 
       endif
    enddo

 case(96)
!
!--triple-integrated M4 (auto-generated by kernels.py)
!
    kernellabel = 'Triple-integrated M_4' 

    radkern = max(radkern,  2.0)
    radkern2 = radkern*radkern
    dq2table = radkern2/real(ikern)
    select case(ndim)
      case(1)
         cnormk = 420./17.
      case(2)
         cnormk = 31680./(511.*pi)
      case(3)
         cnormk = 1575./(31.*pi)
    end select
    do i=0,ikern
       q2 = i*dq2table
       q = sqrt(q2)
       if (q < 1.0) then
          wkern(i)     = -q**9/420. + q**8/128. - q2*q2*q2/48. + 7.*q2*q2/160. - 31.*q2/560. + &
                         127./4032. 
          grwkern(i)   = -3.*q**8/140. + q**7/16. - q2*q2*q/8. + 7.*q2*q/40. - 31.*q/280. 
          grgrwkern(i) = -6.*q**7/35. + 7.*q2*q2*q2/16. - 5.*q2*q2/8. + 21.*q2/40. - 31./280. 
          fsoft(i)     = q*(-1320.*q**9 + 4725.*q**8 - 15400.*q2*q2*q2 + 41580.*q2*q2 - &
                         73656.*q2 + 69850.)/32736. 
          potensoft(i) = -5.*q**11/1364. + 315.*q**10/21824. - 175.*q**8/2976. + &
                         105.*q2*q2*q2/496. - 9.*q2*q2/16. + 3175.*q2/2976. - &
                         17885./10912. 
          dphidh(i)    = 15.*q**11/341. - 315.*q**10/1984. + 525.*q**8/992. - 735.*q2*q2*q2/496. &
                         + 45.*q2*q2/16. - 3175.*q2/992. + 17885./10912. 
       elseif (q < 2.0) then
          wkern(i)     = q**9/1260. - q**8/128. + q**7/35. - q2*q2*q2/24. + q2*q2/20. - &
                         2.*q2/35. + 2./63. 
          grwkern(i)   = q**8/140. - q**7/16. + q2*q2*q2/5. - q2*q2*q/4. + q2*q/5. - 4.*q/35. 
          grgrwkern(i) = 2.*q**7/35. - 7.*q2*q2*q2/16. + 6.*q2*q2*q/5. - 5.*q2*q2/4. + 3.*q2/5. &
                         - 4./35. 
          fsoft(i)     = (440*q**12 - 4725.*q**11 + 19008.*q**10 - 30800.*q**9 + 47520.*q**7 - &
                         76032.*q2*q2*q + 70400.*q2*q - 32.)/(32736.*q2) 
          potensoft(i) = (80*q**12 - 945.*q**11 + 4224.*q**10 - 7700.*q**9 + 15840.*q**7 - &
                         38016.*q2*q2*q + 70400.*q2*q - 107520.*q + 64.)/(65472.*q) 
          dphidh(i)    = -5.*q**11/341. + 315.*q**10/1984. - 20.*q**9/31. + 525.*q**8/496. - &
                         105.*q2*q2*q2/62. + 90.*q2*q2/31. - 100.*q2/31. + 560./341. 
       else
          wkern(i)     = 0.0 
          grwkern(i)   = 0.0 
          grgrwkern(i) = 0.0 
          fsoft(i)     = 1./q2 
          potensoft(i) = -1.0/q 
          dphidh(i)    = 0.0 
       endif
    enddo

  case(98) 
   
!   
!--Optimal kernel OM5 = 4*M_5 - 5*M_4, correctly normalised 
!
    kernellabel = 'OM_5 = optimal kernel of order 5'    

    radkern = max(radkern,  2.5)
    radkern2 = radkern*radkern
    dq2table = radkern2/real(ikern)
    select case(ndim)
      case(1)
         cnormk = 1./66.  !--=1/(7!):OK
      case(2)
         cnormk = 24./(863.*pi)
      case(3)
         cnormk = 1./(60*pi)
    end select  
    do i=0,ikern 
       q2 = i*dq2table
       q  = sqrt(q2)
       if (q.lt.0.5) then
          wkern(i) = 4.*((2.5-q)**4 - 5.*(1.5-q)**4 + 10.*(0.5-q)**4) &
                     -20.*(1. - 1.5*q2 + 0.75*q*q2) 
          grwkern(i) = 4.*(-4.*((2.5-q)**3 - 5.*(1.5-q)**3 + 10.*(0.5-q)**3)) &
                      -20.*(-3.*q+ 2.25*q2)
          grgrwkern(i) = 4.*(12.*((2.5-q)**2 - 5.*(1.5-q)**2 + 10.*(0.5-q)**2)) &
                        -20.*(-3. + 4.5*q) 
       elseif ((q.ge.0.5) .and. (q.lt.1.)) then
          wkern(i) =4.*((2.5-q)**4 - 5.*(1.5-q)**4) &
                      -20.*(1. - 1.5*q2 + 0.75*q*q2)
          grwkern(i) = 4.*(-4.*((2.5-q)**3 - 5.*(1.5-q)**3)) &
                      -20.*(-3.*q+ 2.25*q2)
          grgrwkern(i) =4.*(12.*((2.5-q)**2 - 5.*(1.5-q)**2) ) &
                      -20.*(-3. + 4.5*q)
       elseif ((q.ge.1.) .and. (q.lt.1.5)) then
          wkern(i) =4.*((2.5-q)**4 - 5.*(1.5-q)**4) &
                      -20.*(0.25*(2.-q)**3)
          grwkern(i) = 4.*(-4.*((2.5-q)**3 - 5.*(1.5-q)**3)) &
                      -20.*(-0.75*(2.-q)**2)
          grgrwkern(i) =  4.*(12.*((2.5-q)**2 - 5.*(1.5-q)**2) )   &
                      -20.*(1.5*(2.-q))         
       elseif ((q.ge.1.5).and. (q.lt.2.)) then
          wkern(i)=4.*((2.5-q)**4) &
                      -20.*(0.25*(2.-q)**3)
          grwkern(i) = 4.*(-4.*((2.5-q)**3)) &
                      -20.*(-0.75*(2.-q)**2)
          grgrwkern(i) = 4.*(12.*((2.5-q)**2)) &
                      -20.*(1.5*(2.-q))
       elseif ((q.ge.2.).and. (q.lt.2.5)) then
          wkern(i) =4.*((2.5-q)**4)
          grwkern(i) =4.*(-4.*((2.5-q)**3))
          grgrwkern(i) =4.*(12.*((2.5-q)**2))
       else
          wkern(i) = 0.
          grwkern(i) = 0.
          grgrwkern(i) = 0.
       endif
    enddo
    
  case default
    
!--default is cubic spline (see monaghan 1992; monaghan & lattanzio 1985)
!   
    kernellabel = 'M_4 cubic spline'    
  
    radkern = max(radkern,  2.0)      ! interaction radius of kernel
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

!!----------------------------------------------------------------------
!! function to interpolate linearly from drag kernel tables
!!----------------------------------------------------------------------
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
!--calculate slope for w and interpolate for each
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
!! same but for kernel *and* modified kernel in anticlumping term
!!----------------------------------------------------------------------
subroutine interpolate_kernels2(q2,w,walt,gradw,gradwalt)
 implicit none
 integer :: index,index1
 real, intent(in) :: q2
 real, intent(out) :: w,walt,gradw,gradwalt
 real :: dxx,dwdx,dwaltdx,dgrwdx,dgrwaltdx
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
 walt = wijalt(index)
 dwaltdx =  (wijalt(index1)-walt)*ddq2table
 walt = walt + dwaltdx*dxx

 gradwalt = grwijalt(index)
 dgrwaltdx =  (grwijalt(index1)-gradwalt)*ddq2table
 gradwalt = gradwalt + dgrwaltdx*dxx

end subroutine interpolate_kernels2

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
