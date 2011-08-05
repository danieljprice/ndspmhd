!!-----------------------------------------------------------------
!! Sets up the tables for the kernel
!! Returns kernel, and derivative.
!!
!! Default kernel is the cubic spline, however several alternatives
!! are given, namely:
!!  1) cubic spline
!!  2) cubic spline with constant gradient for r/h < 2/3
!!  3) quintic spline (max r/h = 3)
!!  4) squashed quintic (max r/h = 2)
!!  5) another squashed quintic (max r/h = 2)
!!  6) a general class of quintic splines (max r/h = 2)
!!
!!  Note in the ND case, the normalisation constants are right
!!  only for the cubic and quintic splines in > 1D.
!!
!!-----------------------------------------------------------------

SUBROUTINE setkern  
 USE dimen_mhd
 USE debug
 USE loguns
 USE kernel
 USE kernelextra
 USE options
 USE setup_params   ! for hfact in my kernel
 USE anticlumping
 IMPLICIT NONE         !  define local variables
 INTEGER :: i,j,npower
 REAL :: q,q2,q4,cnormk,cnormkaniso
 REAL :: term1,term2,term3,term4
 REAL :: dterm1,dterm2,dterm3,dterm4
 REAL :: ddterm1,ddterm2,ddterm3,ddterm4
 REAL :: alpha,beta,gamma,A,B,C,wdenom,wint
!
!--allow for tracing flow
!
 IF (trace) WRITE(iprint,*) ' Entering subroutine setkern'
!
!--set choice of kernel (this could be read in as a parameter)
!
 cnormk = 0.0

 SELECT CASE(ikernel)

  CASE(2)
!
!--Quartic spline
!  
    kernelname = 'Quartic spline'    

    radkern = 2.5
    radkern2 = radkern*radkern
    dq2table = radkern2/REAL(ikern)
    SELECT CASE(ndim)
      CASE(1)
         cnormk = 1./24.
      CASE(2)
         cnormk = 96./(1199*pi)
      CASE(3)
         cnormk = 1./(20*pi)
    END SELECT  
    DO i=0,ikern
       q2 = i*dq2table
       q = SQRT(q2)
       IF (q.LT.0.5) THEN
          wij(i) = (2.5-q)**4 - 5.*(1.5-q)**4 + 10.*(0.5-q)**4
          grwij(i) = -4.*((2.5-q)**3 - 5.*(1.5-q)**3 + 10*(0.5-q)**4)
          grgrwij(i) = 12.*((2.5-q)**2 - 5.*(1.5-q)**2 + 10*(0.5-q)**2)
       ELSEIF (q.LT.1.5) THEN
          wij(i) = (2.5-q)**4 - 5.*(1.5-q)**4
          grwij(i) = -4.*((2.5-q)**3 - 5.*(1.5-q)**3)
          grgrwij(i) = 12.*((2.5-q)**2 - 5.*(1.5-q)**2)   
       ELSEIF (q.LT.2.5) THEN
          wij(i) = (2.5-q)**4
          grwij(i) = -4.*((2.5-q)**3)
          grgrwij(i) = 12.*((2.5-q)**2)
       ELSE
          wij(i) = 0.
          grwij(i) = 0.
          grgrwij(i) = 0.
       ENDIF
    ENDDO 

  CASE(3)
!
!--this is the M_6 quintic spline (see e.g. Morris 1996, PhD thesis)
!
    kernelname = 'Quintic spline'  
    radkern = 3.0
    radkern2 = radkern*radkern
    dq2table = radkern2/REAL(ikern)
    SELECT CASE(ndim)
      CASE(1) 
       cnormk = 1./120.
      CASE(2)
       cnormk = 7./(478*pi)
      CASE(3)
       cnormk = 1./(120.*pi)
    END SELECT
    DO i=0,ikern         
       q2 = i*dq2table
       q4 = q2*q2
       q = SQRT(q2)
       term1 = -5.*(3.-q)**4.
       IF (q.LT.1.0) THEN
          wij(i) = 66.-60.*q2 + 30.*q4 - 10.*q4*q
          grwij(i) = term1 + 30*(2.-q)**4. - 75.*(1.-q)**4.
          grgrwij(i) = 20.*(3.-q)**3. - 120.*(2.-q)**3. + 300.*(1.-q)**3.
       ELSEIF ((q.GE.1.0).AND.(q.LT.2.0)) THEN
          wij(i) = (3.-q)**5. - 6.*(2.-q)**5.
          grwij(i) = term1 + 30*(2.-q)**4.
          grgrwij(i) = 20.*(3.-q)**3. - 120.*(2.-q)**3.
       ELSEIF ((q.GE.2.0).AND.(q.LE.3.0)) THEN
          wij(i) = (3.-q)**5.
          grwij(i) = term1
          grgrwij(i) = 20.*(3.-q)**3.
       ELSE
          wij(i) = 0.0
          grwij(i) = 0.0
          grgrwij(i) = 0.
       ENDIF
    ENDDO

  CASE(5,6,7)
!
!--this is the Do-It-Yourself quintic kernel (general class of quintic splines)
!
    radkern = 2.0
    radkern2 = radkern*radkern
    dq2table = radkern2/REAL(ikern)

    gamma = 0.
    IF (ikernel.EQ.5) THEN
       beta = 0.5
       alpha = 1.7   !!1.4
       kernelname = 'New quintic (1)'    
    ELSEIF (ikernel.EQ.6) THEN
       !--very poor on sound waves
       beta = 0.7
       alpha = 1.5
       kernelname = 'New quintic (2)'    
    ELSE
    !--match to cubic spline, ie W''(0) = -2
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
    ENDIF
    C =0.
    A =  (-radkern**4 + (radkern**2 + C*gamma**2)*beta**2)      &
        /(alpha**2*(alpha**2-beta**2))      
    B = -(radkern**4 + A*alpha**4 + C*gamma**4)/(beta**4)
    print*,'matching points = ',beta,alpha
    SELECT CASE(ndim)
      CASE(1)
        cnormk = 3./(A*alpha**6 + B*beta**6 + C*gamma**6 + radkern**6)   ! for radkern = 2 and 1D
        print*,'1D cnormk = ',cnormk,' A,B = ',A,B
      CASE(2)
        cnormk = 42./(2.*pi*(A*alpha**7 + B*beta**7 + C*gamma**7 + radkern**7))
        print*,'2D cnormk = ',cnormk,' A,B = ',A,B,beta,alpha
      CASE DEFAULT
       write(iprint,666)
       stop  
    END SELECT
  
    DO i=0,ikern         
      q2 = i*dq2table
      q = SQRT(q2)
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
      IF (q.LT.gamma) THEN
         wij(i) = term1 + A*term2  + B*term3 + C*term4
         grwij(i) = dterm1 + A*dterm2 + B*dterm3 + C*dterm4
         grgrwij(i) = ddterm1 + A*ddterm2 + B*ddterm3 + C*ddterm4   
      ELSEIF ((q.GE.gamma).AND.(q.LT.beta)) THEN
         wij(i) = term1 + A*term2  + B*term3
         grwij(i) = dterm1 + A*dterm2 + B*dterm3
         grgrwij(i) = ddterm1 + A*ddterm2 + B*ddterm3       
      ELSEIF ((q.GE.beta).AND.(q.LT.alpha)) THEN
         wij(i) = term1 + A*term2
         grwij(i) = dterm1 + A*dterm2
         grgrwij(i) = ddterm1 + A*ddterm2
      ELSEIF ((q.GE.alpha).AND.(q.LT.radkern)) THEN
         wij(i) = term1
         grwij(i) = dterm1
         grgrwij(i) = ddterm1
      ELSE
         wij(i) = 0.0
         grwij(i) = 0.0
         grgrwij(i) = 0.
      ENDIF
    ENDDO

  CASE (8)
!
!--(1-r^2)^3
!   
    npower = 5
    write(kernelname,"(a,i1)") '(1-r^2)^',npower     

    radkern = 2.0
    radkern2 = radkern*radkern
    dq2table = radkern2/REAL(ikern)
    SELECT CASE(ndim)
      CASE(1)
         SELECT CASE(npower)
         CASE(1)
            cnormk = 0.5*3./(2.*radkern**3)
         CASE(2)
            cnormk = 0.5*15./(8.*radkern**5)
         CASE(3)
            cnormk = 0.5*35./(16.*radkern**7)
         CASE(4)
            cnormk = 0.5*315./(128.*radkern**9)
         CASE(5)
            cnormk = 0.5*693./(256.*radkern**11)
         CASE(6)
            cnormk = 0.5*3003./(1024.*radkern**13)
         CASE(7)
            cnormk = 0.5*6435./(2048.*radkern**15)    
         CASE(8)
            cnormk = 0.5*109395./(32768.*radkern**17)    
         CASE DEFAULT
            cnormk = 0.
         END SELECT    
      CASE(2)
         cnormk = 3./(64.*pi)
      CASE(3)
         cnormk = 105./(4096.*pi)
    END SELECT  
    DO i=0,ikern
       q2 = i*dq2table
       q = SQRT(q2)
       IF (q.LT.radkern) THEN
          wij(i) = (radkern**2-q2)**npower
          grwij(i) = -2.*npower*q*(radkern**2-q2)**(npower-1)
          grgrwij(i) = (4.*npower*(npower-1)*q2*(radkern**2-q2)**(npower-2) &
                      - 2.*npower*(radkern**2-q2)**(npower-1))
       ELSE
          wij(i) = 0.
          grwij(i) = 0.
          grgrwij(i) = 0.
       ENDIF
    ENDDO
    
   CASE (9)
!
!--(1-r)^n - a peaked kernel (deriv non-zero at origin) - truly awful
!   
    npower = 4
    write(kernelname,"(a,i1)") '(2-r)^',npower     

    radkern = 2.0
    radkern2 = radkern*radkern
    dq2table = radkern2/REAL(ikern)
    SELECT CASE(ndim)
      CASE(1)
         cnormk = 0.5*(npower+1)/radkern**(npower+1)
      CASE(2,3)
         STOP 'normalisation const not defined in kernel'
    END SELECT  
    DO i=0,ikern
       q2 = i*dq2table
       q = SQRT(q2)
       IF (q.LT.radkern) THEN
          wij(i) = (radkern-q)**npower
          grwij(i) = -npower*(radkern-q)**(npower-1)
          grgrwij(i) = npower*(npower-1)*(radkern-q)**(npower-2)
       ELSE
          wij(i) = 0.
          grwij(i) = 0.
          grgrwij(i) = 0.
       ENDIF
    ENDDO   

  CASE(10)
!
!--Gaussian
!  
    kernelname = 'Gaussian'    

    radkern = 5.0
    radkern2 = radkern*radkern
    dq2table = radkern2/REAL(ikern)
    SELECT CASE(ndim)
      CASE(1)
         cnormk = 1./sqrt(pi)
      CASE(2)
         cnormk = 1./pi
      CASE(3)
         cnormk = 1./(pi*sqrt(pi))
    END SELECT  
    DO i=0,ikern
       q2 = i*dq2table
       q = SQRT(q2)
       IF (q.LT.radkern) THEN
          wij(i) = exp(-q2)
          grwij(i) = -2.*q*wij(i)
          grgrwij(i) = -2.*q*grwij(i) - 2.*wij(i)
       ELSE
          wij(i) = 0.
          grwij(i) = 0.
          grgrwij(i) = 0.
       ENDIF
    ENDDO

  CASE(11)
!
!--this is the usual spline based kernel modified for r/h < 2/3 to
!   prevent particles from clumping (see Thomas & Couchman '92)
!      
    kernelname = 'Thomas & Couchman anti-clumping'
    radkern = 2.0      ! interaction radius of kernel
    radkern2 = radkern*radkern
    dq2table = radkern*radkern/REAL(ikern)
    SELECT CASE(ndim)
      CASE(1)
        cnormk = 0.66666666666   ! normalisation constant
      CASE(2)
        cnormk = 10./(7.*pi)
      CASE(3)
        cnormk = 1./pi
    END SELECT
   
    DO i=0,ikern
       q2 = i*dq2table
       q = SQRT(q2)
       IF (q.LT.1.0) THEN
          wij(i) = 1. - 1.5*q2 + 0.75*q*q2
          IF (q.LT.2./3.) THEN
             grwij(i) = -1.
             grgrwij(i) = 0.
          ELSE
             grwij(i) = -3.*q+ 2.25*q2
             grgrwij(i) = -3. + 4.5*q
          ENDIF
       ELSEIF ((q.GE.1.0).AND.(q.LE.2.0)) THEN
          wij(i) = 0.25*(2.-q)**3.
          grwij(i) = -0.75*(2.-q)**2.
          grgrwij(i) = 1.5*(2.-q)
       ELSE
          wij(i) = 0.0
          grwij(i) = 0.0
          grgrwij(i) = 0.
       ENDIF
    ENDDO
  
  CASE(12)
!
!--this is the squashed quintic spline from Bonet & Kulesegaram
!
    kernelname = 'BK squashed quintic spline'    
   
    radkern = 2.0
    radkern2 = radkern*radkern
    dq2table = radkern2/REAL(ikern)
    SELECT CASE(ndim)
      CASE(1)
       cnormk = 1./16.
      CASE DEFAULT
       write(iprint,666)
       stop
    END SELECT
    DO i=0,ikern
      q2 = i*dq2table
      q = SQRT(q2)
      IF (q.LT.1.0) THEN
         wij(i) = (2.-q)**5 - 16.*(1.-q)**5
         grwij(i) = -5.*(2.-q)**4 + 80.*(1.-q)**4.
         grgrwij(i) = 20.*(2.-q)**3 - 320.*(1.-q)**3.
      ELSEIF ((q.GE.1.0).AND.(q.LE.2.0)) THEN
         wij(i) = (2.-q)**5
         grwij(i) = -5.*(2.-q)**4
         grgrwij(i) = 20.*(2.-q)**3
      ELSE
         wij(i) = 0.0
         grwij(i) = 0.0
         grgrwij(i) = 0.
      ENDIF
    ENDDO  
    
  CASE(13)
!
!--this is the Lucy kernel (to 2h not to 1h)
!
    kernelname = 'Lucy kernel'    
   
    radkern = 2.0
    radkern2 = radkern*radkern
    dq2table = radkern2/REAL(ikern)
    SELECT CASE(ndim)
      CASE(1)
       cnormk = 3./5.
      CASE DEFAULT
       write(iprint,666)
       STOP
    END SELECT
    DO i=0,ikern
      q2 = i*dq2table
      q = SQRT(q2)
      IF (q.LT.2.0) THEN
         wij(i) = (1.+0.5*q)*(1.-0.5*q)**3
         grwij(i) = (0.5*(1.-0.5*q)**3 - 1.5*(1.+0.5*q)*(1.-0.5*q)**2)
         grgrwij(i) = (-1.5*(1.-0.5*q)**2 + 1.5*(1.+0.5*q)*(1.-0.5*q))
      ELSE
         wij(i) = 0.0
         grwij(i) = 0.0
         grgrwij(i) = 0.
      ENDIF
    ENDDO
        
  CASE DEFAULT  
!
!--default is cubic spline (see Monaghan 1992; Monaghan & Lattanzio 1985)
!   
    kernelname = 'Cubic spline'    
  
    radkern = 2.0      ! interaction radius of kernel
    radkern2 = radkern*radkern
    dq2table = radkern*radkern/REAL(ikern)    
    SELECT CASE(ndim)
      CASE(1)
        cnormk = 0.66666666666
      CASE(2)
        cnormk = 10./(7.*pi)
      CASE(3)
        cnormk = 1./pi
    END SELECT
!
!--setup kernel table
!   
    DO i=0,ikern
       q2 = i*dq2table
       q = SQRT(q2)
       IF (q.LT.1.0) THEN
          wij(i) = 1. - 1.5*q2 + 0.75*q*q2
          grwij(i) = -3.*q+ 2.25*q2
          grgrwij(i) = -3. + 4.5*q
       ELSEIF ((q.GE.1.0).AND.(q.LE.2.0)) THEN
          wij(i) = 0.25*(2.-q)**3.
          grwij(i) = -0.75*(2.-q)**2.
          grgrwij(i) = 1.5*(2.-q)
       ELSE
          wij(i) = 0.0
          grwij(i) = 0.0
          grgrwij(i) = 0.
       ENDIF
    ENDDO

 END SELECT

!
!--calculate modified kernel to use on anisotropic forces
!
 IF (ianticlump.eq.1) THEN
!
!--this is Joe's standard anticlumping term
!
    j = NINT((1./hfact)**2/dq2table)
    !!j = NINT((1./1.5)**2/dq2table)
    wdenom = wij(j)

    SELECT CASE(neps) ! integral of W^{n+1} dV
    CASE(3)
       wint = 122559./160160.
    CASE(4)
       wint = 200267./292864.
    CASE(5)
       wint = 825643615./1324331008.
    END SELECT   
    cnormkaniso = cnormk
    !!cnormkaniso = 1./(1./cnormk + eps/(neps+1.)*wint/wdenom**neps)
    
    DO i=0,ikern
       q2 = i*dq2table
       q = SQRT(q2)
       IF (i.EQ.j) WRITE(iprint,*) ' j = ',j,q,q2,wij(j)
       grwijaniso(i) = cnormkaniso*grwij(i)*(1. - 0.5*eps*(wij(i)/wdenom)**neps)
       wijaniso(i) = cnormkaniso*wij(i)*(1. - 0.5*eps/(neps+1.)*(wij(i)/wdenom)**neps) 
       grgrwijaniso(i) = cnormkaniso*(grgrwij(i)*  &
                  (1. - 0.5*eps*(wij(i)/wdenom)**neps) &
                 - 0.5*eps*neps*wij(i)**(neps-1.)/wdenom**neps*grwij(i)*grwij(i))   
    ENDDO

    print*,'cnormk = ',cnormkaniso,eps,neps
    kernelname=TRIM(kernelname)//' & anti-clumping cubic spline' 
 
 ELSEIF (ianticlump.eq.2) THEN
!
!--this is a modified version of the cubic spline
!
    beta = eps
    print*,'enter beta (0.6-0.99) , currently = ',beta
    !!read*,beta
    A = -0.25*(2.-beta)**2/(1.-beta)**2
    if (beta.eq.0.) then
       B = 0.
    else
       B = -(A+1.)/beta**2
    endif
    cnormkaniso = cnormk !!!0.5/(1. + 0.25*(A + B*beta**4))
    DO i=0,ikern
       q2 = i*dq2table
       q = SQRT(q2)
       if (q.lt.beta) then
          wijaniso(i) = cnormkaniso*(0.25*(2.-q)**3 + A*(1.-q)**3 + B*(beta-q)**3)
          grwijaniso(i) = -3.*cnormkaniso*(0.25*(2.-q)**2 + A*(1.-q)**2 + B*(beta-q)**2)
          grgrwijaniso(i) = 6.*cnormkaniso*(0.25*(2.-q) + A*(1.-q) + B*(beta-q))       
       elseif(q.lt.1.) then
          wijaniso(i) = cnormkaniso*(0.25*(2.-q)**3 + A*(1.-q)**3)
          grwijaniso(i) = -3.*cnormkaniso*(0.25*(2.-q)**2 + A*(1.-q)**2)
          grgrwijaniso(i) = 6.*cnormkaniso*(0.25*(2.-q) + A*(1.-q))      	  
       else
          wijaniso(i) = cnormkaniso*wij(i)
          grwijaniso(i) = cnormkaniso*grwij(i)
	  grgrwijaniso(i) = cnormkaniso*grgrwij(i)
       endif

    ENDDO

    print*,'cnormk = ',cnormkaniso,eps,neps
    kernelname=TRIM(kernelname)//' & modified cubic spline'   
 
 ENDIF
!
!--normalise kernel
!
 wij = cnormk*wij
 grwij = cnormk*grwij
 grgrwij = cnormk*grgrwij

666 FORMAT(/,'ERROR!!! Normalisation constant not defined in kernel',/)
!
!--write kernel name to log file
!
 WRITE(iprint,10) TRIM(kernelname)
10 FORMAT(/,' Smoothing kernel = ',a,/)
!
!--the variable ddq2table is used in interpolate_kernel
!
 ddq2table = 1./dq2table 

END SUBROUTINE setkern
