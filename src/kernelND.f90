!!-----------------------------------------------------------------
!! Sets up the tables for the kernel
!! Returns kernel, and derivative.
!!
!! Default kernel is the cubic spline, however several alternatives
!! are given, namely:
!!  1) cubic spline with constant gradient for r/h < 2/3
!!  2) quintic spline (max r/h = 3)
!!  3) squashed quintic (max r/h = 2)
!!  4) another squashed quintic (max r/h = 2)
!!  5) a general class of quintic splines (max r/h = 2)
!!
!!  Note in the ND case, the normalisation constants are right
!!  only for the cubic spline in > 1D.
!!
!! Changes log:
!! 4/8/03 - select case and variable ikernel
!!-----------------------------------------------------------------

SUBROUTINE setkern  
 USE dimen_mhd
 USE debug
 USE loguns
 USE kernel
 USE options
 USE setup_params	! for hfact in my kernel
 USE anticlumping
 IMPLICIT NONE			!  define local variables
 INTEGER :: ikernel
 INTEGER :: i,j,iteration,iC
 REAL :: q,q2,q4,cnormk,term1
 REAL, DIMENSION(0:ikern) :: dqkern,grwijplot,grgrwij 	! only to plot kernel
 REAL :: term1,term2,term3,term4
 REAL :: dterm1,dterm2,dterm3,dterm4
 REAL :: ddterm1,ddterm2,ddterm3,ddterm4
 REAL :: alpha,beta,gamma,A,B,C,aa,bb,cc,dd,ee,W0,grW0
 LOGICAL :: iplot
!
!--allow for tracing flow
!
 IF (trace) WRITE(iprint,*) ' Entering subroutine setkern'
!
!--set choice of kernel (this could be read in as a parameter)
!
 ikernel = 0

 iplot = .false. 	! plot kernel using PGPLOT

 SELECT CASE(ikernel)
  CASE(1)
!
!--this is the usual spline based kernel modified for r/h < 2/3 to
!   prevent particles from clumping (see Thomas & Couchman '92)
!      
    WRITE(iprint,*) 'Using anti-clumping kernel'  

    radkern = 2.0		! interaction radius of kernel
    radkern2 = radkern*radkern
    dq2table = radkern*radkern/REAL(ikern)
    SELECT CASE(ndim)
      CASE(1)
        cnormk = 0.66666666666	! normalisation constant
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
       dqkern(i) = q		! to plot kernel
       IF (q.LT.1.0) THEN
          wij(i) = cnormk*(1. - 1.5*q2 + 0.75*q*q2)
          IF (q.LT.2./3.) THEN
             grwij(i) = -cnormk
	     grgrwij(i) = 0.
	  ELSE
	     grwij(i) = cnormk*(-3.*q+ 2.25*q2)
             grgrwij(i) = cnormk*(-3. + 4.5*q)
	  ENDIF
       ELSEIF ((q.GE.1.0).AND.(q.LE.2.0)) THEN
          wij(i) = cnormk*(0.25*(2.-q)**3.)
          grwij(i) = cnormk*(-0.75*(2.-q)**2.)
          grgrwij(i) = cnormk*(1.5*(2.-q))
       ELSE
          wij(i) = 0.0
          grwij(i) = 0.0
          grgrwij(i) = 0.
       ENDIF
    ENDDO

  CASE(2)
!
!--this is the quintic spline (see e.g. Morris 1996, PhD thesis)
!
    WRITE(iprint,*) 'Using quintic spline kernel'  
  
    radkern = 3.0
    radkern2 = radkern*radkern
    dq2table = radkern2/REAL(ikern)
    SELECT CASE(ndim)
      CASE(1) 
       cnormk = 0.5/60.
      CASE DEFAULT
       STOP ' normalisation constant not defined in kernel '
    END SELECT
    DO i=0,ikern         
       q2 = i*dq2table
       q4 = q2*q2
       q = SQRT(q2)
       dqkern(i) = q		! to plot kernel
       term1 = -5.*(3.-q)**4.
       IF (q.LT.1.0) THEN
          wij(i) = cnormk*(66.-60.*q2 + 30.*q4 - 10.*q4*q)
          grwij(i) = cnormk*(term1 + 30*(2.-q)**4. - 75.*(1.-q)**4.)
          grgrwij(i) = cnormk*(20.*(3.-q)**3. - 120.*(2.-q)**3.	&
     	    			 + 300.*(1.-q)**3.)
       ELSEIF ((q.GE.1.0).AND.(q.LT.2.0)) THEN
          wij(i) = cnormk*((3.-q)**5. - 6.*(2.-q)**5.)
          grwij(i) = cnormk*(term1 + 30*(2.-q)**4.)
          grgrwij(i) = cnormk*(20.*(3.-q)**3. - 120.*(2.-q)**3.)
       ELSEIF ((q.GE.2.0).AND.(q.LE.3.0)) THEN
          wij(i) = cnormk*((3.-q)**5.)
          grwij(i) = cnormk*(term1)
          grgrwij(i) = cnormk*(20.*(3.-q)**3.)
       ELSE
          wij(i) = 0.0
          grwij(i) = 0.0
          grgrwij(i) = 0.
       ENDIF
    ENDDO

  CASE(3)
!
!--this is my squashed quintic spline (quintic rescaled to radius 2)
!
    WRITE(iprint,*) 'Using Dan''s squashed quintic spline kernel'    
  
    radkern = 2.0
    radkern2 = radkern*radkern
    dq2table = radkern2/REAL(ikern)
    SELECT CASE(ndim)
      CASE(1)
       cnormk = 0.5/1280.
      CASE DEFAULT
       STOP ' normalisation constant not defined in kernel '
    END SELECT
    DO i=0,ikern         
       q2 = i*dq2table
       q4 = q2*q2
       q = SQRT(q2)
       dqkern(i) = q		! to plot kernel
       term1 = -15.*(6.-3.*q)**4.
       IF (q.LT.2./3.) THEN
          wij(i) = cnormk*((6.-3.*q)**5. - 6.*(4.-3.*q)**5. + 15.*(2.-3.*q)**5.)
          grwij(i) = cnormk*(term1 + 90*(4.-3.*q)**4. - 225.*(2.-3.*q)**4.)
          grgrwij(i) = cnormk*(180.*(6.-3.*q)**3. - 1080.*(4.-3.*q)**3.	&
     	    			 + 2700.*(2.-3.*q)**3.)
       ELSEIF ((q.GE.2./3.).AND.(q.LT.4./3.)) THEN
          wij(i) = cnormk*((6.-3.*q)**5. - 6.*(4.-3.*q)**5.)
          grwij(i) = cnormk*(term1 + 90*(4.-3.*q)**4.)
          grgrwij(i) = cnormk*(180.*(6.-3.*q)**3. - 1080.*(4.-3.*q)**3.)
       ELSEIF ((q.GE.4./3.).AND.(q.LE.2.0)) THEN
          wij(i) = cnormk*((6.-3.*q)**5.)
          grwij(i) = cnormk*(term1)
          grgrwij(i) = cnormk*(180.*(6.-3.*q)**3.)
       ELSE
          wij(i) = 0.0
          grwij(i) = 0.0
          grgrwij(i) = 0.
       ENDIF
    ENDDO

  CASE(4)
!
!--this is the squashed quintic spline from Bonet & Kulesegaram
!
    WRITE(iprint,*) 'Using BK squashed quintic spline kernel'    
   
    radkern = 2.0
    radkern2 = radkern*radkern
    dq2table = radkern2/REAL(ikern)
    SELECT CASE(ndim)
      CASE(1)
       cnormk = 1./16.
      CASE DEFAULT
       STOP ' normalisation constant not defined in kernel '
    END SELECT
    DO i=0,ikern
      q2 = i*dq2table
      q = SQRT(q2)
      dqkern(i) = q		! to plot kernel
      IF (q.LT.1.0) THEN
         wij(i) = cnormk*((2.-q)**5 - 16.*(1.-q)**5)
         grwij(i) = cnormk*(-5.*(2.-q)**4 + 80.*(1.-q)**4.)
         grgrwij(i) = cnormk*(20.*(2.-q)**3 - 320.*(1.-q)**3.)
      ELSEIF ((q.GE.1.0).AND.(q.LE.2.0)) THEN
         wij(i) = cnormk*(2.-q)**5
         grwij(i) = cnormk*(-5.*(2.-q)**4)
         grgrwij(i) = cnormk*(20.*(2.-q)**3)
      ELSE
         wij(i) = 0.0
         grwij(i) = 0.0
         grgrwij(i) = 0.
      ENDIF
    ENDDO  

  CASE(5)
!
!--this is the Do-It-Yourself quintic kernel (general class of quintic splines)
!
    WRITE(iprint,*) 'Using the Do-It-Yourself quintic kernel'    

    radkern = 2.0
    radkern2 = radkern*radkern
    dq2table = radkern2/REAL(ikern)
    SELECT CASE(ndim)
      CASE(1)
!   gamma = 0.5
!   beta = 1.0	!2./3.
!   alpha = 1.5	!4./3.
!   C =0.
!   A = -6.
!   B = 25.
!   C = 1.    
    A =  (-radkern**4 + (radkern**2 + C*gamma**2)*beta**2)		&
        /(alpha**2*(alpha**2-beta**2))      
    B = -(radkern**4 + A*alpha**4 + C*gamma**4)/(beta**4)
    cnormk = 3./(A*alpha**6 + B*beta**6 + C*gamma**6 + 64.)	! for radkern = 2 and 1D
    print*,'cnormk = ',cnormk,' A,B = ',A,B
      CASE DEFAULT
       STOP ' normalisation constant not defined in kernel '      
    END SELECT
  
    DO i=0,ikern         
      q2 = i*dq2table
      q = SQRT(q2)
      dqkern(i) = q		! to plot kernel
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
         wij(i) = cnormk*(term1 + A*term2  + B*term3 + C*term4)
         grwij(i) = cnormk*(dterm1 + A*dterm2 + B*dterm3 + C*dterm4)
         grgrwij(i) = cnormk*(ddterm1 + A*ddterm2 + B*ddterm3 + C*ddterm4)    
      ELSEIF ((q.GE.gamma).AND.(q.LT.beta)) THEN
         wij(i) = cnormk*(term1 + A*term2  + B*term3)
         grwij(i) = cnormk*(dterm1 + A*dterm2 + B*dterm3)
         grgrwij(i) = cnormk*(ddterm1 + A*ddterm2 + B*ddterm3)          
      ELSEIF ((q.GE.beta).AND.(q.LT.alpha)) THEN
         wij(i) = cnormk*(term1 + A*term2)
         grwij(i) = cnormk*(dterm1 + A*dterm2)
         grgrwij(i) = cnormk*(ddterm1 + A*ddterm2)
      ELSEIF ((q.GE.alpha).AND.(q.LT.radkern)) THEN
         wij(i) = cnormk*term1
         grwij(i) = cnormk*dterm1
         grgrwij(i) = cnormk*ddterm1
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
    WRITE(iprint,*) 'Using the cubic spline kernel'    
  
    radkern = 2.0		! interaction radius of kernel
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
       dqkern(i) = q		! to plot kernel
       IF (q.LT.1.0) THEN
          wij(i) = cnormk*(1. - 1.5*q2 + 0.75*q*q2)
          grwij(i) = cnormk*(-3.*q+ 2.25*q2)
          grgrwij(i) = cnormk*(-3. + 4.5*q)
       ELSEIF ((q.GE.1.0).AND.(q.LE.2.0)) THEN
          wij(i) = cnormk*(0.25*(2.-q)**3.)
          grwij(i) = cnormk*(-0.75*(2.-q)**2.)
          grgrwij(i) = cnormk*(1.5*(2.-q))
       ELSE
          wij(i) = 0.0
          grwij(i) = 0.0
          grgrwij(i) = 0.
       ENDIF
    ENDDO

 END SELECT

!
!--this is to see what effect the anticlumping term has on the kernel
!
 IF (idivBzero.EQ.5) THEN
    j = NINT((1./hfact)**2./dq2table)
    DO i=0,ikern
       q2 = i*dq2table
       q = SQRT(q2)
!       dqkern(i) = q		! to plot kernel
       IF (i.EQ.j) WRITE(iprint,*) ' j = ',j,q,q2,wij(j)
       grwij(i) = grwij(i) + eps*(wij(i)/wij(j))**neps
    ENDDO
    WRITE(iprint,*) ' ...modified with the anti-clumping term' 
 ENDIF	
!
!--the variable ddq2table is used in interpolate_kernel
! 
 ddq2table = 1./dq2table 
      
END SUBROUTINE setkern
