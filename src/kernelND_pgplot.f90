!!-----------------------------------------------------------------
!! Sets up the tables for the kernel
!! Returns kernel, and derivative.
!!
!!  Note for the ND case, the normalisation constants are right
!!  only for the cubic spline in > 1D.
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
 INTEGER :: i,j,iteration,iC
 REAL :: q,q2,q4,cnormk
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

 iplot = .true. 	! plot kernel using PGPLOT

!
!--set kernel related parameters
!      
 radkern = 2.0		! interaction radius of kernel
 radkern2 = radkern*radkern
 IF (ndim.EQ.1) THEN	! normalization constant of kernel
    cnormk = 0.66666666666
 ELSEIF (ndim.EQ.2) THEN
    cnormk = 10./(7.*pi)
 ELSEIF (ndim.EQ.3) THEN
    cnormk = 1./pi
 ENDIF
 dq2table = radkern*radkern/REAL(ikern)
!
!--setup kernel table (this is the usual spline-based kernel)
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
!
!--plot cubic spline kernel
!
 IF (iplot) THEN
    CALL PGBEGIN (0,'?',1,1)
    CALL PGENV (0.0,3.0,-3.5,1.7,0,1)
    CALL PGLABEL ('r/h','  ','SPH cubic spline kernel')  
    CALL PGSCI(2)
    CALL PGSLS(1)
    print*,'plotting cubic spline...'
    print*,'will crash here if compiled in double precision'
    CALL PGLINE(ikern+1,dqkern(0:ikern),wij(0:ikern))
    grwijplot(0:ikern) = grwij(0:ikern)	!*dqkern(0:ikern)
    CALL PGSLS(2)
    CALL PGLINE(ikern+1,dqkern(0:ikern),grwijplot(0:ikern))
    CALL PGSLS(3)
    CALL PGLINE(ikern+1,dqkern(0:ikern),grgrwij(0:ikern))
    CALL PGSCI(1)
 ENDIF

!
!--this is the quintic spline
!
 IF (idivBzero.EQ.6 .OR. idivBzero.EQ.9) THEN
 
 cnormk = 0.5/60.
 radkern = 3.0
 radkern2 = radkern*radkern
 dq2table = radkern2/REAL(ikern)
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

!
!--plot quintic spline
!
  print*,' plotting quintic spline...'
  CALL PGSCI(3)
  CALL PGSLS(1)
  CALL PGLINE(ikern+1,dqkern(0:ikern),wij(0:ikern))
  grwijplot(0:ikern) = grwij(0:ikern)	!*dqkern(0:ikern)
  CALL PGSLS(2)
  CALL PGLINE(ikern+1,dqkern(0:ikern),grwijplot(0:ikern))
  CALL PGSLS(3)
  CALL PGLINE(ikern+1,dqkern(0:ikern),grgrwij(0:ikern))
  CALL PGSCI(1)
 
 ENDIF

!
!--this is the squashed quintic spline
!
 IF (idivBzero.EQ.7) THEN
 
 cnormk = 0.5/1280.
 radkern = 2.0
 radkern2 = radkern*radkern
 dq2table = radkern2/REAL(ikern)
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
!
!--plot squashed quintic spline
!
!  CALL PGSCI(4)
!  CALL PGSLS(1)
!  CALL PGLINE(ikern+1,dqkern(0:ikern),wij(0:ikern))
!  grwijplot(0:ikern) = grwij(0:ikern)	!*dqkern(0:ikern)
!  CALL PGSLS(2)
!  CALL PGLINE(ikern+1,dqkern(0:ikern),grwijplot(0:ikern))
!  CALL PGSLS(3)
!  CALL PGLINE(ikern+1,dqkern(0:ikern),grgrwij(0:ikern))
!  CALL PGSCI(1) 

 ENDIF

!
!--this is the quintic spline from Bonet & Kulesegaram
!
 IF (idivBzero.EQ.8) THEN

  cnormk = 1./16.
  radkern = 2.0
  radkern2 = radkern*radkern
  dq2table = radkern2/REAL(ikern)
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
 
 ENDIF
!
!--this is Dan's quintic kernel (build your own quintic)
!
  DO iteration = 1,3
!  PRINT*,iteration,' Enter beta, alpha'
!  READ*,beta,alpha
  
!  DO iC = 1,3
!      alpha = beta + (radkern-beta)*(iC)/11
!      C = 3*REAL(iC-1)
!      PRINT*,' C = ',C
      C = 0.
      gamma = 0.
      PRINT*,' alpha = ',alpha, ' beta = ',beta

 IF (idivBzero.EQ.9) THEN
  radkern = 2.0
  radkern2 = radkern*radkern
  dq2table = radkern2/REAL(ikern)
!  gamma = 0.5
!  beta = 1.0	!2./3.
!  alpha = 1.5	!4./3.
!  C =0.
!  A = -6.
!  B = 25.
!  C = 1.

  A =  (-radkern**4 + (radkern**2 + C*gamma**2)*beta**2)		&
      /(alpha**2*(alpha**2-beta**2))      
  B = -(radkern**4 + A*alpha**4 + C*gamma**4)/(beta**4)
  cnormk = 3./(A*alpha**6 + B*beta**6 + C*gamma**6 + 64.)	! for radkern = 2 and 1D
!  print*,'cnormk = ',cnormk,' A,B = ',A,B
  
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

 ENDIF

 IF (idivBzero.EQ.10) THEN	! modified quintic
 
 cnormk = 0.5/60.
 radkern = 2.0
 radkern2 = radkern*radkern
 dq2table = radkern2/REAL(ikern)
 alpha = 1.4	!(9. + SQRT(6.))/5.	!4./3.
 W0 = cnormk*((3.-alpha)**5 - 6.*(2.-alpha)**5)
 grW0 = cnormk*(-5.*(3.-alpha)**4 + 30.*(2.-alpha)**4)
 
 A = cnormk*(20.*(3-alpha)**3 - 120*(2-alpha)**3)
 B = cnormk*(-60*(3-alpha)**2 + 360*(2-alpha)**2)
 C = cnormk*(120*(3-alpha) - 720*(2-alpha))
 aa = (B*alpha - 2.*(A + B))/(alpha-2.)**3
 bb = (B - 3.*aa*(alpha**2 - 4.))/(2.*(alpha-2.))
 cc = -12.*aa - 4.*bb
 dd = -8.*aa - 4.*bb - 2.*cc
 ee = grW0 - 0.25*aa*alpha**4 - (bb*alpha**3)/3. -0.5*cc*alpha**2 -dd*alpha
! ff = W0 - 
 print*,'a,b,c,d = ',aa,bb,cc,dd,ee
 
 DO i=0,ikern         
    q2 = i*dq2table
    q4 = q2*q2
    q = SQRT(q2)
    dqkern(i) = q		! to plot kernel
    term1 = -5.*(3.-q)**4.
    IF (q.LT.1.0) THEN
       wij(i) = cnormk*(66.-60.*q2 + 30.*q4 - 10.*q4*q)
       grwij(i) = cnormk*(term1 + 30*(2.-q)**4 - 75.*(1.-q)**4)
       grgrwij(i) = cnormk*(20.*(3.-q)**3 - 120.*(2.-q)**3	&
     	    			 + 300.*(1.-q)**3)
    ELSEIF ((q.GE.1.0).AND.(q.LT.alpha)) THEN
       wij(i) = cnormk*((3.-q)**5 - 6.*(2.-q)**5)
       grwij(i) = cnormk*(term1 + 30.*(2.-q)**4)
       grgrwij(i) = cnormk*(20.*(3.-q)**3 - 120.*(2.-q)**3)
    ELSEIF ((q.GE.alpha).AND.(q.LE.radkern)) THEN
       wij(i) = cnormk*(0.)
       grwij(i) = 0.25*aa*q**4 + (bb*q**3)/3. +0.5*cc*q**2 +dd*q + ee
       grgrwij(i) = aa*q**3 + bb*q**2 + cc*q + dd
       print*,' grgrwij = ',q,grgrwij(i),grwij(i),grW0
!       read*
    ELSE
       wij(i) = 0.0
       grwij(i) = 0.0
       grgrwij(i) = 0.
    ENDIF
 ENDDO 
 
 ENDIF

!
!--plot kernel before applying Joe's correction term
!
 IF (iplot) THEN
   IF (idivBzero.EQ.7) THEN
      CALL PGLABEL ('r/h','  ','SPH squashed quintic kernel')
   ELSEIF (idivBzero.EQ.6) THEN
      CALL PGLABEL ('r/h','  ','SPH quintic kernel') 
   ELSEIF (idivBzero.EQ.8) THEN
      CALL PGLABEL ('r/h','  ','Dan''s quintic kernel')  
   ELSEIF (idivBzero.EQ.9) THEN
      CALL PGLABEL ('r/h','  ','Dan''s DIY-quintic kernel')   
   ELSE
      CALL PGLABEL ('r/h','  ','SPH cubic spline kernel')  
   ENDIF   
   CALL PGSLS(iteration)
   CALL PGLINE(ikern+1,dqkern(0:ikern),wij(0:ikern))
      
   grwijplot(0:ikern) = grwij(0:ikern)	!*dqkern(0:ikern)
!   CALL PGENV (0.0,2.1,minval(grwijplot),0.5,0,1)
!   CALL PGLABEL ('r/h','dW\dab\u/dr','SPH Kernel')
   CALL PGSLS(2)
   CALL PGLINE(ikern+1,dqkern(0:ikern),grwijplot(0:ikern))
   CALL PGSLS(3)
   CALL PGLINE(ikern+1,dqkern(0:ikern),grgrwij(0:ikern))
 ENDIF

!
!--use different kernel for gradient (this one is used by Thomas & Couchman '92
!  to prevent particle clumping)
!
 IF ((imhd.NE.0).AND.(idivBzero.EQ.2)) THEN
    DO i=0,ikern
       q2 = i*dq2table
       q = SQRT(q2)
       IF (q.LT.2./3.) THEN
          grwij(i) = -cnormk
       ELSEIF ((q.GE.2./3.).AND.(q.LT.1.0)) THEN
          grwij(i) = cnormk*(-3.*q+ 2.25*q2)
       ELSEIF ((q.GE.1.0).AND.(q.LE.2.0)) THEN
	  grwij(i) = cnormk*(-0.75*(2.-q)**2.)
       ELSE
	  grwij(i) = 0.0
       ENDIF
    ENDDO
    WRITE(iprint,*) 'setkern: using anti-clumping kernel'
 ELSEIF ((idivBzero.EQ.5)) THEN
    j = NINT((1./hfact)**2./dq2table)
    DO i=0,ikern
       q2 = i*dq2table
       q = SQRT(q2)
!       dqkern(i) = q		! to plot kernel
       IF (i.EQ.j) WRITE(iprint,*) ' j = ',j,q,q2,wij(j)
       grwij(i) = grwij(i)*(1. - eps*(wij(i)/wij(j))**neps)
    ENDDO
    WRITE(iprint,*) 'setkern: using Dan''s anti-clumping kernel' 
 ENDIF	 

 IF (iplot) THEN
  grwijplot(0:ikern) = grwij(0:ikern)	!*dqkern(0:ikern)
  CALL PGSLS(3)
  CALL PGLINE(ikern+1,dqkern(0:ikern),grwijplot(0:ikern))      
 ENDIF

! READ*

 ENDDO	! my parameter loops
! ENDDO

 IF (iplot) CALL PGEND
!
!--the variable ddq2table is used in interpolate_kernel
! 
 ddq2table = 1./dq2table 
      
END SUBROUTINE setkern
