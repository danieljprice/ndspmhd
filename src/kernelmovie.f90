!!-----------------------------------------------------------------
!! Sets up the tables for the kernel
!! Returns kernel, and derivative.
!!
!!  Note for the ND case, the normalisation constants are right
!!  only for the cubic spline in > 1D.
!!
!!  MUST BE COMPILED IN SINGLE PRECISION
!!-----------------------------------------------------------------

subroutine setkern    
  use dimen_mhd
  use debug
  use loguns
  use kernel
  use options
  use setup_params	! for hfact in my kernel
  use anticlumping
  implicit none			!  define local variables
  integer :: i,j,iteration,iC,MM,PP,nc,nalpha
  real :: q,q2,q4,cnormk
  real, dimension(0:ikern) :: dqkern	! only to plot kernel
  real :: term1,term2,term3,term4,dbeta,dalpha
  real :: dterm1,dterm2,dterm3,dterm4
  real :: ddterm1,ddterm2,ddterm3,ddterm4
  real :: alpha,beta,gamma,A,B,C,aa,bb,cc,dd,ee,W0,grW0
  logical :: iplot
  character(LEN=20) :: string
  !
  !--allow for tracing flow
  !
  if (trace) write(iprint,*) ' Entering subroutine setkern'

  iplot = .true. 	! plot kernel using PGPLOT
  !
  !--this is Dan's quintic kernel (build your own quintic)
  !
  radkern = 2.0
  beta = 0.0
  dbeta = (radkern-beta-0.01)/20
  do iteration = 1,20
     print*,iteration
     beta = beta + dbeta
     alpha = beta
     dalpha = dbeta
     nalpha = int((radkern-beta-0.01)/(dbeta))
     print*,'dalpha = ',dalpha,nalpha, 'dbeta = ',dbeta,(radkern-beta-0.01)/(dbeta)
     do iC = 1,nalpha
        alpha = alpha + dalpha
        C = 0.
        gamma = 0.
        
        print*,' alpha = ',alpha, ' beta = ',beta
        radkern = 2.0
        radkern2 = radkern*radkern
        dq2table = radkern2/real(ikern)
        
        A = (-radkern**4 + (radkern**2 + C*gamma**2)*beta**2)		&
             /(alpha**2*(alpha**2-beta**2))      
        B = -(radkern**4 + A*alpha**4 + C*gamma**4)/(beta**4)
        cnormk = 3./(A*alpha**6 + B*beta**6 + C*gamma**6 + 64.)	! for radkern = 2 and 1D
        !  print*,'cnormk = ',cnormk,' A,B = ',A,B
        
        do i=0,ikern         
           q2 = i*dq2table
           q = sqrt(q2)
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
           if (q.lt.gamma) then
              wij(i) = cnormk*(term1 + A*term2  + B*term3 + C*term4)
              grwij(i) = cnormk*(dterm1 + A*dterm2 + B*dterm3 + C*dterm4)
              grgrwij(i) = cnormk*(ddterm1 + A*ddterm2 + B*ddterm3 + C*ddterm4)    
           elseif ((q.ge.gamma).and.(q.lt.beta)) then
              wij(i) = cnormk*(term1 + A*term2  + B*term3)
              grwij(i) = cnormk*(dterm1 + A*dterm2 + B*dterm3)
              grgrwij(i) = cnormk*(ddterm1 + A*ddterm2 + B*ddterm3)          
           elseif ((q.ge.beta).and.(q.lt.alpha)) then
              wij(i) = cnormk*(term1 + A*term2)
              grwij(i) = cnormk*(dterm1 + A*dterm2)
              grgrwij(i) = cnormk*(ddterm1 + A*ddterm2)
           elseif ((q.ge.alpha).and.(q.lt.radkern)) then
              wij(i) = cnormk*term1
              grwij(i) = cnormk*dterm1
              grgrwij(i) = cnormk*ddterm1
           else
              wij(i) = 0.0
              grwij(i) = 0.0
              grgrwij(i) = 0.
           endif
        enddo
        
     endif
     
     !
     !--plot kernel before applying Joe's correction term
     !
     if (iplot) then
        call PGENV (0.0,3.0,-3.5,1.7,0,1)
        call PGLABEL ('r/h','  ','Dan''s DIY-quintic kernel')
        call PGLINE(ikern+1,dqkern(0:ikern),wij(0:ikern))
        call PGPT(1,beta,wij(nint(beta**2/dq2table)),17)
        call PGPT(1,alpha,wij(nint(alpha**2/dq2table)),17)
        
        call PGSLS(2)
        call PGLINE(ikern+1,dqkern(0:ikern),grwij(0:ikern))
        call PGSLS(3)
        call PGLINE(ikern+1,dqkern(0:ikern),grgrwij(0:ikern))
        call PGSLS(1)
     endif
     
     !
     !--the variable ddq2table is used in interpolate_kernel
     ! 
     ddq2table = 1./dq2table 
     !
     !--calculate dispersion relation for this kernel
     !
     kernelname = ' '
     call kernelstability1D
     MM=nint(alpha*100)
     PP=nint(log10(alpha)-log10(alpha*100))
     call pgnumb(MM,PP,1,string,nc)  
     call pgmtxt('T',-3.0,0.95,1.0,'alpha = '//string(1:nc))
     
     MM=nint(beta*100)
     PP=nint(log10(beta)-log10(beta*100))
     call pgnumb(MM,PP,1,string,nc)      
     call pgmtxt('T',-4.5,0.95,1.0,'beta  = '//string(1:nc))
     
  enddo	! my parameter loops
enddo

if (iplot) call PGEND
!
!--the variable ddq2table is used in interpolate_kernel
! 
ddq2table = 1./dq2table 

end subroutine setkern
