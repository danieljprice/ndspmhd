!!------------------------------------------------------------------------!!
!!                                                                        !!
!!  Generic setup for 1D MHD shock tubes - equal smoothing around xcentre !!
!!  this version for use with ND code (references xin(1,j) not xin(j) )   !!
!!                                                                        !!
!!  Gives particles a variable separation so the mass per SPH particle    !!
!!  is constant. Shock is smoothed slightly (equal about xcentre).        !!
!!  Based on a setup used by J.J. Monaghan.                               !!
!!                                                                        !!
!!  Note for all MHD setups, only the magnetic field should be setup      !!
!!  Similarly the thermal energy is setup even if using total energy.     !!
!!                                                                        !!
!!------------------------------------------------------------------------!!

SUBROUTINE setup
 USE dimen_mhd
 USE debug
 USE loguns
 
 USE artvi
 USE bound
 USE eos
 USE options
 USE part
 USE part_in
 USE setup_params
!
!--define local variables
!      
 IMPLICIT NONE
 INTEGER :: i,j,ntot
 REAL :: xcentre,rholeft,rhoright,rhohalf,uuleft,uuright
 REAL :: prleft, prright, exx, delta
 REAL :: massp,Bxinit,Byleft,Byright,Bzleft,Bzright
 REAL :: enleft,enright,B2i
 REAL :: dx0,dx1,dxhalf,xhalf,psepleft
 REAL :: const,gam1,dsmooth
 REAL :: exxal, alphamax
 REAL, DIMENSION(ndimV) :: vleft,vright 
 CHARACTER(LEN=20) :: shkfile
!
!--allow for tracing flow
!
 IF (trace) WRITE(iprint,*) ' Entering subroutine setup'
 IF (ndim.GT.1) STOP 'shock setup only for 1D problems: recompile'
!
!--set boundaries
!            	    
 ibound =1
 nbpts = 6	! must use fixed particles if inflow/outflow at boundaries.
 xmin(1) = -0.5
 xmax(1) = 0.5
 xcentre = (xmax(1) + xmin(1))/2.0
!
!--setup parameters
!
 gam1 = gamma - 1.
 const = SQRT(4.*pi)
 alphamax = 0.5
 dsmooth = 0. 
 rholeft = 1.0
 rhoright = 0.125
 prleft = 1.0
 prright = 0.1
 vleft(1) = 0.
 vright(1) = 0.
 vleft(2) = 0.
 vright(2) = 0.
 vleft(3) = 0.
 vright(3) = 0.
 IF (imhd.NE.0) THEN
    Bxinit = 0.75	!2./const
    Byleft = 1.0	!3.6/const
    Byright = -1.0	!4.0/const
    Bzleft = 0.		!2./const
    Bzright = 0.	!2./const
 ELSE
    Bxinit = 0.
    Byleft = 0.
    Byright = 0.
 ENDIF
!
!-- ...or read shock parameters from the .shk file
!
 shkfile = rootname(1:LEN_TRIM(rootname))//'.shk'
 
 OPEN (UNIT=ireadf, FILE=shkfile,STATUS='old',FORM='formatted')
  READ(ireadf,*) rholeft,rhoright
  READ(ireadf,*) prleft,prright
  READ(ireadf,*) vleft(1),vright(1)
  READ(ireadf,*) vleft(2),vright(2)
  READ(ireadf,*) vleft(3),vright(3)
  READ(ireadf,*) Bxinit
  READ(ireadf,*) Byleft,Byright
  READ(ireadf,*) Bzleft,Bzright
 CLOSE (UNIT=ireadf)
 WRITE(iprint,10) rholeft,rhoright,prleft,prright,vleft(1),vright(1),   &
                  vleft(2),vright(2),vleft(3),vright(3),		&
		  Bxinit,Byleft,Byright,Bzleft,Bzright
 
10 FORMAT( ' 1D shock: rho L: ',f16.9,' R: ',f16.9,/,   &
           '           pr  L: ',f16.9,' R: ',f16.9,/,   &
           '           vx  L: ',f16.9,' R: ',f16.9,/,   &
           '           vy  L: ',f16.9,' R: ',f16.9,/,   &
           '           vz  L: ',f16.9,' R: ',f16.9,/,   &
           '           Bx   : ',f16.9,/,   &	   	   	   
           '           By  L: ',f16.9,' R: ',f16.9,/,   &
           '           Bz  L: ',f16.9,' R: ',f16.9,/)



 IF (ABS(gam1).GT.1.e-3) THEN
    uuleft = prleft/(gam1*rholeft)
    uuright = prright/(gam1*rhoright)
 ELSE
    uuleft = 3.*prleft/(2.*rholeft)
    uuright = 3.*prright/(2.*rhoright)
 ENDIF
 massp = rhoright*psep      
 psepleft = massp/rholeft			
!
! allocate memory to start with
!
 ntot = 2000
 CALL alloc(ntot)
!
!--setup the particles using the parameters given
! 
 xin(1,1) = xmin(1) + 0.5*psepleft
 xin(1,2) = xmin(1) + psep*rhoright/rholeft + 0.5*psepleft
 DO i=1,2
    rhoin(i) = rholeft
    uuin(i) = uuleft
    enin(i) = enleft
    pmass(i) = massp
    velin(:,i) = vleft(:)
    hhin(i) = hfact*(pmass(i)/rhoin(i))**hpower
    alphain(i) = alphamin
    IF (imhd.NE.0) THEN
       Bin(1,i) = Bxinit
       Bin(2,i) = Byleft
       Bin(3,i) = Bzleft
    ENDIF 
 ENDDO
 
 j = 2     
 DO WHILE (xin(1,j).LT.xmax(1))
    dx0 = massp/rhoin(j)
    xhalf = xin(1,j) + 0.5*dx0  !x at the mid point
    delta = xhalf/psep
    j = j + 1 
    IF (delta.lt.-dsmooth) THEN
       rhohalf = rholeft
    ELSEIF (delta.gt.dsmooth) THEN
       rhohalf = rhoright
    ELSE
       exx = exp(delta)
       rhohalf = (rholeft+ rhoright*exx)/(1+exx)
    ENDIF
!---------------------------------------------------	    
!     calculate half steps then final step
!---------------------------------------------------	    
    dxhalf = massp/rhohalf
    dx1 = 2.*dxhalf - dx0
    xin(1,j) = xhalf + 0.5*dx1
    delta = (xin(1,j)-xcentre)/psep
    IF (delta.lt.-dsmooth) THEN
       rhoin(j) = rholeft
       enin(j) = enleft
       uuin(j) = uuleft
       velin(:,j) = vleft(:)
       Bin(2,j) = Byleft
       Bin(3,j) = Bzleft
       alphain(j) = alphamin
    ELSEIF (delta.gt.dsmooth) THEN
       rhoin(j) = rhoright
       enin(j) = enright
       uuin(j) = uuright
       velin(:,j) = vright(:)
       Bin(2,j) = Byright
       Bin(3,j) = Bzright
       alphain(j) = alphamin
    ELSE
       exx = exp(delta)
       rhoin(j) = (rholeft + rhoright*exx)/(1.0 +exx)
       enin(j) = (enleft + enright*exx)/(1.0 + exx)
       uuin(j) = (uuleft + uuright*exx)/(1.0 + exx)
       
       exxal = exp(-(xin(1,j)/(0.5*dsmooth*psep))**2)
       alphain(j) = alphamax*exxal + alphamin
!       IF (prleft.GT.prright) THEN
!          uuin(j) = (prleft + prright*exx)/((1.0 + exx)*gam1*rhoin(i))
!       ELSE
!          uuin(j) = (prright + prleft*exx)/((1.0 + exx)*gam1*rhoin(i))
!       ENDIF
      velin(:,j) = (vleft(:) + vright(:)*exx)/(1.0 + exx)
!       IF (delta.LT.0.) THEN
!          velin(:,j) = vleft(:)
!       ELSE
!          velin(:,j) = vright(:)
!       ENDIF  
       Bin(2,j) = (Byleft + Byright*exx)/(1.0 + exx)
       Bin(3,j) = (Bzleft + Bzright*exx)/(1.0 + exx)
    ENDIF
    pmass(j) = massp
    hhin(j) = hfact*(pmass(j)/rhoin(j))**hpower	!rhohalf
    Bin(1,j) = Bxinit
 ENDDO

 npart = j-1
 ntotal = npart
!
!--if using ghosts adjust outer boundary to lie half a particle spacing away
!
 IF (ibound.GT.1) THEN
    xmax(1) = xin(1,npart) + 0.5*(xin(1,npart)-xin(1,npart-1))
 ENDIF

!
!--allow for tracing flow
!
 IF (trace) WRITE(iprint,*) '  Exiting subroutine setup'
    
 RETURN
END SUBROUTINE setup
