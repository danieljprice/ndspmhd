!!------------------------------------------------------------------------!!
!!                                                                        !!
!!  Generic setup for 2D MHD shock tubes - simple smoothing               !!
!!                                                                        !!
!!  Gives particles a variable separation so the mass per SPH particle    !!
!!  is constant. Shock is smoothed slightly (simple smoothing).           !!
!!  Based on a setup used by J.J. Monaghan.                               !!
!!                                                                        !!
!!  Note for all MHD setups, only the magnetic field should be setup      !!
!!  Similarly the thermal energy is setup even if using total energy.     !!
!!                                                                        !!
!!------------------------------------------------------------------------!!

SUBROUTINE setup
!
!--include relevant global variables
!
 USE dimen_mhd
 USE debug
 USE loguns
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
 INTEGER :: i,j,ntot,ipart,nparty
 REAL :: xcentre,rholeft,rhoright,uuleft,uuright
 REAL :: prleft, prright, exx, delta
 REAL :: massp,Bxinit,Byleft,Byright,Bzleft,Bzright,enleft,enright
 REAL :: v2left,v2right,B2left,B2right
 REAL :: const,gam1,dsmooth
 REAL :: xstart,ystart,psepleft
 REAL, DIMENSION(ndimV) :: vleft,vright
!
!--set boundaries
!            	    
 ibound = 2	! reflecting ghosts
 nbpts = 0	! must use fixed particles if inflow/outflow at boundaries
 xmin(1) = 0.0
 xmax(1) = 1.0
 xmin(2) = 0.0
 xmax(2) = 0.1
!
!--setup parameters
!
 gam1 = gamma - 1.
 const = SQRT(4.*pi)
 dsmooth = 10.      
 rholeft = 1.0
 rhoright = 1.0
 prleft = 1.0
 prright = 1.0
 vleft(1) = 0.	!0.2
 vright(1) = 0.
 vleft(2) = 0.	!0.01
 vright(2) = 0.
 vleft(3) = 0.	!0.5
 vright(3) = 0.
 IF (imhd.NE.0) THEN
    Bxinit = 4./const
    Byleft = 4./const
    Byright = 4.0/const
    Bzleft = 1.0/const
    Bzright = 1.0/const
 ELSE
    Bxinit = 0.
    Byleft = 0.
    Byright = 0.
 ENDIF
 uuleft = prleft/(gam1*rholeft)
 uuright = prright/(gam1*rhoright)
 massp = rhoright*psep**ndim      
 nparty = INT((xmax(2)-xmin(2))/psep)	
 psepleft = (massp/rholeft)**power
 xstart = 0.5*psepleft      ! offset from xmin boundary
 ystart = 0.5*psepleft      ! offset from ymin boundary
!
! allocate memory to start with (need to guess the particle number)
!
 ntot = 1000*nparty
 CALL alloc(ntot)
!
!--setup particles using the parameters given
!      	
 ipart = 0
 
 DO j=1,nparty
    ipart = ipart + 1 	 ! create the first particle on the left
    xcentre = (xmax(1) + xmin(1))/2.0		! this specifies the location of the 
    					! discontinuity in x at this y

    xin(1,ipart) = xmin(1) + xstart    
    xin(1,ipart + 1) = xmin(1) + psep*rhoright/rholeft + xstart

    DO i=ipart,ipart+1
       xin(2,i) = xmin(2) + (j-1)*psep	+ ystart! y position
       rhoin(i) = rholeft
       uuin(i) = uuleft
       enin(i) = enleft
       pmass(i) = massp
       velin(:,i) = vleft(:)	! overwrite vx
       hhin(i) = hfact*(pmass(i)/rholeft)**hpower
       IF (imhd.NE.0) THEN
          Bin(1,i) = Bxinit
          Bin(2,i) = Byleft
          Bin(3,i) = Bzleft
       ENDIF      
    ENDDO

    ipart = ipart + 1   ! already made second particle on the left
    DO WHILE (xin(1,ipart).LT.xmax(1))
       ipart = ipart + 1         
       delta = (xin(1,ipart-1) - xcentre)/psep
       IF (delta.GT.dsmooth) THEN
          rhoin(ipart) = rhoright
          uuin(ipart) = uuright 
          velin(:,ipart) = vright(:) 
          Bin(2,ipart) = Byright
          Bin(3,ipart) = Bzright
       ELSEIF (delta.LT.-dsmooth) THEN
          rhoin(ipart) = rholeft
          uuin(ipart) = uuleft
          velin(:,ipart) = vleft(:)
          Bin(2,ipart) = Byleft
          Bin(3,ipart) = Bzleft
       ELSE
          exx = exp(delta)
          rhoin(ipart) = (rholeft + rhoright*exx)/(1.0 +exx)
!         uuin(ipart) = (uuleft + uuright*exx)/(1.0 + exx)
          uuin(ipart) = (prleft + prright*exx)/((1.0 + exx)*gam1*rhoin(i))
          velin(:,ipart) = (vleft(:) + vright(:)*exx)/(1.0 + exx)
          Bin(2,ipart) = (Byleft + Byright*exx)/(1.0 + exx)
          Bin(3,ipart) = (Bzleft + Bzright*exx)/(1.0 + exx)
       ENDIF
       xin(1,ipart) = xin(1,ipart-2) + 2.*(massp/rhoin(ipart-1))**hpower
       xin(2,ipart) = xmin(2) + (j-1)*psep + ystart
       pmass(ipart) = massp
       hhin(ipart) = hfact*(massp/rhoin(ipart))**hpower
       Bin(1,ipart) = Bxinit
    ENDDO
    ipart = ipart - 1      ! remove the one that is over the xmax boundary

 ENDDO	! over nparty

 npart = ipart
 ntotal = npart
!
!--allow for tracing flow
!
 IF (trace) WRITE(iprint,*) '  Exiting subroutine setup'
            
 RETURN
END
