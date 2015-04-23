!!------------------------------------------------------------------------!!
!!                                                                        !!
!!  Generic setup for 1D MHD shock tubes - simple smoothing               !!
!!                                                                        !!
!!  Gives particles a variable separation so the mass per SPH particle    !!
!!  is constant. Shock is smoothed slightly (simple smoothing).           !!
!!  Setup & smoothing described in Monaghan (1997) JCP 136, 298           !!
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
 USE setup_params
!
!--define local variables
!      
 IMPLICIT NONE
 INTEGER :: i,j,ntot
 REAL :: xcentre,densleft,densright,uuleft,uuright
 REAL :: prleft, prright, exx, delta
 REAL :: massp,Bxinit,Byleft,Byright,Bzleft,Bzright,enleft,enright
 REAL :: v2left,v2right
 REAL :: B2left,B2right,B2densleft,B2densright,B2densi,B2y,B2ydensleft,B2ydensright
 REAL :: const,gam1,dsmooth
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
 ibound = 1	! reflecting ghosts
 nbpts = 6	! must use fixed particles if inflow/outflow at boundaries
 xmin(1) = -0.5
 xmax(1) = 0.5
 xcentre = (xmax(1) + xmin(1))/2.0
!
!--setup parameters here...
!
 gam1 = gamma - 1.
 const = 1./SQRT(4.*pi)
 dsmooth = 20.      
 densleft = 1.0
 densright = 1.0
 prleft = 1.0
 prright = 1.0
 vleft(1) = 36.87	!0.2
 vright(1) = -36.87
 vleft(2) = -0.155	!0.01
 vright(2) = 0.
 vleft(3) = -0.0386	!0.5
 vright(3) = 0.
 IF (imhd.NE.0) THEN
    Bxinit = 4.*const
    Byleft = 4.*const
    Byright = 4.0*const
    Bzleft = 1.0*const
    Bzright = 1.0*const
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
  READ(ireadf,*) densleft,densright
  READ(ireadf,*) prleft,prright
  READ(ireadf,*) vleft(1),vright(1)
  READ(ireadf,*) vleft(2),vright(2)
  READ(ireadf,*) vleft(3),vright(3)
  READ(ireadf,*) Bxinit
  READ(ireadf,*) Byleft,Byright
  READ(ireadf,*) Bzleft,Bzright
 CLOSE (UNIT=ireadf)
 WRITE(iprint,10) densleft,densright,prleft,prright,vleft(1),vright(1),   &
                  vleft(2),vright(2),vleft(3),vright(3),		&
		  Bxinit,Byleft,Byright,Bzleft,Bzright
 
10 FORMAT( '1D shock: dens L: ',f8.3,' R: ',f8.3,/,   &
           '          pr  L: ',f8.3,' R: ',f8.3,/,   &
           '          vx  L: ',f8.3,' R: ',f8.3,/,   &
           '          vy  L: ',f8.3,' R: ',f8.3,/,   &
           '          vz  L: ',f8.3,' R: ',f8.3,/,   &
           '          Bx   : ',f8.3,/,   &	   	   	   
           '          By  L: ',f8.3,' R: ',f8.3,/,   &
           '          Bz  L: ',f8.3,' R: ',f8.3,/)
  
 IF (ABS(gam1).GT.1.e-3) THEN
    uuleft = prleft/(gam1*densleft)
    uuright = prright/(gam1*densright)
 ELSE
    uuleft = 3.*prleft/(2.*densleft)
    uuright = 3.*prright/(2.*densright)
 ENDIF
 massp = densright*psep      
 
 B2left = Bxinit**2 + Byleft**2 + Bzleft**2
 B2right = Bxinit**2 + Byright**2 + Bzright**2
 enleft = uuleft + 0.5*DOT_PRODUCT(vleft,vleft) + 0.5*B2left/densleft
 enright = uuright + 0.5*DOT_PRODUCT(vright,vright) + 0.5*B2right/densright
 B2densleft = (Bxinit**2 + Byleft**2)/densleft
 B2densright = (Bxinit**2 + Byright**2)/densright
 B2ydensleft = (Byleft**2)/densleft
 B2ydensright = (Byright**2)/densright


!
! allocate memory to start with (need to guess the particle number)
!
 ntot = 2000
 CALL alloc(ntot)
!
!--setup particles using the parameters given
!      	
 x(1,1) = xmin(1) + 0.5*massp/densleft
 x(1,2) = xmin(1) + psep*densright/densleft + 0.5*massp/densleft
 DO i=1,2
    dens(i) = densleft
    uu(i) = uuleft
    enin(i) = enleft
    pmass(i) = massp
    vel(:,i) = vleft(:)	! overwrite vx
    IF (imhd.NE.0) THEN
       Bfield(1,i) = Bxinit
       Bfield(2,i) = Byleft
       Bfield(3,i) = Bzleft
    ENDIF      
 ENDDO

 i = 2
 DO WHILE (x(1,i).LT.xmax(1))
    i = i + 1         
    delta = 2.*(x(1,i-1) - xcentre)/psep
    Bfield(1,i) = Bxinit    
    IF (delta.GT.dsmooth) THEN
       dens(i) = densright
       uu(i) = uuright 
       vel(:,i) = vright(:) 
       Bfield(2,i) = Byright
       Bfield(3,i) = Bzright
       enin(i) = enright
    ELSEIF (delta.LT.-dsmooth) THEN
       dens(i) = densleft
       uu(i) = uuleft
       vel(:,i) = vleft(:)
       Bfield(2,i) = Byleft
       Bfield(3,i) = Bzleft
       enin(i) = enleft
    ELSE
       exx = exp(delta)       
       dens(i) = (densleft + densright*exx)/(1.0 +exx)
!       uu(i) = (uuleft + uuright*exx)/(1.0 + exx)
       uu(i) = (prleft + prright*exx)/((1.0 + exx)*gam1*dens(i))
       IF (delta.GT.0.) THEN
          vel(:,i) = vright(:)
       ELSE
          vel(:,i) = vleft(:)
       ENDIF
       Bfield(2,i) = (Byleft + Byright*exx)/(1.0 + exx)
       Bfield(3,i) = (Bzleft + Bzright*exx)/(1.0 + exx)
    ENDIF
    x(1,i) = x(1,i-2) + 2.*massp/dens(i-1)
    pmass(i) = massp
 ENDDO
 
 npart = i-1
 ntotal = npart
!
!--allow for tracing flow
!
 IF (trace) WRITE(iprint,*) '  Exiting subroutine setup'
            
 RETURN
END

!
! use this routine to modify the dump upon code restart
!
subroutine modify_dump()
 implicit none

end subroutine modify_dump
