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

!!------------------------------------------------------------------------!!
!!                                                                        !!
!!  Generic setup for 1D MHD shock tubes - equal smoothing around xcentre !!
!!  this version for use with ND code (references x(1,j) not x(j) )   !!
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
 USE setup_params
 USE mem_allocation, only:alloc
!
!--define local variables
!      
 IMPLICIT NONE
 INTEGER :: i,j,ntot,npartleft,npartright
 REAL :: xcentre,densleft,densright,denshalf,uuleft,uuright
 REAL :: prleft, prright, exx, delta
 REAL :: massp,Bxinit,Byleft,Byright,Bzleft,Bzright
 REAL :: dx0,dx1,dxhalf,xhalf,psepleft
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
 ibound =1
 nbpts = 6      ! must use fixed particles if inflow/outflow at boundaries.
 xmin(1) = -0.5
 xmax(1) = 0.5
 xcentre = (xmax(1) + xmin(1))/2.0
!
!--setup parameters
!
 gam1 = gamma - 1.
 const = SQRT(4.*pi)
 dsmooth = 0. !!20. 
 densleft = 1.0
 densright = 1.0
 prleft = 1000.0
 prright = 0.01
 vleft(1) = 0.
 vright(1) = 0.
 vleft(2) = 0.
 vright(2) = 0.
 vleft(3) = 0.
 vright(3) = 0.
 IF (imhd.NE.0) THEN
    Bxinit = 0.75       !2./const
    Byleft = 1.0        !3.6/const
    Byright = -1.0      !4.0/const
    Bzleft = 0.         !2./const
    Bzright = 0.        !2./const
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
 WRITE(iprint,10) densleft,densright,prleft,prright,vleft(1),vright(1), &
                  vleft(2),vright(2),vleft(3),vright(3),                &
                  Bxinit,Byleft,Byright,Bzleft,Bzright
 
10 FORMAT( ' 1D shock: dens L: ',f16.9,' R: ',f16.9,/,   &
           '           pr  L: ',f16.9,' R: ',f16.9,/,   &
           '           vx  L: ',f16.9,' R: ',f16.9,/,   &
           '           vy  L: ',f16.9,' R: ',f16.9,/,   &
           '           vz  L: ',f16.9,' R: ',f16.9,/,   &
           '           Bx   : ',f16.9,/,                &
           '           By  L: ',f16.9,' R: ',f16.9,/,   &
           '           Bz  L: ',f16.9,' R: ',f16.9,/)



 IF (ABS(gam1).GT.1.e-3) THEN
    uuleft = prleft/(gam1*densleft)
    uuright = prright/(gam1*densright)
 ELSE
    uuleft = 3.*prleft/(2.*densleft)
    uuright = 3.*prright/(2.*densright)
 ENDIF
 massp = densright*psep      
 psepleft = massp/densleft
!
! allocate memory to start with
!
 ntot = (xcentre - xmin(1))/psepleft + (xmax(1) - xcentre)/psep
 npartleft = int((xcentre - xmin(1))/psepleft)
 npartright = int((xmax(1) - xcentre)/psep)
 print*,'npartleft = ',npartleft,' npartright = ',npartright
 CALL alloc(ntot)
!
!--setup the particles using the parameters given
! 
 x(1,1) = xmin(1) + 0.5*psepleft
 x(1,2) = xmin(1) + 1.5*psepleft
 DO i=1,2
    dens(i) = densleft
    uu(i) = uuleft
    pmass(i) = massp
    vel(:,i) = vleft(:)
    IF (imhd.NE.0) THEN
       Bfield(1,i) = Bxinit
       Bfield(2,i) = Byleft
       Bfield(3,i) = Bzleft
    ENDIF 
 ENDDO
 
 j = 2     
 npart = npartleft + npartright
 DO WHILE (j.LT.npart)
! DO j=1,npart
!    if (j.le.npartleft) then
!       x(1,j) = xcentre - (npartleft-j-1)*psepleft - 0.5*psepleft
!    else
!       x(1,j) = xcentre + (j-npartleft-1)*psep + 1.5*psep
!    endif
    dx0 = massp/dens(j)
    xhalf = x(1,j) + 0.5*dx0  !x at the mid point
    delta = xhalf/psep
    j = j + 1 
    IF (delta.lt.-dsmooth) THEN
       denshalf = densleft
    ELSEIF (delta.gt.dsmooth) THEN
       denshalf = densright
    ELSE
       exx = exp(delta)
       denshalf = (densleft+ densright*exx)/(1+exx)
    ENDIF
!---------------------------------------------------	    
!     calculate half steps then final step
!---------------------------------------------------	    
    dxhalf = massp/denshalf
    dx1 = 2.*dxhalf - dx0
    x(1,j) = xhalf + 0.5*dx1

    delta = (x(1,j)-xcentre)/psep
    IF (delta.lt.-dsmooth) THEN
       dens(j) = densleft
       uu(j) = uuleft
       vel(:,j) = vleft(:)
       Bfield(2,j) = Byleft
       Bfield(3,j) = Bzleft
    ELSEIF (delta.gt.dsmooth) THEN
       dens(j) = densright
       uu(j) = uuright
       vel(:,j) = vright(:)
       Bfield(2,j) = Byright
       Bfield(3,j) = Bzright
    ELSE
       exx = exp(delta)
       dens(j) = (densleft + densright*exx)/(1.0 +exx)
       uu(j) = (uuleft + uuright*exx)/(1.0 + exx)
!       IF (prleft.GT.prright) THEN
!          uu(j) = (prleft + prright*exx)/((1.0 + exx)*gam1*dens(i))
!       ELSE
!          uu(j) = (prright + prleft*exx)/((1.0 + exx)*gam1*dens(i))
!       ENDIF
!      vel(:,j) = (vleft(:) + vright(:)*exx)/(1.0 + exx)
       IF (delta.LT.0.) THEN
          vel(:,j) = vleft(:)
       ELSE
          vel(:,j) = vright(:)
       ENDIF  
       Bfield(2,j) = (Byleft + Byright*exx)/(1.0 + exx)
       Bfield(3,j) = (Bzleft + Bzright*exx)/(1.0 + exx)
    ENDIF
    pmass(j) = massp
    Bfield(1,j) = Bxinit
 ENDDO

!
!--constant component of force for subtraction 
!
 Bconst(:) = 0.
 Bconst(1) = Bxinit
!
!--setup vector potential if necessary
!
 if (imhd.lt.0) then
    do i=1,npart
       Bevol(:,i) = 0.
       Bevol(3,i) = -Bfield(2,i)*x(1,i)
       Bevol(2,i) = Bfield(3,i)*x(1,i)
    enddo
 endif
 !npart = j-1
 ntotal = npart
 print*,'end of setup, npart = ',npart
!
!--if using ghosts adjust outer boundary to lie half a particle spacing away
!
 IF (ibound(1).GT.1) THEN
    xmax(1) = x(1,npart) + 0.5*(x(1,npart)-x(1,npart-1))
 ENDIF

!
!--allow for tracing flow
!
 IF (trace) WRITE(iprint,*) '  Exiting subroutine setup'
    
 RETURN
END SUBROUTINE setup

!
! use this routine to modify the dump upon code restart
!
subroutine modify_dump()
 implicit none

end subroutine modify_dump
