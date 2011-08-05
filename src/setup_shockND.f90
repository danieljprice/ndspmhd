!!------------------------------------------------------------------------!!
!!  Generic setup for ND shock tubes in rectangular boxes                 !!
!!  This version does equal mass particles                                !!
!!------------------------------------------------------------------------!!
subroutine setup
 use dimen_mhd
 use debug
 use loguns
 use bound
 use eos
 use options
 use part
 use part_in
 use setup_params
 implicit none
 integer :: i,j
 real :: rholeft,rhoright,prleft,prright
 real :: uuleft, uuright
 real :: dsmooth, exx, delta
 real :: massp,Bxinit,Byleft,Byright,Bzleft,Bzright
 real, dimension(ndimV) :: vleft, vright
 real, dimension(ndim) :: sidelength,xminleft,xminright,xmaxleft,xmaxright
 real :: boxlength, xshock, gam1, psepleft, psepright
 real :: total_mass, volume
 character(len=20) :: shkfile
!
!--allow for tracing flow
!
 IF (trace) WRITE(iprint,*) ' Entering subroutine setup'
!
!--set boundaries
!            	    
 ibound(1) = 2		! reflecting in x
 ibound(2:3) = 3	! periodic in yz
 nbpts = 0		! must use fixed particles if inflow/outflow at boundaries
 boxlength = 1.0
 sidelength(1) = 120.	! relative dimensions of boundaries
 sidelength(2) = 6.
 sidelength(3) = 6.
 xmin(1) = -0.5
 xmax(1) = xmin(1) + boxlength
 xmin(2:3) = -0.5*boxlength*sidelength(2:3)/sidelength(1)
 xmax(2:3) = abs(xmin(2:3))
 
 xshock = (xmax(1) + xmin(1))/2.0
!
!--set default values
!
 dsmooth = 20.    
 rholeft = 4.0
 rhoright = 1.0
 prleft = 1.0
 prright = 0.1795
 vleft = 0.
 vright = 0.
 if (imhd.ne.0) then
    Bxinit = 0.
    Byleft = 0.
    Byright = 0.
    Bzleft = 0.
    Bzright = 0.
 endif
!
!--read shock parameters from the .shk file
!
 shkfile = rootname(1:LEN_TRIM(rootname))//'.shk'
 
 open(UNIT=ireadf,FILE=shkfile,STATUS='old',FORM='formatted',ERR=666)
   read(ireadf,*,ERR=667,END=667) rholeft,rhoright
   read(ireadf,*,ERR=667,END=667) prleft,prright
   read(ireadf,*,ERR=667,END=667) vleft(1),vright(1)
   read(ireadf,*,ERR=667,END=667) vleft(2),vright(2)
   read(ireadf,*,ERR=667,END=667) vleft(3),vright(3)
   if (imhd.ne.0) then
      read(ireadf,*,ERR=667,END=667) Bxinit
      read(ireadf,*,ERR=667,END=667) Byleft,Byright
      read(ireadf,*,ERR=667,END=667) Bzleft,Bzright
   endif
 close(UNIT=ireadf)
 goto 668
666 write(iprint,*) 'shock parameters file not found, using defaults...'
    goto 668
667 write(iprint,*) 'error in shock parameters file, using some defaults...'  
668 continue
!
!--print setup parameters to the log file
!
 write(iprint,10) ndim,rholeft,rhoright,prleft,prright,vleft(1),vright(1),   &
                  vleft(2),vright(2),vleft(3),vright(3)
 if (imhd.ne.0) write(iprint,20) Bxinit,Byleft,Byright,Bzleft,Bzright

10 FORMAT( 1x,i1,'D shock: rho L: ',f8.3,' R: ',f8.3,/,   &
           '           pr  L: ',f8.3,' R: ',f8.3,/,   &
           '           vx  L: ',f8.3,' R: ',f8.3,/,   &
           '           vy  L: ',f8.3,' R: ',f8.3,/,   &
           '           vz  L: ',f8.3,' R: ',f8.3)
20 FORMAT( '           Bx   : ',f8.3,/,   &	   	   	   
           '           By  L: ',f8.3,' R: ',f8.3,/,   &
           '           Bz  L: ',f8.3,' R: ',f8.3,/)
!
!--now setup the shock
! 
!--set boundaries of regions to initially cover the whole domain
 xminleft(:) = xmin(:)
 xmaxleft(:) = xmax(:)
 xminright(:) = xmin(:)
 xmaxright(:) = xmax(:)
!
!--then divide the x axis into two halves at xshock
!
 xmaxleft(1) = xshock
 xminright(1) = xshock
 
 print*,' left half  ',xminleft,' to ',xmaxleft
 print*,' right half ',xminright,' to ',xmaxright
! massp = (psep**ndim)*rhoright
 psepleft = psep
 psepright = psep*(rholeft/rhoright)**(1./ndim)
 
 psep = psepleft
 call set_uniform_cartesian(1,xminleft,xmaxleft,.false.)  ! set left half
 volume = PRODUCT(xmaxleft-xminleft)
 total_mass = volume*rholeft
 massp = total_mass/npart
 
 if (rholeft.ne.rhoright) then 
    psep = psepright
 endif
 call set_uniform_cartesian(1,xminright,xmaxright,.false.) ! set right half

 print*,'npart = ',npart
!
!--now set particle properties
!
 gam1 = gamma - 1.
 if (ABS(gam1).GT.1.e-3) then
    uuleft = prleft/(gam1*rholeft)
    uuright = prright/(gam1*rhoright)
 else
    uuleft = 3.*prleft/(2.*rholeft)
    uuright = 3.*prright/(2.*rhoright)
 endif

 do i=1,npart
    delta = (x(1,i) - xshock)/psep
    if (delta.GT.dsmooth) then
       rho(i) = rhoright
       uu(i) = uuright 
       vel(:,i) = vright(:) 
       Bfield(2,i) = Byright
       Bfield(3,i) = Bzright
    elseif (delta.LT.-dsmooth) then
       rho(i) = rholeft
       uu(i) = uuleft
       vel(:,i) = vleft(:)
       Bfield(2,i) = Byleft
       Bfield(3,i) = Bzleft
    else
       exx = exp(delta)       
       rho(i) = (rholeft + rhoright*exx)/(1.0 +exx)
!       uu(i) = (uuleft + uuright*exx)/(1.0 + exx)
       uu(i) = (prleft + prright*exx)/((1.0 + exx)*gam1*rho(i))
       if (delta.GT.0.) THEN
          vel(:,i) = vright(:)
       else
          vel(:,i) = vleft(:)
       endif
       Bfield(2,i) = (Byleft + Byright*exx)/(1.0 + exx)
       Bfield(3,i) = (Bzleft + Bzright*exx)/(1.0 + exx)      
    endif       
    pmass(i) = massp    
    hh(i) = hfact*(pmass(i)/rho(i))**hpower    
 enddo
 
 return
end subroutine setup
