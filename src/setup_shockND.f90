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
 use setup_params
 implicit none
 integer :: i,j
 real :: densleft,densright,prleft,prright
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
 densleft = 4.0
 densright = 1.0
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
   read(ireadf,*,ERR=667,END=667) densleft,densright
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
 write(iprint,10) ndim,densleft,densright,prleft,prright,vleft(1),vright(1),   &
                  vleft(2),vright(2),vleft(3),vright(3)
 if (imhd.ne.0) write(iprint,20) Bxinit,Byleft,Byright,Bzleft,Bzright

10 FORMAT( 1x,i1,'D shock: dens L: ',f8.3,' R: ',f8.3,/,   &
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
! massp = (psep**ndim)*densright
 psepleft = psep
 psepright = psep*(densleft/densright)**(1./ndim)
 
 call set_uniform_cartesian(1,psepleft,xminleft,xmaxleft,.false.)  ! set left half
 volume = PRODUCT(xmaxleft-xminleft)
 total_mass = volume*densleft
 massp = total_mass/npart
 
 call set_uniform_cartesian(1,psepright,xminright,xmaxright,.false.) ! set right half

 print*,'npart = ',npart
!
!--now set particle properties
!
 gam1 = gamma - 1.
 if (ABS(gam1).GT.1.e-3) then
    uuleft = prleft/(gam1*densleft)
    uuright = prright/(gam1*densright)
 else
    uuleft = 3.*prleft/(2.*densleft)
    uuright = 3.*prright/(2.*densright)
 endif

 do i=1,npart
    delta = (x(1,i) - xshock)/psep
    if (delta.GT.dsmooth) then
       dens(i) = densright
       uu(i) = uuright 
       vel(:,i) = vright(:) 
       Bfield(2,i) = Byright
       Bfield(3,i) = Bzright
    elseif (delta.LT.-dsmooth) then
       dens(i) = densleft
       uu(i) = uuleft
       vel(:,i) = vleft(:)
       Bfield(2,i) = Byleft
       Bfield(3,i) = Bzleft
    else
       exx = exp(delta)       
       dens(i) = (densleft + densright*exx)/(1.0 +exx)
!       uu(i) = (uuleft + uuright*exx)/(1.0 + exx)
       uu(i) = (prleft + prright*exx)/((1.0 + exx)*gam1*dens(i))
       if (delta.GT.0.) THEN
          vel(:,i) = vright(:)
       else
          vel(:,i) = vleft(:)
       endif
       Bfield(2,i) = (Byleft + Byright*exx)/(1.0 + exx)
       Bfield(3,i) = (Bzleft + Bzright*exx)/(1.0 + exx)      
    endif       
    pmass(i) = massp    
    hh(i) = hfact*(pmass(i)/dens(i))**hpower    
 enddo
 
 return
end subroutine setup
