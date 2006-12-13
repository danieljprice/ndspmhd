!!------------------------------------------------------------------------!!
!!  Generic setup for ND shock tubes in rectangular boxes                 !!
!!  This version does equal mass particles so that the particle           !!
!!  spacing changes if densleft .ne. densright. Does this by setting up   !!
!!  a uniform cartesian distribution in each half of the tube. This means !!
!!  that the density is *not* smoothed at the interface and furthermore   !!
!!  that the density will be *wrong* unless                               !!
!!  psepright = (densleft/densright)**(1./ndim)                           !!
!!  is an exact division of the box width                                 !!
!!                                                                        !!
!!  NOTE: If this setup is used in 1D must use boundaryND.f90             !!
!!        NOT boundaryND_1D.f90 (ie. no inflow/outflow)                   !!
!!------------------------------------------------------------------------!!
subroutine setup
 use dimen_mhd, only:ndim,ndimV
 use debug, only:trace
 use loguns, only:rootname,ireadf,iprint
 use bound, only:xmin,xmax
 use eos
 use options
 use part
 use setup_params
 use timestep, only:tmax
 
 use uniform_distributions
 implicit none
 integer :: i,npartold,nparty
 real :: densleft,densright,prleft,prright
 real :: uuleft, uuright
 real :: dsmooth, exx, delta, const
 real :: massp,masspleft,masspright,Bxinit,Byleft,Byright,Bzleft,Bzright
 real :: vxleft, vxright, vyleft, vyright, vzleft, vzright
 real, dimension(ndim) :: xminleft,xminright,xmaxleft,xmaxright
 real :: boxlength, xshock, gam1, psepleft, psepright, psepleftx, pseprightx
 real :: total_mass, volume, cs_L,cs2_L,cs2_R,mach_R,mach_L,vjump,gamm1
 character(len=20) :: shkfile
 logical :: equalmass, stretchx
!
!--allow for tracing flow
!
 if (trace) write(iprint,*) ' Entering subroutine setup'
!! if (ndim.lt.2) stop 'ndim needs to be > 1D for this setup'
 if (ndimV.lt.3) write(iprint,*) 'Warning: ndimV < 3: Bz not used'
!
!--set default values
!
 dsmooth = 0.   
 equalmass = .true.   ! use equal mass particles??
 stretchx = .false.    ! stretch in x-direction only to give density contrast?
 const = sqrt(4.*pi)
 gamm1 = gamma - 1.
 
 mach_R = 5.
 mach_L = sqrt((2. + gamm1*mach_R**2)/(2.*gamma*mach_R**2 - gamm1))
 vjump = (2. + (gamm1)*mach_R**2)/((gamma + 1.)*mach_R**2)
 cs_L = 1.0
 cs2_L = cs_L**2
 cs2_R = cs2_L*((gamma+1.)**2*mach_R**2)/ &
        ((2.*gamma*mach_R**2 - gamm1)*(2. + gamm1*mach_R**2))
 densleft = 1.0
 densright = vjump*densleft
 prleft = cs2_L*densleft/gamma
 prright = cs2_R*densright/gamma
 vxright = -sqrt(mach_R**2*cs2_R)
 vxleft = vjump*vxright
 vyleft = 0.01
 vyright = 0.
 vzleft = 0.5
 vzright = 0.
 if (imhd.ne.0) then
    Bxinit = 2./const
    Byleft = 3.6/const
    Byright = 4./const
    Bzleft = 2./const
    Bzright = 2./const
 endif 
!
!--read shock parameters from the .shk file
!
 shkfile = rootname(1:LEN_TRIM(rootname))//'.shk'
 
 open(UNIT=ireadf,FILE=shkfile,STATUS='old',FORM='formatted',ERR=666)
   read(ireadf,*,ERR=667,END=667) densleft,densright
   read(ireadf,*,ERR=667,END=667) prleft,prright
   read(ireadf,*,ERR=667,END=667) vxleft,vxright
   read(ireadf,*,ERR=667,END=667) vyleft,vyright
   read(ireadf,*,ERR=667,END=667) vzleft,vzright
   if (imhd.ne.0) then
      read(ireadf,*,ERR=667,END=667) Bxinit
      read(ireadf,*,ERR=667,END=667) Byleft,Byright
      read(ireadf,*,ERR=667,END=667) Bzleft,Bzright
   endif
 close(UNIT=ireadf)
 goto 668
666 write(iprint,*) 'shock parameters file not found, using defaults...'
    write(iprint,*) 'Writing ',shkfile,' with initial left/right states'
    open(UNIT=ireadf,FILE=shkfile,STATUS='replace',FORM='formatted')
       write(ireadf,*) densleft,densright
       write(ireadf,*) prleft,prright
       write(ireadf,*) vxleft,vxright
       write(ireadf,*) vyleft,vyright
       write(ireadf,*) vzleft,vzright
       write(ireadf,*) Bxinit
       write(ireadf,*) Byleft,Byright
       write(ireadf,*) Bzleft,Bzright       
    close(UNIT=ireadf)
    

    goto 668
667 write(iprint,*) 'error in shock parameters file, using some defaults...'  
    close(UNIT=ireadf)
668 continue
!
!--print setup parameters to the log file
!
 write(iprint,10) ndim,densleft,densright,prleft,prright,vxleft,vxright,   &
                  vyleft,vyright,vzleft,vzright
 if (imhd.ne.0) write(iprint,20) Bxinit,Byleft,Byright,Bzleft,Bzright

10 FORMAT(/,1x,i1,'D shock: dens L: ',f8.3,' R: ',f8.3,/,   &
           '           pr  L: ',f8.3,' R: ',f8.3,/,   &
           '           vx  L: ',f8.3,' R: ',f8.3,/,   &
           '           vy  L: ',f8.3,' R: ',f8.3,/,   &
           '           vz  L: ',f8.3,' R: ',f8.3)
20 FORMAT( '           Bx   : ',f8.3,/,   &                                 
           '           By  L: ',f8.3,' R: ',f8.3,/,   &
           '           Bz  L: ',f8.3,' R: ',f8.3,/)
!
!--set boundaries
!                        
 if ((abs(vxleft).gt.tiny(vxleft)).or.(abs(vxright).gt.tiny(vxright))) then
    ibound(1) = 1               ! fixed x particles
 else
    ibound(1) = 1                ! reflecting in x
 endif
 if (ndim.ge.2) ibound(2:ndim) = 3        ! periodic in yz
 nbpts = 0                ! must use fixed particles if inflow/outflow at boundaries
 boxlength = 1.0
 xmin(1) = -0.5
 xmax(1) = xmin(1) + boxlength
 nparty = 6
 if (ndim.ge.2) then
    xmin(2:ndim) = 0.0
    xmax(2:ndim) = xmin(2:ndim) + nparty*psep !!!abs(xmin(2:ndim))
 endif
!
!--extend boundaries if inflow
!
 if (vxleft.gt.tiny(vxleft)) then
    xmin(1) = xmin(1) - vxleft*tmax - 6.*psep
 else
    xmin(1) = xmin(1) - 6.*psep
 endif
 if (vxright.lt.-tiny(vxright)) then
    xmax(1) = xmax(1) - vxright*tmax + 6.*psep
 else
    xmax(1) = xmax(1) + 6.*psep
 endif
 
 xshock = 0.0 !!(xmax(1) + xmin(1))/2.0

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
 psepleft = psep
 psepright = psep*(densleft/densright)**(1./ndim)

 if (abs(densleft-densright).gt.1.e-6 .and. equalmass) then
    if (stretchx .and. ndim.ge.2) then
       !--this makes the density jump by stretching in the x direction only
       psepleftx = psep
       pseprightx = psep*(densleft/densright)
       print*,' left half  ',xminleft,' to ',xmaxleft,' psepleft = ',psepleft
       print*,' right half ',xminright,' to ',xmaxright,' psepright = ',psepright
       call set_uniform_cartesian(1,psep,xminleft,xmaxleft,psepx=psepleftx)  ! set left half
       xmin = xminleft
       xminleft(1) = xshock - npart/nparty*psep !!! = xminleft
       xmin(1) = xshock - npart/nparty*psep !!! = xminleft
       volume = PRODUCT(xmaxleft-xminleft)
       total_mass = volume*densleft
       massp = total_mass/npart
       print*,'particle mass = ',massp
       masspleft = massp
       masspright = massp

       npartold = npart
       call set_uniform_cartesian(1,psep,xminright,xmaxright,psepx=pseprightx) ! set right half
       xmax = xmaxright
       xmax(1) = xshock + (npart-npartold)/nparty*pseprightx 
    else
       !--this sets up a square lattice    
       print*,' left half  ',xminleft,' to ',xmaxleft,' psepleft = ',psepleft
       print*,' right half ',xminright,' to ',xmaxright,' psepright = ',psepright
   !!    massp = (psep**ndim)*densright
       call set_uniform_cartesian(1,psepleft,xminleft,xmaxleft,fill=.true.)  ! set left half
       xmin = xminleft
       volume = PRODUCT(xmaxleft-xminleft)
       total_mass = volume*densleft
       massp = total_mass/npart
       masspleft = massp
       masspright = massp

       call set_uniform_cartesian(1,psepright,xminright,xmaxright,fill=.true.) ! set right half
       xmax = xmaxright
    endif
 else  ! set all of volume if densities are equal
    call set_uniform_cartesian(2,psep,xmin,xmax,adjustbound=.true.)
    volume = PRODUCT(xmax-xmin)
!    vol_left = PRODUCT(xmaxleft-xminleft)
    masspleft = densleft*volume/REAL(npart)
!    vol_right = PRODUCT(xmaxright-xminright)
    masspright = densright*volume/REAL(npart)
 endif

 print*,'npart = ',npart
!
!--if using moving boundaries, fix the particles near the boundaries
!
 nbpts = 0
 if (ibound(1).eq.1) then
    do i=1,npart
       if ((x(1,i).lt.(xmin(1) + 2.*hfact*psepleft)).or. &
           (x(1,i).gt.(xmax(1) - 2.*hfact*psepright))) then
          itype(i) = 1
          nbpts = nbpts + 1
       endif
    enddo
 endif   
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
       vel(1,i) = vxright
       pmass(i) = masspright
       if (ndimV.ge.2) vel(2,i) = vyright
       if (ndimV.ge.3) vel(3,i) = vzright
       if (ndimV.ge.2) Bfield(2,i) = Byright
       if (ndimV.ge.3) Bfield(3,i) = Bzright
    elseif (delta.LT.-dsmooth) then
       dens(i) = densleft
       pmass(i) = masspleft
       uu(i) = uuleft
       vel(1,i) = vxleft
       if (ndimV.ge.2) vel(2,i) = vyleft
       if (ndimV.ge.3) vel(3,i) = vzleft
       if (ndimV.ge.2) Bfield(2,i) = Byleft
       if (ndimV.ge.3) Bfield(3,i) = Bzleft
    else
       exx = exp(delta)       
       dens(i) = (densleft + densright*exx)/(1.0 +exx)
       pmass(i) = (masspleft + masspright*exx)/(1.0 + exx)
!       uu(i) = (uuleft + uuright*exx)/(1.0 + exx)
       uu(i) = (prleft + prright*exx)/((1.0 + exx)*gam1*dens(i))
!       vel(1,i) = (vxleft + vxright*exx)/(1.0 + exx)
!       if (ndimV.ge.2) vel(2,i) = (vyleft + vyright*exx)/(1.0 + exx)
!       if (ndimV.ge.3) vel(3,i) = (vzleft + vzright*exx)/(1.0 + exx)
       if (delta.GT.0.) THEN
          vel(1,i) = vxright
          if (ndimV.ge.2) vel(2,i) = vyright
          if (ndimV.ge.3) vel(3,i) = vzright 
       else
          vel(1,i) = vxleft
          if (ndimV.ge.2) vel(2,i) = vyleft
          if (ndimV.ge.3) vel(3,i) = vzleft       
       endif
       if (ndimV.ge.2) Bfield(2,i) = (Byleft + Byright*exx)/(1.0 + exx)
       if (ndimV.ge.3) Bfield(3,i) = (Bzleft + Bzright*exx)/(1.0 + exx)      
    endif           
    Bfield(1,i) = Bxinit
 enddo
!
!--setup const component of mag field which can be subtracted
!
 Bconst(:) = 0.
 Bconst(1) = Bxinit
! if (abs(Byleft - Byright).lt.tiny(Byleft)) then
!    Bconst(2) = Byleft
!    write(iprint,*) 'treating constant y-field as external'
! endif
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
!
!--if no smoothing applied, get rho from a sum and then set u to give a
!  smooth pressure jump (no spikes)
!
! if (abs(dsmooth).lt.tiny(dsmooth) .and. abs(gam1).gt.1.e-3) then
 if (.false.) then
    write(iprint,*) 'calling density to make smooth pressure jump...'
    call primitive2conservative
    do i=1,npart
       if (x(1,i).le.xshock) then
          uu(i) = prleft/(gam1*dens(i))
       else
          uu(i) = prright/(gam1*dens(i))
       endif
    enddo
 endif 
 
 return
end subroutine setup
