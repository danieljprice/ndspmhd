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
 use dimen_mhd
 use debug
 use loguns
 use bound
 use eos
 use options
 use part
 use setup_params
 use timestep  ! uses tmax for moving boundaries
 
 use uniform_distributions
 implicit none
 integer :: i,j
 real :: densleft,densright,prleft,prright
 real :: uuleft, uuright
 real :: dx, const
 real :: massp,masspleft,masspright,Bxinit,Byleft,Byright,Bzleft,Bzright
 real :: vxleft, vxright, vyleft, vyright, vzleft, vzright
 real, dimension(ndim) :: sidelength,xminleft,xminright,xmaxleft,xmaxright
 real, dimension(ndim) :: xminleftleft,xmaxleftleft
 real :: boxlength, xshock, gam1, psepleft, psepright, psepleftleft
 real :: total_mass, volume, cs_L,cs2_L,cs2_R,mach_R,mach_L,vjump,gamm1
 real :: Alayer,mach_ratio,cs
 real :: densleftleft,vxleftleft,vyleftleft,vzleftleft
 real :: Bxleftleft,Byleftleft,Bzleftleft
 real :: vxi,densi,masspleftleft
 character(len=20) :: shkfile
 logical :: equalmass
!
!--allow for tracing flow
!
 if (trace) write(iprint,*) ' Entering subroutine setup'
!! if (ndim.lt.2) stop 'ndim needs to be > 1D for this setup'
 if (ndimV.lt.3) write(iprint,*) 'Warning: ndimV < 3: Bz not used'
!
!--set default values
!
 equalmass = .true.   ! use equal mass particles??
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
 vyleft = 0.0
 vyright = 0.
 vzleft = 0.
 vzright = 0.
 if (imhd.ne.0) then
    Bxinit = 2./const
    Byleft = 3.6/const
    Byright = 4./const
    Bzleft = 2./const
    Bzright = 2./const
 endif
!
!--now setup parameters specific to the acoustic-shock interaction problem
!
 xlayer = 0.225
 dwidthlayer = 1./0.05
 mach_ratio = 0.5

 Alayer = 0.5*(1./gamm1 + 0.5*mach_L**2 &
       - (1./gamm1 + 0.5*mach_R**2)*(mach_ratio)**(-2.*gamm1/(gamma + 1.)))
 Alayercs = Alayer*cs_L**2
!
!--get density after layer
!
 dx = (xmin(1)-xlayer)*dwidthlayer
 call get_profile(dx,vxi,cs,cs_L,mach_L,Alayer)
 densleftleft = vxleft*densleft/vxi
 print*,'Alayer = ',Alayer
 print*,'mach_L = ',mach_L
 print*,'machleftleft = ',mach_L*mach_ratio
 print*,'densleftleft = ',densleftleft
 print*,'csleftleft = ',cs
 print*,'vxleftleft = ',vxi
 print*,'machleftleft = ',vxi/cs
!
!--print setup parameters to the log file
!
 write(iprint,10) ndim,densleft,densright,prleft,prright,vxleft,vxright,   &
                  vyleft,vyright,vzleft,vzright
 if (imhd.ne.0) write(iprint,20) Bxinit,Byleft,Byright,Bzleft,Bzright

10 FORMAT(/,1x,i1,'D acoustic-shock: dens L: ',f8.3,' R: ',f8.3,/,   &
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
 if ((abs(vxleft).gt.1.e-4).or.(abs(vxright).gt.1.e-4)) then
    ibound(1) = 1               ! fixed x particles
 else
    ibound(1) = 2                ! reflecting in x
 endif
 if (ndim.ge.2) ibound(2:ndim) = 3        ! periodic in yz
 nbpts = 0                ! must use fixed particles if inflow/outflow at boundaries
 boxlength = 1.0
 sidelength(1) = 512.        ! relative dimensions of boundaries
 if (ndim.ge.2) sidelength(2:ndim) = 6.
 xmin(1) = -0.5
 xmax(1) = xmin(1) + boxlength
 if (ndim.ge.2) then
    xmin(2:ndim) = 0.0   !-0.5*boxlength*sidelength(2:ndim)/sidelength(1)
    xmax(2:ndim) = xmin(2:ndim) + 4.*psep !!!abs(xmin(2:ndim))
 endif
!
!--extend boundaries if inflow
!
 if (vxleft.gt.0.) then
    xmin(1) = xmin(1) - vxleft*tmax - 6.*psep
 endif
! if (vxright.lt.0.) then
!    xmax(1) = xmax(1) - vxright*tmax + 6.*psep
! endif
 
 xshock = 0.8 !!(xmax(1) + xmin(1))/2.0

!
!--now setup the shock
! 
!--set boundaries of regions to initially cover the whole domain
 xminleft(:) = xmin(:)
 xmaxleft(:) = xmax(:)
 xminright(:) = xmin(:)
 xmaxright(:) = xmax(:)
 xminleftleft(:) = xmin(:)
 xmaxleftleft(:) = xmax(:)
!
!--then divide the x axis into three regions halves at xshock
!
 xmaxleftleft(1) = xlayer
 xminleft(1) = xlayer
 xmaxleft(1) = xshock
 xminright(1) = xshock
!
!--particle separation is relative to middle region
!
 psepleft = psep
 psepright = psep*(densleft/densright)**(1./ndim)
 psepleftleft = psep*(densleft/densleftleft)**(1./ndim)

 if (abs(densleft-densright).gt.1.e-6 .and. equalmass) then
    print*,' left region  ',xminleftleft,' to ',xmaxleftleft,' psepleft = ',psepleftleft
    print*,' middle region  ',xminleft,' to ',xmaxleft,' psepleft = ',psepleft
    print*,' right region ',xminright,' to ',xmaxright,' psepright = ',psepright
!!    massp = (psep**ndim)*densright

    call set_uniform_cartesian(2,psepleftleft,xminleftleft,xmaxleftleft,.false.)  ! set left half
    xmin = xminleftleft
!
!--particle volume is relative to middle region
!
    call set_uniform_cartesian(2,psepleft,xminleft,xmaxleft,.false.)  ! set left half
    volume = PRODUCT(xmaxleft-xminleft)
    total_mass = volume*densleft
    massp = total_mass/npart
    masspleftleft = massp
    masspleft = massp
    masspright = massp
 
    call set_uniform_cartesian(2,psepright,xminright,xmaxright,.false.) ! set right half
    xmax = xmaxright
 else  ! set all of volume if densities are equal
    call set_uniform_cartesian(1,psep,xmin,xmax,.false.)
    volume = PRODUCT(xmax-xmin)
!    vol_left = PRODUCT(xmaxleft-xminleft)
    masspleft = densleft*volume/REAL(npart)
!    vol_right = PRODUCT(xmaxright-xminright)
    masspright = densright*volume/REAL(npart)
    masspleftleft = densleftleft*volume/REAL(npart)
 endif

 print*,'npart = ',npart
!
!--if using moving boundaries, fix the particles near the boundaries
!
 nbpts = 0
 if (ibound(1).eq.1) then
    do i=1,npart
       if ((x(1,i).lt.(xmin(1) + 4.*psepleft)).or. &
           (x(1,i).gt.(xmax(1) - 4.*psepright))) then
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
    if (x(1,i).GT.xminright(1)) then
       dens(i) = densright
       uu(i) = uuright 
       vel(1,i) = vxright
       pmass(i) = masspright
       if (ndimV.ge.2) vel(2,i) = vyright
       if (ndimV.ge.3) vel(3,i) = vzright 
       if (ndimV.ge.2) Bfield(2,i) = Byright
       if (ndimV.ge.3) Bfield(3,i) = Bzright
    elseif (x(1,i).gt.xminleft(1)) then
       dens(i) = densleft
       pmass(i) = masspleft
       uu(i) = uuleft
       vel(1,i) = vxleft
       if (ndimV.ge.2) vel(2,i) = vyleft
       if (ndimV.ge.3) vel(3,i) = vzleft
       if (ndimV.ge.2) Bfield(2,i) = Byleft
       if (ndimV.ge.3) Bfield(3,i) = Bzleft
    elseif (x(1,i).gt.xminleftleft(1)) then
       dx = (x(1,i) - xlayer)*dwidthlayer
       call get_profile(dx,vxi,cs,cs_L,mach_L,Alayer)
       densi = vxleft*densleft/vxi
       dens(i) = densi
       pmass(i) = masspleftleft
       uu(i) = cs**2/(gamma*gamm1)
       vel(1,i) = vxleftleft
       if (ndimV.ge.2) vel(2,i) = vyleftleft
       if (ndimV.ge.3) vel(3,i) = vzleftleft
       if (ndimV.ge.2) Bfield(2,i) = Byleftleft
       if (ndimV.ge.3) Bfield(3,i) = Bzleftleft    
    endif           
    Bfield(1,i) = Bxinit
 enddo
!
!--setup const component of mag field which can be subtracted
!
 Bconst(:) = 0.
 Bconst(1) = Bxinit

 return
 
contains

subroutine get_profile(dx,velx,cs,cs_L,mach_L,Alayer)
 implicit none
 real, intent(in) :: dx,mach_L,cs_L,Alayer
 real, intent(out) :: velx,cs
 integer, parameter :: maxits = 1000
 integer :: its
 real, parameter :: tol = 1.e-10
 real :: machno,machmin,machmax,test
 
 machmin = 0.
 machmax = 1.0
 its = 0
!
!--iterate using bisection method
!
 do while (abs(1.-machmin/machmax).gt.tol .and. its.le.maxits)
    machno = 0.5*(machmin + machmax)
    test = (machno/mach_L)**(2.*(gamma-1.)/(gamma+1.))* &
           (Alayer*(TANH(dx)-1.) + 1./(gamma-1.) + 0.5*mach_L**2) &
           -0.5*machno**2 - 1./(gamma-1.)
    if (test.lt.0.) then
       machmin = machno
    else
       machmax = machno
    endif
    its = its + 1
 enddo
 
 if (its.eq.maxits) STOP 'ERROR: mach number not converged in setup'
 
 print*,'machno = ',machno
 machno = 0.5*(machmin + machmax)
 cs = cs_L*(mach_L/machno)**((gamma-1.)/(gamma + 1.))
 velx = -machno*cs
 
end subroutine get_profile

end subroutine setup
