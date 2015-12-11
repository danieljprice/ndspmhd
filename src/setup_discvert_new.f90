!----------------------------------------------------------------
!     Set up a periodic box (no shear, radial motion neglected)
!     with a vertical linear force which mimics the vertical
!     component of the star's gravity (2D problem)
!     fixed density, variable cs
!----------------------------------------------------------------

subroutine setup
!
!--include relevant global variables
!
 use dimen_mhd
 use debug
 use loguns
 use bound
 use eos
 use options
 use part
 use setup_params
 use mem_allocation, only:alloc
 use uniform_distributions 
!
!--define local variables
!            
 implicit none
 integer :: i
 integer :: ntypes
 integer :: ngas,npdust,jtype
 integer :: iwhichpot,idust_toponly
 real, dimension(ndimV) :: Bzero
 real :: massp,masspdust,massp_without
 real :: spsoundi
 real :: denszero,uuzero,przero,denszerodust
 real :: dust_to_gas_ratio
 real :: x_min,x_max,z_min,z_max,z2_max,Lx,Hz,HH2,voltot
 real :: cs2,cs2max,zz2_part
 
 if (ndim.ne.2) stop 'error: this is a 2D problem and ndim.ne.2'
 
!--idust_toponly
!--idust_toponly = 0 ; dust everywhere
!--idust_toponly = 1 ; dust only in the top half midplane
 idust_toponly = 0
 if (idust_toponly.ne.0 .and. idust_toponly.ne.1) then
    print*,"discvert setup: this value for idust_toponly does not exists",idust_toponly
    stop
 endif
 if (idust_toponly.eq.1 .and. idust.eq.2) print*,"WARNING: dust in a half space implemented for a single fluid only"

!-setup physical quatities
 dust_to_gas_ratio = 0.01
 Bzero    = 0.
 denszero = 1.0
 Lx       = 0.5
 Hz       = 1.0
 cs2max   = 5.0
 HH2      = Hz*Hz
 if (idust.eq.1) then
    ntypes = 1
 else
    ntypes = 2
 endif

! 
!--set external forces
!--iwhichpot = 1 ; quadratic potential (linear force)
!--iwhichpot = 2 ; quartic potential (cubic force)
!--iwhichpot = 3 ; potential in z**3/2 (square root force)
! 
 iwhichpot = 1
 select case(iwhichpot)
 case(1)
    iexternal_force = 12
 case(2)
    iexternal_force = 13
 case(3)
    iexternal_force = 14 
 case default 
    print*,"discvert setup: this potential does not exist"
    stop
 end select
 
!
!--set boundaries
!

 ibound(1) = 3 ! periodic in x
 ibound(2) = 3 ! periodic in y (fixed sucks ; use the symmetry of the pb instead)
 nbpts  = 0
 x_min   = -Lx
 x_max   =  Lx
 z_min   = -2.0*Hz
 z_max   =  2.0*Hz
 z2_max  =  z_max*z_max
 xmin(1) =  x_min
 xmax(1) =  x_max
 xmin(2) =  z_min
 xmax(2) =  z_max
 
!
!--initially set up a uniform density grid (also determines npart)
!  (the call to set_uniform_cartesian means this works in 1,2 and 3D)
!

 ngas  = 0
 npdust = 0
 do jtype=1,ntypes
    call set_uniform_cartesian(1,psep,xmin,xmax,adjustbound=.true.)
    if (jtype.eq.1) then
       ngas = npart
       itype(1:ngas) = itypegas
    elseif (jtype.eq.2) then
       npdust = npart - ngas
       itype(ngas+1:ngas+npdust) = itypedust
    endif
 enddo
 
 voltot = (z_max - z_min)*(x_max - x_min)
 massp  = voltot/FLOAT(ngas)        ! average particle mass
 massp_without = massp
 
 if (idust.eq.1) then ! one fluid dust
    massp = massp*(1. + dust_to_gas_ratio)
 endif

 ! two fluid dust (if npdust > 0)
 masspdust = 0.
 if (npdust.gt.0) masspdust = dust_to_gas_ratio*voltot/FLOAT(npdust) ! average particle mass
 denszerodust = dust_to_gas_ratio*denszero
 if (ntypes.gt.1) print*,' ngas = ',ngas,' npdust = ',npdust
!
!--allocate memory here
!
 call alloc(npart)
!
!--setup uniform density grid of particles with variable cs
!

 do i=1,npart
    vel(:,i) = 0.
    zz2_part = x(2,i)*x(2,i)
    if (itype(i).eq.itypedust) then
       dens(i)  = denszerodust
       pmass(i) = masspdust
       uu(i)    = 0.
    else
       dens(i)  = denszero
       pmass(i) = massp        
       if (idust.eq.1) then
          dustfrac(1,i) = dust_to_gas_ratio/(1. + dust_to_gas_ratio)
          deltav(1,i)  = 0.
          if (idust_toponly.eq.1) then
             if (x(2,i).lt.0.) then
                pmass(i)     = massp_without
                dustfrac(1,i) = 0.
             endif
          endif      
       endif 
       select case(iwhichpot)
       case(1)
          cs2 = cs2max + gamma*0.5*(z2_max - zz2_part)
       case(2)
          cs2 = cs2max + gamma*0.25*(z2_max**2 - zz2_part**2)
       case(3)
          cs2 = cs2max + gamma*(2./3.)*(z2_max**0.75 - zz2_part**0.75)
       end select
       if (gamma.gt.1.) then
          uu(i) = cs2/(gamma - 1.)/gamma
       else
          uu(i) = 1.5*cs2
       endif
    endif
    if (imhd.GT.0) then 
       Bfield(:,i) = Bzero
    else
       Bfield(:,i) = 0.
    endif 
 ENDDO

 ntotal = npart
!
!--get sound speed from equation of state (want average sound speed, so
!  before the density is perturbed)
!
 call equation_of_state1(przero,spsoundi,uuzero,denszero)
 print*,' gamma = ',gamma
 print*,' pr = ',przero,' cs = ',spsoundi,' u = ',uu(1),' dens = ',dens(1)


 return
end subroutine setup

!
! use this routine to modify the dump upon code restart
!
subroutine modify_dump()
 implicit none

end subroutine modify_dump
