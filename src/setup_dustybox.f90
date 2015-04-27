!----------------------------------------------------------------
!     Set up a box for a dustybox test
!     should work in 1, 2 and 3 dimensions
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
 use uniform_distributions, only:set_uniform_cartesian
 use mem_allocation, only:alloc
!
!--define local variables
!            
 implicit none
 integer :: i
 integer :: ntypes
 integer :: ngas,ndust,jtype
 real, dimension(ndimV) :: Bzero
 real :: massp,masspdust
 real :: spsoundi
 real :: denszero,uuzero,przero,denszerodust
 real :: dust_to_gas_ratio
!
!--setup parameters (could read in from a file)
!
 Bzero    = 0.
 denszero = 1.0
 uuzero   = 1.0
 if (idust.eq.1) then
    ntypes = 1
 else
    ntypes = 2
 endif
!
!--set boundaries
!
 ibound = 3        ! periodic boundaries
 nbpts = 0                ! use ghosts not fixed
 xmin(:) = 0.   ! set position of boundaries
 xmax(1) = 1.0 
 if (ndim.GE.2) then
    xmax(2:ndim) = 11.*psep ! would need to adjust this depending on grid setup
 endif
!
!--initially set up a uniform density grid (also determines npart)
!  (the call to set_uniform_cartesian means this works in 1,2 and 3D)
!
 ngas  = 0
 ndust = 0
 do jtype=1,ntypes
    call set_uniform_cartesian(1,psep,xmin,xmax,adjustbound=.true.)
    if (jtype.eq.1) then
       ngas = npart
       itype(1:ngas) = itypegas
    elseif (jtype.eq.2) then
       ndust = npart - ngas
       itype(ngas+1:ngas+ndust) = itypedust
    endif
 enddo

 massp = 1.0/FLOAT(ngas)        ! average particle mass
 
 dust_to_gas_ratio = 1.
 if (idust.eq.1) then ! one fluid dust
    massp = massp*(1. + dust_to_gas_ratio)
 endif

 ! two fluid dust (if ndust > 0)
 masspdust = 0.
 if (ndust.gt.0) masspdust = dust_to_gas_ratio*1.0/FLOAT(ndust) ! average particle mass
 denszerodust = dust_to_gas_ratio*denszero
 if (ntypes.gt.1) print*,' ngas = ',ngas,' ndust = ',ndust
!
!--allocate memory here
!
 call alloc(npart)
!
!--setup uniform density grid of particles
! 
 do i=1,npart
    vel(:,i) = 0.
    if (itype(i).eq.itypedust) then
       dens(i)  = denszerodust
       pmass(i) = masspdust
       uu(i)    = 0. 
    else
       if (idust.eq.1) then
          dustfrac(i) = dust_to_gas_ratio/(1. + dust_to_gas_ratio)
          deltav(1,i)  = -1. ! deltav = vdust - vgas
          vel(1,i) = 0.5     ! v = vg + rhod/rho*deltav
       else
          vel(1,i) = 1. !--gas velocity       
       endif
       dens(i)  = denszero
       pmass(i) = massp
       uu(i)    = uuzero
    endif
    if (imhd.GT.0) then 
       Bfield(:,i) = Bzero
    else
       Bfield(:,i) = 0.
    endif 
 enddo

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
