!----------------------------------------------------------------
!     Set up a Kelvin-Helmholtz instability test as
!     in Jim Stone et al. (Athena test suite)
!----------------------------------------------------------------

subroutine setup
!
!--include relevant global variables
!
 use dimen_mhd
 use debug
 use loguns
 use bound
 use options
 use part
 use setup_params, only:psep,pi
 use eos, only:gamma
 use mem_allocation, only:alloc
 
 use uniform_distributions
 use cons2prim, only:primitive2conservative
!
!--define local variables
!            
 implicit none
 integer :: i,iseed,ipart
 real :: massp,volume,totmass,ran1
 real :: denszero,densmedium,przero,psepmedium
 real, dimension(ndim) :: xminregion,xmaxregion
 logical, parameter :: equalmass = .true.
!
!--allow for tracing flow
!
 if (trace) write(iprint,*) ' entering subroutine setup(unifdis)'
 if (ndim.ne.2) stop 'error: need ndim=2 for this problem'
!
!--set boundaries
! 	    
 nbpts = 0	! use ghosts not fixed
 xmin(:) = -0.5	! set position of boundaries
 xmax(:) = 0.5
!
!--set up the uniform density grid
!
 denszero = 1.0
 densmedium = 2.0
 przero = 2.5
 if (.not.equalmass) then
!
!--unequal masses setup whole grid
! 
    call set_uniform_cartesian(2,psep,xmin,xmax,fill=.true.)
    volume = product(xmax(:)-xmin(:))
    totmass = denszero*volume
    massp = totmass/float(npart)
 else
!
!--for equal mass particles setup each region separately
!
    xminregion(1) = xmin(1)
    xmaxregion(1) = xmax(1)
    xminregion(2) = -0.5
    xmaxregion(2) = -0.25
    call set_uniform_cartesian(1,psep,xminregion,xmaxregion,fill=.true.)
!
!--determine particle mass from npart in first region
!
    volume = product(xmaxregion(:)-xminregion(:))
    totmass = denszero*volume
    massp = totmass/float(npart) ! average particle mass

!    xminregion(2) = 0.25
!    xmaxregion(2) = 0.5
!    call set_uniform_cartesian(1,psep,xminregion,xmaxregion,fill=.true.)
!
!--setup -0.25 < y < 0.25
! 
    xminregion(2) = -0.25
    xmaxregion(2) = 0.
    psepmedium = psep*(denszero/densmedium)**(1./ndim)
    call set_uniform_cartesian(1,psepmedium,xminregion,xmaxregion,fill=.true.)
!
!--reallocate memory to new size of list
!
    call alloc(2*npart)
!
!--reflect particles above and below the y axis
!
    ipart = npart
    do i=1,npart
       ipart = ipart + 1
       x(1,ipart) = x(1,i)
       x(2,ipart) = -x(2,i)
    enddo
    npart = ipart
    ntotal = npart
 endif
 npart = ntotal
 print*,'npart =',npart
!
!--now declare periodic boundaries
!
 ibound = 3	! boundaries
 iseed = -23864
!
!--now assign particle properties
! 
 do i=1,ntotal
    vel(:,i) = 0.
    if (abs(x(2,i)).lt.0.25) then
       vel(1,i) = 0.5 
       dens(i) = densmedium
       if (equalmass) then
          pmass(i) = massp
       else
          pmass(i) = massp*(densmedium/denszero)
       endif
    else
       vel(1,i) = -0.5
       dens(i) = denszero
       pmass(i) = massp
    endif
    uu(i) = przero/((gamma-1.)*dens(i))
    Bfield(:,i) = 0.
    Bfield(1,i) = 0.5
!
!--add random velocity perturbation
!
    !!vel(1,i) = vel(1,i)*(1.0 + 0.01*(ran1(iseed)-0.5))
    !!vel(2,i) = 0.01*(ran1(iseed)-0.5)
    if (abs(x(2,i)-0.25).lt.0.025) then
       vel(2,i) = 0.025*sin(-2.*pi*(x(1,i)+0.5)*6.)
    elseif (abs(x(2,i)+0.25).lt.0.025) then
       vel(2,i) = 0.025*sin(2.*pi*(x(1,i)+0.5)*6.)
    endif
 enddo
 
 !!omegakh = sqrt(densmedium/denszero)*1.0/(denszero + densmedium)
 print*,' tau_kh = ',sqrt(densmedium/denszero)*(1./6.)
!
!--get rho from a sum and then set u to give a
!  smooth pressure
!
 write(iprint,*) 'calling density to make smooth pressure jump...'
 call primitive2conservative
 do i=1,npart
    uu(i) = przero/((gamma - 1.)*rho(i))
 enddo
!
!--allow for tracing flow
!
 if (trace) write(iprint,*) '  exiting subroutine setup'
  
 return
end

subroutine modify_dump
 use loguns, only:iprint
 use part
 use options, only:imhd
 use timestep, only:time
 use setup_params, only:pi
 implicit none
 integer :: i
!
!--now assign particle properties
!
 write(iprint,*) 'modifying dump with velocities for Kelvin-Helmholtz run'
 do i=1,ntotal
    vel(:,i) = 0.
    if (abs(x(2,i)).lt.0.25) then
       vel(1,i) = 0.5 
    else
       vel(1,i) = -0.5
    endif
!
!--add random velocity perturbation
!
    !!vel(1,i) = vel(1,i)*(1.0 + 0.01*(ran1(iseed)-0.5))
    !!vel(2,i) = 0.01*(ran1(iseed)-0.5)
    if (abs(x(2,i)-0.25).lt.0.025 .or. abs(x(2,i)+0.25).lt.0.025) then
       vel(2,i) = 0.025*sin(2.*pi*(x(1,i)+0.5)*6.)
    endif
 enddo

 time = 0.
 
end subroutine modify_dump
