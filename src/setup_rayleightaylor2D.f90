!----------------------------------------------------------------
!     Set up a Rayleigh-Taylor instability test as
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
 use random,    only:ran1
!
!--define local variables
!            
 implicit none
 integer :: i,iseed,ipart
 real :: massp,volume,totmass,massmedium,Rfunc2,z
 real :: denszero,densmedium,przero,psepmedium,pri
 real, dimension(ndim) :: xminregion,xmaxregion
 logical, parameter :: equalmass = .true.
 logical, parameter :: perturbinterface = .true.

 Rfunc2(z) = -0.01*cos(6.*pi*z)  ! interface definition
!
!--allow for tracing flow
!
 if (trace) write(iprint,*) ' entering subroutine setup(unifdis)'
 if (ndim.ne.2) stop 'error: need ndim=2 for this problem'
!
!--set boundaries
! 	    
 ibound(1) = 2 ! periodic in x
 ibound(2) = 1 ! fixed in y
 nbpts = 0     ! use ghosts not fixed
 if (perturbinterface) then
    xmin(1) = 0. ! set position of boundaries
    xmax(1) = 1./6.
    xmin(2) = -0.5 !- 6.*psep	! set position of boundaries
    xmax(2) = 0.5 !+ 6.*psep 
 else
    xmin(1) = -0.25 ! set position of boundaries
    xmax(1) = 0.25
    xmin(2) = -0.75 !- 6.*psep	! set position of boundaries
    xmax(2) = 0.75 !+ 6.*psep
 endif
 if (iexternal_force.ne.8) stop 'need iexternal force = 8 for this problem'
!
!--set up the uniform density grid
!
 denszero = 1.0
 densmedium = 2.0
 przero = 2.5
 psepmedium = psep*(denszero/densmedium)**(1./ndim)
 massmedium = massp*(densmedium/denszero)
 volume = product(xmax(:)-xmin(:))
 totmass = denszero*volume

 if (.not.equalmass) then
!
!--unequal masses setup whole grid
!
    call set_uniform_cartesian(2,psep,xmin,xmax,fill=.true.)
    massp = totmass/float(npart)

 elseif (perturbinterface) then
!
!--setup whole domain to get particle mass
!
    call set_uniform_cartesian(2,psep,xminregion,xmaxregion,fill=.true.)
    massp = totmass/float(npart)
    npart = 0
!
!--for equal mass particles setup each region separately
!  (this uses a mask to get a perturbed interface)
!
    call set_uniform_cartesian(1,psep,xmin,xmax,fill=.true.,mask=-5)
    call set_uniform_cartesian(1,psepmedium,xmin,xmax,fill=.true.,mask=5)
    ntotal = npart

 else
!
!--for equal mass particles setup each region separately
!
    xminregion(1) = xmin(1)
    xmaxregion(1) = xmax(1) ! 0.
    xminregion(2) = xmin(2)
    xmaxregion(2) = 0.
    call set_uniform_cartesian(2,psep,xminregion,xmaxregion,fill=.true.)
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
    xminregion(2) = 0.
    xmaxregion(2) = xmax(2)
    call set_uniform_cartesian(2,psepmedium,xminregion,xmaxregion,fill=.true.)
!
!--reallocate memory to new size of list
!
!    call alloc(2*npart)
!
!--reflect particles above and below the y axis
!
    ipart = npart
!    do i=1,npart
!       ipart = ipart + 1
!       x(1,ipart) = -x(1,i)
!       x(2,ipart) = x(2,i)
!    enddo
    npart = ipart
    ntotal = npart
 endif
 npart = ntotal
 print*,'npart =',npart
!
!--now declare periodic boundaries
!
 xmin(2) = -0.75
 xmax(2) = 0.75
!
!--now assign particle properties
!
 do i=1,ntotal
    vel(:,i) = 0.
    if (perturbinterface) then
       if (x(2,i) .gt. Rfunc2(x(1,i))) then
          dens(i) = densmedium
          if (equalmass) then
             pmass(i) = massp
          else
             pmass(i) = massmedium
          endif
       else
          dens(i) = denszero
          pmass(i) = massp
       endif
    else
       if (abs(x(2,i)).lt.1./3.) then
          vel(2,i) = 0.01*(1. + cos(4.*pi*x(1,i)))*(1. + cos(3.*pi*x(2,i)))/4.
       else
          vel(2,i) = 0.
       endif
       if (x(2,i).gt.0.) then
          dens(i) = densmedium
          if (equalmass) then
             pmass(i) = massp
          else
             pmass(i) = massmedium
          endif
       else
          dens(i) = denszero
          pmass(i) = massp
       endif
    endif
    pri = przero - 0.1*dens(i)*x(2,i)
    uu(i) = pri/((gamma-1.)*dens(i))
    !--set fixed particles in y
    if (x(2,i).gt.xmax(2) .or. x(2,i).lt.xmin(2)) then
       nbpts = nbpts + 1
       itype(i) = 1
    endif
!    Bfield(:,i) = 0.
!    Bfield(1,i) = 0.5
 enddo
!
!--get rho from a sum and then set u to give a
!  smooth pressure
!
 write(iprint,*) 'calling density to make smooth pressure jump...'
 call primitive2conservative
 do i=1,npart
    pri = przero - 0.1*rho(i)*x(2,i)
    uu(i) = pri/((gamma-1.)*dens(i))
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
 write(iprint,*) 'modifying dump with velocities for Rayleigh-Taylor run'
 do i=1,ntotal
    vel(1,i) = 0.
    if (abs(x(2,i)).lt.1./3.) then
       vel(2,i) = 0.5*(1. + cos(4.*pi*x(1,i)))*(1. + cos(3.*pi*x(2,i)))/4.
    else
       vel(2,i) = 0.
    endif
 enddo
 time = 0.
 
end subroutine modify_dump
